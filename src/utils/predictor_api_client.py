"""
Predictor API Client
====================

Multi-source API client for fetching in silico predictor scores from
external databases. Implements source priority fallback and graceful
degradation when sources are unavailable.

Implements a strict, validated caching layer:
- Cache → API → (optional user fallback) → no data
- All cached entries are validated before use
- Invalid/corrupted cache entries are rejected

Key Design Principles:
- All factual predictor data is fetched from external APIs
- No hardcoded predictor scores - only thresholds/weights are local
- Multiple sources supported with configurable priority
- Graceful degradation when predictors are unavailable

Supported Sources:
- dbNSFP (via myvariant.info): Most comprehensive source
- AlphaMissense API: Official Google DeepMind API
- CADD API: Official CADD scoring service
- VEP (Ensembl): Variant Effect Predictor annotations

Author: Can Sevilmiş
License: MIT License
"""

import requests
from typing import Optional, Any, Union
from dataclasses import dataclass

try:
    from colorama import Fore, Style
except ImportError:
    class Fore:
        RED = YELLOW = GREEN = CYAN = MAGENTA = ''
    class Style:
        RESET_ALL = ''

from config.predictors import (
    PredictorScore,
    PopulationStats,
    PREDICTOR_SOURCE_PRIORITY,
    POPULATION_SOURCE_PRIORITY,
    INVERTED_PREDICTORS,
    INSILICO_WEIGHTS,
    validate_predictor_score,
    validate_population_stats,
    validate_cached_predictor_data,
    validate_cached_population_data,
)
from config.constants import API_SETTINGS

# Import validated cache (optional - falls back to dict if unavailable)
try:
    from utils.cache import (
        ResultCache, CacheKey,
        build_predictor_cache_key, build_population_cache_key
    )
    CACHE_AVAILABLE = True
except ImportError:
    CACHE_AVAILABLE = False
    ResultCache = None


# =============================================================================
# API Endpoints Configuration
# =============================================================================

PREDICTOR_API_ENDPOINTS = {
    # myvariant.info aggregates dbNSFP and other sources
    'myvariant': 'https://myvariant.info/v1/variant',
    
    # AlphaMissense Google Cloud API (requires setup)
    'alphamissense': 'https://alphamissense.hegelab.org/api/variant',
    
    # CADD scoring service
    'cadd': 'https://cadd.gs.washington.edu/api/v1.0/scores',
    
    # Ensembl VEP REST API
    'vep': 'https://rest.ensembl.org/vep/human/hgvs',
}


class PredictorAPIClient:
    """
    Multi-source API client for in silico predictor scores.
    
    This client fetches predictor scores from multiple external databases
    with configurable priority and automatic fallback.
    
    Implements validated caching:
    - All cached entries are validated before use
    - Invalid entries are rejected and re-fetched from API
    - Corrupted cache files are handled gracefully
    
    Usage:
        cache = ResultCache(cache_dir='./api_cache')
        client = PredictorAPIClient(api_enabled=True, result_cache=cache)
        scores = client.get_predictor_scores(
            chrom='17', pos=43092919, ref='G', alt='A'
        )
        # Returns: {'revel': PredictorScore(...), 'cadd_phred': PredictorScore(...), ...}
    """
    
    def __init__(
        self,
        api_enabled: bool = True,
        timeout: int = 15,
        cache: Optional[dict] = None,
        test_mode: bool = False,
        result_cache: Optional['ResultCache'] = None
    ):
        """
        Initialize the predictor API client.
        
        Args:
            api_enabled: Whether to make actual API calls
            timeout: Request timeout in seconds
            cache: Optional cache dictionary for API responses (legacy)
            test_mode: If True, return mock data instead of API calls
            result_cache: Optional ResultCache instance for validated caching
        """
        self.api_enabled = api_enabled
        self.timeout = timeout
        self.test_mode = test_mode
        
        # Use ResultCache if provided, otherwise fall back to simple dict
        if result_cache is not None and CACHE_AVAILABLE:
            self.result_cache = result_cache
            self.cache = {}  # Keep simple cache for compatibility
            self._use_validated_cache = True
        else:
            self.result_cache = None
            self.cache = cache if cache is not None else {}
            self._use_validated_cache = False
        
        # Track which sources are available/healthy
        self._source_health = {}
    
    def _get_cached_predictor_data(
        self,
        source: str,
        chrom: str,
        pos: int,
        ref: str,
        alt: str
    ) -> Optional[dict]:
        """
        Get validated cached predictor data.
        
        Returns None if:
        - Cache miss
        - Cached data is invalid (rejected and invalidated)
        - Cache is disabled
        """
        if not self._use_validated_cache or self.result_cache is None:
            return None
        
        cache_key = build_predictor_cache_key(source, chrom, pos, ref, alt)
        cached_data = self.result_cache.get(cache_key)
        
        if cached_data is None:
            return None
        
        # Validate cached data
        is_valid, errors = validate_cached_predictor_data(cached_data)
        if not is_valid:
            print(f"{Fore.YELLOW}⚠️  Invalid cached {source} data - invalidating: {errors}{Style.RESET_ALL}")
            self.result_cache.invalidate(cache_key)
            return None
        
        return cached_data
    
    def _set_cached_predictor_data(
        self,
        source: str,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        data: dict
    ) -> None:
        """
        Set validated predictor data in cache.
        
        Only caches valid data; invalid data is rejected.
        """
        if not self._use_validated_cache or self.result_cache is None:
            return
        
        # Validate before caching
        is_valid, errors = validate_cached_predictor_data(data)
        if not is_valid:
            print(f"{Fore.YELLOW}⚠️  Not caching invalid {source} data: {errors}{Style.RESET_ALL}")
            return
        
        cache_key = build_predictor_cache_key(source, chrom, pos, ref, alt)
        self.result_cache.set(cache_key, data, validated=True)
    
    def get_predictor_scores(
        self,
        chrom: Optional[str] = None,
        pos: Optional[int] = None,
        ref: Optional[str] = None,
        alt: Optional[str] = None,
        gene: Optional[str] = None,
        hgvs_c: Optional[str] = None,
        transcript: Optional[str] = None
    ) -> dict[str, PredictorScore]:
        """
        Fetch all available predictor scores for a variant.
        
        Args:
            chrom: Chromosome (e.g., '17', 'X')
            pos: Genomic position (GRCh38)
            ref: Reference allele
            alt: Alternate allele
            gene: Gene symbol (optional, for context)
            hgvs_c: HGVS coding notation (optional alternative to coordinates)
            transcript: Transcript ID (optional, for HGVS)
            
        Returns:
            Dictionary mapping predictor names to PredictorScore objects.
            Missing predictors will have value=None.
        """
        if not self.api_enabled:
            return self._get_empty_scores()
        
        # Initialize result with empty scores for all configured predictors
        results = {}
        for predictor in INSILICO_WEIGHTS.keys():
            results[predictor] = PredictorScore(
                predictor=predictor,
                value=None,
                source=None,
                is_inverted=(predictor in INVERTED_PREDICTORS)
            )
        
        if self.test_mode:
            return self._get_mock_scores()
        
        # Primary source: myvariant.info (aggregates dbNSFP)
        if chrom and pos and ref and alt:
            myvariant_scores = self._fetch_from_myvariant(chrom, pos, ref, alt)
            self._merge_scores(results, myvariant_scores)
        
        # Secondary source: AlphaMissense API for alphamissense specifically
        if 'alphamissense' in results and results['alphamissense'].value is None:
            if chrom and pos and ref and alt:
                am_score = self._fetch_alphamissense(chrom, pos, ref, alt)
                if am_score:
                    results['alphamissense'] = am_score
        
        # Tertiary source: CADD API for cadd_phred specifically
        if 'cadd_phred' in results and results['cadd_phred'].value is None:
            if chrom and pos and ref and alt:
                cadd_score = self._fetch_cadd(chrom, pos, ref, alt)
                if cadd_score:
                    results['cadd_phred'] = cadd_score
        
        return results
    
    def _fetch_from_myvariant(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str
    ) -> dict[str, PredictorScore]:
        """
        Fetch predictor scores from myvariant.info (aggregates dbNSFP).
        
        myvariant.info provides access to dbNSFP which contains:
        - REVEL, CADD, AlphaMissense, SIFT, PolyPhen-2, MetaSVM, VEST4, FATHMM, etc.
        """
        scores = {}
        
        # Build variant ID in HGVS-like format for myvariant.info
        # Format: chr{chrom}:g.{pos}{ref}>{alt}
        variant_id = f"chr{chrom}:g.{pos}{ref}>{alt}"
        legacy_cache_key = f"myvariant_{variant_id}"
        
        # Check validated cache first (if available)
        if self._use_validated_cache:
            cached_data = self._get_cached_predictor_data('dbNSFP', chrom, pos, ref, alt)
            if cached_data:
                return self._parse_myvariant_response(cached_data)
        elif legacy_cache_key in self.cache:
            # Legacy simple dict cache
            return self._parse_myvariant_response(self.cache[legacy_cache_key])
        
        try:
            print(f"{Fore.YELLOW}🔍 Querying myvariant.info for predictor scores: {variant_id}...{Style.RESET_ALL}")
            
            # Request dbNSFP fields
            params = {
                'fields': 'dbnsfp.revel.score,dbnsfp.cadd.phred,dbnsfp.alphamissense.am_pathogenicity,'
                         'dbnsfp.sift.score,dbnsfp.polyphen2.hdiv.score,dbnsfp.metasvm.score,'
                         'dbnsfp.vest4.score,dbnsfp.fathmm.score,dbnsfp.bayesdel.addaf_score,'
                         'dbnsfp.primateai.score,dbnsfp.mpc.score',
            }
            
            response = requests.get(
                f"{PREDICTOR_API_ENDPOINTS['myvariant']}/{variant_id}",
                params=params,
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                data = response.json()
                
                # Parse scores first to extract just the score values for caching
                scores = self._parse_myvariant_response(data)
                
                # Cache validated score values
                if self._use_validated_cache:
                    cache_data = {
                        name: score.value
                        for name, score in scores.items()
                        if score.value is not None
                    }
                    self._set_cached_predictor_data('dbNSFP', chrom, pos, ref, alt, cache_data)
                else:
                    self.cache[legacy_cache_key] = data
                
                available = sum(1 for s in scores.values() if s.value is not None)
                print(f"{Fore.GREEN}✅ myvariant.info: Retrieved {available} predictor scores{Style.RESET_ALL}")
            elif response.status_code == 404:
                print(f"{Fore.YELLOW}⚠️  Variant not found in myvariant.info{Style.RESET_ALL}")
            else:
                print(f"{Fore.RED}❌ myvariant.info error: HTTP {response.status_code}{Style.RESET_ALL}")
                
        except requests.exceptions.Timeout:
            print(f"{Fore.RED}❌ myvariant.info timeout{Style.RESET_ALL}")
        except Exception as e:
            print(f"{Fore.RED}❌ myvariant.info error: {str(e)}{Style.RESET_ALL}")
        
        return scores
    
    def _parse_myvariant_response(self, data: dict) -> dict[str, PredictorScore]:
        """
        Parse myvariant.info response into PredictorScore objects.
        
        Also handles cached score dictionaries (where keys are predictor names
        and values are floats directly).
        """
        scores = {}
        
        # Check if this is a simplified cached format (predictor -> value)
        # vs raw API response (with nested dbnsfp structure)
        if 'dbnsfp' not in data:
            # This is a cached simplified format
            for predictor, value in data.items():
                if value is not None and validate_predictor_score(predictor, value):
                    scores[predictor] = PredictorScore(
                        predictor=predictor,
                        value=float(value),
                        source='dbNSFP',
                        version='4.x_cached',
                        is_inverted=(predictor in INVERTED_PREDICTORS)
                    )
            return scores
        
        dbnsfp = data.get('dbnsfp', {})
        if not dbnsfp:
            return scores
        
        # Extract scores from dbNSFP nested structure
        score_mappings = {
            'revel': ('revel', 'score'),
            'cadd_phred': ('cadd', 'phred'),
            'alphamissense': ('alphamissense', 'am_pathogenicity'),
            'sift': ('sift', 'score'),
            'polyphen2': ('polyphen2', 'hdiv', 'score'),
            'metasvm': ('metasvm', 'score'),
            'vest4': ('vest4', 'score'),
            'fathmm': ('fathmm', 'score'),
            'bayesdel': ('bayesdel', 'addaf_score'),
            'primateai': ('primateai', 'score'),
            'mpc': ('mpc', 'score'),
        }
        
        for predictor, path in score_mappings.items():
            value = self._extract_nested_value(dbnsfp, path)
            if value is not None:
                # Handle list values (take mean or first)
                if isinstance(value, list):
                    value = value[0] if len(value) == 1 else sum(value) / len(value)
                
                try:
                    value = float(value)
                    # Validate before adding
                    if validate_predictor_score(predictor, value):
                        scores[predictor] = PredictorScore(
                            predictor=predictor,
                            value=value,
                            source='dbNSFP',
                            version='4.x',
                            raw=dbnsfp,
                            is_inverted=(predictor in INVERTED_PREDICTORS)
                        )
                except (ValueError, TypeError):
                    pass
        
        return scores
    
    def _extract_nested_value(self, data: dict, path: tuple) -> Optional[Any]:
        """Extract a value from nested dictionary using a path tuple."""
        current = data
        for key in path:
            if isinstance(current, dict):
                current = current.get(key)
            else:
                return None
            if current is None:
                return None
        return current
    
    def _fetch_alphamissense(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str
    ) -> Optional[PredictorScore]:
        """
        Fetch AlphaMissense score from dedicated API.
        
        AlphaMissense is a deep learning model from Google DeepMind.
        This uses a community API endpoint since Google's official API
        requires special access.
        """
        legacy_cache_key = f"alphamissense_{chrom}_{pos}_{ref}_{alt}"
        
        # Check validated cache first
        if self._use_validated_cache:
            cached_data = self._get_cached_predictor_data('AlphaMissense_API', chrom, pos, ref, alt)
            if cached_data:
                score = cached_data.get('alphamissense')
                if score is not None and validate_predictor_score('alphamissense', score):
                    return PredictorScore(
                        predictor='alphamissense',
                        value=float(score),
                        source='AlphaMissense_API',
                        version='1.0_cached',
                        is_inverted=False
                    )
        elif legacy_cache_key in self.cache:
            cached = self.cache[legacy_cache_key]
            if cached:
                return PredictorScore(
                    predictor='alphamissense',
                    value=cached.get('score'),
                    source='AlphaMissense_API',
                    version='1.0',
                    raw=cached,
                    is_inverted=False
                )
            return None
        
        try:
            # Use HegeLab's AlphaMissense API
            response = requests.get(
                f"{PREDICTOR_API_ENDPOINTS['alphamissense']}/{chrom}/{pos}/{ref}/{alt}",
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                data = response.json()
                
                score = data.get('am_pathogenicity') or data.get('score')
                if score is not None:
                    score = float(score)
                    
                    # Validate before caching
                    if validate_predictor_score('alphamissense', score):
                        if self._use_validated_cache:
                            self._set_cached_predictor_data(
                                'AlphaMissense_API', chrom, pos, ref, alt,
                                {'alphamissense': score}
                            )
                        else:
                            self.cache[legacy_cache_key] = data
                        
                        print(f"{Fore.GREEN}✅ AlphaMissense API: score={score:.4f}{Style.RESET_ALL}")
                        return PredictorScore(
                            predictor='alphamissense',
                            value=score,
                            source='AlphaMissense_API',
                            version='1.0',
                            raw=data,
                            is_inverted=False
                        )
            elif response.status_code == 404:
                self.cache[legacy_cache_key] = None
                print(f"{Fore.YELLOW}⚠️  Variant not found in AlphaMissense{Style.RESET_ALL}")
                
        except Exception as e:
            print(f"{Fore.YELLOW}⚠️  AlphaMissense API error: {str(e)}{Style.RESET_ALL}")
        
        return None
    
    def _fetch_cadd(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str
    ) -> Optional[PredictorScore]:
        """
        Fetch CADD score from official CADD API.
        
        CADD (Combined Annotation Dependent Depletion) is one of the
        most widely used variant impact predictors.
        """
        legacy_cache_key = f"cadd_{chrom}_{pos}_{ref}_{alt}"
        
        # Check validated cache first
        if self._use_validated_cache:
            cached_data = self._get_cached_predictor_data('CADD_API', chrom, pos, ref, alt)
            if cached_data:
                phred = cached_data.get('cadd_phred')
                if phred is not None and validate_predictor_score('cadd_phred', phred):
                    return PredictorScore(
                        predictor='cadd_phred',
                        value=float(phred),
                        source='CADD_API',
                        version='1.6_cached',
                        is_inverted=False
                    )
        elif legacy_cache_key in self.cache:
            cached = self.cache[legacy_cache_key]
            if cached:
                return PredictorScore(
                    predictor='cadd_phred',
                    value=cached.get('phred'),
                    source='CADD_API',
                    version='1.6',
                    raw=cached,
                    is_inverted=False
                )
            return None
        
        try:
            # CADD API uses a different variant format
            # Format: chrom-pos-ref-alt
            variant_str = f"{chrom}-{pos}-{ref}-{alt}"
            
            response = requests.get(
                f"{PREDICTOR_API_ENDPOINTS['cadd']}/{variant_str}",
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                data = response.json()
                
                # CADD returns both raw and phred-scaled scores
                phred_score = data.get('PHRED') or data.get('phred')
                if phred_score is not None:
                    phred_score = float(phred_score)
                    
                    # Validate before caching
                    if validate_predictor_score('cadd_phred', phred_score):
                        if self._use_validated_cache:
                            self._set_cached_predictor_data(
                                'CADD_API', chrom, pos, ref, alt,
                                {'cadd_phred': phred_score}
                            )
                        else:
                            self.cache[legacy_cache_key] = data
                        
                        print(f"{Fore.GREEN}✅ CADD API: PHRED={phred_score:.2f}{Style.RESET_ALL}")
                        return PredictorScore(
                            predictor='cadd_phred',
                            value=phred_score,
                            source='CADD_API',
                            version='1.6',
                            raw=data,
                            is_inverted=False
                        )
            elif response.status_code == 404:
                self.cache[legacy_cache_key] = None
                
        except Exception as e:
            print(f"{Fore.YELLOW}⚠️  CADD API error: {str(e)}{Style.RESET_ALL}")
        
        return None
    
    def _merge_scores(
        self,
        results: dict[str, PredictorScore],
        new_scores: dict[str, PredictorScore]
    ) -> None:
        """Merge new scores into results, only if result is currently None."""
        for predictor, score in new_scores.items():
            if predictor in results and results[predictor].value is None:
                if score.value is not None:
                    results[predictor] = score
    
    def _get_empty_scores(self) -> dict[str, PredictorScore]:
        """Return empty PredictorScore objects for all configured predictors."""
        return {
            predictor: PredictorScore(
                predictor=predictor,
                value=None,
                source=None,
                is_inverted=(predictor in INVERTED_PREDICTORS)
            )
            for predictor in INSILICO_WEIGHTS.keys()
        }
    
    def _get_mock_scores(self) -> dict[str, PredictorScore]:
        """Return mock predictor scores for testing."""
        mock_values = {
            'revel': 0.75,
            'cadd_phred': 25.0,
            'alphamissense': 0.8,
            'sift': 0.01,  # Inverted: low = deleterious
            'polyphen2': 0.95,
            'metasvm': 0.9,
            'vest4': 0.8,
            'fathmm': -2.5,  # Inverted: negative = deleterious
        }
        
        return {
            predictor: PredictorScore(
                predictor=predictor,
                value=mock_values.get(predictor),
                source='mock',
                version='test',
                is_inverted=(predictor in INVERTED_PREDICTORS)
            )
            for predictor in INSILICO_WEIGHTS.keys()
        }
    
    def get_available_predictor_count(
        self,
        scores: dict[str, PredictorScore]
    ) -> int:
        """Count how many predictors have valid scores."""
        return sum(1 for s in scores.values() if s.is_available())
    
    def get_weighted_predictors(
        self,
        scores: dict[str, PredictorScore]
    ) -> list[tuple[str, PredictorScore, float]]:
        """
        Get available predictors with their weights, sorted by weight descending.
        
        Returns:
            List of (predictor_name, PredictorScore, weight) tuples.
        """
        available = [
            (name, score, INSILICO_WEIGHTS.get(name, 0.0))
            for name, score in scores.items()
            if score.is_available()
        ]
        return sorted(available, key=lambda x: x[2], reverse=True)


class PopulationAPIClient:
    """
    Multi-source API client for population allele frequency data.
    
    This client fetches population frequency data from gnomAD and other
    population databases with automatic fallback.
    
    Implements validated caching:
    - All cached entries are validated before use
    - Invalid AF/AC/AN values are rejected
    - Corrupted cache entries are invalidated
    
    Usage:
        cache = ResultCache(cache_dir='./api_cache')
        client = PopulationAPIClient(api_enabled=True, result_cache=cache)
        stats = client.get_population_stats(
            chrom='17', pos=43092919, ref='G', alt='A'
        )
        # Returns: {'gnomad_v4': PopulationStats(...), 'exac': PopulationStats(...)}
    """
    
    def __init__(
        self,
        api_enabled: bool = True,
        timeout: int = 15,
        cache: Optional[dict] = None,
        test_mode: bool = False,
        result_cache: Optional['ResultCache'] = None
    ):
        """
        Initialize the population API client.
        
        Args:
            api_enabled: Whether to make actual API calls
            timeout: Request timeout in seconds
            cache: Optional cache dictionary for API responses (legacy)
            test_mode: If True, return mock data instead of API calls
            result_cache: Optional ResultCache instance for validated caching
        """
        self.api_enabled = api_enabled
        self.timeout = timeout
        self.test_mode = test_mode
        
        # Use ResultCache if provided, otherwise fall back to simple dict
        if result_cache is not None and CACHE_AVAILABLE:
            self.result_cache = result_cache
            self.cache = {}  # Keep simple cache for compatibility
            self._use_validated_cache = True
        else:
            self.result_cache = None
            self.cache = cache if cache is not None else {}
            self._use_validated_cache = False
    
    def _get_cached_population_data(
        self,
        source: str,
        chrom: str,
        pos: int,
        ref: str,
        alt: str
    ) -> Optional[dict]:
        """
        Get validated cached population data.
        
        Returns None if:
        - Cache miss
        - Cached data is invalid (rejected and invalidated)
        - Cache is disabled
        """
        if not self._use_validated_cache or self.result_cache is None:
            return None
        
        cache_key = build_population_cache_key(source, chrom, pos, ref, alt)
        cached_data = self.result_cache.get(cache_key)
        
        if cached_data is None:
            return None
        
        # Validate cached data
        is_valid, errors = validate_cached_population_data(cached_data)
        if not is_valid:
            print(f"{Fore.YELLOW}⚠️  Invalid cached {source} population data - invalidating: {errors}{Style.RESET_ALL}")
            self.result_cache.invalidate(cache_key)
            return None
        
        return cached_data
    
    def _set_cached_population_data(
        self,
        source: str,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        data: dict
    ) -> None:
        """
        Set validated population data in cache.
        
        Only caches valid data; invalid data is rejected.
        """
        if not self._use_validated_cache or self.result_cache is None:
            return
        
        # Validate before caching
        is_valid, errors = validate_cached_population_data(data)
        if not is_valid:
            print(f"{Fore.YELLOW}⚠️  Not caching invalid {source} population data: {errors}{Style.RESET_ALL}")
            return
        
        cache_key = build_population_cache_key(source, chrom, pos, ref, alt)
        self.result_cache.set(cache_key, data, validated=True)
    
    def get_population_stats(
        self,
        chrom: Optional[str] = None,
        pos: Optional[int] = None,
        ref: Optional[str] = None,
        alt: Optional[str] = None
    ) -> dict[str, PopulationStats]:
        """
        Fetch population frequency data for a variant from all sources.
        
        Args:
            chrom: Chromosome (e.g., '17', 'X')
            pos: Genomic position (GRCh38)
            ref: Reference allele
            alt: Alternate allele
            
        Returns:
            Dictionary mapping population names to PopulationStats objects.
        """
        if not self.api_enabled:
            return {}
        
        if self.test_mode:
            return self._get_mock_population_stats()
        
        results = {}
        
        # Primary source: gnomAD v4 via GraphQL
        if chrom and pos and ref and alt:
            gnomad_stats = self._fetch_gnomad(chrom, pos, ref, alt)
            if gnomad_stats:
                results['gnomad_v4'] = gnomad_stats
        
        # Additional sources could be added here (ExAC, TOPMed, 1000G)
        
        return results
    
    def _fetch_gnomad(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str
    ) -> Optional[PopulationStats]:
        """
        Fetch population frequency from gnomAD v4 GraphQL API.
        
        This is the same query used by the main APIClient but packaged
        into a PopulationStats dataclass for consistent handling.
        """
        from config.api_config import API_ENDPOINTS
        
        legacy_cache_key = f"gnomad_pop_{chrom}_{pos}_{ref}_{alt}"
        
        # Check validated cache first
        if self._use_validated_cache:
            cached_data = self._get_cached_population_data('gnomAD_GraphQL', chrom, pos, ref, alt)
            if cached_data:
                # Reconstruct PopulationStats from cached data
                if validate_population_stats(
                    af=cached_data.get('af'),
                    an=cached_data.get('an'),
                    ac=cached_data.get('ac')
                ):
                    return PopulationStats(
                        population='gnomad_v4',
                        af=cached_data.get('af'),
                        an=cached_data.get('an'),
                        ac=cached_data.get('ac'),
                        homozygote_count=cached_data.get('homozygote_count'),
                        hemizygote_count=cached_data.get('hemizygote_count'),
                        popmax_af=cached_data.get('popmax_af'),
                        popmax_population=cached_data.get('popmax_population'),
                        source='gnomAD_GraphQL',
                        version='v4.0_cached'
                    )
        elif legacy_cache_key in self.cache:
            return self.cache[legacy_cache_key]
        
        graphql_query = """
        query VariantFrequency($variantId: String!, $datasetId: DatasetId!) {
          variant(variantId: $variantId, dataset: $datasetId) {
            variant_id
            genome {
              ac
              an
              af
              ac_hom
              ac_hemi
              filters
              faf95 {
                popmax
                popmax_population
              }
              populations {
                id
                ac
                an
                af
              }
            }
          }
        }
        """
        
        gnomad_variant_id = f"{chrom}-{pos}-{ref}-{alt}"
        variables = {
            'variantId': gnomad_variant_id,
            'datasetId': 'gnomad_r4'
        }
        
        try:
            print(f"{Fore.YELLOW}🔍 Querying gnomAD for population data: {chrom}:{pos}{ref}>{alt}...{Style.RESET_ALL}")
            
            response = requests.post(
                API_ENDPOINTS.get('gnomad_graphql', 'https://gnomad.broadinstitute.org/api'),
                json={'query': graphql_query, 'variables': variables},
                timeout=self.timeout,
                headers={'Content-Type': 'application/json'}
            )
            
            if response.status_code == 200:
                data = response.json()
                
                # Check for variant not found
                variant = data.get('data', {}).get('variant')
                if not variant:
                    # Variant not found - return zero frequency
                    stats = PopulationStats(
                        population='gnomad_v4',
                        af=0.0,
                        an=0,
                        ac=0,
                        source='gnomAD_GraphQL',
                        version='v4.0'
                    )
                    
                    # Cache absent variant data
                    if self._use_validated_cache:
                        self._set_cached_population_data(
                            'gnomAD_GraphQL', chrom, pos, ref, alt,
                            {'af': 0.0, 'an': 0, 'ac': 0}
                        )
                    else:
                        self.cache[legacy_cache_key] = stats
                        
                    print(f"{Fore.YELLOW}⚠️  Variant not found in gnomAD v4 (may support PM2){Style.RESET_ALL}")
                    return stats
                
                # Parse genome-wide data
                genome = variant.get('genome', {})
                
                # Parse sub-populations
                subpop = {}
                for pop in genome.get('populations', []):
                    pop_id = pop.get('id', 'unknown')
                    subpop[pop_id] = {
                        'af': pop.get('af'),
                        'ac': pop.get('ac'),
                        'an': pop.get('an')
                    }
                
                faf95 = genome.get('faf95', {}) or {}
                
                stats = PopulationStats(
                    population='gnomad_v4',
                    af=genome.get('af'),
                    an=genome.get('an'),
                    ac=genome.get('ac'),
                    homozygote_count=genome.get('ac_hom'),
                    hemizygote_count=genome.get('ac_hemi'),
                    subpop=subpop if subpop else None,
                    popmax_af=faf95.get('popmax'),
                    popmax_population=faf95.get('popmax_population'),
                    filters=genome.get('filters'),
                    source='gnomAD_GraphQL',
                    version='v4.0',
                    raw=variant
                )
                
                # Cache validated population data
                if self._use_validated_cache:
                    cache_data = {
                        'af': stats.af,
                        'an': stats.an,
                        'ac': stats.ac,
                        'homozygote_count': stats.homozygote_count,
                        'hemizygote_count': stats.hemizygote_count,
                        'popmax_af': stats.popmax_af,
                        'popmax_population': stats.popmax_population,
                    }
                    self._set_cached_population_data(
                        'gnomAD_GraphQL', chrom, pos, ref, alt, cache_data
                    )
                else:
                    self.cache[legacy_cache_key] = stats
                
                af = stats.af
                if af:
                    print(f"{Fore.GREEN}✅ gnomAD: AF={af:.6f} (AC={stats.ac}, AN={stats.an}){Style.RESET_ALL}")
                else:
                    print(f"{Fore.GREEN}✅ gnomAD: Variant found, AF=0{Style.RESET_ALL}")
                
                return stats
                
        except Exception as e:
            print(f"{Fore.RED}❌ gnomAD API error: {str(e)}{Style.RESET_ALL}")
        
        return None
    
    def _get_mock_population_stats(self) -> dict[str, PopulationStats]:
        """Return mock population stats for testing."""
        return {
            'gnomad_v4': PopulationStats(
                population='gnomad_v4',
                af=0.00001,
                an=152000,
                ac=2,
                homozygote_count=0,
                hemizygote_count=0,
                subpop={
                    'nfe': {'af': 0.00001, 'ac': 1, 'an': 100000},
                    'afr': {'af': 0.00002, 'ac': 1, 'an': 50000}
                },
                source='mock',
                version='test'
            )
        }
    
    def get_max_frequency(
        self,
        stats: dict[str, PopulationStats]
    ) -> Optional[float]:
        """Get maximum allele frequency across all population sources."""
        afs = [s.af for s in stats.values() if s.af is not None]
        return max(afs) if afs else None
    
    def is_absent_from_all(
        self,
        stats: dict[str, PopulationStats]
    ) -> bool:
        """Check if variant is absent from all population databases."""
        if not stats:
            return True  # No data means effectively absent
        return all(s.is_absent() for s in stats.values())
