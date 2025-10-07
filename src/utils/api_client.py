import time

class APICache:
    def __init__(self):
        self.cache = {}
        self.ttl = 3600  # seconds
        self.timestamps = {}
    def get(self, key):
        if key in self.cache and (time.time() - self.timestamps[key]) < self.ttl:
            return self.cache[key]
        return None
    def set(self, key, value):
        self.cache[key] = value
        self.timestamps[key] = time.time()

class RateLimiter:
    def __init__(self, max_calls=10, period=1):
        self.max_calls = max_calls
        self.period = period
        self.calls = []
    def allow(self):
        now = time.time()
        self.calls = [t for t in self.calls if now - t < self.period]
        if len(self.calls) < self.max_calls:
            self.calls.append(now)
            return True
        return False

class CachedAPIClient:
    def __init__(self):
        self.cache = APICache()
        self.rate_limiter = RateLimiter()
    def get_variant_annotations(self, variant_key):
        cached = self.cache.get(variant_key)
        if cached:
            return cached
        if not self.rate_limiter.allow():
            raise Exception("Rate limit exceeded")
        result = self._fetch_variant_annotations(variant_key)
        self.cache.set(variant_key, result)
        return result
    def _fetch_variant_annotations(self, variant_key):
        # Placeholder: Simulate annotation
        return {"variant": variant_key, "annotation": "simulated_result"}
    def batch_annotate_variants(self, variant_list):
        cached_results = {}
        uncached_variants = []
        for variant in variant_list:
            cached_result = self.cache.get(variant)
            if cached_result:
                cached_results[variant] = cached_result
            else:
                uncached_variants.append(variant)
        if uncached_variants:
            batch_results = {v: self._fetch_variant_annotations(v) for v in uncached_variants}
            for variant, result in batch_results.items():
                self.cache.set(variant, result)
                cached_results[variant] = result
        return cached_results
"""
API Client Module
================

This module handles external API calls for enriching variant data
with information from Ensembl, ClinVar, and other databases.
"""

import requests
import json
import time
from typing import Dict, Optional, Any, Tuple
from datetime import datetime, timedelta
import os
from colorama import Fore, Style, init
from config.constants import API_ENDPOINTS, OUTPUT_SETTINGS, COLORAMA_COLORS
from utils.api_error_handler import get_error_handler

# Initialize colorama
init()


class APIClient:
    """
    Client for external API calls and data enrichment.
    
    Handles requests to Ensembl, ClinVar, and other external databases
    with caching and error handling.
    """
    
    def __init__(self, cache_enabled: bool = True):
        """
        Initialize the API client.
        
        Args:
            cache_enabled (bool): Whether to enable response caching
        """
        self.cache_enabled = cache_enabled
        self.cache = {}
        self.cache_file = OUTPUT_SETTINGS['cache_filename']
        self.max_cache_age = timedelta(hours=OUTPUT_SETTINGS['max_cache_age_hours'])
        
        # Initialize error handler
        self.error_handler = get_error_handler()
        
        # Load existing cache
        self._load_cache()
    
    def _load_cache(self):
        """Load cache from file if it exists."""
        if self.cache_enabled and os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, 'r') as f:
                    self.cache = json.load(f)
                # Clean expired entries
                self._clean_expired_cache()
            except (json.JSONDecodeError, IOError):
                self.cache = {}
    
    def _save_cache(self):
        """Save cache to file."""
        if self.cache_enabled:
            try:
                with open(self.cache_file, 'w') as f:
                    json.dump(self.cache, f, indent=2)
            except IOError:
                print("Warning: Could not save cache to file")
    
    def _clean_expired_cache(self):
        """Remove expired cache entries."""
        current_time = datetime.now()
        expired_keys = []
        
        for key, entry in self.cache.items():
            if 'timestamp' in entry:
                cached_time = datetime.fromisoformat(entry['timestamp'])
                if current_time - cached_time > self.max_cache_age:
                    expired_keys.append(key)
        
        for key in expired_keys:
            del self.cache[key]
    
    def _get_cached_response(self, cache_key: str) -> Optional[Dict[str, Any]]:
        """Get cached response if available and not expired."""
        if not self.cache_enabled:
            return None
        
        if cache_key in self.cache:
            entry = self.cache[cache_key]
            if 'timestamp' in entry:
                cached_time = datetime.fromisoformat(entry['timestamp'])
                if datetime.now() - cached_time <= self.max_cache_age:
                    return entry.get('data')
        
        return None
    
    def _cache_response(self, cache_key: str, data: Dict[str, Any]):
        """Cache API response."""
        if self.cache_enabled:
            self.cache[cache_key] = {
                'data': data,
                'timestamp': datetime.now().isoformat()
            }
            self._save_cache()
    
    def _api_call_with_retry(
        self, 
        api_name: str, 
        func, 
        fallback_value: Optional[Dict[str, Any]] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Execute API call with error handling and retry logic.
        
        Args:
            api_name (str): Name of the API
            func (callable): Function to execute
            fallback_value (dict, optional): Value to return on failure
            **kwargs: Arguments to pass to func
        
        Returns:
            Dict: API response or fallback value
        """
        result, error = self.error_handler.handle_api_call(
            api_name=api_name,
            func=func,
            max_retries=3,
            retry_delay=1.0,
            exponential_backoff=True,
            fallback_return=fallback_value,
            **kwargs
        )
        
        if error:
            # Error already logged by handler
            if result is None and fallback_value is not None:
                return fallback_value
            elif result is None:
                return {'error': error, 'source': api_name}
        
        return result
    
    def get_error_statistics(self) -> Dict[str, Any]:
        """
        Get API error statistics from error handler.
        
        Returns:
            Dict with error counts and details
        """
        return self.error_handler.get_error_statistics()
    
    def print_error_summary(self):
        """Print API error summary to console."""
        self.error_handler.print_error_summary()
    
    def _aa_code_to_full(self, aa_code: str) -> str:
        """Convert single-letter amino acid code to three-letter code."""
        aa_map = {
            'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
            'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
            'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
            'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
        }
        return aa_map.get(aa_code, aa_code)
    
    def _aa_full_to_code(self, aa_full: str) -> str:
        """Convert three-letter amino acid code to single-letter code."""
        aa_map = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
            'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
            'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
            'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
        }
        return aa_map.get(aa_full, aa_full)
    
    def get_chromosome_from_ensembl(self, gene_symbol: str) -> Optional[str]:
        """
        Get chromosome information for a gene from Ensembl with MyGene.info fallback.
        
        Args:
            gene_symbol (str): Gene symbol
            
        Returns:
            Optional[str]: Chromosome or None if not found
        """
        if not gene_symbol or gene_symbol.upper() == 'NOT SPECIFIED':
            return None
        
        cache_key = f"ensembl_chr_{gene_symbol.upper()}"
        cached_result = self._get_cached_response(cache_key)
        
        if cached_result is not None:
            print(f"{COLORAMA_COLORS['CYAN']}💾 Using cached chromosome data for {gene_symbol}{COLORAMA_COLORS['RESET']}")
            return cached_result.get('chromosome')
        
        # Try Ensembl first
        try:
            print(f"{COLORAMA_COLORS['BLUE']}🔍 Fetching chromosome info for {gene_symbol} from Ensembl...{COLORAMA_COLORS['RESET']}")
            server = API_ENDPOINTS['ensembl_rest']
            ext = f"/lookup/symbol/homo_sapiens/{gene_symbol}?expand=0"
            
            response = requests.get(
                server + ext,
                headers={"Content-Type": "application/json"},
                timeout=10
            )
            
            if response.status_code == 200:
                data = response.json()
                chromosome = data.get("seq_region_name")
                
                if chromosome:
                    print(f"{COLORAMA_COLORS['GREEN']}✅ Found chromosome {chromosome} for gene {gene_symbol}{COLORAMA_COLORS['RESET']}")
                    # Cache the result
                    self._cache_response(cache_key, {'chromosome': chromosome})
                    return chromosome
            elif response.status_code == 404:
                print(f"{COLORAMA_COLORS['YELLOW']}⚠️ Gene {gene_symbol} not found in Ensembl, trying MyGene.info...{COLORAMA_COLORS['RESET']}")
            else:
                print(f"{COLORAMA_COLORS['YELLOW']}⚠️ Ensembl API returned status {response.status_code}, trying MyGene.info...{COLORAMA_COLORS['RESET']}")
                
        except requests.Timeout:
            print(f"{COLORAMA_COLORS['YELLOW']}⏱️ Timeout from Ensembl, trying MyGene.info...{COLORAMA_COLORS['RESET']}")
        except requests.ConnectionError:
            print(f"{COLORAMA_COLORS['YELLOW']}🌐 Connection error with Ensembl, trying MyGene.info...{COLORAMA_COLORS['RESET']}")
        except requests.RequestException as e:
            print(f"{COLORAMA_COLORS['YELLOW']}⚠️ Error from Ensembl ({e}), trying MyGene.info...{COLORAMA_COLORS['RESET']}")
        except Exception as e:
            print(f"{COLORAMA_COLORS['YELLOW']}⚠️ Unexpected error with Ensembl ({e}), trying MyGene.info...{COLORAMA_COLORS['RESET']}")
        
        # Fallback to MyGene.info
        try:
            print(f"{COLORAMA_COLORS['BLUE']}🔍 Fetching chromosome info from MyGene.info...{COLORAMA_COLORS['RESET']}")
            response = requests.get(
                f"https://mygene.info/v3/query",
                params={
                    'q': f'symbol:{gene_symbol}',
                    'species': 'human',
                    'fields': 'genomic_pos'
                },
                timeout=10
            )
            
            if response.status_code == 200:
                data = response.json()
                hits = data.get('hits', [])
                
                if hits:
                    genomic_pos = hits[0].get('genomic_pos')
                    if genomic_pos:
                        # Handle both single dict and list of dicts
                        if isinstance(genomic_pos, list):
                            chromosome = genomic_pos[0].get('chr')
                        else:
                            chromosome = genomic_pos.get('chr')
                        
                        if chromosome:
                            print(f"{COLORAMA_COLORS['GREEN']}✅ Found chromosome {chromosome} for gene {gene_symbol} (MyGene.info){COLORAMA_COLORS['RESET']}")
                            # Cache the result
                            self._cache_response(cache_key, {'chromosome': chromosome})
                            return chromosome
            
            print(f"{COLORAMA_COLORS['RED']}❌ Could not find chromosome information for {gene_symbol}{COLORAMA_COLORS['RESET']}")
            return None
            
        except Exception as e:
            print(f"{COLORAMA_COLORS['RED']}❌ Error fetching from MyGene.info: {e}{COLORAMA_COLORS['RESET']}")
            return None
    
    def get_clinvar_status(self, chromosome: str, position: int, 
                          ref_allele: str, alt_allele: str) -> Dict[str, Any]:
        """
        Get ClinVar status for a variant.
        
        Args:
            chromosome (str): Chromosome
            position (int): Genomic position
            ref_allele (str): Reference allele
            alt_allele (str): Alternate allele
            
        Returns:
            Dict[str, Any]: ClinVar information
        """
        cache_key = f"clinvar_{chromosome}_{position}_{ref_allele}_{alt_allele}"
        cached_result = self._get_cached_response(cache_key)
        
        if cached_result is not None:
            return cached_result
        
        try:
            # Search ClinVar using E-utilities
            search_term = f"{chromosome}[chr] AND {position}[chrpos]"
            
            # First, search for the variant
            search_url = f"{API_ENDPOINTS['clinvar']}/esearch.fcgi"
            search_params = {
                'db': 'clinvar',
                'term': search_term,
                'retmode': 'json',
                'retmax': 10
            }
            
            search_response = requests.get(search_url, params=search_params, timeout=10)
            
            if search_response.status_code != 200:
                result = {'status': 'not_found', 'significance': None, 'review_status': None}
                self._cache_response(cache_key, result)
                return result
            
            search_data = search_response.json()
            
            if not search_data.get('esearchresult', {}).get('idlist'):
                result = {'status': 'not_found', 'significance': None, 'review_status': None}
                self._cache_response(cache_key, result)
                return result
            
            # Get detailed information for the first result
            clinvar_id = search_data['esearchresult']['idlist'][0]
            
            summary_url = f"{API_ENDPOINTS['clinvar']}/esummary.fcgi"
            summary_params = {
                'db': 'clinvar',
                'id': clinvar_id,
                'retmode': 'json'
            }
            
            summary_response = requests.get(summary_url, params=summary_params, timeout=10)
            
            if summary_response.status_code == 200:
                summary_data = summary_response.json()
                
                # Extract relevant information
                result = self._parse_clinvar_summary(summary_data, clinvar_id)
                self._cache_response(cache_key, result)
                return result
            else:
                result = {'status': 'error', 'significance': None, 'review_status': None}
                self._cache_response(cache_key, result)
                return result
                
        except requests.RequestException as e:
            print(f"Error fetching ClinVar data: {e}")
            result = {'status': 'error', 'significance': None, 'review_status': None}
            return result
    
    def _parse_clinvar_summary(self, summary_data: Dict[str, Any], clinvar_id: str) -> Dict[str, Any]:
        """
        Parse ClinVar summary data.
        
        Args:
            summary_data (Dict[str, Any]): Raw ClinVar summary data
            clinvar_id (str): ClinVar ID
            
        Returns:
            Dict[str, Any]: Parsed ClinVar information
        """
        try:
            result = summary_data.get('result', {})
            variant_info = result.get(clinvar_id, {})
            
            # Extract clinical significance from germline_classification (new API format)
            germline_class = variant_info.get('germline_classification', {})
            clinical_significance = germline_class.get('description', variant_info.get('clinical_significance', ''))
            review_status = germline_class.get('review_status', variant_info.get('review_status', ''))
            last_evaluated = germline_class.get('last_evaluated', variant_info.get('last_evaluated', ''))
            
            # Extract variant title
            title = variant_info.get('title', '')
            
            # Map review status to star rating
            star_rating_map = {
                'practice guideline': 4,
                'reviewed by expert panel': 3,
                'criteria provided, multiple submitters, no conflicts': 2,
                'criteria provided, single submitter': 1,
                'criteria provided, conflicting interpretations': 1,
                'no assertion provided': 0,
                'no assertion criteria provided': 0
            }
            star_rating = star_rating_map.get(review_status.lower(), 0)
            
            return {
                'status': 'found',
                'clinvar_id': clinvar_id,
                'significance': clinical_significance,
                'review_status': review_status,
                'star_rating': star_rating,
                'last_evaluated': last_evaluated,
                'title': title,
                'url': f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_id}/"
            }
            
        except (KeyError, TypeError):
            return {
                'status': 'found',
                'clinvar_id': clinvar_id,
                'significance': 'unknown',
                'review_status': 'unknown',
                'last_evaluated': '',
                'title': '',
                'url': f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_id}/"
            }
    
    def generate_varsome_url(self, chromosome: str, position: int, 
                           ref_allele: str, alt_allele: str) -> str:
        """
        Generate Varsome URL for the variant.
        
        Args:
            chromosome (str): Chromosome
            position (int): Genomic position
            ref_allele (str): Reference allele
            alt_allele (str): Alternate allele
            
        Returns:
            str: Varsome URL
        """
        # Format chromosome (remove 'chr' prefix if present)
        chr_formatted = chromosome.replace('chr', '').replace('Chr', '')
        
        # Generate the URL
        return f"{API_ENDPOINTS['varsome']}/{chr_formatted}-{position}-{ref_allele}-{alt_allele}"
    
    def validate_variant_coordinates(self, chromosome: str, position: int, 
                                   ref_allele: str, alt_allele: str) -> Dict[str, Any]:
        """
        Validate variant coordinates using external APIs.
        
        Args:
            chromosome (str): Chromosome
            position (int): Genomic position
            ref_allele (str): Reference allele
            alt_allele (str): Alternate allele
            
        Returns:
            Dict[str, Any]: Validation results
        """
        validation_result = {
            'valid': True,
            'warnings': [],
            'errors': []
        }
        
        # Basic format validation
        if not chromosome or not str(position).isdigit():
            validation_result['valid'] = False
            validation_result['errors'].append("Invalid chromosome or position format")
        
        if not ref_allele or not alt_allele:
            validation_result['valid'] = False
            validation_result['errors'].append("Reference and alternate alleles are required")
        
        # Check if alleles contain only valid DNA bases
        valid_bases = set('ATCG')
        if not set(ref_allele.upper()).issubset(valid_bases):
            validation_result['valid'] = False
            validation_result['errors'].append(f"Invalid reference allele: {ref_allele}")
        
        if not set(alt_allele.upper()).issubset(valid_bases):
            validation_result['valid'] = False
            validation_result['errors'].append(f"Invalid alternate allele: {alt_allele}")
        
        # Check for common issues
        if ref_allele.upper() == alt_allele.upper():
            validation_result['valid'] = False
            validation_result['errors'].append("Reference and alternate alleles are identical")
        
        return validation_result
    
    def get_gene_info(self, gene_symbol: str) -> Dict[str, Any]:
        """
        Get comprehensive gene information from Ensembl with MyGene.info fallback.
        
        Args:
            gene_symbol (str): Gene symbol
            
        Returns:
            Dict[str, Any]: Gene information
        """
        cache_key = f"ensembl_gene_{gene_symbol.upper()}"
        cached_result = self._get_cached_response(cache_key)
        
        if cached_result is not None:
            return cached_result
        
        # Try Ensembl first
        try:
            server = API_ENDPOINTS['ensembl_rest']
            ext = f"/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1"
            
            response = requests.get(
                server + ext,
                headers={"Content-Type": "application/json"},
                timeout=10
            )
            
            if response.status_code == 200:
                data = response.json()
                
                gene_info = {
                    'gene_id': data.get('id'),
                    'gene_name': data.get('display_name'),
                    'description': data.get('description'),
                    'chromosome': data.get('seq_region_name'),
                    'start': data.get('start'),
                    'end': data.get('end'),
                    'strand': data.get('strand'),
                    'biotype': data.get('biotype'),
                    'canonical_transcript': data.get('canonical_transcript')
                }
                
                self._cache_response(cache_key, gene_info)
                return gene_info
                
        except requests.RequestException as e:
            print(f"{COLORAMA_COLORS['YELLOW']}⚠️ Error from Ensembl ({e}), trying MyGene.info...{COLORAMA_COLORS['RESET']}")
        
        # Fallback to MyGene.info
        try:
            response = requests.get(
                f"https://mygene.info/v3/query",
                params={
                    'q': f'symbol:{gene_symbol}',
                    'species': 'human',
                    'fields': 'genomic_pos,name,summary,type_of_gene'
                },
                timeout=10
            )
            
            if response.status_code == 200:
                data = response.json()
                hits = data.get('hits', [])
                
                if hits:
                    hit = hits[0]
                    genomic_pos = hit.get('genomic_pos')
                    
                    gene_info = {
                        'gene_id': hit.get('_id'),
                        'gene_name': hit.get('symbol', gene_symbol),
                        'description': hit.get('name', ''),
                        'chromosome': None,
                        'start': None,
                        'end': None,
                        'strand': None,
                        'biotype': hit.get('type_of_gene'),
                        'canonical_transcript': None
                    }
                    
                    # Extract genomic position
                    if genomic_pos:
                        if isinstance(genomic_pos, list):
                            pos = genomic_pos[0]
                        else:
                            pos = genomic_pos
                        
                        gene_info['chromosome'] = pos.get('chr')
                        gene_info['start'] = pos.get('start')
                        gene_info['end'] = pos.get('end')
                        gene_info['strand'] = pos.get('strand')
                    
                    self._cache_response(cache_key, gene_info)
                    print(f"{COLORAMA_COLORS['GREEN']}✅ Found gene info for {gene_symbol} (MyGene.info){COLORAMA_COLORS['RESET']}")
                    return gene_info
            
            return {'error': f"Gene {gene_symbol} not found"}
                
        except requests.RequestException as e:
            print(f"{COLORAMA_COLORS['RED']}❌ Error fetching gene information: {e}{COLORAMA_COLORS['RESET']}")
            return {'error': str(e)}
    
    def clear_cache(self):
        """Clear all cached data."""
        self.cache = {}
        if os.path.exists(self.cache_file):
            try:
                os.remove(self.cache_file)
            except OSError:
                print("Warning: Could not remove cache file")
    
    def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache statistics."""
        total_entries = len(self.cache)
        ensembl_entries = sum(1 for key in self.cache.keys() if key.startswith('ensembl'))
        clinvar_entries = sum(1 for key in self.cache.keys() if key.startswith('clinvar'))
        gnomad_entries = sum(1 for key in self.cache.keys() if key.startswith('gnomad'))
        
        return {
            'total_entries': total_entries,
            'ensembl_entries': ensembl_entries,
            'clinvar_entries': clinvar_entries,
            'gnomad_entries': gnomad_entries,
            'cache_file': self.cache_file,
            'max_age_hours': OUTPUT_SETTINGS['max_cache_age_hours']
        }
    
    def get_gene_constraint(self, gene_symbol: str) -> Dict[str, Any]:
        """
        Query gnomAD for gene constraint metrics (pLI, LOEUF, oe_lof).
        
        This method retrieves loss-of-function (LOF) constraint metrics from gnomAD v4
        to determine if a gene is intolerant to loss-of-function variants.
        
        Args:
            gene_symbol (str): Gene symbol (e.g., 'BRCA1', 'TP53')
            
        Returns:
            Dict[str, Any]: Dictionary containing:
                - pLI (float): Probability of loss-of-function intolerance (0-1)
                - LOEUF (float): Loss-of-function observed/expected upper bound fraction
                - oe_lof (float): Observed/expected LOF ratio
                - oe_lof_upper (float): Upper bound of oe_lof confidence interval
                - exp_lof (float): Expected number of LOF variants
                - obs_lof (int): Observed number of LOF variants
                - is_lof_intolerant (bool): True if gene is LOF intolerant
                - confidence (str): Confidence level ('high', 'medium', 'low')
                - classification (str): 'LOF_intolerant', 'LOF_tolerant', or 'uncertain'
                - source (str): Data source version
                - error (str, optional): Error message if query failed
                
        Thresholds for LOF intolerance (from CONSTRAINT_THRESHOLDS):
            - pLI ≥ 0.9 → LOF intolerant
            - LOEUF ≤ 0.35 → LOF intolerant
            - oe_lof_upper > 0.6 → LOF tolerant
        """
        from config.constants import API_SETTINGS, CONSTRAINT_THRESHOLDS
        
        # Check if API is enabled
        if not API_SETTINGS.get('enabled', True):
            return {
                'error': 'API integration is disabled',
                'is_lof_intolerant': None,
                'classification': 'unknown'
            }
        
        # Normalize gene symbol
        gene_symbol = gene_symbol.strip().upper()
        
        # Check cache first
        cache_key = f"gnomad_constraint_{gene_symbol}"
        cached = self._get_cached_response(cache_key)
        if cached:
            print(f"{Fore.CYAN}ℹ️  Using cached constraint data for {gene_symbol}{Style.RESET_ALL}")
            return cached
        
        # gnomAD GraphQL query for gene constraint
        graphql_query = """
        query GeneConstraint($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!) {
          gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
            gene_id
            gene_version
            symbol
            gnomad_constraint {
              exp_lof
              obs_lof
              oe_lof
              oe_lof_lower
              oe_lof_upper
              lof_z
              pLI
              exp_mis
              obs_mis
              oe_mis
              oe_mis_lower
              oe_mis_upper
              mis_z
              flags
            }
          }
        }
        """
        
        variables = {
            'geneSymbol': gene_symbol,
            'referenceGenome': 'GRCh38'
        }
        
        try:
            print(f"{Fore.YELLOW}🔍 Querying gnomAD for {gene_symbol} constraint metrics...{Style.RESET_ALL}")
            
            response = requests.post(
                API_ENDPOINTS['gnomad_graphql'],
                json={'query': graphql_query, 'variables': variables},
                timeout=API_SETTINGS.get('timeout', 15),
                headers={'Content-Type': 'application/json'}
            )
            
            if response.status_code == 200:
                data = response.json()
                
                # Check for errors in GraphQL response
                if 'errors' in data:
                    error_msg = data['errors'][0].get('message', 'Unknown GraphQL error')
                    print(f"{Fore.RED}❌ gnomAD API error: {error_msg}{Style.RESET_ALL}")
                    return {
                        'error': error_msg,
                        'is_lof_intolerant': None,
                        'classification': 'unknown'
                    }
                
                # Extract gene constraint data
                gene_data = data.get('data', {}).get('gene')
                if not gene_data:
                    return {
                        'error': f'Gene {gene_symbol} not found in gnomAD',
                        'is_lof_intolerant': None,
                        'classification': 'unknown'
                    }
                
                constraint = gene_data.get('gnomad_constraint', {})
                if not constraint:
                    return {
                        'error': f'No constraint data available for {gene_symbol}',
                        'is_lof_intolerant': None,
                        'classification': 'unknown'
                    }
                
                # Extract metrics
                pLI = constraint.get('pLI')
                oe_lof = constraint.get('oe_lof')
                oe_lof_upper = constraint.get('oe_lof_upper')
                oe_lof_lower = constraint.get('oe_lof_lower')
                exp_lof = constraint.get('exp_lof')
                obs_lof = constraint.get('obs_lof')
                lof_z = constraint.get('lof_z')
                
                # Extract missense constraint metrics
                oe_mis = constraint.get('oe_mis')
                oe_mis_upper = constraint.get('oe_mis_upper')
                oe_mis_lower = constraint.get('oe_mis_lower')
                exp_mis = constraint.get('exp_mis')
                obs_mis = constraint.get('obs_mis')
                mis_z = constraint.get('mis_z')
                
                flags = constraint.get('flags', [])
                
                # Calculate LOEUF (synonymous with oe_lof_upper in gnomAD v4)
                LOEUF = oe_lof_upper
                
                # Determine LOF intolerance based on thresholds
                pLI_threshold = CONSTRAINT_THRESHOLDS['pLI_intolerant']
                LOEUF_threshold = CONSTRAINT_THRESHOLDS['LOEUF_intolerant']
                oe_upper_threshold = CONSTRAINT_THRESHOLDS['oe_lof_upper_tolerant']
                
                is_lof_intolerant = False
                confidence = 'low'
                classification = 'uncertain'
                
                # Classification logic for LOF
                if pLI is not None and pLI >= pLI_threshold:
                    is_lof_intolerant = True
                    confidence = 'high'
                    classification = 'LOF_intolerant'
                elif LOEUF is not None and LOEUF <= LOEUF_threshold:
                    is_lof_intolerant = True
                    confidence = 'high'
                    classification = 'LOF_intolerant'
                elif oe_lof_upper is not None and oe_lof_upper > oe_upper_threshold:
                    is_lof_intolerant = False
                    confidence = 'high'
                    classification = 'LOF_tolerant'
                elif pLI is not None and LOEUF is not None:
                    # Both available but don't meet thresholds
                    if pLI >= 0.5 or LOEUF <= 0.5:
                        confidence = 'medium'
                        classification = 'uncertain_intolerant'
                    else:
                        confidence = 'medium'
                        classification = 'uncertain_tolerant'
                
                # Classification logic for missense constraint
                # Based on gnomAD recommendations:
                # - mis_z > 3.09: Significantly constrained against missense variation
                # - oe_mis_upper < 0.6: Missense constrained
                # - oe_mis > 1.0: Missense tolerant
                is_mis_constrained = False
                mis_classification = 'uncertain'
                mis_confidence = 'low'
                
                if mis_z is not None and mis_z > 3.09:
                    is_mis_constrained = True
                    mis_classification = 'missense_constrained'
                    mis_confidence = 'high'
                elif oe_mis_upper is not None and oe_mis_upper < 0.6:
                    is_mis_constrained = True
                    mis_classification = 'missense_constrained'
                    mis_confidence = 'high'
                elif oe_mis is not None and oe_mis > 1.0:
                    is_mis_constrained = False
                    mis_classification = 'missense_tolerant'
                    mis_confidence = 'medium'
                elif oe_mis is not None:
                    # oe_mis between 0.6 and 1.0
                    if oe_mis < 0.8:
                        mis_classification = 'uncertain_constrained'
                        mis_confidence = 'low'
                    else:
                        mis_classification = 'uncertain_tolerant'
                        mis_confidence = 'low'
                
                result = {
                    'gene_symbol': gene_symbol,
                    # LOF metrics
                    'pLI': pLI,
                    'LOEUF': LOEUF,
                    'oe_lof': oe_lof,
                    'oe_lof_upper': oe_lof_upper,
                    'oe_lof_lower': oe_lof_lower,
                    'exp_lof': exp_lof,
                    'obs_lof': obs_lof,
                    'lof_z': lof_z,
                    # Missense metrics
                    'oe_mis': oe_mis,
                    'oe_mis_upper': oe_mis_upper,
                    'oe_mis_lower': oe_mis_lower,
                    'exp_mis': exp_mis,
                    'obs_mis': obs_mis,
                    'mis_z': mis_z,
                    # Flags and classifications
                    'flags': flags,
                    'is_lof_intolerant': is_lof_intolerant,
                    'confidence': confidence,
                    'classification': classification,
                    'is_mis_constrained': is_mis_constrained,
                    'mis_classification': mis_classification,
                    'mis_confidence': mis_confidence,
                    'source': 'gnomAD v4',
                    'thresholds_used': {
                        'pLI': pLI_threshold,
                        'LOEUF': LOEUF_threshold,
                        'oe_lof_upper_tolerant': oe_upper_threshold,
                        'mis_z_constrained': 3.09,
                        'oe_mis_upper_constrained': 0.6
                    }
                }
                
                # Cache the result
                self._cache_response(cache_key, result)
                
                # Print summary
                print(f"{Fore.GREEN}✅ gnomAD constraint data retrieved for {gene_symbol}{Style.RESET_ALL}")
                print(f"   LOF: pLI={pLI:.3f} | LOEUF={LOEUF:.3f} | {classification}")
                if mis_z is not None or oe_mis is not None:
                    mis_info = f"mis_z={mis_z:.3f}" if mis_z else f"oe_mis={oe_mis:.3f}"
                    print(f"   Missense: {mis_info} | {mis_classification}")
                
                return result
                
            else:
                error_msg = f"HTTP {response.status_code}: {response.text[:100]}"
                print(f"{Fore.RED}❌ gnomAD API request failed: {error_msg}{Style.RESET_ALL}")
                return {
                    'error': error_msg,
                    'is_lof_intolerant': None,
                    'classification': 'unknown'
                }
                
        except requests.Timeout:
            print(f"{Fore.RED}❌ gnomAD API request timed out{Style.RESET_ALL}")
            return {
                'error': 'Request timeout',
                'is_lof_intolerant': None,
                'classification': 'unknown'
            }
        except requests.RequestException as e:
            print(f"{Fore.RED}❌ gnomAD API request failed: {str(e)}{Style.RESET_ALL}")
            return {
                'error': str(e),
                'is_lof_intolerant': None,
                'classification': 'unknown'
            }
        except Exception as e:
            print(f"{Fore.RED}❌ Unexpected error querying gnomAD: {str(e)}{Style.RESET_ALL}")
            return {
                'error': f'Unexpected error: {str(e)}',
                'is_lof_intolerant': None,
                'classification': 'unknown'
            }
    
    def get_clingen_gene_validity(self, gene_symbol: str, disease: Optional[str] = None) -> Dict[str, Any]:
        """
        Query ClinGen for gene-disease validity and disease mechanism.
        
        ClinGen curates gene-disease relationships with evidence-based validity classifications.
        This is crucial for understanding if LOF variants are pathogenic in specific disease contexts.
        
        Args:
            gene_symbol (str): Gene symbol (e.g., 'BRCA1')
            disease (str, optional): Disease name or MONDO ID (e.g., 'Breast-ovarian cancer')
            
        Returns:
            Dict[str, Any]: Dictionary containing:
                - gene_symbol (str): Gene symbol
                - curations (list): List of gene-disease curations with:
                    * disease_label (str): Disease name
                    * disease_id (str): MONDO/OMIM ID
                    * classification (str): 'Definitive', 'Strong', 'Moderate', 'Limited', 'Disputed', 'Refuted'
                    * mechanism (str): 'loss of function', 'gain of function', 'dominant negative', etc.
                    * moi (str): Mode of inheritance (AD, AR, XL, etc.)
                    * pmids (list): Supporting publications
                - primary_mechanism (str): Most common disease mechanism
                - supports_lof_pathogenicity (bool): True if LOF mechanism established
                - confidence (str): Overall confidence level
                - source (str): 'ClinGen'
                - error (str, optional): Error message if query failed
        
        Note: This helps resolve cases like BRCA1 where population constraint (LOF tolerant)
              conflicts with disease-specific pathogenicity (LOF causes cancer).
        """
        from config.constants import API_SETTINGS
        
        # Check if API is enabled
        if not API_SETTINGS.get('enabled', True):
            return {
                'error': 'API integration is disabled',
                'supports_lof_pathogenicity': None,
                'source': 'disabled'
            }
        
        # Normalize gene symbol
        gene_symbol = gene_symbol.strip().upper()
        
        # Check cache first
        cache_key = f"clingen_erepo_{gene_symbol}"
        if disease:
            cache_key += f"_{disease.lower().replace(' ', '_')}"
        
        cached = self._get_cached_response(cache_key)
        if cached:
            print(f"{Fore.CYAN}ℹ️  Using cached ClinGen data for {gene_symbol}{Style.RESET_ALL}")
            return cached
        
        try:
            print(f"{Fore.YELLOW}🔍 Querying ClinGen eRepo for {gene_symbol} variant interpretations...{Style.RESET_ALL}")
            
            # ClinGen Evidence Repository (eRepo) API endpoint
            base_url = "https://erepo.genome.network/evrepo/api"
            
            # Get variant interpretations for gene
            response = requests.get(
                f"{base_url}/interpretations",
                params={
                    'gene': gene_symbol,
                    'matchMode': 'keyword',
                    'matchLimit': '100'  # Get up to 100 interpretations
                },
                timeout=API_SETTINGS.get('timeout', 15),
                headers={'Accept': 'application/json'}
            )
            
            if response.status_code == 200:
                data = response.json()
                interpretations = data.get('variantInterpretations', [])
                
                if not interpretations:
                    result = {
                        'gene_symbol': gene_symbol,
                        'total_interpretations': 0,
                        'supports_lof_pathogenicity': None,
                        'confidence': 'no_data',
                        'source': 'ClinGen eRepo',
                        'message': f'No interpretations found for {gene_symbol}'
                    }
                    self._cache_response(cache_key, result)
                    return result
                
                # OPTIMIZATION: Extract evidence directly from summary instead of detail requests!
                # This eliminates 30+ slow API calls (1.5s each = 45s total)
                pathogenic_count = 0
                benign_count = 0
                pvs1_count = 0
                lof_diseases = set()
                evidence_code_counts = {}
                
                # Analyze interpretations from SUMMARY data (no detail requests needed!)
                for interp_summary in interpretations[:50]:  # Analyze more since we're faster now
                    try:
                        # Get condition/disease from summary
                        condition = interp_summary.get('condition', {})
                        disease_label = condition.get('label', '')
                        if disease_label:
                            lof_diseases.add(disease_label)
                        
                        # Extract evidence codes from guidelines in SUMMARY
                        guidelines = interp_summary.get('guidelines', [])
                        has_pvs1_met = False
                        
                        for guideline in guidelines:
                            if not isinstance(guideline, dict):
                                continue
                            
                            agents = guideline.get('agents', [])
                            for agent in agents:
                                if not isinstance(agent, dict):
                                    continue
                                
                                evidence_codes = agent.get('evidenceCodes', [])
                                for ev_code in evidence_codes:
                                    if not isinstance(ev_code, dict):
                                        continue
                                    
                                    label = ev_code.get('label', '')
                                    status = ev_code.get('status', '')
                                    
                                    # Track evidence codes
                                    if label:
                                        if label not in evidence_code_counts:
                                            evidence_code_counts[label] = {'met': 0, 'not_met': 0}
                                        
                                        if status == 'Met':
                                            evidence_code_counts[label]['met'] += 1
                                            
                                            # PVS1 met = LOF evidence
                                            if label == 'PVS1':
                                                has_pvs1_met = True
                                                pathogenic_count += 1
                                        else:
                                            evidence_code_counts[label]['not_met'] += 1
                        
                        if has_pvs1_met:
                            pvs1_count += 1
                    
                    except Exception as e:
                        # Skip failed interpretations
                        continue
                
                # Determine if LOF mechanism is supported
                # Strong evidence: >5 pathogenic interpretations with PVS1
                # Moderate evidence: >3 pathogenic interpretations with PVS1
                # Weak evidence: >0 pathogenic interpretations with PVS1
                supports_lof = pvs1_count > 0
                
                if pvs1_count >= 5:
                    confidence = 'high'
                elif pvs1_count >= 3:
                    confidence = 'medium'
                elif pvs1_count > 0:
                    confidence = 'low'
                else:
                    confidence = 'no_data'
                
                result = {
                    'gene_symbol': gene_symbol,
                    'total_interpretations': len(interpretations),
                    'analyzed_interpretations': min(50, len(interpretations)),
                    'pathogenic_count': pathogenic_count,
                    'benign_count': benign_count,
                    'pvs1_lof_count': pvs1_count,
                    'lof_diseases': list(lof_diseases),
                    'evidence_codes': evidence_code_counts,
                    'supports_lof_pathogenicity': supports_lof,
                    'confidence': confidence,
                    'source': 'ClinGen eRepo',
                    'primary_mechanism': 'loss_of_function' if supports_lof else 'unknown',
                    'optimization_note': 'Fast mode: extracted evidence from summary data'
                }
                
                # Cache the result
                self._cache_response(cache_key, result)
                
                # Print summary
                print(f"{Fore.GREEN}✅ ClinGen eRepo data retrieved for {gene_symbol}{Style.RESET_ALL}")
                print(f"   Interpretations: {len(interpretations)} (analyzed: {min(50, len(interpretations))})")
                if supports_lof:
                    print(f"   ✅ LOF mechanism: {pvs1_count} variants with PVS1 evidence")
                    if lof_diseases:
                        print(f"   Diseases: {', '.join(list(lof_diseases)[:2])}")
                else:
                    print(f"   ⚠️  No PVS1 LOF evidence found")
                print(f"   ⚡ Optimization: Summary extraction (no detail requests)")
                
                return result
            
            elif response.status_code == 404:
                # Gene not found in ClinGen eRepo
                result = {
                    'gene_symbol': gene_symbol,
                    'total_interpretations': 0,
                    'supports_lof_pathogenicity': None,
                    'confidence': 'no_data',
                    'source': 'ClinGen',
                    'confidence': 'no_data',
                    'source': 'ClinGen eRepo',
                    'note': f'Gene {gene_symbol} not found in ClinGen eRepo'
                }
                
                self._cache_response(cache_key, result)
                print(f"{Fore.YELLOW}⚠️  {gene_symbol} not found in ClinGen eRepo{Style.RESET_ALL}")
                return result
            
            else:
                error_msg = f"HTTP {response.status_code}: {response.text[:100]}"
                print(f"{Fore.RED}❌ ClinGen eRepo API request failed: {error_msg}{Style.RESET_ALL}")
                return {
                    'error': error_msg,
                    'supports_lof_pathogenicity': None,
                    'source': 'ClinGen eRepo'
                }
        
        except requests.Timeout:
            print(f"{Fore.RED}❌ ClinGen eRepo API request timed out{Style.RESET_ALL}")
            return {
                'error': 'Request timeout',
                'supports_lof_pathogenicity': None,
                'source': 'ClinGen eRepo'
            }
        except requests.RequestException as e:
            print(f"{Fore.RED}❌ ClinGen eRepo API request failed: {str(e)}{Style.RESET_ALL}")
            return {
                'error': str(e),
                'supports_lof_pathogenicity': None,
                'source': 'ClinGen eRepo'
            }
        except Exception as e:
            print(f"{Fore.RED}❌ Unexpected error querying ClinGen: {str(e)}{Style.RESET_ALL}")
            return {
                'error': f'Unexpected error: {str(e)}',
                'supports_lof_pathogenicity': None,
                'source': 'ClinGen'
            }
    
    def get_clinvar_classification(self, variant_id: Optional[str] = None, gene: Optional[str] = None, 
                                   hgvs: Optional[str] = None) -> Dict[str, Any]:
        """
        Query ClinVar for variant classification and review status.
        
        Critical for PP5 (ClinVar Pathogenic) and BP6 (ClinVar Benign) criteria.
        Validates review status (star rating) for evidence strength.
        
        Args:
            variant_id (str, optional): ClinVar Variation ID (e.g., 'VCV000012345')
            gene (str, optional): Gene symbol for search
            hgvs (str, optional): HGVS notation for search
        
        Returns:
            Dict with:
                - variation_id (str): ClinVar Variation ID
                - classification (str): Clinical significance
                - review_status (str): Review status
                - star_rating (int): 0-4 stars
                - assertion_criteria (str): Guidelines used
                - date_last_evaluated (str): Last review date
                - submitters (list): List of submitting organizations
                - evidence_count (int): Number of submissions
                - conflicts (bool): True if conflicting interpretations
                - source (str): 'ClinVar'
        """
        from config.constants import API_SETTINGS
        
        if not API_SETTINGS.get('enabled', True):
            return {'error': 'API disabled', 'source': 'ClinVar'}
        
        # Check cache - include all parameters to prevent collisions
        cache_key = f"clinvar_classification_{variant_id}_{gene}_{hgvs}"
        cached = self._get_cached_response(cache_key)
        if cached:
            print(f"{Fore.CYAN}ℹ️  Using cached ClinVar data{Style.RESET_ALL}")
            return cached
        
        try:
            print(f"{Fore.YELLOW}🔍 Querying ClinVar...{Style.RESET_ALL}")
            
            # NCBI E-utilities API for ClinVar
            base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
            
            # Build search term based on input
            if variant_id:
                # If VCV ID provided, convert to UID (remove VCV prefix)
                if variant_id.startswith('VCV'):
                    search_term = variant_id.replace('VCV', '').lstrip('0')  # VCV000055406 -> 55406
                else:
                    search_term = variant_id
            elif gene and hgvs:
                search_term = f"{gene}[gene] AND {hgvs}[variant name]"
            else:
                return {'error': 'Need variant_id or (gene + hgvs)', 'source': 'ClinVar'}
            
            # ESearch to find ClinVar ID
            search_url = f"{base_url}/esearch.fcgi"
            search_params = {
                'db': 'clinvar',
                'term': search_term,
                'retmode': 'json',
                'retmax': 1
            }
            
            search_response = requests.get(search_url, params=search_params, timeout=15)
            
            if search_response.status_code == 200:
                search_data = search_response.json()
                id_list = search_data.get('esearchresult', {}).get('idlist', [])
                
                if not id_list:
                    result = {
                        'classification': 'not_found',
                        'review_status': 'not_found',
                        'star_rating': 0,
                        'source': 'ClinVar',
                        'message': f'No ClinVar entry found for {search_term}'
                    }
                    self._cache_response(cache_key, result)
                    return result
                
                clinvar_id = id_list[0]
                
                # ESummary to get details
                summary_url = f"{base_url}/esummary.fcgi"
                summary_params = {
                    'db': 'clinvar',
                    'id': clinvar_id,
                    'retmode': 'json'
                }
                
                summary_response = requests.get(summary_url, params=summary_params, timeout=15)
                
                if summary_response.status_code == 200:
                    summary_data = summary_response.json()
                    record = summary_data.get('result', {}).get(clinvar_id, {})
                    
                    # Extract classification from germline_classification field
                    germline_class = record.get('germline_classification', {})
                    classification = germline_class.get('description', 'Unknown')
                    review_status = germline_class.get('review_status', 'no assertion provided')
                    
                    # Map review status to star rating
                    star_rating_map = {
                        'practice guideline': 4,
                        'reviewed by expert panel': 3,
                        'criteria provided, multiple submitters, no conflicts': 2,
                        'criteria provided, single submitter': 1,
                        'criteria provided, conflicting interpretations': 1,
                        'no assertion provided': 0,
                        'no assertion criteria provided': 0
                    }
                    
                    star_rating = star_rating_map.get(review_status.lower(), 0)
                    
                    # Check for conflicts
                    conflicts = 'conflict' in review_status.lower()
                    
                    result = {
                        'variation_id': clinvar_id,
                        'classification': classification,
                        'review_status': review_status,
                        'star_rating': star_rating,
                        'conflicts': conflicts,
                        'date_last_evaluated': germline_class.get('last_evaluated', 'Unknown'),
                        'source': 'ClinVar',
                        'gene': record.get('genes', [{}])[0].get('symbol', 'Unknown') if record.get('genes') else 'Unknown'
                    }
                    
                    self._cache_response(cache_key, result)
                    
                    
                    print(f"{Fore.GREEN}✅ ClinVar: {classification} ({star_rating}★){Style.RESET_ALL}")
                    return result
            
            return {'error': f'HTTP {search_response.status_code}', 'source': 'ClinVar'}
        
        except Exception as e:
            print(f"{Fore.RED}❌ ClinVar API error: {str(e)}{Style.RESET_ALL}")
            return {'error': str(e), 'source': 'ClinVar'}
    
    def search_clinvar_variants_at_position(self, gene: str, hgvs_p: str) -> Dict[str, Any]:
        """
        Search ClinVar for variants at the same amino acid position (for PS1/PM5).
        
        **QUICK WIN #3**: Enables PS1 (same AA change) and PM5 (different AA at same position).
        
        Args:
            gene (str): Gene symbol (e.g., 'TP53')
            hgvs_p (str): Protein HGVS (e.g., 'p.R273H')
        
        Returns:
            Dict with:
                - same_aa_pathogenic (bool): Same amino acid change is pathogenic
                - different_aa_pathogenic (list): Different AA changes at same position
                - variants (list): List of variants found
                - source (str): 'ClinVar'
        """
        from config.constants import API_SETTINGS
        
        if not API_SETTINGS.get('enabled', True):
            return {'error': 'API disabled', 'source': 'ClinVar'}
        
        # Extract position from HGVS
        # Supports: p.R273H (1-letter), p.Arg273His (3-letter)
        import re
        
        # Try 3-letter code first (e.g., p.Arg273His)
        position_match = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|=)', hgvs_p)
        if not position_match:
            # Try 1-letter code (e.g., p.R273H)
            position_match = re.search(r'p\.([A-Z])(\d+)([A-Z]|\*|=)', hgvs_p)
        
        if not position_match:
            return {'error': f'Invalid HGVS format: {hgvs_p}', 'source': 'ClinVar'}
        
        ref_aa = position_match.group(1)
        position = position_match.group(2)
        alt_aa = position_match.group(3)
        
        # Build cache key
        cache_key = f"clinvar_position_{gene}_{position}"
        cached = self._get_cached_response(cache_key)
        if cached:
            print(f"{Fore.CYAN}ℹ️  Using cached ClinVar position data for {gene} position {position}{Style.RESET_ALL}")
            return cached
        
        try:
            print(f"{Fore.YELLOW}🔍 Searching ClinVar for {gene} variants at position {position}...{Style.RESET_ALL}")
            
            # NCBI E-utilities API
            base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
            
            # Search for all variants at this gene position
            # Query: gene[gene] AND position (e.g., TP53[gene] AND p.273)
            search_term = f"{gene}[gene] AND p.{position}"
            
            search_url = f"{base_url}/esearch.fcgi"
            search_params = {
                'db': 'clinvar',
                'term': search_term,
                'retmode': 'json',
                'retmax': 50  # Get up to 50 variants at this position
            }
            
            search_response = requests.get(search_url, params=search_params, timeout=15)
            
            if search_response.status_code != 200:
                return {'error': f'HTTP {search_response.status_code}', 'source': 'ClinVar'}
            
            search_data = search_response.json()
            id_list = search_data.get('esearchresult', {}).get('idlist', [])
            
            if not id_list:
                result = {
                    'same_aa_pathogenic': False,
                    'different_aa_pathogenic': [],
                    'variants': [],
                    'source': 'ClinVar',
                    'message': f'No variants found at {gene} p.{position}'
                }
                self._cache_response(cache_key, result)
                return result
            
            # Get summaries for all variants
            summary_url = f"{base_url}/esummary.fcgi"
            summary_params = {
                'db': 'clinvar',
                'id': ','.join(id_list),
                'retmode': 'json'
            }
            
            summary_response = requests.get(summary_url, params=summary_params, timeout=15)
            
            if summary_response.status_code != 200:
                return {'error': f'HTTP {summary_response.status_code}', 'source': 'ClinVar'}
            
            summary_data = summary_response.json()
            
            # Parse variants
            same_aa_pathogenic = False
            different_aa_pathogenic = []
            variants = []
            
            for clinvar_id in id_list:
                record = summary_data.get('result', {}).get(clinvar_id, {})
                
                # Get variant name (protein change)
                title = record.get('title', '')
                protein_change = record.get('protein_change', '')
                
                # Get classification
                germline_class = record.get('germline_classification', {})
                classification = germline_class.get('description', 'Unknown')
                
                # Check if pathogenic
                is_pathogenic = 'pathogenic' in classification.lower() and 'benign' not in classification.lower()
                
                variant_info = {
                    'clinvar_id': clinvar_id,
                    'protein_change': protein_change,
                    'classification': classification,
                    'is_pathogenic': is_pathogenic,
                    'title': title
                }
                variants.append(variant_info)
                
                if is_pathogenic and protein_change:
                    # Handle multi-variant format (e.g., "R141L, R273L, R114L")
                    protein_changes = [pc.strip() for pc in protein_change.split(',')]
                    
                    for pc in protein_changes:
                        # Extract position from protein change (e.g., R273L -> 273)
                        pc_match = re.search(r'([A-Z])(\d+)([A-Z])', pc)
                        if not pc_match:
                            continue
                        
                        pc_ref_aa = pc_match.group(1)
                        pc_position = pc_match.group(2)
                        pc_alt_aa = pc_match.group(3)
                        
                        # Check if this change is at the query position
                        if pc_position == position:
                            # Convert query ref_aa to single letter for comparison
                            query_ref_code = self._aa_full_to_code(ref_aa)
                            query_alt_code = self._aa_full_to_code(alt_aa)
                            
                            # Check ref amino acid matches (same position verification)
                            if pc_ref_aa != query_ref_code:
                                continue
                            
                            # Convert to full format for comparison (e.g., R273H -> Arg273His)
                            full_pc = f"p.{self._aa_code_to_full(pc_ref_aa)}{pc_position}{self._aa_code_to_full(pc_alt_aa)}"
                            
                            # Check if same amino acid change (PS1)
                            if full_pc == hgvs_p or pc_alt_aa == query_alt_code:
                                same_aa_pathogenic = True
                            # Different amino acid at same position (PM5)
                            else:
                                different_aa_pathogenic.append({
                                    **variant_info,
                                    'specific_change': pc,
                                    'full_hgvs': full_pc
                                })
            
            result = {
                'same_aa_pathogenic': same_aa_pathogenic,
                'different_aa_pathogenic': different_aa_pathogenic,
                'variants': variants,
                'position': position,
                'gene': gene,
                'query_hgvs': hgvs_p,
                'source': 'ClinVar'
            }
            
            self._cache_response(cache_key, result)
            
            print(f"{Fore.GREEN}✅ Found {len(variants)} variants at {gene} p.{position}{Style.RESET_ALL}")
            if same_aa_pathogenic:
                print(f"   ✓ Same AA change is pathogenic (PS1 applies)")
            if different_aa_pathogenic:
                print(f"   ✓ {len(different_aa_pathogenic)} different pathogenic AA changes (PM5 may apply)")
            
            return result
        
        except Exception as e:
            print(f"{Fore.RED}❌ ClinVar search error: {str(e)}{Style.RESET_ALL}")
            return {'error': str(e), 'source': 'ClinVar'}

    def get_variant_frequency(self, variant_id: Optional[str] = None, 
                               chrom: Optional[str] = None, 
                               pos: Optional[int] = None, 
                               ref: Optional[str] = None, 
                               alt: Optional[str] = None) -> Dict[str, Any]:
        """
        Get population allele frequency data for a variant from gnomAD v4.
        
        Args:
            variant_id (str, optional): Variant ID (rsID or VCV)
            chrom (str, optional): Chromosome (e.g., '17', 'X')
            pos (int, optional): Position (GRCh38)
            ref (str, optional): Reference allele
            alt (str, optional): Alternate allele
            
        Returns:
            Dict[str, Any]: Dictionary containing:
                - allele_count (int): Number of alternate alleles observed
                - allele_number (int): Total number of alleles genotyped
                - allele_frequency (float): Allele frequency (AC/AN)
                - homozygote_count (int): Number of homozygous individuals
                - hemizygote_count (int): Number of hemizygous individuals (X/Y chr)
                - popmax_af (float): Maximum AF across populations
                - popmax_population (str): Population with maximum AF
                - filters (list): Quality filters
                - source (str): 'gnomAD v4'
                - error (str, optional): Error message if query failed
                
        Usage for BA1/BS1/PM2:
            - BA1: allele_frequency > 0.05 (5%)
            - BS1: allele_frequency > expected_max_af for disease
            - PM2: allele_frequency is None or very rare (<0.0001)
        """
        from config.constants import API_SETTINGS
        
        # Check if API is enabled
        if not API_SETTINGS.get('enabled', True):
            return {'error': 'API integration is disabled', 'source': 'gnomAD v4'}
        
        # Validate input
        if not ((chrom and pos and ref and alt) or variant_id):
            return {'error': 'Must provide either (chrom, pos, ref, alt) or variant_id', 'source': 'gnomAD v4'}
        
        # Build cache key
        if chrom and pos and ref and alt:
            cache_key = f"gnomad_freq_{chrom}_{pos}_{ref}_{alt}"
            query_descriptor = f"{chrom}:{pos}{ref}>{alt}"
        else:
            cache_key = f"gnomad_freq_{variant_id}"
            query_descriptor = variant_id
        
        # Check cache
        cached = self._get_cached_response(cache_key)
        if cached:
            print(f"{Fore.CYAN}ℹ️  Using cached frequency data{Style.RESET_ALL}")
            return cached
        
        # gnomAD GraphQL query for variant frequency
        graphql_query = """
        query VariantFrequency($variantId: String!, $datasetId: DatasetId!) {
          variant(variantId: $variantId, dataset: $datasetId) {
            variant_id
            reference_genome
            chrom
            pos
            ref
            alt
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
            }
          }
        }
        """
        
        # Build variant ID for gnomAD (format: chrom-pos-ref-alt)
        if chrom and pos and ref and alt:
            gnomad_variant_id = f"{chrom}-{pos}-{ref}-{alt}"
        else:
            # If only variant_id provided, we need to resolve it first
            # This is complex - would need ClinVar → genomic coordinates → gnomAD
            return {
                'error': 'Variant ID lookup not yet implemented. Please provide genomic coordinates.',
                'source': 'gnomAD v4'
            }
        
        variables = {
            'variantId': gnomad_variant_id,
            'datasetId': 'gnomad_r4'
        }
        
        try:
            print(f"{Fore.YELLOW}🔍 Querying gnomAD for variant frequency: {query_descriptor}...{Style.RESET_ALL}")
            
            response = requests.post(
                API_ENDPOINTS['gnomad_graphql'],
                json={'query': graphql_query, 'variables': variables},
                timeout=API_SETTINGS.get('timeout', 15),
                headers={'Content-Type': 'application/json'}
            )
            
            if response.status_code == 200:
                data = response.json()
                
                # Check for errors
                if 'errors' in data:
                    error_msg = data['errors'][0].get('message', 'Unknown GraphQL error')
                    
                    # "Variant not found" is not an error - it's useful info for PM2
                    if 'not found' in error_msg.lower():
                        print(f"{Fore.YELLOW}⚠️  Variant not found in gnomAD v4{Style.RESET_ALL}")
                        result = {
                            'allele_count': 0,
                            'allele_number': 0,
                            'allele_frequency': None,
                            'homozygote_count': 0,
                            'hemizygote_count': 0,
                            'popmax_af': None,
                            'popmax_population': None,
                            'filters': [],
                            'source': 'gnomAD v4',
                            'note': 'Variant not found (may support PM2)'
                        }
                        self._cache_response(cache_key, result)
                        return result
                    else:
                        print(f"{Fore.RED}❌ gnomAD API error: {error_msg}{Style.RESET_ALL}")
                        return {'error': error_msg, 'source': 'gnomAD v4'}
                
                # Extract variant data
                variant = data.get('data', {}).get('variant')
                if not variant:
                    # Variant not found - this is useful for PM2 (absent in databases)
                    print(f"{Fore.YELLOW}⚠️  Variant not found in gnomAD v4{Style.RESET_ALL}")
                    result = {
                        'allele_count': 0,
                        'allele_number': 0,
                        'allele_frequency': None,
                        'homozygote_count': 0,
                        'hemizygote_count': 0,
                        'popmax_af': None,
                        'popmax_population': None,
                        'filters': [],
                        'source': 'gnomAD v4',
                        'note': 'Variant not found (may support PM2)'
                    }
                    self._cache_response(cache_key, result)
                    return result
                
                # Extract genome-wide frequency data
                genome = variant.get('genome', {})
                ac = genome.get('ac', 0)
                an = genome.get('an', 0)
                af = genome.get('af')
                homozygote_count = genome.get('ac_hom', 0)
                hemizygote_count = genome.get('ac_hemi', 0)
                filters = genome.get('filters', [])
                
                # Extract filtering allele frequency (faf95) for popmax
                faf95 = genome.get('faf95', {})
                popmax_af = faf95.get('popmax')
                popmax_population = faf95.get('popmax_population')
                
                result = {
                    'allele_count': ac,
                    'allele_number': an,
                    'allele_frequency': af,
                    'homozygote_count': homozygote_count,
                    'hemizygote_count': hemizygote_count,
                    'popmax_af': popmax_af,
                    'popmax_population': popmax_population,
                    'filters': filters,
                    'source': 'gnomAD v4'
                }
                
                self._cache_response(cache_key, result)
                
                if af:
                    print(f"{Fore.GREEN}✅ gnomAD: AF={af:.6f} (AC={ac}, AN={an}){Style.RESET_ALL}")
                else:
                    print(f"{Fore.GREEN}✅ gnomAD: Variant found but AF=0{Style.RESET_ALL}")
                
                return result
            
            return {'error': f'HTTP {response.status_code}', 'source': 'gnomAD v4'}
        
        except Exception as e:
            print(f"{Fore.RED}❌ gnomAD API error: {str(e)}{Style.RESET_ALL}")
            return {'error': str(e), 'source': 'gnomAD v4'}
    
    def get_clingen_dosage_sensitivity(self, gene_symbol: str) -> Dict[str, Any]:
        """
        Query ClinGen Dosage Sensitivity Map for haploinsufficiency/triplosensitivity scores.
        
        HYBRID APPROACH:
        1. First tries local TSV file (fast, no network)
        2. If not found, downloads from ClinGen FTP (reliable, current)
        3. Caches result for future queries
        
        Data Source: https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
        
        Args:
            gene_symbol (str): Gene symbol (e.g., 'BRCA1', 'TP53', 'DMD')
            
        Returns:
            Dict[str, Any]: Dictionary containing:
                - haploinsufficiency_score (int): 0-3, 30, or 40
                    * 3 = Sufficient evidence for dosage pathogenicity
                    * 2 = Some evidence for dosage pathogenicity
                    * 1 = Little evidence for dosage pathogenicity
                    * 0 = No evidence available
                    * 30 = Gene associated with autosomal recessive phenotype
                    * 40 = Dosage sensitivity unlikely
                - triplosensitivity_score (int): 0-3, 30, or 40
                - haploinsufficiency_description (str): Text description
                - triplosensitivity_description (str): Text description
                - pvs1_recommendation (str): PVS1 strength recommendation
                - confidence (str): Confidence level
                - pmids (list): Supporting PMIDs
                - source (str): Data source
                - last_evaluated (str): Last curation date
                
        PVS1 Strength Recommendations:
            - Score 3: Keep PVS1 (Very Strong)
            - Score 2: Downgrade to PS1 (Strong)
            - Score 1: Downgrade to PM2 (Moderate)
            - Score 0/30/40: Consider not applying PVS1
        """
        from config.constants import API_SETTINGS, API_ENDPOINTS
        import os
        
        # Check if API is enabled
        if not API_SETTINGS.get('enabled', True):
            return {
                'error': 'API integration is disabled',
                'haploinsufficiency_score': None,
                'triplosensitivity_score': None,
                'pvs1_recommendation': 'unknown'
            }
        
        # Normalize gene symbol
        gene_symbol = gene_symbol.strip().upper()
        
        # Check cache first
        cache_key = f"clingen_dosage_{gene_symbol}"
        cached = self._get_cached_response(cache_key)
        if cached:
            print(f"{Fore.CYAN}ℹ️  Using cached dosage sensitivity data for {gene_symbol}{Style.RESET_ALL}")
            return cached
        
        # HYBRID APPROACH: Try local file first, then FTP
        local_tsv_path = 'ClinGen_gene_curation_list_GRCh38.tsv'
        tsv_content = None
        data_source = None
        
        # Step 1: Try local TSV file
        if os.path.exists(local_tsv_path):
            try:
                print(f"{Fore.CYAN}� Reading local ClinGen TSV file...{Style.RESET_ALL}")
                with open(local_tsv_path, 'r', encoding='utf-8') as f:
                    tsv_content = f.read()
                data_source = 'Local ClinGen TSV (GRCh38)'
            except Exception as e:
                print(f"{Fore.YELLOW}⚠️  Could not read local file: {e}{Style.RESET_ALL}")
        
        # Step 2: If local file not available, download from FTP
        if tsv_content is None:
            try:
                print(f"{Fore.YELLOW}🌐 Downloading ClinGen TSV from FTP...{Style.RESET_ALL}")
                url = API_ENDPOINTS['clingen_dosage_tsv']
                
                response = requests.get(
                    url,
                    timeout=API_SETTINGS.get('timeout', 30),
                    headers={'Accept': 'text/tab-separated-values'}
                )
                
                if response.status_code == 200:
                    tsv_content = response.text
                    data_source = 'ClinGen FTP (GRCh38)'
                    
                    # Save to local file for future use
                    try:
                        with open(local_tsv_path, 'w', encoding='utf-8') as f:
                            f.write(tsv_content)
                        print(f"{Fore.GREEN}✅ Saved TSV file locally for future queries{Style.RESET_ALL}")
                    except Exception as e:
                        print(f"{Fore.YELLOW}⚠️  Could not save local file: {e}{Style.RESET_ALL}")
                else:
                    return {
                        'error': f'HTTP {response.status_code} - Could not download ClinGen TSV',
                        'haploinsufficiency_score': None,
                        'triplosensitivity_score': None,
                        'pvs1_recommendation': 'unknown',
                        'source': 'ClinGen Dosage Sensitivity Map'
                    }
            
            except requests.RequestException as e:
                return {
                    'error': f'Network error: {str(e)}',
                    'haploinsufficiency_score': None,
                    'triplosensitivity_score': None,
                    'pvs1_recommendation': 'unknown',
                    'source': 'ClinGen Dosage Sensitivity Map'
                }
        
        # Step 3: Parse TSV content and find gene
        try:
            print(f"{Fore.YELLOW}🔍 Searching for {gene_symbol} in ClinGen data...{Style.RESET_ALL}")
            
            lines = tsv_content.split('\n')
            
            for line in lines:
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.split('\t')
                if len(fields) < 21:
                    continue
                
                current_gene = fields[0].strip().upper()
                if current_gene == gene_symbol:
                    # Extract fields (TSV format)
                    gene_id = fields[1]
                    cytoband = fields[2]
                    genomic_location = fields[3]
                    
                    # Haploinsufficiency data
                    hi_score_str = fields[4].strip()
                    hi_description = fields[5].strip()
                    hi_pmids = [fields[i].strip() for i in range(6, 12) if i < len(fields) and fields[i].strip()]
                    
                    # Triplosensitivity data
                    ts_score_str = fields[12].strip() if len(fields) > 12 else '0'
                    ts_description = fields[13].strip() if len(fields) > 13 else 'Not evaluated'
                    
                    # Last evaluated date
                    last_evaluated = fields[20].strip() if len(fields) > 20 else 'Unknown'
                    
                    # Parse scores
                    try:
                        hi_score = int(hi_score_str) if hi_score_str and hi_score_str != 'Not yet evaluated' else None
                    except ValueError:
                        hi_score = None
                    
                    try:
                        ts_score = int(ts_score_str) if ts_score_str and ts_score_str != 'Not yet evaluated' else None
                    except ValueError:
                        ts_score = None
                    
                    # Map score to PVS1 recommendation
                    pvs1_rec = self._determine_pvs1_strength(hi_score)
                    
                    result = {
                        'haploinsufficiency_score': hi_score,
                        'triplosensitivity_score': ts_score,
                        'haploinsufficiency_description': hi_description,
                        'triplosensitivity_description': ts_description,
                        'pvs1_recommendation': pvs1_rec,
                        'confidence': self._dosage_confidence_level(hi_score),
                        'pmids': hi_pmids,
                        'gene_id': gene_id,
                        'cytoband': cytoband,
                        'genomic_location': genomic_location,
                        'source': data_source,
                        'last_evaluated': last_evaluated,
                        'clingen_gene_url': f"https://search.clinicalgenome.org/kb/gene-dosage/{gene_symbol}",
                        'gene_symbol': gene_symbol
                    }
                    
                    self._cache_response(cache_key, result)
                    
                    # Color-coded output
                    if hi_score == 3:
                        color = Fore.GREEN
                        msg = f"✅ ClinGen: {gene_symbol} HI Score={hi_score} (Sufficient) → Keep PVS1 Very Strong"
                    elif hi_score == 2:
                        color = Fore.YELLOW
                        msg = f"⚠️  ClinGen: {gene_symbol} HI Score={hi_score} (Some) → Downgrade PVS1 to PS1 Strong"
                    elif hi_score == 1:
                        color = Fore.YELLOW
                        msg = f"⚠️  ClinGen: {gene_symbol} HI Score={hi_score} (Little) → Downgrade PVS1 to PM2 Moderate"
                    elif hi_score == 30:
                        color = Fore.MAGENTA
                        msg = f"ℹ️  ClinGen: {gene_symbol} HI Score={hi_score} (AR phenotype) → Consider not applying PVS1"
                    elif hi_score == 40:
                        color = Fore.RED
                        msg = f"❌ ClinGen: {gene_symbol} HI Score={hi_score} (Unlikely) → Do not apply PVS1"
                    elif hi_score == 0:
                        color = Fore.RED
                        msg = f"❌ ClinGen: {gene_symbol} HI Score={hi_score} (No evidence) → Consider not applying PVS1"
                    else:
                        color = Fore.CYAN
                        msg = f"ℹ️  ClinGen: {gene_symbol} HI Score={hi_score} → Review manually"
                    
                    print(f"{color}{msg}{Style.RESET_ALL}")
                    return result
            
            # Gene not found
            result = {
                'error': f'Gene {gene_symbol} not found in ClinGen Dosage Sensitivity Map',
                'haploinsufficiency_score': None,
                'triplosensitivity_score': None,
                'pvs1_recommendation': 'insufficient_data',
                'confidence': 'none',
                'source': data_source,
                'gene_symbol': gene_symbol
            }
            self._cache_response(cache_key, result)
            print(f"{Fore.YELLOW}⚠️  {gene_symbol} not found in ClinGen dosage map{Style.RESET_ALL}")
            return result
        
        except Exception as e:
            result = {
                'error': f'Unexpected error: {str(e)}',
                'haploinsufficiency_score': None,
                'triplosensitivity_score': None,
                'pvs1_recommendation': 'unknown',
                'source': 'ClinGen Dosage Sensitivity Map'
            }
            print(f"{Fore.RED}❌ Unexpected error: {str(e)}{Style.RESET_ALL}")
            return result
    
    def _determine_pvs1_strength(self, hi_score: Optional[int]) -> str:
        """
        Determine PVS1 strength recommendation based on haploinsufficiency score.
        
        Args:
            hi_score: Haploinsufficiency score (0-3, 30, or 40)
            
        Returns:
            str: PVS1 strength recommendation
        """
        if hi_score is None:
            return 'insufficient_data'
        elif hi_score == 3:
            return 'pvs1_very_strong'  # Keep PVS1
        elif hi_score == 2:
            return 'ps1_strong'  # Downgrade to PS1
        elif hi_score == 1:
            return 'pm2_moderate'  # Downgrade to PM2
        elif hi_score == 0 or hi_score == 30 or hi_score == 40:
            return 'consider_not_applying'  # No evidence / AR phenotype / Unlikely
        else:
            return 'unknown'
    
    def _dosage_confidence_level(self, hi_score: Optional[int]) -> str:
        """
        Determine confidence level based on haploinsufficiency score.
        
        Args:
            hi_score: Haploinsufficiency score (0-3, 30, or 40)
            
        Returns:
            str: Confidence level
        """
        if hi_score is None:
            return 'none'
        elif hi_score == 3:
            return 'high'
        elif hi_score == 2:
            return 'medium'
        elif hi_score == 1:
            return 'low'
        elif hi_score == 0 or hi_score == 30 or hi_score == 40:
            return 'very_low'
        else:
            return 'none'






