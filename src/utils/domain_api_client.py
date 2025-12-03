"""
Domain and Hotspot API Client
=============================

Provides protein domain and mutational hotspot information from remote API sources.

DESIGN PRINCIPLE:
    This module MUST NOT contain hardcoded biological data (gene names, hotspot
    residues, domain coordinates, etc.). All such data MUST be fetched from 
    remote APIs at runtime. Local configuration is limited to:
    - API endpoints and timeouts
    - Confidence thresholds for interpretation
    - Cache settings

Data Sources:
    1. CancerHotspots.org - Cancer mutation hotspot positions
    2. UniProt REST API - Protein domain annotations
    3. (Future) ClinVar API - Pathogenic variant clustering

Author: Can Sevilmiş
License: MIT License
"""

import requests
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field
from datetime import datetime, timedelta
import json
import os


@dataclass
class HotspotAnnotation:
    """
    Result of hotspot/domain annotation lookup from remote APIs.
    
    This is the canonical return type for PM1 evidence evaluation.
    All fields are populated from remote API responses only.
    
    Attributes:
        is_hotspot: True if position is a documented mutational hotspot
        in_critical_domain: True if position is within a functional domain
        confidence: Float 0.0-1.0 indicating evidence strength
        source: API source(s) that provided the data
        details: Human-readable description of findings
        domain_name: Name of the domain if in_critical_domain is True
        hotspot_count: Number of mutations at this position (from CancerHotspots)
        tumor_count: Number of tumors with mutations at this position
        raw_api_response: Original API response for debugging
    """
    is_hotspot: bool = False
    in_critical_domain: bool = False
    confidence: float = 0.0
    source: str = "none"
    details: str = ""
    domain_name: Optional[str] = None
    hotspot_count: int = 0
    tumor_count: int = 0
    raw_api_response: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            'is_hotspot': self.is_hotspot,
            'in_critical_domain': self.in_critical_domain,
            'confidence': self.confidence,
            'source': self.source,
            'details': self.details,
            'domain_name': self.domain_name,
            'hotspot_count': self.hotspot_count,
            'tumor_count': self.tumor_count,
        }


class DomainAPIClient:
    """
    Client for fetching protein domain and hotspot information from remote APIs.
    
    IMPORTANT: This class does NOT contain any hardcoded biological data.
    All hotspot/domain information is fetched from external APIs:
    - CancerHotspots.org for mutation hotspots
    - UniProt REST API for protein domains
    
    If APIs are unavailable, methods return empty/negative results rather
    than falling back to hardcoded data. This ensures data freshness and
    avoids embedding biological facts that may become outdated.
    
    Usage:
        client = DomainAPIClient()
        annotation = client.get_hotspot_annotation("TP53", position=248)
        if annotation.is_hotspot:
            print(f"Hotspot found: {annotation.details}")
    """
    
    # API configuration (not biological data)
    CANCER_HOTSPOTS_BASE_URL = "https://www.cancerhotspots.org/api"
    UNIPROT_BASE_URL = "https://rest.uniprot.org"
    
    # Confidence thresholds for interpreting API responses
    # These are interpretive rules, not biological facts
    CONFIDENCE_THRESHOLDS = {
        'high_hotspot_tumor_count': 10,      # >=10 tumors = high confidence hotspot
        'moderate_hotspot_tumor_count': 3,   # >=3 tumors = moderate confidence
        'critical_domain_types': {'Domain', 'Active site', 'Binding site'},
        'functional_region_types': {'Region', 'Motif'},
    }
    
    def __init__(self, cache_enabled: bool = True, timeout: int = 10):
        """
        Initialize domain API client with caching.
        
        Args:
            cache_enabled: Whether to cache API responses locally
            timeout: HTTP request timeout in seconds
        """
        self.cache_enabled = cache_enabled
        self.timeout = timeout
        self.cache: Dict[str, Any] = {}
        
        # Cache file in temp directory to avoid permission issues
        import tempfile
        cache_dir = os.path.join(tempfile.gettempdir(), 'acmg_assistant')
        os.makedirs(cache_dir, exist_ok=True)
        self.cache_file = os.path.join(cache_dir, 'domain_api_cache.json')
        self.max_cache_age = timedelta(days=7)  # Shorter cache for API data
        
        if cache_enabled:
            self._load_cache()
    
    def _load_cache(self):
        """Load cache from file."""
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, 'r') as f:
                    self.cache = json.load(f)
                self._clean_expired_cache()
            except (json.JSONDecodeError, IOError):
                self.cache = {}
    
    def _save_cache(self):
        """Save cache to file."""
        if self.cache_enabled:
            try:
                with open(self.cache_file, 'w') as f:
                    json.dump(self.cache, f, indent=2)
            except IOError as e:
                # Don't fail, but inform user cache couldn't be saved
                try:
                    from config.constants import COLORAMA_COLORS
                    print(f"{COLORAMA_COLORS['YELLOW']}⚠ Cache could not be saved: {e}{COLORAMA_COLORS['RESET']}")
                except ImportError:
                    print(f"Warning: Cache could not be saved: {e}")
    
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
        """Get cached response if available."""
        if not self.cache_enabled or cache_key not in self.cache:
            return None
        
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
    
    # =========================================================================
    # PRIMARY API: get_hotspot_annotation
    # =========================================================================
    
    def get_hotspot_annotation(
        self, 
        gene: str, 
        position: Optional[int] = None,
        hgvs_p: Optional[str] = None
    ) -> HotspotAnnotation:
        """
        Get hotspot/domain annotation for a variant from remote APIs only.
        
        This is the primary method for PM1 evidence evaluation. It queries
        remote APIs and returns a structured annotation. NO hardcoded 
        biological data is used.
        
        Args:
            gene: Gene symbol (e.g., 'TP53')
            position: Amino acid position (e.g., 248)
            hgvs_p: Protein HGVS notation (e.g., 'p.Arg248Gln') - used to
                    extract position if not provided
        
        Returns:
            HotspotAnnotation with is_hotspot, in_critical_domain, confidence,
            and source information. Returns empty annotation (all False/0) if
            APIs are unavailable or return no data.
        
        Example:
            >>> client = DomainAPIClient()
            >>> ann = client.get_hotspot_annotation("TP53", position=248)
            >>> if ann.is_hotspot and ann.confidence >= 0.9:
            ...     print("PM1 applies")
        """
        # Extract position from hgvs_p if not provided
        if position is None and hgvs_p:
            position = self._extract_position_from_hgvs_p(hgvs_p)
        
        # Check cache first
        cache_key = f"annotation_{gene.upper()}_{position or 'none'}"
        cached = self._get_cached_response(cache_key)
        if cached:
            return HotspotAnnotation(**cached)
        
        # Initialize empty result
        annotation = HotspotAnnotation()
        sources: List[str] = []
        
        # Query CancerHotspots.org if position is available
        if position:
            hotspot_result = self._query_cancer_hotspots(gene, position)
            if hotspot_result:
                annotation.is_hotspot = True
                annotation.hotspot_count = hotspot_result.get('mutation_count', 0)
                annotation.tumor_count = hotspot_result.get('tumor_count', 0)
                annotation.raw_api_response['cancer_hotspots'] = hotspot_result
                sources.append('CancerHotspots.org')
                
                # Calculate confidence based on tumor count
                if annotation.tumor_count >= self.CONFIDENCE_THRESHOLDS['high_hotspot_tumor_count']:
                    annotation.confidence = 0.95
                elif annotation.tumor_count >= self.CONFIDENCE_THRESHOLDS['moderate_hotspot_tumor_count']:
                    annotation.confidence = 0.75
                else:
                    annotation.confidence = 0.5
        
        # Query UniProt for domain information
        domain_result = self._query_uniprot_domains(gene, position)
        if domain_result:
            annotation.raw_api_response['uniprot'] = domain_result
            sources.append('UniProt')
            
            if domain_result.get('in_domain'):
                annotation.in_critical_domain = True
                annotation.domain_name = domain_result.get('domain_name')
                domain_type = domain_result.get('domain_type', '')
                
                # Boost confidence if in critical domain type
                if domain_type in self.CONFIDENCE_THRESHOLDS['critical_domain_types']:
                    annotation.confidence = max(annotation.confidence, 0.8)
                elif domain_type in self.CONFIDENCE_THRESHOLDS['functional_region_types']:
                    annotation.confidence = max(annotation.confidence, 0.6)
        
        # Build source and details strings
        annotation.source = ' + '.join(sources) if sources else 'none'
        annotation.details = self._build_annotation_details(annotation, gene, position)
        
        # Cache the result
        self._cache_response(cache_key, annotation.to_dict())
        
        return annotation
    
    def _extract_position_from_hgvs_p(self, hgvs_p: str) -> Optional[int]:
        """Extract amino acid position from protein HGVS notation."""
        if not hgvs_p:
            return None
        import re
        # Match patterns like:
        # - p.Arg248Gln (three-letter code)
        # - p.R248Q (single-letter code)
        # - p.Arg248* (nonsense)
        # - p.248 (position only, rare)
        match = re.search(r'p\.(?:[A-Za-z]{1,3})?(\d+)', hgvs_p)
        if match:
            return int(match.group(1))
        return None
    
    def _build_annotation_details(
        self, 
        annotation: HotspotAnnotation, 
        gene: str, 
        position: Optional[int]
    ) -> str:
        """Build human-readable details string from annotation data."""
        parts = []
        
        if annotation.is_hotspot:
            msg = f"Position {position} in {gene} is a documented cancer hotspot"
            if annotation.tumor_count > 0:
                msg += f" ({annotation.tumor_count} tumors, {annotation.hotspot_count} mutations)"
            parts.append(msg)
        
        if annotation.in_critical_domain:
            msg = f"Position is within functional domain"
            if annotation.domain_name:
                msg += f": {annotation.domain_name}"
            parts.append(msg)
        
        if not parts:
            if annotation.source == 'none':
                return "No hotspot/domain data available from remote APIs"
            return f"No hotspot/domain evidence found (checked: {annotation.source})"
        
        return "; ".join(parts)
    
    # =========================================================================
    # LEGACY API: get_hotspot_info (backward compatibility)
    # =========================================================================
    
    def get_hotspot_info(self, gene: str, position: Optional[int] = None) -> Dict[str, Any]:
        """
        Get hotspot information for a gene/position.
        
        DEPRECATED: Use get_hotspot_annotation() instead for new code.
        This method is maintained for backward compatibility.
        
        Args:
            gene: Gene symbol (e.g., 'TP53')
            position: Optional amino acid position
        
        Returns:
            Dict with hotspot information (legacy format)
        """
        # Use new annotation method internally
        annotation = self.get_hotspot_annotation(gene, position)
        
        # Convert to legacy format
        result = {
            'is_hotspot_gene': annotation.is_hotspot or annotation.in_critical_domain,
            'hotspot_residues': [position] if annotation.is_hotspot and position else [],
            'domains': [annotation.domain_name] if annotation.domain_name else [],
            'position_is_hotspot': annotation.is_hotspot,
            'source': annotation.source,
        }
        
        if annotation.is_hotspot:
            result['hotspot_details'] = {
                'tumor_count': annotation.tumor_count,
                'mutation_count': annotation.hotspot_count,
            }
        
        return result
    
    # =========================================================================
    # API Query Methods (Internal)
    # =========================================================================
    
    def _query_cancer_hotspots(self, gene: str, position: int) -> Optional[Dict[str, Any]]:
        """
        Query CancerHotspots.org API for mutation hotspots.
        
        Args:
            gene: Gene symbol
            position: Amino acid position
        
        Returns:
            Dict with hotspot data or None if not a hotspot/error
        """
        try:
            url = f"{self.CANCER_HOTSPOTS_BASE_URL}/hotspots/single/{gene.upper()}/{position}"
            response = requests.get(url, timeout=self.timeout)
            
            if response.status_code == 200:
                data = response.json()
                if data and len(data) > 0:
                    return {
                        'tumor_count': data[0].get('tumorCount', 0),
                        'mutation_count': data[0].get('count', 0),
                        'residue': position,
                    }
        except (requests.RequestException, KeyError, IndexError, ValueError):
            # API unavailable or error - return None (no fallback to hardcoded data)
            pass
        
        return None
    
    def _query_uniprot_domains(
        self, 
        gene: str, 
        position: Optional[int] = None
    ) -> Optional[Dict[str, Any]]:
        """
        Query UniProt REST API for protein domain information.
        
        Args:
            gene: Gene symbol
            position: Optional amino acid position to check domain membership
        
        Returns:
            Dict with domain data or None if error
        """
        try:
            # Search for UniProt accession using gene name
            search_url = (
                f"{self.UNIPROT_BASE_URL}/uniprotkb/search"
                f"?query=gene:{gene}+AND+organism_id:9606&format=json&size=1"
            )
            response = requests.get(search_url, timeout=self.timeout)
            
            if response.status_code != 200:
                return None
            
            search_data = response.json()
            if not search_data.get('results'):
                return None
            
            # Get primary accession
            accession = search_data['results'][0].get('primaryAccession')
            if not accession:
                return None
            
            # Fetch detailed entry
            entry_url = f"{self.UNIPROT_BASE_URL}/uniprotkb/{accession}.json"
            response = requests.get(entry_url, timeout=self.timeout)
            
            if response.status_code != 200:
                return None
            
            entry_data = response.json()
            
            # Extract domain information
            domain_types = (
                self.CONFIDENCE_THRESHOLDS['critical_domain_types'] |
                self.CONFIDENCE_THRESHOLDS['functional_region_types']
            )
            
            domains = []
            features = entry_data.get('features', [])
            
            for feature in features:
                feature_type = feature.get('type', '')
                if feature_type in domain_types:
                    description = feature.get('description', '')
                    start = feature.get('location', {}).get('start', {}).get('value')
                    end = feature.get('location', {}).get('end', {}).get('value')
                    
                    if description and start and end:
                        domains.append({
                            'type': feature_type,
                            'description': description,
                            'start': start,
                            'end': end,
                        })
            
            result = {
                'accession': accession,
                'domains': domains,
                'in_domain': False,
                'domain_name': None,
                'domain_type': None,
            }
            
            # Check if position is in any domain
            if position:
                for domain in domains:
                    if domain['start'] <= position <= domain['end']:
                        result['in_domain'] = True
                        result['domain_name'] = domain['description']
                        result['domain_type'] = domain['type']
                        break
            
            return result
        
        except (requests.RequestException, KeyError, ValueError, TypeError):
            # API unavailable or error - return None
            pass
        
        return None
    
    # =========================================================================
    # Utility Methods
    # =========================================================================
    
    def check_position_in_domain(self, gene: str, position: int) -> Dict[str, Any]:
        """
        Check if a specific position is within a functional domain.
        
        Args:
            gene: Gene symbol
            position: Amino acid position
        
        Returns:
            Dict with domain membership information
        """
        result = {
            'in_domain': False,
            'domain_name': None,
            'domain_type': None
        }
        
        # Use internal query method
        uniprot_data = self._query_uniprot_domains(gene, position)
        if uniprot_data:
            result['in_domain'] = uniprot_data.get('in_domain', False)
            result['domain_name'] = uniprot_data.get('domain_name')
            result['domain_type'] = uniprot_data.get('domain_type')
        
        return result
    
    def is_api_available(self) -> Dict[str, bool]:
        """
        Check availability of remote APIs.
        
        Returns:
            Dict with API name -> availability status
        """
        status = {
            'cancer_hotspots': False,
            'uniprot': False,
        }
        
        try:
            # Check CancerHotspots.org
            response = requests.get(
                f"{self.CANCER_HOTSPOTS_BASE_URL}/hotspots/single/TP53/248",
                timeout=5
            )
            status['cancer_hotspots'] = response.status_code == 200
        except requests.RequestException:
            pass
        
        try:
            # Check UniProt
            response = requests.get(
                f"{self.UNIPROT_BASE_URL}/uniprotkb/search?query=gene:TP53&format=json&size=1",
                timeout=5
            )
            status['uniprot'] = response.status_code == 200
        except requests.RequestException:
            pass
        
        return status
