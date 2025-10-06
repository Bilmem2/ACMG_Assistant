"""
Domain and Hotspot API Client
Provides protein domain and mutational hotspot information from multiple sources.
"""

import requests
from typing import Dict, Any, List, Optional
from datetime import datetime, timedelta
import json
import os

class DomainAPIClient:
    """
    Client for fetching protein domain and hotspot information.
    
    Uses multiple data sources:
    1. UniProt - Protein domain annotations
    2. CancerHotspots.org - Cancer mutation hotspots
    3. ClinVar - Pathogenic variant clustering
    4. Fallback - Hardcoded known hotspot genes
    """
    
    # Known hotspot genes (fallback)
    KNOWN_HOTSPOTS = {
        'TP53': {
            'domains': ['DNA_binding_domain', 'tetramerization_domain'],
            'hotspot_residues': [175, 245, 248, 249, 273, 282],
            'source': 'Literature'
        },
        'KRAS': {
            'domains': ['GTPase_domain'],
            'hotspot_residues': [12, 13, 61, 146],
            'source': 'Literature'
        },
        'PIK3CA': {
            'domains': ['helical_domain', 'kinase_domain'],
            'hotspot_residues': [542, 545, 546, 1047],
            'source': 'Literature'
        },
        'BRAF': {
            'domains': ['kinase_domain'],
            'hotspot_residues': [600],
            'source': 'Literature'
        },
        'EGFR': {
            'domains': ['kinase_domain', 'tyrosine_kinase_domain'],
            'hotspot_residues': [719, 746, 790, 858],
            'source': 'Literature'
        },
        'PTEN': {
            'domains': ['phosphatase_domain', 'C2_domain'],
            'hotspot_residues': [130, 233],
            'source': 'Literature'
        },
        'BRCA1': {
            'domains': ['RING_domain', 'BRCT_domain'],
            'hotspot_residues': [],
            'source': 'Literature'
        },
        'BRCA2': {
            'domains': ['DNA_binding_domain'],
            'hotspot_residues': [],
            'source': 'Literature'
        }
    }
    
    def __init__(self, cache_enabled: bool = True):
        """Initialize domain API client with caching."""
        self.cache_enabled = cache_enabled
        self.cache = {}
        
        # Cache file in user's home directory or temp directory to avoid permission issues
        import tempfile
        cache_dir = os.path.join(tempfile.gettempdir(), 'acmg_assistant')
        os.makedirs(cache_dir, exist_ok=True)
        self.cache_file = os.path.join(cache_dir, 'domain_api_cache.json')
        self.max_cache_age = timedelta(days=30)  # Domain data rarely changes
        
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
    
    def get_hotspot_info(self, gene: str, position: Optional[int] = None) -> Dict[str, Any]:
        """
        Get comprehensive hotspot information for a gene/position.
        
        Args:
            gene: Gene symbol (e.g., 'TP53')
            position: Optional amino acid position
        
        Returns:
            Dict with hotspot information:
            {
                'is_hotspot_gene': bool,
                'hotspot_residues': List[int],
                'domains': List[str],
                'position_is_hotspot': bool (if position provided),
                'source': str
            }
        """
        cache_key = f"hotspot_{gene.upper()}_{position or 'all'}"
        cached = self._get_cached_response(cache_key)
        if cached:
            return cached
        
        result = {
            'is_hotspot_gene': False,
            'hotspot_residues': [],
            'domains': [],
            'position_is_hotspot': False,
            'source': 'none'
        }
        
        # Try multiple sources in order
        # 1. CancerHotspots.org (if position provided)
        if position:
            cancer_hotspot = self._check_cancer_hotspots(gene, position)
            if cancer_hotspot:
                result.update(cancer_hotspot)
                self._cache_response(cache_key, result)
                return result
        
        # 2. UniProt API for domains
        uniprot_data = self._get_uniprot_domains(gene)
        if uniprot_data:
            result['domains'] = uniprot_data.get('domains', [])
            result['is_hotspot_gene'] = len(result['domains']) > 0
            result['source'] = 'UniProt'
        
        # 3. Fallback to known hotspots
        if gene.upper() in self.KNOWN_HOTSPOTS:
            known = self.KNOWN_HOTSPOTS[gene.upper()]
            result['is_hotspot_gene'] = True
            result['hotspot_residues'] = known['hotspot_residues']
            if not result['domains']:
                result['domains'] = known['domains']
            result['source'] = f"{result.get('source', '')} + {known['source']}"
            
            if position and position in known['hotspot_residues']:
                result['position_is_hotspot'] = True
        
        self._cache_response(cache_key, result)
        return result
    
    def _check_cancer_hotspots(self, gene: str, position: int) -> Optional[Dict[str, Any]]:
        """
        Check CancerHotspots.org for mutation hotspots.
        
        Args:
            gene: Gene symbol
            position: Amino acid position
        
        Returns:
            Dict with hotspot info or None
        """
        try:
            url = f"https://www.cancerhotspots.org/api/hotspots/single/{gene.upper()}/{position}"
            response = requests.get(url, timeout=5)
            
            if response.status_code == 200:
                data = response.json()
                if data and len(data) > 0:
                    return {
                        'is_hotspot_gene': True,
                        'position_is_hotspot': True,
                        'hotspot_residues': [position],
                        'domains': [],
                        'source': 'CancerHotspots.org',
                        'hotspot_details': {
                            'tumor_count': data[0].get('tumorCount', 0),
                            'mutation_count': data[0].get('count', 0)
                        }
                    }
        except (requests.RequestException, KeyError, IndexError):
            pass
        
        return None
    
    def _get_uniprot_domains(self, gene: str) -> Optional[Dict[str, Any]]:
        """
        Fetch protein domain information from UniProt.
        
        Args:
            gene: Gene symbol
        
        Returns:
            Dict with domain information or None
        """
        try:
            # Search for UniProt accession using gene name
            search_url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene}+AND+organism_id:9606&format=json&size=1"
            response = requests.get(search_url, timeout=10)
            
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
            entry_url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
            response = requests.get(entry_url, timeout=10)
            
            if response.status_code != 200:
                return None
            
            entry_data = response.json()
            
            # Extract domain information
            domains = []
            features = entry_data.get('features', [])
            
            for feature in features:
                feature_type = feature.get('type', '')
                if feature_type in ['Domain', 'Region', 'Motif', 'Active site', 'Binding site']:
                    description = feature.get('description', '')
                    if description:
                        domains.append({
                            'type': feature_type,
                            'description': description,
                            'start': feature.get('location', {}).get('start', {}).get('value'),
                            'end': feature.get('location', {}).get('end', {}).get('value')
                        })
            
            return {
                'domains': [d['description'] for d in domains],
                'domain_details': domains,
                'accession': accession
            }
        
        except (requests.RequestException, KeyError, ValueError):
            return None
    
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
        
        uniprot_data = self._get_uniprot_domains(gene)
        if not uniprot_data or 'domain_details' not in uniprot_data:
            return result
        
        for domain in uniprot_data['domain_details']:
            start = domain.get('start')
            end = domain.get('end')
            
            if start and end and start <= position <= end:
                result['in_domain'] = True
                result['domain_name'] = domain.get('description')
                result['domain_type'] = domain.get('type')
                break
        
        return result
    
    def get_hotspot_residues_from_clinvar(self, gene: str, threshold: int = 5) -> List[int]:
        """
        Identify hotspot residues from ClinVar pathogenic variant clustering.
        
        Args:
            gene: Gene symbol
            threshold: Minimum number of pathogenic variants at residue
        
        Returns:
            List of hotspot residue positions
        """
        # This would require ClinVar API integration
        # For now, return empty list (can be implemented later)
        return []
