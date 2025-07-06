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
    
    def get_chromosome_from_ensembl(self, gene_symbol: str) -> Optional[str]:
        """
        Get chromosome information for a gene from Ensembl.
        
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
        
        try:
            print(f"{COLORAMA_COLORS['BLUE']}🔍 Fetching chromosome info for {gene_symbol} from Ensembl...{COLORAMA_COLORS['RESET']}")
            server = API_ENDPOINTS['ensembl']
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
                else:
                    print(f"{COLORAMA_COLORS['YELLOW']}⚠️ No chromosome information found for {gene_symbol}{COLORAMA_COLORS['RESET']}")
                    return None
            elif response.status_code == 404:
                print(f"{COLORAMA_COLORS['YELLOW']}⚠️ Gene {gene_symbol} not found in Ensembl database{COLORAMA_COLORS['RESET']}")
                return None
            else:
                print(f"{COLORAMA_COLORS['RED']}❌ Ensembl API returned status {response.status_code}{COLORAMA_COLORS['RESET']}")
                return None
                
        except requests.Timeout:
            print(f"{COLORAMA_COLORS['RED']}⏱️ Timeout while fetching data from Ensembl{COLORAMA_COLORS['RESET']}")
            return None
        except requests.ConnectionError:
            print(f"{COLORAMA_COLORS['RED']}🌐 Connection error - please check your internet connection{COLORAMA_COLORS['RESET']}")
            return None
        except requests.RequestException as e:
            print(f"{COLORAMA_COLORS['RED']}❌ Error fetching chromosome from Ensembl: {e}{COLORAMA_COLORS['RESET']}")
            return None
        except Exception as e:
            print(f"{COLORAMA_COLORS['RED']}❌ Unexpected error: {e}{COLORAMA_COLORS['RESET']}")
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
            
            # Extract clinical significance
            clinical_significance = variant_info.get('clinical_significance', '')
            
            # Extract review status
            review_status = variant_info.get('review_status', '')
            
            # Extract last evaluated date
            last_evaluated = variant_info.get('last_evaluated', '')
            
            # Extract variant title
            title = variant_info.get('title', '')
            
            return {
                'status': 'found',
                'clinvar_id': clinvar_id,
                'significance': clinical_significance,
                'review_status': review_status,
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
        Get comprehensive gene information from Ensembl.
        
        Args:
            gene_symbol (str): Gene symbol
            
        Returns:
            Dict[str, Any]: Gene information
        """
        cache_key = f"ensembl_gene_{gene_symbol.upper()}"
        cached_result = self._get_cached_response(cache_key)
        
        if cached_result is not None:
            return cached_result
        
        try:
            server = API_ENDPOINTS['ensembl']
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
            else:
                return {'error': f"Gene {gene_symbol} not found"}
                
        except requests.RequestException as e:
            print(f"Error fetching gene information: {e}")
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
        
        return {
            'total_entries': total_entries,
            'ensembl_entries': ensembl_entries,
            'clinvar_entries': clinvar_entries,
            'cache_file': self.cache_file,
            'max_age_hours': OUTPUT_SETTINGS['max_cache_age_hours']
        }
