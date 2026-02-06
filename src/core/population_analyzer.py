# =============================================================================
# Population Analyzer Module
# =============================================================================
# Changes in this documentation pass (Dec 2024):
# - Added module docstring
# - Added class and method docstrings with type hints
# - Documented the ethnicity-aware analysis approach
#
# Changes (Jan 2025):
# - Refactored to use pre-fetched PopulationStats from multi-source API
# - Removed local JSON file dependency for population frequencies
# - Added backward compatibility with legacy population_data dict
# - Kept local ethnicity thresholds (configuration data, not variant data)
# =============================================================================
"""
Population Analyzer Module
==========================

Provides population frequency analysis for ACMG criteria evaluation,
including ethnicity-aware threshold adjustments for BA1, BS1, and PM2.

Multi-Source Data Flow (Jan 2025):
- Population data is pre-fetched by EvidenceEvaluator._fetch_external_data()
- This module acts as a "pure interpreter" of pre-fetched PopulationStats
- Falls back to legacy population_data dict for backward compatibility
- Ethnicity thresholds remain local (configuration, not variant data)

Classes:
    GnomADClient: Client for retrieving population data (now uses API data)
    EthnicityAwarePopulationAnalyzer: Main analyzer with ethnicity adjustments
"""

import json
from pathlib import Path
from typing import Dict, Optional, Any


class GnomADClient:
    """
    Client for retrieving population frequency data.
    
    Multi-Source Data Flow:
    - Primary: Uses pre-fetched PopulationStats from variant_data.population_stats
    - Fallback: Legacy population_data dict for backward compatibility
    - Last resort: Local JSON cache (deprecated)
    
    Note: For live gnomAD API queries, see utils/api_client.py and
    utils/predictor_api_client.py.
    """
    
    def get_population_frequencies(self, variant_data) -> Dict[str, float]:
        """
        Get population frequencies for a variant.
        
        This method now primarily uses pre-fetched PopulationStats from the
        multi-source API, with fallbacks for backward compatibility.
        
        Args:
            variant_data: VariantData object with population_stats, population_data,
                         or gene/hgvs_c properties
            
        Returns:
            Dict mapping population codes to allele frequencies
        """
        frequencies = {}
        
        # Prefer manual population_data when provided by user
        if not frequencies:
            population_data = getattr(variant_data, 'population_data', {}) or {}
            if self._has_manual_population_data(population_data):
                # Map common keys to frequency dict
                if population_data.get('gnomad_af') is not None:
                    frequencies['ALL'] = population_data['gnomad_af']
                if population_data.get('gnomad_af_popmax') is not None:
                    frequencies['POPMAX'] = population_data['gnomad_af_popmax']
                
                # Check for ethnicity-specific frequencies
                ethnicity_keys = {
                    'gnomad_af_nfe': 'NFE',
                    'gnomad_af_afr': 'AFR',
                    'gnomad_af_eas': 'EAS',
                    'gnomad_af_amr': 'AMR',
                    'gnomad_af_sas': 'SAS',
                    'gnomad_af_fin': 'FIN',
                    'gnomad_af_asj': 'ASJ',
                }
                for key, code in ethnicity_keys.items():
                    if population_data.get(key) is not None:
                        frequencies[code] = population_data[key]

        # Use pre-fetched PopulationStats (typed API data) if no manual override
        if not frequencies:
            population_stats = getattr(variant_data, 'population_stats', None)
            if population_stats:
                for source, stats in population_stats.items():
                    # Handle PopulationStats dataclass
                    if hasattr(stats, 'af') and stats.af is not None:
                        frequencies['ALL'] = stats.af
                        
                        # Extract sub-population frequencies if available
                        subpop = getattr(stats, 'subpop', None)
                        if subpop:
                            for pop_code, pop_data in subpop.items():
                                if isinstance(pop_data, dict):
                                    pop_af = pop_data.get('af')
                                else:
                                    pop_af = getattr(pop_data, 'af', None)
                                if pop_af is not None:
                                    frequencies[pop_code.upper()] = pop_af
                        
                        # Use popmax if available
                        popmax_af = getattr(stats, 'popmax_af', None)
                        if popmax_af is not None:
                            frequencies['POPMAX'] = popmax_af
                        
                        break  # Use first available source
                    
                    # Handle dict format
                    elif isinstance(stats, dict):
                        if stats.get('af') is not None:
                            frequencies['ALL'] = stats['af']
                            
                            # Extract sub-populations
                            subpop = stats.get('subpop', {})
                            for pop_code, pop_data in subpop.items():
                                if isinstance(pop_data, dict) and pop_data.get('af') is not None:
                                    frequencies[pop_code.upper()] = pop_data['af']
                            
                            if stats.get('popmax_af') is not None:
                                frequencies['POPMAX'] = stats['popmax_af']
                            
                            break
        
        # Last resort: Try local JSON cache (deprecated, for legacy compatibility)
        if not frequencies:
            gene = getattr(variant_data, 'gene', None)
            hgvs_c = getattr(variant_data, 'hgvs_c', None)
            
            if gene and hgvs_c:
                key = f"{gene}:{hgvs_c}"
                try:
                    fallback_path = Path(__file__).resolve().parent.parent.parent / 'data' / 'population_data' / 'population_frequencies.json'
                    with open(fallback_path, 'r') as f:
                        freq_db = json.load(f)
                    frequencies = freq_db.get(key, {})
                except Exception:
                    pass
        
        return frequencies

    def _has_manual_population_data(self, population_data: Dict[str, Any]) -> bool:
        """Check if the user provided any population data overrides."""
        manual_keys = [
            'gnomad_af', 'gnomad_af_popmax', 'exac_af',
            'gnomad_af_nfe', 'gnomad_af_afr', 'gnomad_af_eas',
            'gnomad_af_amr', 'gnomad_af_sas', 'gnomad_af_fin',
            'gnomad_af_asj'
        ]
        return any(population_data.get(key) is not None for key in manual_keys)


class EthnicityAwarePopulationAnalyzer:
    """
    Population frequency analyzer with ethnicity-specific thresholds.
    
    Applies different frequency thresholds based on patient ethnicity
    to account for population-specific allele frequency variations.
    
    Attributes:
        gnomad_client: GnomADClient instance for frequency lookups
        ethnicity_thresholds: Dict mapping ethnicity codes to BA1/BS1 thresholds
    """
    
    def __init__(self):
        """Initialize analyzer with gnomAD client and ethnicity thresholds."""
        self.gnomad_client = GnomADClient()
        self.ethnicity_thresholds = self._load_ethnicity_thresholds()
    
    def _load_ethnicity_thresholds(self) -> Dict[str, float]:
        """
        Load ethnicity-specific frequency thresholds.
        
        Returns:
            Dict mapping ethnicity codes to threshold values
        """
        try:
            thresholds_path = Path(__file__).resolve().parent.parent.parent / 'data' / 'population_data' / 'ethnicity_thresholds.json'
            with open(thresholds_path, 'r') as f:
                return json.load(f)
        except Exception:
            # Default thresholds if file not found
            return {'ALL': 0.01, 'EUR': 0.01, 'AFR': 0.01}

    
    def analyze_population_frequency(self, variant_data, 
                                     patient_ethnicity: Optional[str] = None) -> Optional[str]:
        """
        Analyze population frequency with ethnicity-aware thresholds.
        
        Evaluates variant frequency against population-specific thresholds
        to determine applicable ACMG criteria (BA1, BS1, or PM2).
        
        Args:
            variant_data: VariantData object containing variant information
            patient_ethnicity: Optional ethnicity code (e.g., 'EUR', 'AFR', 'EAS')
            
        Returns:
            str or None: 'BA1' if frequency very high (benign stand-alone),
                        'BS1' if frequency high (benign strong),
                        'PM2' if frequency very low/absent (pathogenic moderate),
                        None if frequency is intermediate
        """
        frequencies = self.gnomad_client.get_population_frequencies(variant_data)
        
        if patient_ethnicity:
            relevant_freq = frequencies.get(patient_ethnicity, frequencies.get('ALL'))
            threshold = self.ethnicity_thresholds.get(patient_ethnicity, 0.01)
        else:
            relevant_freq = frequencies.get('ALL')
            threshold = 0.005
        if relevant_freq is None:
            return None
        if relevant_freq >= threshold:
            return 'BA1'
        if relevant_freq >= threshold * 0.1:
            return 'BS1'
        if relevant_freq <= 0.0001:
            return 'PM2'
        return None
