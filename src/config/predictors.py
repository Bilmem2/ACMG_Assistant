"""
In Silico Predictor Configuration
================================

Contains weights, thresholds, and data structures for computational 
pathogenicity predictors and population frequency data.

Multi-source data fetching is supported with configurable source priority.
All factual variant-level data is fetched via external APIs - no hardcoding.

Author: Can Sevilmiş
License: MIT License
"""

from dataclasses import dataclass, field
from typing import Optional, Any


# =============================================================================
# Typed Data Structures for Multi-Source API Responses
# =============================================================================

@dataclass
class PredictorScore:
    """
    Typed container for in silico predictor scores from external APIs.
    
    This dataclass ensures consistent handling of predictor data fetched
    from multiple sources (dbNSFP, AlphaMissense API, etc.).
    
    Attributes:
        predictor: Name of the predictor (e.g., 'revel', 'cadd_phred')
        value: Numeric score (None if unavailable from all sources)
        source: API source that provided this score (e.g., 'dbNSFP', 'AlphaMissense_API')
        version: Version of the source database/API
        raw: Raw API response for debugging/auditing
        is_inverted: True for predictors where lower = pathogenic (SIFT, FATHMM)
    """
    predictor: str
    value: Optional[float] = None
    source: Optional[str] = None
    version: Optional[str] = None
    raw: Optional[dict] = None
    is_inverted: bool = False
    
    def is_available(self) -> bool:
        """Check if a valid score is available."""
        return self.value is not None
    
    def get_normalized_value(self) -> Optional[float]:
        """
        Get normalized score where higher = more pathogenic.
        Handles inverted predictors like SIFT and FATHMM.
        """
        if self.value is None:
            return None
        if self.is_inverted:
            # For SIFT: 0 = deleterious, 1 = tolerated → invert
            # For FATHMM: negative = deleterious → shift and normalize
            if self.predictor == 'sift':
                return 1.0 - self.value
            elif self.predictor == 'fathmm':
                # FATHMM ranges roughly from -10 to +10, normalize to 0-1
                # More negative = more deleterious
                return max(0.0, min(1.0, (5.0 - self.value) / 10.0))
        return self.value


@dataclass
class PopulationStats:
    """
    Typed container for population allele frequency data from external APIs.
    
    This dataclass ensures consistent handling of population frequency data
    fetched from gnomAD and other population databases.
    
    Attributes:
        population: Population database name (e.g., 'gnomad_v4', 'exac')
        af: Overall allele frequency
        an: Total allele number (sample size)
        ac: Allele count
        homozygote_count: Number of homozygous individuals
        hemizygote_count: Number of hemizygous individuals (X/Y chr)
        subpop: Sub-population frequencies (e.g., {'afr': 0.001, 'eas': 0.002})
        popmax_af: Maximum AF across sub-populations
        popmax_population: Sub-population with maximum AF
        filters: Quality filters applied
        source: API source (e.g., 'gnomAD_GraphQL')
        version: Database version (e.g., 'v4.0')
        raw: Raw API response for debugging/auditing
    """
    population: str
    af: Optional[float] = None
    an: Optional[int] = None
    ac: Optional[int] = None
    homozygote_count: Optional[int] = None
    hemizygote_count: Optional[int] = None
    subpop: Optional[dict[str, dict]] = None
    popmax_af: Optional[float] = None
    popmax_population: Optional[str] = None
    filters: Optional[list[str]] = None
    source: Optional[str] = None
    version: Optional[str] = None
    raw: Optional[dict] = None
    
    def is_available(self) -> bool:
        """Check if frequency data is available."""
        return self.af is not None or self.ac is not None
    
    def is_absent(self) -> bool:
        """Check if variant is absent from this population (AF=0 or not found)."""
        if self.af is not None:
            return self.af == 0.0
        if self.ac is not None:
            return self.ac == 0
        return True  # No data means effectively absent
    
    def get_max_subpop_af(self) -> Optional[float]:
        """Get maximum allele frequency across sub-populations."""
        if self.popmax_af is not None:
            return self.popmax_af
        if self.subpop:
            afs = [sp.get('af') for sp in self.subpop.values() if sp.get('af') is not None]
            return max(afs) if afs else None
        return self.af


# =============================================================================
# Multi-Source Predictor Configuration
# =============================================================================

# Priority order for predictor data sources per predictor type.
# First available source with valid data is used.
# Format: predictor_name -> list of source names in priority order
PREDICTOR_SOURCE_PRIORITY = {
    # AlphaMissense has a dedicated API, use it first
    'alphamissense': ['AlphaMissense_API', 'dbNSFP', 'VEP'],
    
    # REVEL and CADD are widely available
    'revel': ['dbNSFP', 'VEP', 'UCSC'],
    'cadd_phred': ['CADD_API', 'dbNSFP', 'VEP', 'UCSC'],
    
    # Classic predictors
    'sift': ['dbNSFP', 'VEP', 'Ensembl'],
    'polyphen2': ['dbNSFP', 'VEP', 'Ensembl'],
    
    # Meta-predictors
    'metasvm': ['dbNSFP'],
    'vest4': ['dbNSFP'],
    'fathmm': ['dbNSFP', 'VEP'],
    
    # Additional predictors for future expansion
    'bayesdel': ['dbNSFP'],
    'primateai': ['dbNSFP'],
    'mpc': ['dbNSFP'],
}

# Population frequency source priority
POPULATION_SOURCE_PRIORITY = {
    'gnomad': ['gnomAD_GraphQL', 'gnomAD_REST', 'Ensembl'],
    'exac': ['ExAC_REST', 'dbNSFP'],
    'topmed': ['TOPMed_API', 'dbNSFP'],
    '1000g': ['Ensembl', 'dbNSFP'],
}

# Predictors that use inverted scoring (lower = more pathogenic)
INVERTED_PREDICTORS = {'sift', 'fathmm'}

# Minimum number of predictors required for reliable composite score
MIN_PREDICTORS_FOR_COMPOSITE = 3


# =============================================================================
# In Silico Predictor Weights and Thresholds
# =============================================================================

# In silico predictor weights for Computational Metascore
INSILICO_WEIGHTS = {
    'revel': 0.25,
    'cadd_phred': 0.20,
    'alphamissense': 0.15,
    'sift': 0.10,
    'polyphen2': 0.10,
    'metasvm': 0.10,
    'vest4': 0.05,
    'fathmm': 0.05
}

# In silico predictor thresholds
INSILICO_THRESHOLDS = {
    'revel': {'pathogenic': 0.75, 'benign': 0.25},
    'cadd_phred': {'pathogenic': 25, 'benign': 10},
    'alphamissense': {'pathogenic': 0.564, 'benign': 0.34},
    'sift': {'pathogenic': 0.05, 'benign': 0.95},  # SIFT is inverted
    'polyphen2': {'pathogenic': 0.85, 'benign': 0.15},
    'metasvm': {'pathogenic': 0.83, 'benign': 0.17},
    'vest4': {'pathogenic': 0.7, 'benign': 0.3},
    'fathmm': {'pathogenic': -1.5, 'benign': 1.5}  # FATHMM is inverted
}

# VAMPP score thresholds for enhanced metascore
VAMPP_SCORE_THRESHOLDS = {
    'pathogenic_strong': 0.85,
    'pathogenic_moderate': 0.7,
    'benign_moderate': 0.3,
    'benign_strong': 0.15
}

# =============================================================================
# Missense Composite Score Thresholds
# =============================================================================
# These thresholds map the composite missense score (0.0-1.0) to ACMG evidence
# categories for PP3 (pathogenic computational) and BP4 (benign computational).
#
# IMPORTANT: These are research/educational approximations, NOT clinically
# validated thresholds. They should be calibrated against known pathogenic
# and benign variant datasets before clinical use.
#
# Score interpretation:
#   - Higher scores (closer to 1.0) = more likely damaging/pathogenic
#   - Lower scores (closer to 0.0) = more likely benign/tolerated
# =============================================================================
MISSENSE_COMPOSITE_THRESHOLDS = {
    # Pathogenic evidence thresholds (PP3)
    'PP3_strong': 0.90,       # Very high confidence damaging
    'PP3_moderate': 0.75,     # High confidence damaging
    'PP3_supporting': 0.60,   # Moderate confidence damaging
    
    # Neutral zone - no clear evidence either way
    'neutral_upper': 0.60,    # Below this, not enough for PP3
    'neutral_lower': 0.40,    # Above this, not enough for BP4
    
    # Benign evidence thresholds (BP4)
    'BP4_supporting': 0.40,   # Moderate confidence benign
    'BP4_moderate': 0.25,     # High confidence benign
    'BP4_strong': 0.10,       # Very high confidence benign
}

# Weights for combining sub-scores into composite missense score
# These weights reflect the relative importance of each evidence type
# based on literature and clinical practice guidelines.
MISSENSE_SCORE_WEIGHTS = {
    'conservation': 0.25,     # Evolutionary conservation (phyloP, GERP, phastCons)
    'functional': 0.35,       # In silico functional predictors (REVEL, CADD, etc.)
    'structural': 0.15,       # Structural impact (protein stability, interactions)
    'domain': 0.15,           # Functional domain context
    'population': 0.10,       # Population frequency context
}


# =============================================================================
# Validation Functions for Cached Data
# =============================================================================

# Valid ranges for predictor scores
PREDICTOR_VALUE_RANGES = {
    # Standard 0-1 range predictors
    'revel': (0.0, 1.0),
    'alphamissense': (0.0, 1.0),
    'sift': (0.0, 1.0),
    'polyphen2': (0.0, 1.0),
    'metasvm': (0.0, 1.0),
    'vest4': (0.0, 1.0),
    'bayesdel': (0.0, 1.0),
    'primateai': (0.0, 1.0),
    
    # CADD PHRED scores: 0 to ~60
    'cadd_phred': (0.0, 60.0),
    
    # FATHMM: typically -20 to +20
    'fathmm': (-20.0, 20.0),
    
    # MPC: 0 to ~5
    'mpc': (0.0, 5.0),
}


def validate_predictor_score(
    predictor: str,
    value: Optional[float],
    strict: bool = False
) -> bool:
    """
    Validate a predictor score value.
    
    Args:
        predictor: Name of the predictor (e.g., 'revel', 'cadd_phred')
        value: Numeric score value (can be None)
        strict: If True, requires value to be present; if False, None is valid
        
    Returns:
        True if value is valid, False otherwise
        
    Examples:
        >>> validate_predictor_score('revel', 0.85)
        True
        >>> validate_predictor_score('revel', 1.5)  # Out of range
        False
        >>> validate_predictor_score('revel', None)  # None is valid
        True
        >>> validate_predictor_score('revel', None, strict=True)  # strict mode
        False
    """
    # None is valid unless strict mode
    if value is None:
        return not strict
    
    # Check type
    if not isinstance(value, (int, float)):
        return False
    
    # Check for NaN or Inf
    import math
    if math.isnan(value) or math.isinf(value):
        return False
    
    # Normalize predictor name
    predictor_lower = predictor.lower().replace('-', '_').replace(' ', '_')
    
    # Get valid range
    valid_range = PREDICTOR_VALUE_RANGES.get(predictor_lower)
    
    if valid_range:
        min_val, max_val = valid_range
        if value < min_val or value > max_val:
            return False
    else:
        # Unknown predictor - allow any finite number
        pass
    
    return True


def validate_predictor_score_object(score: 'PredictorScore') -> bool:
    """
    Validate a PredictorScore object.
    
    Args:
        score: PredictorScore object to validate
        
    Returns:
        True if the score object is valid, False otherwise
    """
    if not isinstance(score, PredictorScore):
        return False
    
    # Predictor name must be present
    if not score.predictor or not isinstance(score.predictor, str):
        return False
    
    # Validate the score value
    return validate_predictor_score(score.predictor, score.value)


def validate_population_stats(
    af: Optional[float] = None,
    an: Optional[int] = None,
    ac: Optional[int] = None,
    strict: bool = False
) -> bool:
    """
    Validate population allele frequency statistics.
    
    Args:
        af: Allele frequency (0.0 to 1.0)
        an: Allele number (positive integer)
        ac: Allele count (non-negative integer)
        strict: If True, requires at least af or (ac, an) to be present
        
    Returns:
        True if values are valid, False otherwise
        
    Examples:
        >>> validate_population_stats(af=0.001)
        True
        >>> validate_population_stats(af=1.5)  # AF > 1
        False
        >>> validate_population_stats(ac=5, an=1000)
        True
        >>> validate_population_stats(ac=5, an=0)  # AN must be > 0
        False
        >>> validate_population_stats(ac=100, an=50)  # AC > AN
        False
    """
    import math
    
    # Check for required data in strict mode
    if strict:
        if af is None and (ac is None or an is None):
            return False
    
    # Validate allele frequency
    if af is not None:
        if not isinstance(af, (int, float)):
            return False
        if math.isnan(af) or math.isinf(af):
            return False
        if af < 0.0 or af > 1.0:
            return False
    
    # Validate allele number
    if an is not None:
        if not isinstance(an, int):
            # Allow float if it's a whole number
            if isinstance(an, float) and an.is_integer():
                an = int(an)
            else:
                return False
        if an < 0:
            return False
        # AN should not be 0 if we have AC
        if ac is not None and an == 0:
            return False
    
    # Validate allele count
    if ac is not None:
        if not isinstance(ac, int):
            # Allow float if it's a whole number
            if isinstance(ac, float) and ac.is_integer():
                ac = int(ac)
            else:
                return False
        if ac < 0:
            return False
    
    # Cross-validation: AC cannot exceed AN
    if ac is not None and an is not None:
        if ac > an:
            return False
    
    return True


def validate_population_stats_object(stats: 'PopulationStats') -> bool:
    """
    Validate a PopulationStats object.
    
    Args:
        stats: PopulationStats object to validate
        
    Returns:
        True if the stats object is valid, False otherwise
    """
    if not isinstance(stats, PopulationStats):
        return False
    
    # Population name must be present
    if not stats.population or not isinstance(stats.population, str):
        return False
    
    # Validate main frequency data
    if not validate_population_stats(
        af=stats.af,
        an=stats.an,
        ac=stats.ac
    ):
        return False
    
    # Validate popmax_af if present
    if stats.popmax_af is not None:
        import math
        if not isinstance(stats.popmax_af, (int, float)):
            return False
        if math.isnan(stats.popmax_af) or math.isinf(stats.popmax_af):
            return False
        if stats.popmax_af < 0.0 or stats.popmax_af > 1.0:
            return False
    
    # Validate homozygote count if present
    if stats.homozygote_count is not None:
        if not isinstance(stats.homozygote_count, int):
            return False
        if stats.homozygote_count < 0:
            return False
    
    return True


def validate_cached_predictor_data(data: dict) -> tuple[bool, list[str]]:
    """
    Validate a dictionary of cached predictor data.
    
    Args:
        data: Dictionary with predictor names as keys and scores as values
        
    Returns:
        Tuple of (is_valid, list of error messages)
        
    Example:
        >>> data = {'revel': 0.85, 'cadd_phred': 25.0, 'sift': 0.01}
        >>> valid, errors = validate_cached_predictor_data(data)
        >>> valid
        True
    """
    errors = []
    
    if not isinstance(data, dict):
        return False, ['Data must be a dictionary']
    
    for predictor, value in data.items():
        if not isinstance(predictor, str):
            errors.append(f'Invalid predictor key type: {type(predictor)}')
            continue
        
        if not validate_predictor_score(predictor, value):
            errors.append(f'Invalid {predictor} score: {value}')
    
    return len(errors) == 0, errors


def validate_cached_population_data(data: dict) -> tuple[bool, list[str]]:
    """
    Validate a dictionary of cached population data.
    
    Args:
        data: Dictionary with population stats (af, an, ac, etc.)
        
    Returns:
        Tuple of (is_valid, list of error messages)
        
    Example:
        >>> data = {'af': 0.001, 'an': 10000, 'ac': 10}
        >>> valid, errors = validate_cached_population_data(data)
        >>> valid
        True
    """
    errors = []
    
    if not isinstance(data, dict):
        return False, ['Data must be a dictionary']
    
    af = data.get('af')
    an = data.get('an')
    ac = data.get('ac')
    
    if not validate_population_stats(af=af, an=an, ac=ac):
        errors.append(f'Invalid population stats: af={af}, an={an}, ac={ac}')
    
    # Validate popmax_af if present
    popmax_af = data.get('popmax_af')
    if popmax_af is not None:
        if not validate_population_stats(af=popmax_af):
            errors.append(f'Invalid popmax_af: {popmax_af}')
    
    return len(errors) == 0, errors