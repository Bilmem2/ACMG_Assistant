"""
ACMG Configuration Module
========================

Centralized configuration imports for ACMG Variant Classification Assistant.

Author: Can Sevilmiş
License: MIT License
"""

# Version information
from .version import VERSION_INFO

# ACMG criteria and rules
from .acmg_criteria import (
    EVIDENCE_WEIGHTS,
    CLASSIFICATION_RULES,
    STATISTICAL_THRESHOLDS_2015,
    STATISTICAL_THRESHOLDS_2023,
    STATISTICAL_THRESHOLDS,
    DOSAGE_SENSITIVITY_THRESHOLDS,
    CONFIDENCE_LEVELS,
    REPUTABLE_SOURCE_REQUIREMENTS
)

# Gene-specific configurations
from .gene_rules import (
    GENE_SPECIFIC_THRESHOLDS,
    LOF_INTOLERANT_GENES,
    LOF_TOLERANT_GENES
)

# In silico predictor configurations
from .predictors import (
    INSILICO_WEIGHTS,
    INSILICO_THRESHOLDS,
    VAMPP_SCORE_THRESHOLDS
)

# UI and messaging
from .ui_config import (
    COLORAMA_COLORS,
    SECTION_HEADERS,
    CLASSIFICATION_COLORS,
    ERROR_MESSAGES,
    SUCCESS_MESSAGES,
    WARNING_MESSAGES,
    INFO_MESSAGES
)

# Validation patterns and aliases
from .validation import (
    VALIDATION_PATTERNS,
    VARIANT_CONSEQUENCES,
    ALL_VARIANT_CONSEQUENCES,
    ALIASES
)

# Test scenarios
from .test_scenarios import (
    TEST_SCENARIOS,
    TEST_MODE_DATA
)

# API configuration
from .api_config import (
    API_ENDPOINTS,
    API_SETTINGS
)

# Backward compatibility - export all constants as they were
__all__ = [
    # Version
    'VERSION_INFO',
    
    # ACMG Criteria
    'EVIDENCE_WEIGHTS',
    'CLASSIFICATION_RULES', 
    'STATISTICAL_THRESHOLDS_2015',
    'STATISTICAL_THRESHOLDS_2023',
    'STATISTICAL_THRESHOLDS',
    'DOSAGE_SENSITIVITY_THRESHOLDS',
    'CONFIDENCE_LEVELS',
    'REPUTABLE_SOURCE_REQUIREMENTS',
    
    # Gene Rules
    'GENE_SPECIFIC_THRESHOLDS',
    'LOF_INTOLERANT_GENES',
    'LOF_TOLERANT_GENES',
    
    # Predictors
    'INSILICO_WEIGHTS',
    'INSILICO_THRESHOLDS', 
    'VAMPP_SCORE_THRESHOLDS',
    
    # UI Config
    'COLORAMA_COLORS',
    'SECTION_HEADERS',
    'CLASSIFICATION_COLORS',
    'ERROR_MESSAGES',
    'SUCCESS_MESSAGES',
    'WARNING_MESSAGES',
    'INFO_MESSAGES',
    
    # Validation
    'VALIDATION_PATTERNS',
    'VARIANT_CONSEQUENCES',
    'ALL_VARIANT_CONSEQUENCES',
    'ALIASES',
    
    # Test Data
    'TEST_SCENARIOS',
    'TEST_MODE_DATA',
    
    # API Config
    'API_ENDPOINTS',
    'API_SETTINGS'
]