"""
Configuration constants for ACMG Variant Classification Assistant
===============================================================

This file provides backward compatibility by importing from the new modular structure.
For new code, import directly from specific modules (e.g., from config.acmg_criteria import EVIDENCE_WEIGHTS)

Author: Can Sevilmiş
License: MIT License
Version: 3.0.0
Last Updated: July 10, 2025
"""

# Import all constants from the new modular structure for backward compatibility
from .version import VERSION_INFO
from .acmg_criteria import (
    EVIDENCE_WEIGHTS, CLASSIFICATION_RULES, STATISTICAL_THRESHOLDS_2015,
    STATISTICAL_THRESHOLDS_2023, STATISTICAL_THRESHOLDS, DOSAGE_SENSITIVITY_THRESHOLDS,
    CONFIDENCE_LEVELS, REPUTABLE_SOURCE_REQUIREMENTS
)
from .gene_rules import GENE_SPECIFIC_THRESHOLDS, LOF_INTOLERANT_GENES, LOF_TOLERANT_GENES
from .predictors import INSILICO_WEIGHTS, INSILICO_THRESHOLDS, VAMPP_SCORE_THRESHOLDS
from .ui_config import (
    COLORAMA_COLORS, SECTION_HEADERS, CLASSIFICATION_COLORS,
    ERROR_MESSAGES, SUCCESS_MESSAGES, WARNING_MESSAGES, INFO_MESSAGES
)
from .validation import VALIDATION_PATTERNS, VARIANT_CONSEQUENCES, ALL_VARIANT_CONSEQUENCES, ALIASES
from .test_scenarios import TEST_SCENARIOS, TEST_MODE_DATA
from .api_config import API_ENDPOINTS, API_SETTINGS

# =============================================================================
# PM1 Hotspot/Domain Evidence Thresholds
# =============================================================================
# Interpretive thresholds for converting remote API confidence scores to
# PM1 evidence levels. These are NOT biological facts - they define how
# to interpret the confidence returned by DomainAPIClient.
#
# DESIGN PRINCIPLE: All actual hotspot/domain data comes from remote APIs
# (CancerHotspots.org, UniProt). These thresholds only control interpretation.

PM1_CONFIDENCE_THRESHOLDS = {
    'PM1': 0.85,            # High confidence -> full PM1 evidence
    'PM1_supporting': 0.60,  # Moderate confidence -> PM1_supporting
}

# Minimum tumor count from CancerHotspots.org to consider a position
# as a high-confidence mutational hotspot.
PM1_HOTSPOT_TUMOR_THRESHOLDS = {
    'high': 10,     # >=10 tumors = high confidence
    'moderate': 3,  # >=3 tumors = moderate confidence
}

# Domain types from UniProt that qualify as "critical functional domains"
# for PM1 purposes. These are UniProt feature type names, not gene names.
PM1_CRITICAL_DOMAIN_TYPES = {
    'Domain',        # Named protein domain
    'Active site',   # Catalytic site
    'Binding site',  # Substrate/cofactor binding
}

PM1_FUNCTIONAL_REGION_TYPES = {
    'Region',        # Named functional region
    'Motif',         # Conserved motif
}

# =============================================================================
# Phenotype Matching Configuration
# =============================================================================
# Thresholds for PP4/BP5 evidence based on phenotype-genotype similarity.
# These control when phenotype matching triggers ACMG evidence codes.

PHENOTYPE_SIMILARITY_THRESHOLDS = {
    'PP4': 0.8,           # High match -> PP4 (phenotype highly specific for gene-disease)
    'PP4_SUPPORTING': 0.5, # Moderate match -> PP4_supporting  
    'BP5': 0.2,           # Low match -> BP5 (phenotype inconsistent with gene-disease)
}

# Low-information HPO terms that should receive reduced weight in similarity calculations.
# These terms are overly generic and can cause false-positive PP4 evidence when used alone.
# 
# WHY DOWN-WEIGHTED: Terms like "Neoplasm" (HP:0002664) or "Phenotypic abnormality" 
# (HP:0000118) match many genes non-specifically. Weighting them at 1.0 would inflate
# similarity scores for patients who present with "cancer" even when the variant is
# in a non-cancer gene. Down-weighting prevents false PP4 due to highly generic terms.
LOW_INFORMATION_HPO = {
    "HP:0002664",   # Neoplasm - very generic cancer term
    "HP:0000118",   # Phenotypic abnormality - root catch-all term
    "HP:0000001",   # All - HPO root term
    "HP:0012823",   # Clinical modifier - not a real phenotype
    "HP:0040279",   # Frequency - metadata, not phenotype
    "HP:0003674",   # Onset - temporal modifier, not phenotype
}

# Weight assigned to low-information HPO terms in weighted Jaccard calculation.
# Normal terms = 1.0, low-information terms = this value (0.3 recommended).
LOW_INFORMATION_HPO_WEIGHT = 0.3

# Minimum number of total terms (patient ∪ gene) required for reliable similarity.
# If fewer terms are available, no phenotype-based evidence (PP4/BP5) is assigned.
# This prevents spurious evidence from very sparse phenotype data.
MIN_TERMS_FOR_PHENOTYPE_EVIDENCE = 3

# Legacy constants that might be missing - add them here if needed
OUTPUT_SETTINGS = {
    'enabled': True,
    'format': 'text',
    'cache_filename': 'api_cache.json',
    'error_log': 'api_errors.log',
    'report_filename': 'variant_classification_report.txt',
    'log_filename': 'classification.log',
    'max_cache_age_hours': 24
}

CONSTRAINT_THRESHOLDS = {
    'lof_intolerant': 0.35,
    'lof_tolerant': 0.6,
    'pLI_intolerant': 0.9,
    'LOEUF_intolerant': 0.35,
    'oe_lof_upper_tolerant': 0.6
}

def get_colored_message(message: str, color: str = 'WHITE') -> str:
    """Get a colored message for terminal output."""
    return f"{COLORAMA_COLORS.get(color, COLORAMA_COLORS['WHITE'])}{message}{COLORAMA_COLORS['RESET']}"