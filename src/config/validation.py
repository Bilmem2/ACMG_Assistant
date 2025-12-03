"""
Validation Patterns and Input Aliases
====================================

Contains regex patterns for input validation and user convenience aliases.

Author: Can Sevilmiş
License: MIT License
"""

# Validation patterns
VALIDATION_PATTERNS = {
    'gene': r'^[A-Z][A-Z0-9-]*$',
    'gene_symbol': r'^[A-Z][A-Z0-9-]*$',  # Alias for gene
    'chromosome': r'^(chr)?(([1-9]|1[0-9]|2[0-2])|[XYM])$',
    'position': r'^\d+$',
    'allele': r'^[ATCG]+$',
    # HGVS cDNA patterns - supports both simple and full format
    # Full format: NM_000546.6:c.1528C>T
    # Simple format: c.1528C>T or 1528C>T
    'hgvs_cdna_full': r'^(NM_\d+\.\d+):c\.[\d\+\-\*]+[A-Z]>[A-Z]',
    'hgvs_cdna_simple': r'^c\.[\d\+\-\*]+[A-Z]>[A-Z]',
    'hgvs_cdna': r'^(NM_\d+\.\d+:)?c\.[\d\+\-\*]+[A-Z]>[A-Z]',
    # Also support position-only format like 1528C>T
    'position_variant': r'^[\d\+\-\*]+[A-Z]>[A-Z]$'
}

# Variant consequences
VARIANT_CONSEQUENCES = [
    'missense_variant', 'nonsense', 'frameshift_variant', 'splice_donor_variant',
    'splice_acceptor_variant', 'start_lost', 'stop_gained', 'synonymous_variant',
    'intronic_variant', 'regulatory_region_variant'
]

ALL_VARIANT_CONSEQUENCES = VARIANT_CONSEQUENCES

# Input aliases for user convenience
ALIASES = {
    'miss': 'missense',
    'fs': 'frameshift',
    'syn': 'synonymous',
    'het': 'heterozygous',
    'hom': 'homozygous',
    'path': 'pathogenic',
    'ben': 'benign',
    'vus': 'uncertain_significance'
}