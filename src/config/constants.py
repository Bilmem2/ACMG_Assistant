"""
Configuration constants for ACMG Variant Classification Assistant
===============================================================

This file contains all the constants, thresholds, and configuration
values used throughout the application.

Author: Can Sevilmiş
License: MIT License
Version: 3.0.0
Last Updated: July 10, 2025
"""

# Version and metadata information
VERSION_INFO = {
    'version': '3.4.0',
    'author': 'Can Sevilmiş',
    'license': 'MIT License',
    'last_updated': 'October 6, 2025',
    'guidelines': 'ACMG/AMP 2015 & 2023',
    'description': 'ACMG Variant Classification Assistant - Dual Guidelines Support',
    'major_features': [
        'Dual-mode ACMG 2015/2023 guidelines support',
        'Complete 28 ACMG/AMP criteria implementation',
        'Interactive evidence evaluation',
        'Enhanced computational metascore',
        'User-driven criteria assignment',
        'Dynamic evidence weighting',
        'Gene-specific thresholds for BA1/BS1',
        'LOF intolerant/tolerant gene classification',
        'Ensembl API integration for chromosome lookup',
        'Improved input validation and error handling',
        '✨ Confidence & provenance tracking',
        '✨ Automated Fisher\'s exact test (PS4)',
        '✨ Automated LOD scoring (PP1/BS4)',
        '✨ Strict PS2 2023 upgrade rules',
        '✨ Stricter PP5/BP6 source validation',
        '✨ HGVS format support with auto-extraction'
    ]
}

# ACMG Evidence Criteria Weights
EVIDENCE_WEIGHTS = {
    'PVS1': 8, 'PS1': 4, 'PS2': 4, 'PS3': 4, 'PS4': 4,
    'PM1': 2, 'PM2': 2, 'PM3': 2, 'PM4': 2, 'PM5': 2, 'PM6': 2,
    'PP1': 1, 'PP2': 1, 'PP3': 1, 'PP4': 1, 'PP5': 1,
    'BA1': -8, 'BS1': -4, 'BS2': -4, 'BS3': -4, 'BS4': -4,
    'BP1': -1, 'BP2': -1, 'BP3': -1, 'BP4': -1, 'BP5': -1, 'BP6': -1, 'BP7': -1
}

# Classification rules for ACMG/AMP guidelines
CLASSIFICATION_RULES = {
    'Pathogenic': {
        'rules': [
            {'very_strong': 1, 'strong': 0, 'moderate': 0, 'supporting': 0},
            {'very_strong': 0, 'strong': 2, 'moderate': 0, 'supporting': 0},
            {'very_strong': 0, 'strong': 1, 'moderate': 3, 'supporting': 0},
            {'very_strong': 0, 'strong': 1, 'moderate': 2, 'supporting': 2},
            {'very_strong': 0, 'strong': 1, 'moderate': 1, 'supporting': 4}
        ]
    },
    'Likely Pathogenic': {
        'rules': [
            {'very_strong': 0, 'strong': 1, 'moderate': 1, 'supporting': 0},
            {'very_strong': 0, 'strong': 1, 'moderate': 0, 'supporting': 2},
            {'very_strong': 0, 'strong': 0, 'moderate': 3, 'supporting': 0},
            {'very_strong': 0, 'strong': 0, 'moderate': 2, 'supporting': 2},
            {'very_strong': 0, 'strong': 0, 'moderate': 1, 'supporting': 4}
        ]
    },
    'Benign': {
        'rules': [
            {'stand_alone': 1},
            {'strong': 2},
            {'strong': 1, 'supporting': 2}  # ACMG 2015: 1 BS + ≥2 BP required for Benign
        ]
    },
    'Likely Benign': {
        'rules': [
            {'strong': 1, 'supporting': 0},
            {'strong': 0, 'supporting': 2}
        ]
    }
}

# Gene-specific thresholds
GENE_SPECIFIC_THRESHOLDS = {
    'BRCA1': {'BA1': 0.05, 'BS1': 0.01, 'PM2': 0.0001},
    'BRCA2': {'BA1': 0.05, 'BS1': 0.01, 'PM2': 0.0001},
    'TTN': {'BA1': 0.1, 'BS1': 0.05, 'PM2': 0.001},
    'MUC16': {'BA1': 0.1, 'BS1': 0.05, 'PM2': 0.001},
    'OBSCN': {'BA1': 0.1, 'BS1': 0.05, 'PM2': 0.001},
    'default': {'BA1': 0.05, 'BS1': 0.01, 'PM2': 0.0001}
}

# LOF mechanism gene classification for PVS1 evaluation
LOF_INTOLERANT_GENES = {
    'BRCA1', 'BRCA2', 'TP53', 'RB1', 'APC', 'VHL', 'NF1', 'NF2',
    'CDKN2A', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'ATM', 'CHEK2',
    'PALB2', 'MYH7', 'MYBPC3', 'SCN1A', 'SCN2A', 'MECP2', 'PAH', 'CFTR'
}

LOF_TOLERANT_GENES = {
    'TTN', 'MUC16', 'OBSCN', 'PCLO', 'RYR1', 'SYNE1', 'SYNE2', 'USH2A', 'FLG'
}

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

# Statistical thresholds for various analyses
# DEFAULT: ACMG 2015 Guidelines
STATISTICAL_THRESHOLDS_2015 = {
    'fisher_exact_p_value': 0.05,
    'case_control_odds_ratio': 2.0,
    'case_control_min_cases': 5,  # Minimum affected probands
    'case_control_min_controls': 1000,  # Minimum control individuals
    'segregation_lod_score': 3.0,  # LOD ≥ 3.0 for strong support
    'segregation_lod_supporting': 1.5,  # LOD ≥ 1.5 for PP1
    'segregation_families_min': 3,  # Minimum families for PP1/BS4
    'splice_ai_threshold': 0.5,
    'conservation_threshold': 2.0,
    'pm5_min_pathogenic_variants': 1  # ≥1 pathogenic variant at same codon
}

# ACMG 2023 Updates (ClinGen SVI Working Group)
STATISTICAL_THRESHOLDS_2023 = {
    'fisher_exact_p_value': 0.05,
    'case_control_odds_ratio': 5.0,  # Increased from 2.0 to 5.0 for PS4
    'case_control_min_cases': 10,  # Increased minimum affected probands
    'case_control_min_controls': 2000,  # Increased minimum controls
    'segregation_lod_score': 3.0,  # LOD ≥ 3.0 for PP1_Moderate
    'segregation_lod_supporting': 1.5,  # LOD 1.5-2.99 for PP1_Supporting
    'segregation_lod_moderate': 3.0,  # LOD 3.0-4.99 for PP1_Moderate
    'segregation_lod_strong': 5.0,  # LOD ≥ 5.0 for PP1_Strong
    'segregation_lod_bs4': -2.0,  # LOD ≤ -2.0 for BS4
    'segregation_families_min': 3,  # Minimum informative meioses
    'splice_ai_threshold': 0.5,
    'conservation_threshold': 2.0,
    'pm5_min_pathogenic_variants': 2  # ≥2 pathogenic variants at same codon
}

# Default to ACMG 2015
STATISTICAL_THRESHOLDS = STATISTICAL_THRESHOLDS_2015

# Confidence levels for evidence criteria
CONFIDENCE_LEVELS = {
    'high': 'High confidence - deterministic or well-validated',
    'medium': 'Medium confidence - data-driven with assumptions',
    'low': 'Low confidence - user input or incomplete data',
    'very_low': 'Very low confidence - placeholder or uncertain'
}

# Reputable source requirements for PP5/BP6
REPUTABLE_SOURCE_REQUIREMENTS = {
    'expert_panels': ['ClinGen', 'ENIGMA', 'InSiGHT', 'ACMG', 'ClinVar Expert Panel'],
    'certified_labs': True,  # Require certified clinical lab
    'min_stars': 2,  # Minimum ClinVar review stars
    'max_age_years': 5  # Maximum age of classification in years
}

# Colorama color codes for terminal output
COLORAMA_COLORS = {
    'RED': '\033[91m',
    'GREEN': '\033[92m',
    'YELLOW': '\033[93m',
    'BLUE': '\033[94m',
    'MAGENTA': '\033[95m',
    'CYAN': '\033[96m',
    'WHITE': '\033[97m',
    'ENDC': '\033[0m',  # End color
    'RESET': '\033[0m',  # Reset color (alias for ENDC)
    'BOLD': '\033[1m',
    'UNDERLINE': '\033[4m'
}

# Section headers for organized output
SECTION_HEADERS = {
    'variant_info': f"{COLORAMA_COLORS['CYAN']}🧬 VARIANT INFORMATION{COLORAMA_COLORS['ENDC']}",
    'pathogenic_evidence': f"{COLORAMA_COLORS['RED']}🔴 PATHOGENIC EVIDENCE{COLORAMA_COLORS['ENDC']}",
    'benign_evidence': f"{COLORAMA_COLORS['GREEN']}🟢 BENIGN EVIDENCE{COLORAMA_COLORS['ENDC']}",
    'final_classification': f"{COLORAMA_COLORS['BOLD']}📊 FINAL CLASSIFICATION{COLORAMA_COLORS['ENDC']}",
    'recommendations': f"{COLORAMA_COLORS['YELLOW']}💡 RECOMMENDATIONS{COLORAMA_COLORS['ENDC']}",
    'report_generation': f"{COLORAMA_COLORS['BLUE']}📋 REPORT GENERATION{COLORAMA_COLORS['ENDC']}"
}

# Classification-specific colors
CLASSIFICATION_COLORS = {
    'Pathogenic': f"{COLORAMA_COLORS['RED']}{COLORAMA_COLORS['BOLD']}",
    'Likely Pathogenic': f"{COLORAMA_COLORS['RED']}",
    'Uncertain Significance': f"{COLORAMA_COLORS['YELLOW']}",
    'Likely Benign': f"{COLORAMA_COLORS['GREEN']}",
    'Benign': f"{COLORAMA_COLORS['GREEN']}{COLORAMA_COLORS['BOLD']}"
}

# Error and message constants
ERROR_MESSAGES = {
    'invalid_gene': 'Invalid gene symbol provided',
    'invalid_chromosome': 'Invalid chromosome format',
    'invalid_position': 'Invalid genomic position',
    'invalid_allele': 'Invalid allele format',
    'missing_data': 'Required data is missing',
    'api_error': 'API request failed',
    'file_not_found': 'File not found',
    'invalid_format': 'Invalid input format'
}

SUCCESS_MESSAGES = {
    'data_loaded': 'Data loaded successfully',
    'analysis_complete': 'Analysis completed successfully',
    'report_generated': 'Report generated successfully'
}

WARNING_MESSAGES = {
    'missing_optional_data': 'Optional data missing, proceeding with available data',
    'low_confidence': 'Low confidence result due to limited data'
}

INFO_MESSAGES = {
    'interactive_mode': 'Interactive mode enabled - user input required',
    'test_mode': 'Test mode enabled - using default responses'
}

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

# Test mode data
TEST_MODE_DATA = {
    'default_responses': {
        'PS1': False,
        'PM5': False,
        'PP1': False,
        'PP4': False,
        'BS4': False
    },
    'basic_info': {
        'gene': 'BRCA1',
        'chromosome': '17',
        'position': '43093464',
        'ref_allele': 'C',
        'alt_allele': 'T',
        'amino_acid_change': 'p.Gln356Ter',
        'variant_type': 'nonsense',
        'consequence': 'stop_gained',
        'transcript': 'NM_007294.4',
        'hgvs_cdna': 'c.1066C>T',
        'hgvs_protein': 'p.Gln356Ter'
    },
    'population_data': {
        'gnomad_af': 0.0,
        'gnomad_homozygous': 0,
        'exac_af': 0.0,
        'thousand_genomes_af': 0.0,
        'esp6500_af': None
    },
    'insilico_data': {
        'revel_score': None,  # REVEL doesn't score nonsense variants
        'cadd_phred': 35.0,   # High CADD score for nonsense
        'sift_score': None,   # SIFT doesn't score nonsense variants
        'polyphen2_score': None,  # PolyPhen doesn't score nonsense variants
        'vest4_score': None,
        'alphamissense_score': None  # AlphaMissense only for missense
    },
    'functional_data': {
        'case_control': 'not_available',
        'segregation': 'not_available',
        'functional_studies': 'not_available'
    },
    'genetic_data': {
        'de_novo_status': 'not_reported',
        'parental_confirmation': 'not_confirmed',
        'inheritance_pattern': 'autosomal_dominant',
        'zygosity': 'heterozygous'
    }
}

# API endpoints
API_ENDPOINTS = {
    'clinvar': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/',
    'gnomad': 'https://gnomad.broadinstitute.org/api/',
    'varsome': 'https://api.varsome.com/lookup/',
    'ensembl_rest': 'https://rest.ensembl.org'
}

# Output settings
OUTPUT_SETTINGS = {
    'report_format': 'txt',
    'include_details': True,
    'color_output': True,
    'cache_filename': 'api_cache.json',
    'max_cache_age_hours': 24,
    'output_directory': 'reports',
    'log_level': 'INFO',
    'report_filename': 'variant_classification_report.txt',
    'log_filename': 'acmg_analysis.log'
}

def get_colored_message(message: str, color_key: str) -> str:
    """Get a colored message for terminal output."""
    color = COLORAMA_COLORS.get(color_key, '')
    return f"{color}{message}{COLORAMA_COLORS['ENDC']}"
