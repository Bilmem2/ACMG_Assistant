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
    'version': '3.5.0',
    'author': 'Can Sevilmiş',
    'license': 'MIT License',
    'last_updated': 'October 6, 2025',
    'guidelines': 'ACMG/AMP 2015 & 2023',
    'description': 'ACMG Variant Classification Assistant - Comprehensive API Integration',
    'major_features': [
        'Dual-mode ACMG 2015/2023 guidelines support',
        'Complete 28 ACMG/AMP criteria implementation',
        'Interactive evidence evaluation',
        'Enhanced computational metascore',
        'User-driven criteria assignment',
        'Dynamic evidence weighting',
        '🔬 gnomAD v4 constraint API (LOF intolerance)',
        '🔬 ClinGen eRepo API (disease-specific LOF mechanism)',
        '🔬 Multi-source LOF validation (population + clinical)',
        '🔬 Transparent data provenance tracking',
        '✨ Gene-specific BA1/BS1 thresholds',
        '✨ Automated Fisher\'s exact test (PS4)',
        '✨ Automated LOD scoring (PP1/BS4)',
        '✨ Strict PS2 2023 upgrade rules',
        '✨ Enhanced PP5/BP6 source validation',
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
            {'very_strong': 0, 'strong': 0, 'moderate': 1, 'supporting': 4},
            {'very_strong': 0, 'strong': 0, 'moderate': 1, 'supporting': 2}  # PM + 2 PP (e.g., PM2+PP3+PP5)
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
    'PALB2', 'MYH7', 'MYBPC3', 'SCN1A', 'SCN2A', 'MECP2', 'PAH', 'CFTR', 'RUNX1'
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

# ClinGen Dosage Sensitivity thresholds for PVS1 modulation
DOSAGE_SENSITIVITY_THRESHOLDS = {
    'haploinsufficiency': {
        'sufficient': 3,      # Sufficient evidence → Keep PVS1 Very Strong
        'some': 2,           # Some evidence → Downgrade to PS1 Strong
        'little': 1,         # Little evidence → Downgrade to PM2 Moderate
        'no_evidence': 0,    # No evidence → Consider not applying
        'unlikely': 40       # Unlikely HI → Consider not applying
    },
    'triplosensitivity': {
        'sufficient': 3,
        'some': 2,
        'little': 1,
        'no_evidence': 0,
        'unlikely': 40
    }
}

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
TEST_SCENARIOS = {
    '1_pvs1_nonsense': {
        'name': 'PVS1 - Nonsense Variant',
        'description': 'Classic loss-of-function nonsense mutation in critical gene (BRCA1)',
        'expected_classification': 'Pathogenic',
        'expected_criteria': ['PVS1', 'PM2'],
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
            'revel_score': None,
            'cadd_phred': 35.0,
            'sift_score': None,
            'polyphen2_score': None,
            'vest4_score': None,
            'alphamissense_score': None
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
    },
    '2_ps2_denovo': {
        'name': 'PS2 - De Novo Variant',
        'description': 'De novo missense in autism-associated gene (SCN2A) with confirmed paternity/maternity',
        'expected_classification': 'Likely Pathogenic',
        'expected_criteria': ['PS2', 'PM1', 'PM2', 'PP2', 'PP3'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': True,  # Patient phenotype matches
            'BS4': False
        },
        'basic_info': {
            'gene': 'SCN2A',
            'chromosome': '2',
            'position': '165310411',
            'ref_allele': 'C',
            'alt_allele': 'T',
            'amino_acid_change': 'p.Arg853Gln',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_001040142.2',
            'hgvs_cdna': 'c.2558G>A',
            'hgvs_protein': 'p.Arg853Gln'
        },
        'population_data': {
            'gnomad_af': 0.0,
            'gnomad_homozygous': 0,
            'exac_af': 0.0,
            'thousand_genomes_af': 0.0,
            'esp6500_af': None
        },
        'insilico_data': {
            'revel_score': 0.82,
            'cadd_phred': 28.5,
            'sift_score': 0.01,
            'polyphen2_score': 0.98,
            'metasvm_score': 0.91,
            'alphamissense_score': 0.85
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'not_available',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'confirmed',
            'parental_confirmation': 'both_confirmed',
            'maternity_confirmed': True,
            'paternity_confirmed': True,
            'inheritance_pattern': 'de_novo',
            'zygosity': 'heterozygous'
        },
        'patient_phenotypes': 'epilepsy, developmental delay, autism spectrum disorder'
    },
    '3_pm1_hotspot': {
        'name': 'PM1 - Critical Domain Hotspot',
        'description': 'Missense in well-established functional domain (TP53 DNA-binding)',
        'expected_classification': 'Likely Pathogenic',
        'expected_criteria': ['PM1', 'PM2', 'PM5', 'PP2', 'PP3'],
        'default_responses': {
            'PS1': False,
            'PM5': True,  # Same AA different nt change is pathogenic
            'PP1': False,
            'PP4': True,
            'BS4': False
        },
        'basic_info': {
            'gene': 'TP53',
            'chromosome': '17',
            'position': '7674220',
            'ref_allele': 'C',
            'alt_allele': 'A',
            'amino_acid_change': 'p.Arg248Trp',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_000546.6',
            'hgvs_cdna': 'c.742C>T',
            'hgvs_protein': 'p.Arg248Trp'
        },
        'population_data': {
            'gnomad_af': 0.000002,
            'gnomad_homozygous': 0,
            'exac_af': 0.000001,
            'thousand_genomes_af': 0.0,
            'esp6500_af': None
        },
        'insilico_data': {
            'revel_score': 0.95,
            'cadd_phred': 32.0,
            'sift_score': 0.0,
            'polyphen2_score': 1.0,
            'metasvm_score': 0.98,
            'alphamissense_score': 0.92
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
        },
        'patient_phenotypes': 'Li-Fraumeni syndrome, multiple primary cancers'
    },
    '4_ba1_common': {
        'name': 'BA1 - Common Benign Variant',
        'description': 'High frequency missense variant (>5% in population)',
        'expected_classification': 'Benign',
        'expected_criteria': ['BA1'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': False,
            'BS4': False
        },
        'basic_info': {
            'gene': 'MTHFR',
            'chromosome': '1',
            'position': '11796321',
            'ref_allele': 'G',
            'alt_allele': 'A',
            'amino_acid_change': 'p.Ala222Val',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_005957.5',
            'hgvs_cdna': 'c.665C>T',
            'hgvs_protein': 'p.Ala222Val'
        },
        'population_data': {
            'gnomad_af': 0.32,
            'gnomad_homozygous': 15000,
            'exac_af': 0.31,
            'thousand_genomes_af': 0.30,
            'esp6500_af': 0.29
        },
        'insilico_data': {
            'revel_score': 0.15,
            'cadd_phred': 15.2,
            'sift_score': 0.45,
            'polyphen2_score': 0.10,
            'metasvm_score': 0.08,
            'alphamissense_score': 0.12
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'not_available',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'autosomal_recessive',
            'zygosity': 'heterozygous'
        }
    },
    '5_bp7_synonymous': {
        'name': 'BP7 - Synonymous Variant',
        'description': 'Silent mutation with no predicted splice impact',
        'expected_classification': 'Likely Benign',
        'expected_criteria': ['BP7', 'BP4'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': False,
            'BS4': False
        },
        'basic_info': {
            'gene': 'CFTR',
            'chromosome': '7',
            'position': '117559590',
            'ref_allele': 'C',
            'alt_allele': 'T',
            'amino_acid_change': 'p.Thr854=',
            'variant_type': 'synonymous',
            'consequence': 'synonymous_variant',
            'transcript': 'NM_000492.4',
            'hgvs_cdna': 'c.2562G>A',
            'hgvs_protein': 'p.Thr854='
        },
        'population_data': {
            'gnomad_af': 0.0012,
            'gnomad_homozygous': 2,
            'exac_af': 0.0011,
            'thousand_genomes_af': 0.0010,
            'esp6500_af': None
        },
        'insilico_data': {
            'revel_score': None,
            'cadd_phred': 8.5,
            'sift_score': None,
            'polyphen2_score': None,
            'spliceai_ag_score': 0.01,
            'spliceai_al_score': 0.02,
            'spliceai_dg_score': 0.01,
            'spliceai_dl_score': 0.01
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'not_available',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'autosomal_recessive',
            'zygosity': 'heterozygous'
        }
    },
    '6_vus_novel': {
        'name': 'VUS - Novel Missense',
        'description': 'Novel missense with conflicting predictions, no functional data',
        'expected_classification': 'Uncertain Significance',
        'expected_criteria': ['PM2', 'PP3'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': False,
            'BS4': False
        },
        'basic_info': {
            'gene': 'BRCA2',
            'chromosome': '13',
            'position': '32339832',
            'ref_allele': 'G',
            'alt_allele': 'A',
            'amino_acid_change': 'p.Ala1570Thr',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_000059.4',
            'hgvs_cdna': 'c.4708G>A',
            'hgvs_protein': 'p.Ala1570Thr'
        },
        'population_data': {
            'gnomad_af': 0.0,
            'gnomad_homozygous': 0,
            'exac_af': 0.0,
            'thousand_genomes_af': 0.0,
            'esp6500_af': None
        },
        'insilico_data': {
            'revel_score': 0.52,
            'cadd_phred': 22.5,
            'sift_score': 0.08,
            'polyphen2_score': 0.62,
            'metasvm_score': 0.48,
            'alphamissense_score': 0.55
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
    },
    '7_ps3_functional': {
        'name': 'PS3 - Functional Study Proven',
        'description': 'Missense with well-established functional evidence of pathogenicity',
        'expected_classification': 'Likely Pathogenic',
        'expected_criteria': ['PS3', 'PM1', 'PM2', 'PP2', 'PP3'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': True,
            'PP4': True,
            'BS4': False
        },
        'basic_info': {
            'gene': 'KCNQ1',
            'chromosome': '11',
            'position': '2478327',
            'ref_allele': 'C',
            'alt_allele': 'T',
            'amino_acid_change': 'p.Arg190Gln',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_000218.3',
            'hgvs_cdna': 'c.569G>A',
            'hgvs_protein': 'p.Arg190Gln'
        },
        'population_data': {
            'gnomad_af': 0.000001,
            'gnomad_homozygous': 0,
            'exac_af': 0.0,
            'thousand_genomes_af': 0.0,
            'esp6500_af': None
        },
        'insilico_data': {
            'revel_score': 0.88,
            'cadd_phred': 29.5,
            'sift_score': 0.0,
            'polyphen2_score': 0.99,
            'metasvm_score': 0.92,
            'alphamissense_score': 0.89
        },
        'functional_data': {
            'case_control': 'pathogenic',
            'segregation': 'cosegregates',
            'functional_studies': 'damaging',
            'functional_study_type': 'electrophysiology',
            'functional_details': 'Patch-clamp studies show complete loss of K+ channel function'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'autosomal_dominant',
            'zygosity': 'heterozygous',
            'segregation_families': 5,
            'segregation_affected': 8,
            'segregation_unaffected': 0
        },
        'patient_phenotypes': 'Long QT syndrome, syncope, cardiac arrhythmia'
    },
    '8_bs3_functional_benign': {
        'name': 'BS3 - Functional Study Benign',
        'description': 'Missense with strong functional evidence of benign impact',
        'expected_classification': 'Likely Benign',
        'expected_criteria': ['BS3', 'BP4'],
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
            'position': '43106487',
            'ref_allele': 'C',
            'alt_allele': 'T',
            'amino_acid_change': 'p.Pro871Leu',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_007294.4',
            'hgvs_cdna': 'c.2612C>T',
            'hgvs_protein': 'p.Pro871Leu'
        },
        'population_data': {
            'gnomad_af': 0.0003,
            'gnomad_homozygous': 0,
            'exac_af': 0.0002,
            'thousand_genomes_af': 0.0001,
            'esp6500_af': None
        },
        'insilico_data': {
            'revel_score': 0.25,
            'cadd_phred': 18.5,
            'sift_score': 0.35,
            'polyphen2_score': 0.15,
            'metasvm_score': 0.18,
            'alphamissense_score': 0.22
        },
        'functional_data': {
            'case_control': 'benign',
            'segregation': 'not_available',
            'functional_studies': 'benign',
            'functional_study_type': 'homologous_recombination',
            'functional_details': 'HR assay shows normal BRCA1 function (>90% of wild-type)'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'autosomal_dominant',
            'zygosity': 'heterozygous'
        }
    },
    '9_pvs1_frameshift': {
        'name': 'PVS1 - Frameshift Variant',
        'description': 'Frameshift deletion causing loss of function in critical gene',
        'expected_classification': 'Pathogenic',
        'expected_criteria': ['PVS1', 'PM2', 'PP3'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': True,
            'BS4': False
        },
        'basic_info': {
            'gene': 'BRCA2',
            'chromosome': '13',
            'position': '32914137',
            'ref_allele': 'GT',
            'alt_allele': 'G',
            'amino_acid_change': 'p.Val1283fs',
            'variant_type': 'frameshift',
            'consequence': 'frameshift_variant',
            'transcript': 'NM_000059.4',
            'hgvs_cdna': 'c.3847delT',
            'hgvs_protein': 'p.Val1283fs'
        },
        'population_data': {
            'gnomad_af': 0.0,
            'gnomad_homozygous': 0,
            'exac_af': 0.0,
            'thousand_genomes_af': 0.0
        },
        'insilico_data': {
            'cadd_phred': 33.0,
            'revel_score': None
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
        },
        'patient_phenotypes': 'breast cancer, family history'
    },
    '10_pvs1_splice': {
        'name': 'PVS1 - Canonical Splice Site',
        'description': 'Variant at canonical splice site (±1 or 2 positions)',
        'expected_classification': 'Pathogenic',
        'expected_criteria': ['PVS1', 'PM2'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': False,
            'BS4': False
        },
        'basic_info': {
            'gene': 'DMD',
            'chromosome': 'X',
            'position': '31160983',
            'ref_allele': 'G',
            'alt_allele': 'A',
            'amino_acid_change': 'p.?',
            'variant_type': 'splice_donor',
            'consequence': 'splice_donor_variant',
            'transcript': 'NM_004006.3',
            'hgvs_cdna': 'c.9164+1G>A',
            'hgvs_protein': 'p.?'
        },
        'population_data': {
            'gnomad_af': 0.0,
            'gnomad_homozygous': 0,
            'exac_af': 0.0
        },
        'insilico_data': {
            'spliceai_dg_score': 0.95,
            'spliceai_dl_score': 0.02,
            'cadd_phred': 34.0
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'not_available',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'x_linked',
            'zygosity': 'hemizygous'
        },
        'patient_phenotypes': 'Duchenne muscular dystrophy'
    },
    '11_ps1_same_aa': {
        'name': 'PS1 - Same Amino Acid Change',
        'description': 'Same AA change as known pathogenic (different nucleotide)',
        'expected_classification': 'Likely Pathogenic',
        'expected_criteria': ['PS1', 'PM2', 'PP2', 'PP3'],
        'default_responses': {
            'PS1': True,
            'PM5': False,
            'PP1': False,
            'PP4': True,
            'BS4': False
        },
        'basic_info': {
            'gene': 'RET',
            'chromosome': '10',
            'position': '43113112',
            'ref_allele': 'T',
            'alt_allele': 'A',
            'amino_acid_change': 'p.Cys634Trp',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_020975.6',
            'hgvs_cdna': 'c.1902T>A',
            'hgvs_protein': 'p.Cys634Trp'
        },
        'population_data': {
            'gnomad_af': 0.0,
            'gnomad_homozygous': 0,
            'exac_af': 0.0
        },
        'insilico_data': {
            'revel_score': 0.91,
            'cadd_phred': 31.0,
            'sift_score': 0.0,
            'polyphen2_score': 1.0,
            'alphamissense_score': 0.94
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
        },
        'patient_phenotypes': 'MEN2A, medullary thyroid carcinoma'
    },
    '12_ps4_case_control': {
        'name': 'PS4 - Case-Control Evidence',
        'description': 'Significantly enriched in cases vs controls (ACMG 2023 criteria)',
        'expected_classification': 'Likely Pathogenic',
        'expected_criteria': ['PS4', 'PM2', 'PP2', 'PP3'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': True,
            'BS4': False
        },
        'basic_info': {
            'gene': 'APOE',
            'chromosome': '19',
            'position': '44908684',
            'ref_allele': 'C',
            'alt_allele': 'T',
            'amino_acid_change': 'p.Arg158Cys',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_000041.4',
            'hgvs_cdna': 'c.472C>T',
            'hgvs_protein': 'p.Arg158Cys'
        },
        'population_data': {
            'gnomad_af': 0.015,
            'gnomad_homozygous': 5,
            'exac_af': 0.014
        },
        'insilico_data': {
            'revel_score': 0.75,
            'cadd_phred': 26.5,
            'sift_score': 0.02,
            'polyphen2_score': 0.95
        },
        'functional_data': {
            'case_control': 'pathogenic',
            'case_control_or': 3.8,
            'case_control_pvalue': 0.00001,
            'case_count': 250,
            'control_count': 5000,
            'segregation': 'not_available',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'complex',
            'zygosity': 'heterozygous'
        },
        'patient_phenotypes': 'Alzheimer disease, early onset'
    },
    '13_pm3_recessive': {
        'name': 'PM3 - Recessive Trans Configuration',
        'description': 'Detected in trans with pathogenic variant (recessive disorder)',
        'expected_classification': 'Likely Pathogenic',
        'expected_criteria': ['PM3', 'PM2', 'PP3'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': True,
            'BS4': False
        },
        'basic_info': {
            'gene': 'CFTR',
            'chromosome': '7',
            'position': '117534438',
            'ref_allele': 'G',
            'alt_allele': 'A',
            'amino_acid_change': 'p.Gly551Asp',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_000492.4',
            'hgvs_cdna': 'c.1652G>A',
            'hgvs_protein': 'p.Gly551Asp'
        },
        'population_data': {
            'gnomad_af': 0.00003,
            'gnomad_homozygous': 0,
            'exac_af': 0.00002
        },
        'insilico_data': {
            'revel_score': 0.87,
            'cadd_phred': 29.0,
            'sift_score': 0.0,
            'polyphen2_score': 0.99
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'not_available',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'autosomal_recessive',
            'zygosity': 'compound_heterozygous',
            'trans_pathogenic': 'c.1521_1523delCTT (p.Phe508del)',
            'second_variant_classification': 'Pathogenic'
        },
        'patient_phenotypes': 'cystic fibrosis, pancreatic insufficiency'
    },
    '14_pm4_inframe': {
        'name': 'PM4 - In-frame Indel in Non-repeat',
        'description': 'In-frame deletion in critical region without repeat',
        'expected_classification': 'Likely Pathogenic',
        'expected_criteria': ['PM4', 'PM1', 'PM2'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': True,
            'BS4': False
        },
        'basic_info': {
            'gene': 'PTEN',
            'chromosome': '10',
            'position': '87933147',
            'ref_allele': 'AGGCGC',
            'alt_allele': 'A',
            'amino_acid_change': 'p.Arg130_Ala131del',
            'variant_type': 'inframe_deletion',
            'consequence': 'inframe_deletion',
            'transcript': 'NM_000314.8',
            'hgvs_cdna': 'c.388_393delGGCGCA',
            'hgvs_protein': 'p.Arg130_Ala131del'
        },
        'population_data': {
            'gnomad_af': 0.0,
            'gnomad_homozygous': 0,
            'exac_af': 0.0
        },
        'insilico_data': {
            'cadd_phred': 27.5,
            'provean_score': -8.5
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
        },
        'patient_phenotypes': 'Cowden syndrome, macrocephaly'
    },
    '15_pm6_denovo_no_paternity': {
        'name': 'PM6 - De Novo (No Paternity Test)',
        'description': 'Assumed de novo without paternity/maternity confirmation',
        'expected_classification': 'Likely Pathogenic',
        'expected_criteria': ['PM6', 'PM1', 'PM2', 'PP3'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': True,
            'BS4': False
        },
        'basic_info': {
            'gene': 'MECP2',
            'chromosome': 'X',
            'position': '154030892',
            'ref_allele': 'C',
            'alt_allele': 'T',
            'amino_acid_change': 'p.Arg306Cys',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_004992.4',
            'hgvs_cdna': 'c.916C>T',
            'hgvs_protein': 'p.Arg306Cys'
        },
        'population_data': {
            'gnomad_af': 0.0,
            'gnomad_homozygous': 0,
            'exac_af': 0.0
        },
        'insilico_data': {
            'revel_score': 0.92,
            'cadd_phred': 32.0,
            'sift_score': 0.0,
            'polyphen2_score': 1.0
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'not_available',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'assumed',
            'parental_confirmation': 'not_confirmed',
            'maternity_confirmed': False,
            'paternity_confirmed': False,
            'inheritance_pattern': 'x_linked_dominant',
            'zygosity': 'heterozygous'
        },
        'patient_phenotypes': 'Rett syndrome, developmental regression'
    },
    '16_pp1_cosegregation': {
        'name': 'PP1 - Cosegregation in Family',
        'description': 'Variant cosegregates with disease in multiple affected family members',
        'expected_classification': 'Likely Pathogenic',
        'expected_criteria': ['PP1', 'PM2', 'PP3'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': True,
            'PP4': True,
            'BS4': False
        },
        'basic_info': {
            'gene': 'PKD1',
            'chromosome': '16',
            'position': '2138710',
            'ref_allele': 'G',
            'alt_allele': 'A',
            'amino_acid_change': 'p.Arg2220Gln',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_001009944.3',
            'hgvs_cdna': 'c.6659G>A',
            'hgvs_protein': 'p.Arg2220Gln'
        },
        'population_data': {
            'gnomad_af': 0.0,
            'gnomad_homozygous': 0,
            'exac_af': 0.0
        },
        'insilico_data': {
            'revel_score': 0.68,
            'cadd_phred': 25.5,
            'sift_score': 0.01,
            'polyphen2_score': 0.92
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'cosegregates',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'autosomal_dominant',
            'zygosity': 'heterozygous',
            'segregation_families': 1,
            'segregation_affected': 5,
            'segregation_unaffected': 0
        },
        'patient_phenotypes': 'polycystic kidney disease, renal cysts'
    },
    '17_pp5_clinvar_path': {
        'name': 'PP5 - ClinVar Pathogenic',
        'description': 'Reported as pathogenic in reputable source (ClinVar)',
        'expected_classification': 'Likely Pathogenic',
        'expected_criteria': ['PP5', 'PM2', 'PP3'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': False,
            'BS4': False
        },
        'basic_info': {
            'gene': 'BRCA2',
            'chromosome': '13',
            'position': '32340300',
            'ref_allele': 'C',
            'alt_allele': 'T',
            'amino_acid_change': 'p.Arg1699Trp',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_000059.4',
            'hgvs_cdna': 'c.5095C>T',
            'hgvs_protein': 'p.Arg1699Trp'
        },
        'population_data': {
            'gnomad_af': 0.00001,
            'gnomad_homozygous': 0,
            'exac_af': 0.0
        },
        'insilico_data': {
            'revel_score': 0.75,
            'cadd_phred': 26.0,
            'sift_score': 0.02,
            'polyphen2_score': 0.90
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
        },
        'clinvar_data': {
            'clinical_significance': 'Pathogenic',
            'review_status': 'criteria provided, multiple submitters',
            'submitter': 'Multiple ClinVar submitters',
            'submitter_count': 3,
            'last_evaluated': '2024-05-15'
        }
    },
    '18_bs1_high_maf': {
        'name': 'BS1 - High MAF for Disorder',
        'description': 'MAF too high for rare disorder prevalence',
        'expected_classification': 'Likely Benign',
        'expected_criteria': ['BS1', 'BP4'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': False,
            'BS4': False
        },
        'basic_info': {
            'gene': 'HFE',
            'chromosome': '6',
            'position': '26093141',
            'ref_allele': 'G',
            'alt_allele': 'A',
            'amino_acid_change': 'p.Cys282Tyr',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_000410.4',
            'hgvs_cdna': 'c.845G>A',
            'hgvs_protein': 'p.Cys282Tyr'
        },
        'population_data': {
            'gnomad_af': 0.065,
            'gnomad_homozygous': 850,
            'exac_af': 0.064,
            'thousand_genomes_af': 0.048
        },
        'insilico_data': {
            'revel_score': 0.45,
            'cadd_phred': 22.0,
            'sift_score': 0.05,
            'polyphen2_score': 0.78
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'not_available',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'autosomal_recessive',
            'zygosity': 'heterozygous'
        }
    },
    '19_bs2_healthy_homozygotes': {
        'name': 'BS2 - Healthy Adult Homozygotes',
        'description': 'Found in homozygous state in healthy adults (recessive disorder)',
        'expected_classification': 'Likely Benign',
        'expected_criteria': ['BS2', 'BP4'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': False,
            'BS4': False
        },
        'basic_info': {
            'gene': 'GJB2',
            'chromosome': '13',
            'position': '20189544',
            'ref_allele': 'G',
            'alt_allele': 'A',
            'amino_acid_change': 'p.Val153Ile',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_004004.6',
            'hgvs_cdna': 'c.457G>A',
            'hgvs_protein': 'p.Val153Ile'
        },
        'population_data': {
            'gnomad_af': 0.018,
            'gnomad_homozygous': 45,
            'exac_af': 0.017
        },
        'insilico_data': {
            'revel_score': 0.32,
            'cadd_phred': 19.5,
            'sift_score': 0.12,
            'polyphen2_score': 0.42
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'not_available',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'autosomal_recessive',
            'zygosity': 'heterozygous',
            'healthy_homozygotes': 45
        }
    },
    '20_bs4_nonsegregation': {
        'name': 'BS4 - Non-segregation',
        'description': 'Does not segregate with disease in affected family members',
        'expected_classification': 'Likely Benign',
        'expected_criteria': ['BS4', 'BP4'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': False,
            'BS4': True
        },
        'basic_info': {
            'gene': 'TTN',
            'chromosome': '2',
            'position': '179390733',
            'ref_allele': 'G',
            'alt_allele': 'A',
            'amino_acid_change': 'p.Arg15440His',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_001267550.2',
            'hgvs_cdna': 'c.46319G>A',
            'hgvs_protein': 'p.Arg15440His'
        },
        'population_data': {
            'gnomad_af': 0.0008,
            'gnomad_homozygous': 0,
            'exac_af': 0.0007
        },
        'insilico_data': {
            'revel_score': 0.38,
            'cadd_phred': 20.5,
            'sift_score': 0.15,
            'polyphen2_score': 0.55
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'does_not_segregate',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'autosomal_dominant',
            'zygosity': 'heterozygous',
            'segregation_families': 1,
            'segregation_affected': 3,
            'segregation_unaffected': 2,
            'nonsegregation_details': '2 affected without variant, 2 unaffected with variant'
        }
    },
    '21_bp1_missense_truncating_gene': {
        'name': 'BP1 - Missense in Truncating Mechanism Gene',
        'description': 'Missense change in gene where LOF is mechanism (truncating expected)',
        'expected_classification': 'Likely Benign',
        'expected_criteria': ['BP1', 'BP4'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': False,
            'BS4': False
        },
        'basic_info': {
            'gene': 'TSC1',
            'chromosome': '9',
            'position': '132910394',
            'ref_allele': 'G',
            'alt_allele': 'A',
            'amino_acid_change': 'p.Ala456Thr',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_000368.5',
            'hgvs_cdna': 'c.1366G>A',
            'hgvs_protein': 'p.Ala456Thr'
        },
        'population_data': {
            'gnomad_af': 0.0002,
            'gnomad_homozygous': 0,
            'exac_af': 0.00015
        },
        'insilico_data': {
            'revel_score': 0.28,
            'cadd_phred': 18.0,
            'sift_score': 0.25,
            'polyphen2_score': 0.35
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
    },
    '22_bp2_trans_with_pathogenic': {
        'name': 'BP2 - Trans with Pathogenic (Dominant Gene)',
        'description': 'Observed in trans with pathogenic variant in dominant gene',
        'expected_classification': 'Likely Benign',
        'expected_criteria': ['BP2', 'BP4'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': False,
            'BS4': False
        },
        'basic_info': {
            'gene': 'NF1',
            'chromosome': '17',
            'position': '31226810',
            'ref_allele': 'C',
            'alt_allele': 'T',
            'amino_acid_change': 'p.Pro1427Leu',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_000267.3',
            'hgvs_cdna': 'c.4280C>T',
            'hgvs_protein': 'p.Pro1427Leu'
        },
        'population_data': {
            'gnomad_af': 0.0001,
            'gnomad_homozygous': 0,
            'exac_af': 0.00008
        },
        'insilico_data': {
            'revel_score': 0.35,
            'cadd_phred': 19.0,
            'sift_score': 0.18,
            'polyphen2_score': 0.48
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
            'zygosity': 'compound_heterozygous',
            'trans_pathogenic': 'c.2970_2972delAAT (p.Met990del)',
            'second_variant_classification': 'Pathogenic',
            'individual_phenotype': 'unaffected'
        }
    },
    '23_bp3_inframe_repeat': {
        'name': 'BP3 - In-frame Indel in Repeat',
        'description': 'In-frame deletion in repetitive region without known function',
        'expected_classification': 'Likely Benign',
        'expected_criteria': ['BP3', 'BP4'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': False,
            'BS4': False
        },
        'basic_info': {
            'gene': 'AR',
            'chromosome': 'X',
            'position': '67545316',
            'ref_allele': 'GCAGCAGCA',
            'alt_allele': 'G',
            'amino_acid_change': 'p.Gln18_Gln20del',
            'variant_type': 'inframe_deletion',
            'consequence': 'inframe_deletion',
            'transcript': 'NM_000044.6',
            'hgvs_cdna': 'c.52_60delCAGCAGCAG',
            'hgvs_protein': 'p.Gln18_Gln20del'
        },
        'population_data': {
            'gnomad_af': 0.008,
            'gnomad_homozygous': 0,
            'exac_af': 0.007
        },
        'insilico_data': {
            'cadd_phred': 12.5
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'not_available',
            'functional_studies': 'benign'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'x_linked',
            'zygosity': 'hemizygous',
            'repeat_region': 'CAG polyglutamine tract'
        }
    },
    '24_bp6_clinvar_benign': {
        'name': 'BP6 - ClinVar Benign',
        'description': 'Reported as benign in reputable source (ClinVar)',
        'expected_classification': 'Likely Benign',
        'expected_criteria': ['BP6', 'BP4'],
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
            'position': '43104121',
            'ref_allele': 'A',
            'alt_allele': 'C',
            'amino_acid_change': 'p.Glu1038Gly',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_007294.4',
            'hgvs_cdna': 'c.3113A>G',
            'hgvs_protein': 'p.Glu1038Gly'
        },
        'population_data': {
            'gnomad_af': 0.0005,
            'gnomad_homozygous': 0,
            'exac_af': 0.0004
        },
        'insilico_data': {
            'revel_score': 0.22,
            'cadd_phred': 17.0,
            'sift_score': 0.35,
            'polyphen2_score': 0.28
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
        },
        'clinvar_data': {
            'clinical_significance': 'Benign',
            'review_status': 'reviewed by expert panel',
            'submitter': 'ClinGen Expert Panel',
            'submitter_count': 5,
            'last_evaluated': '2024-08-20'
        }
    },
    '25_pvs1_splice_acceptor': {
        'name': 'PVS1 - Splice Acceptor Variant',
        'description': 'Canonical splice acceptor site variant (±1 or 2 positions)',
        'expected_classification': 'Pathogenic',
        'expected_criteria': ['PVS1', 'PM2'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': True,
            'BS4': False
        },
        'basic_info': {
            'gene': 'PAH',
            'chromosome': '12',
            'position': '103232816',
            'ref_allele': 'G',
            'alt_allele': 'A',
            'amino_acid_change': 'p.?',
            'variant_type': 'splice_acceptor',
            'consequence': 'splice_acceptor_variant',
            'transcript': 'NM_000277.3',
            'hgvs_cdna': 'c.1066-2A>G',
            'hgvs_protein': 'p.?'
        },
        'population_data': {
            'gnomad_af': 0.0,
            'gnomad_homozygous': 0,
            'exac_af': 0.0
        },
        'insilico_data': {
            'spliceai_ag_score': 0.98,
            'spliceai_al_score': 0.03,
            'cadd_phred': 33.5
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'not_available',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'autosomal_recessive',
            'zygosity': 'compound_heterozygous',
            'trans_pathogenic': 'c.1222C>T (p.Arg408Trp)',
            'second_variant_classification': 'Pathogenic'
        },
        'patient_phenotypes': 'phenylketonuria, intellectual disability'
    },
    '26_pvs1_start_loss': {
        'name': 'PVS1 - Start Loss Variant',
        'description': 'Loss of initiator methionine (start codon)',
        'expected_classification': 'Pathogenic',
        'expected_criteria': ['PVS1', 'PM2'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': False,
            'BS4': False
        },
        'basic_info': {
            'gene': 'HBB',
            'chromosome': '11',
            'position': '5227002',
            'ref_allele': 'A',
            'alt_allele': 'C',
            'amino_acid_change': 'p.Met1?',
            'variant_type': 'start_lost',
            'consequence': 'start_lost',
            'transcript': 'NM_000518.5',
            'hgvs_cdna': 'c.1A>C',
            'hgvs_protein': 'p.Met1?'
        },
        'population_data': {
            'gnomad_af': 0.0,
            'gnomad_homozygous': 0,
            'exac_af': 0.0
        },
        'insilico_data': {
            'cadd_phred': 35.0
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'not_available',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'autosomal_recessive',
            'zygosity': 'homozygous'
        },
        'patient_phenotypes': 'beta-thalassemia major'
    },
    '27_pvs1_stop_loss': {
        'name': 'PVS1 - Stop Loss Variant',
        'description': 'Loss of stop codon leading to read-through',
        'expected_classification': 'Pathogenic',
        'expected_criteria': ['PVS1', 'PM2'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': False,
            'BS4': False
        },
        'basic_info': {
            'gene': 'RUNX1',
            'chromosome': '21',
            'position': '36231780',
            'ref_allele': 'TGA',
            'alt_allele': 'TTA',
            'amino_acid_change': 'p.Ter481Leuext*?',
            'variant_type': 'stop_lost',
            'consequence': 'stop_lost',
            'transcript': 'NM_001754.5',
            'hgvs_cdna': 'c.1441_1443delTGAinsTTA',
            'hgvs_protein': 'p.Ter481Leuext*?'
        },
        'population_data': {
            'gnomad_af': 0.0,
            'gnomad_homozygous': 0,
            'exac_af': 0.0
        },
        'insilico_data': {
            'cadd_phred': 31.0,
            'revel_score': 0.82
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
        },
        'patient_phenotypes': 'familial platelet disorder, AML predisposition'
    },
    '28_pm4_inframe_insertion': {
        'name': 'PM4 - In-frame Insertion',
        'description': 'In-frame insertion in critical domain',
        'expected_classification': 'Likely Pathogenic',
        'expected_criteria': ['PM4', 'PM1', 'PM2', 'PP3'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': False,
            'BS4': False
        },
        'basic_info': {
            'gene': 'FLT3',
            'chromosome': '13',
            'position': '28608258',
            'ref_allele': 'A',
            'alt_allele': 'AGCCAATATATGAACAGAAACTCTTT',
            'amino_acid_change': 'p.Tyr597_Glu598insTyrGluValAsnIlePhe',
            'variant_type': 'inframe_insertion',
            'consequence': 'inframe_insertion',
            'transcript': 'NM_004119.3',
            'hgvs_cdna': 'c.1791_1792insTTCAATATATGGAACAAAAGTCTTT',
            'hgvs_protein': 'p.Tyr597_Glu598insTyrGluValAsnIlePhe'
        },
        'population_data': {
            'gnomad_af': 0.0,
            'gnomad_homozygous': 0,
            'exac_af': 0.0
        },
        'insilico_data': {
            'cadd_phred': 28.0
        },
        'functional_data': {
            'case_control': 'pathogenic',
            'segregation': 'not_available',
            'functional_studies': 'damaging'
        },
        'genetic_data': {
            'de_novo_status': 'somatic',
            'parental_confirmation': 'not_applicable',
            'inheritance_pattern': 'somatic',
            'zygosity': 'heterozygous'
        },
        'patient_phenotypes': 'acute myeloid leukemia'
    },
    '29_pp1_strong_multi_family': {
        'name': 'PP1_Strong - Multi-family Cosegregation',
        'description': 'Strong cosegregation evidence from multiple families',
        'expected_classification': 'Likely Pathogenic',
        'expected_criteria': ['PP1_Strong', 'PM2', 'PP3'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': True,
            'PP4': True,
            'BS4': False
        },
        'basic_info': {
            'gene': 'KCNH2',
            'chromosome': '7',
            'position': '150645213',
            'ref_allele': 'G',
            'alt_allele': 'A',
            'amino_acid_change': 'p.Arg534His',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_000238.4',
            'hgvs_cdna': 'c.1601G>A',
            'hgvs_protein': 'p.Arg534His'
        },
        'population_data': {
            'gnomad_af': 0.0,
            'gnomad_homozygous': 0,
            'exac_af': 0.0
        },
        'insilico_data': {
            'revel_score': 0.79,
            'cadd_phred': 27.5,
            'sift_score': 0.0,
            'polyphen2_score': 0.95
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'cosegregates',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'autosomal_dominant',
            'zygosity': 'heterozygous',
            'segregation_families': 4,
            'segregation_affected': 12,
            'segregation_unaffected': 0,
            'lod_score': 3.5
        },
        'patient_phenotypes': 'Long QT syndrome type 2, sudden cardiac death'
    },
    '30_pm3_homozygous_recessive': {
        'name': 'PM3 - Homozygous Recessive',
        'description': 'Homozygous variant in recessive disorder',
        'expected_classification': 'Likely Pathogenic',
        'expected_criteria': ['PM3_Strong', 'PM2', 'PP3'],
        'default_responses': {
            'PS1': False,
            'PM5': False,
            'PP1': False,
            'PP4': True,
            'BS4': False
        },
        'basic_info': {
            'gene': 'ABCA4',
            'chromosome': '1',
            'position': '94543994',
            'ref_allele': 'C',
            'alt_allele': 'T',
            'amino_acid_change': 'p.Gly1961Glu',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'transcript': 'NM_000350.3',
            'hgvs_cdna': 'c.5882G>A',
            'hgvs_protein': 'p.Gly1961Glu'
        },
        'population_data': {
            'gnomad_af': 0.00005,
            'gnomad_homozygous': 0,
            'exac_af': 0.00004
        },
        'insilico_data': {
            'revel_score': 0.83,
            'cadd_phred': 28.5,
            'sift_score': 0.0,
            'polyphen2_score': 0.98
        },
        'functional_data': {
            'case_control': 'not_available',
            'segregation': 'not_available',
            'functional_studies': 'not_available'
        },
        'genetic_data': {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'both_heterozygous',
            'inheritance_pattern': 'autosomal_recessive',
            'zygosity': 'homozygous',
            'parental_genotypes': 'both heterozygous carriers'
        },
        'patient_phenotypes': 'Stargardt disease, macular degeneration'
    }
}

# Legacy TEST_MODE_DATA for backward compatibility (points to first scenario)
TEST_MODE_DATA = TEST_SCENARIOS['1_pvs1_nonsense']

# API endpoints
API_ENDPOINTS = {
    # Clinical databases
    'clinvar': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/',
    'clinvar_api': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi',
    'clingen_erepo': 'https://erepo.genome.network/evrepo/api',
    'clingen_search': 'https://search.clinicalgenome.org/kb/gene-validity/',
    'clingen_dosage_tsv': 'https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv',
    
    # Population frequency databases
    'gnomad_graphql': 'https://gnomad.broadinstitute.org/api/',
    'gnomad_v4': 'https://gnomad.broadinstitute.org/api/v4/',
    'gnomad_browser': 'https://gnomad.broadinstitute.org/api',
    
    # In silico prediction tools
    'dbnsfp': 'https://dbnsfp.s3.amazonaws.com/dbNSFP4.4a/',  # Public S3 bucket
    'alphamissense': 'https://storage.googleapis.com/alphamissense/',  # Google DeepMind
    'cravat': 'https://run.opencravat.org/submit/annotate',
    'varsome': 'https://api.varsome.com/lookup/',
    
    # Genome annotation
    'ensembl_rest': 'https://rest.ensembl.org',
    'ensembl_vep': 'https://rest.ensembl.org/vep/human/hgvs/',
    'ucsc_genome': 'https://api.genome.ucsc.edu/'
}

# API Settings for gnomAD and other external services
API_SETTINGS = {
    'enabled': True,  # Master switch for all API integrations
    'timeout': 15,  # Request timeout in seconds
    'max_retries': 3,  # Maximum retry attempts
    'cache_ttl': 86400,  # Cache time-to-live in seconds (24 hours)
    'rate_limit': {
        'calls': 10,  # Maximum calls
        'period': 1  # Per period in seconds
    },
    'fallback_to_manual': True  # Use manual input if API fails
}

# gnomAD Constraint Thresholds for LOF Intolerance
CONSTRAINT_THRESHOLDS = {
    'pLI_intolerant': 0.9,  # pLI ≥ 0.9 indicates LOF intolerant
    'LOEUF_intolerant': 0.35,  # LOEUF ≤ 0.35 indicates LOF intolerant
    'oe_lof_upper_tolerant': 0.6,  # oe_lof_upper > 0.6 indicates LOF tolerant
    'confidence_threshold': 'high'  # Minimum confidence for classification
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
