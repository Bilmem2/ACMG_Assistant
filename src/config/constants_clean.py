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
    'version': '3.0.0',
    'author': 'Can Sevilmiş',
    'license': 'MIT License',
    'last_updated': 'July 10, 2025',
    'guidelines': 'ACMG/AMP 2015 & 2023',
    'description': 'ACMG Variant Classification Assistant with Complete Criteria Implementation',
    'major_features': [
        'Complete 28 ACMG/AMP criteria implementation',
        'Interactive evidence evaluation',
        'Enhanced computational metascore',
        'User-driven criteria assignment',
        'Dynamic evidence weighting'
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
STATISTICAL_THRESHOLDS = {
    'fisher_exact_p_value': 0.05,
    'case_control_odds_ratio': 2.0,
    'segregation_lod_score': 3.0,
    'splice_ai_threshold': 0.5,
    'conservation_threshold': 2.0
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
    'BOLD': '\033[1m',
    'UNDERLINE': '\033[4m'
}

# Section headers for organized output
SECTION_HEADERS = {
    'variant_info': f"{COLORAMA_COLORS['CYAN']}🧬 VARIANT INFORMATION{COLORAMA_COLORS['ENDC']}",
    'pathogenic_evidence': f"{COLORAMA_COLORS['RED']}🔴 PATHOGENIC EVIDENCE{COLORAMA_COLORS['ENDC']}",
    'benign_evidence': f"{COLORAMA_COLORS['GREEN']}🟢 BENIGN EVIDENCE{COLORAMA_COLORS['ENDC']}",
    'final_classification': f"{COLORAMA_COLORS['BOLD']}📊 FINAL CLASSIFICATION{COLORAMA_COLORS['ENDC']}",
    'recommendations': f"{COLORAMA_COLORS['YELLOW']}💡 RECOMMENDATIONS{COLORAMA_COLORS['ENDC']}"
}

# Classification-specific colors
CLASSIFICATION_COLORS = {
    'Pathogenic': f"{COLORAMA_COLORS['RED']}{COLORAMA_COLORS['BOLD']}",
    'Likely Pathogenic': f"{COLORAMA_COLORS['RED']}",
    'Uncertain Significance': f"{COLORAMA_COLORS['YELLOW']}",
    'Likely Benign': f"{COLORAMA_COLORS['GREEN']}",
    'Benign': f"{COLORAMA_COLORS['GREEN']}{COLORAMA_COLORS['BOLD']}"
}

def get_colored_message(message: str, color_key: str) -> str:
    """Get a colored message for terminal output."""
    color = COLORAMA_COLORS.get(color_key, '')
    return f"{color}{message}{COLORAMA_COLORS['ENDC']}"
