"""
ACMG/AMP Criteria Configuration
==============================

Contains ACMG evidence weights, classification rules, and statistical thresholds
for both 2015 and 2023 guidelines.

Author: Can Sevilmiş
License: MIT License
"""

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