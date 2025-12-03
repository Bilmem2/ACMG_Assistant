"""
Gene-Specific Rules and Thresholds
==================================

Contains gene-specific thresholds and LOF gene classifications for ACMG criteria.

Author: Can Sevilmiş
License: MIT License
"""

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