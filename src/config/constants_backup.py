"""
Configuration constants for ACMG Variant Classification Assistant
===============================================================

This file contains all the constants, thresholds, and configuration
values used throughout the application.

Author: Can Sevilmiş
License: MIT License
Version: 2.0.0
Last Updated: July 8, 2025
"""

# Version and metadata information
VERSION_INFO = {
    'version': '3.0.0',
    'author': 'Can Sevilmiş',
    'license': 'MIT License',
    'last_updated': 'July 9, 2025',
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

# ACMG Evidence Criteria Weights and Descriptions
EVIDENCE_WEIGHTS = {
    # Pathogenic criteria
    'PVS1': {'weight': 'Very Strong', 'category': 'Pathogenic', 'description': 'Null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multiexon deletion) in a gene where LOF is a known mechanism of disease'},
    'PS1': {'weight': 'Strong', 'category': 'Pathogenic', 'description': 'Same amino acid change as a previously established pathogenic variant regardless of nucleotide change'},
    'PS2': {'weight': 'Strong', 'category': 'Pathogenic', 'description': 'De novo (both maternity and paternity confirmed) in a patient with the disease and no family history'},
    'PS2_Very_Strong': {'weight': 'Very Strong', 'category': 'Pathogenic', 'description': 'De novo (confirmed) in a patient with the disease and unaffected parents (2023 guidelines)'},
    'PS3': {'weight': 'Strong', 'category': 'Pathogenic', 'description': 'Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product'},
    'PS4': {'weight': 'Strong', 'category': 'Pathogenic', 'description': 'The prevalence of the variant in affected individuals is significantly increased compared to the prevalence in controls'},
    'PM1': {'weight': 'Moderate', 'category': 'Pathogenic', 'description': 'Located in a mutational hot spot and/or critical and well-established functional domain without benign variation'},
    'PM2': {'weight': 'Moderate', 'category': 'Pathogenic', 'description': 'Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium'},
    'PM3': {'weight': 'Moderate', 'category': 'Pathogenic', 'description': 'For recessive disorders, detected in trans with a pathogenic variant'},
    'PM4': {'weight': 'Moderate', 'category': 'Pathogenic', 'description': 'Protein length changes as a result of in-frame deletions/insertions in a non-repeat region or stop-loss variants'},
    'PM5': {'weight': 'Moderate', 'category': 'Pathogenic', 'description': 'Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before'},
    'PM6': {'weight': 'Moderate', 'category': 'Pathogenic', 'description': 'Assumed de novo, but without confirmation of paternity and maternity'},
    'PP1': {'weight': 'Supporting', 'category': 'Pathogenic', 'description': 'Cosegregation with disease in multiple affected family members in a gene definitively known to cause the disease'},
    'PP2': {'weight': 'Supporting', 'category': 'Pathogenic', 'description': 'Missense variant in a gene that has a low rate of benign missense variation and where missense variants are a common mechanism of disease'},
    'PP3': {'weight': 'Supporting', 'category': 'Pathogenic', 'description': 'Multiple lines of computational evidence support a deleterious effect on the gene or gene product'},
    'PP4': {'weight': 'Supporting', 'category': 'Pathogenic', 'description': 'Patient\'s phenotype or family history is highly specific for a disease with a single genetic etiology'},
    'PP5': {'weight': 'Supporting', 'category': 'Pathogenic', 'description': 'Reputable source recently reports variant as pathogenic, but the evidence is not available to the laboratory to perform an independent evaluation'},
    
    # Benign criteria
    'BA1': {'weight': 'Stand-alone', 'category': 'Benign', 'description': 'Allele frequency is >5% in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium'},
    'BS1': {'weight': 'Strong', 'category': 'Benign', 'description': 'Allele frequency is greater than expected for disorder'},
    'BS2': {'weight': 'Strong', 'category': 'Benign', 'description': 'Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked (hemizygous) disorder, with full penetrance expected at an early age'},
    'BS3': {'weight': 'Strong', 'category': 'Benign', 'description': 'Well-established in vitro or in vivo functional studies show no damaging effect on protein function or splicing'},
    'BS4': {'weight': 'Strong', 'category': 'Benign', 'description': 'Lack of segregation in affected members of a family'},
    'BP1': {'weight': 'Supporting', 'category': 'Benign', 'description': 'Missense variant in a gene for which primarily truncating variants are known to cause disease'},
    'BP2': {'weight': 'Supporting', 'category': 'Benign', 'description': 'Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder or observed in cis with a pathogenic variant in any inheritance pattern'},
    'BP3': {'weight': 'Supporting', 'category': 'Benign', 'description': 'In-frame deletions/insertions in a repetitive region without a known function'},
    'BP4': {'weight': 'Supporting', 'category': 'Benign', 'description': 'Multiple lines of computational evidence suggest no impact on gene or gene product'},
    'BP5': {'weight': 'Supporting', 'category': 'Benign', 'description': 'Variant found in a case with an alternate molecular basis for disease'},
    'BP6': {'weight': 'Supporting', 'category': 'Benign', 'description': 'Reputable source recently reports variant as benign, but the evidence is not available to the laboratory to perform an independent evaluation'},
    'BP7': {'weight': 'Supporting', 'category': 'Benign', 'description': 'A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the splice consensus sequence nor the creation of a new splice site AND the nucleotide is not highly conserved'}
}

# Classification rules
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
            {'strong': 1, 'supporting': 1}
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

# LOF mechanism gene classification for PVS1 evaluation
# Based on ExAC/gnomAD constraint scores and clinical literature
LOF_INTOLERANT_GENES = {
    # High-confidence LOF intolerant genes (pLI > 0.9)
    'BRCA1', 'BRCA2', 'TP53', 'RB1', 'APC', 'VHL', 'NF1', 'NF2',
    'CDKN2A', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'ATM', 'CHEK2',
    'PALB2', 'MYH7', 'MYBPC3', 'SCN1A', 'SCN2A', 'MECP2', 'PAH', 'CFTR'
}

LOF_TOLERANT_GENES = {
    # Genes with high tolerance to LOF variants (pLI < 0.1)
    'TTN', 'MUC16', 'OBSCN', 'PCLO', 'RYR1', 'SYNE1', 'SYNE2', 'USH2A', 'FLG'
}

# Genes where LOF mechanism is disease-dependent
LOF_CONTEXT_DEPENDENT = {
    'APOE': 'Alzheimer disease context',
    'HFE': 'Hemochromatosis context', 
    'F5': 'Factor V Leiden context',
    'SERPINA1': 'Alpha-1 antitrypsin deficiency context'
}

# In silico predictor weights for Computational Metascore
INSILICO_WEIGHTS = {
    'revel': 0.25,
    'cadd_phred': 0.20,
    'alphamissense': 0.15,
    'metarnn': 0.08,
    'metarnn_ranked': 0.08,
    'clinpred': 0.10,
    'bayesdel_addaf': 0.06,
    'bayesdel_noaf': 0.06,
    'sift': 0.05,
    'polyphen2': 0.05,
    'mutationtaster': 0.03,
    'mutationtaster_ranked': 0.03,
    'fathmm': 0.02,
    'fathmm_ranked': 0.02,
    # High priority additions
    'vest4': 0.12,           # Cancer-specific predictor
    'primateai': 0.15,       # State-of-the-art missense prediction
    'esm1b': 0.13,           # Protein language model
    'provean': 0.04,         # Protein variation effect
    'gerp_pp': 0.02,
    'phylop_vert': 0.015,
    'phylop_vert_ranked': 0.015,
    'phylop_mamm': 0.015,
    'phylop_mamm_ranked': 0.015,
    'phylop_prim': 0.01,
    'phylop_prim_ranked': 0.01,
    # SpliceAI scores (important for intronic and splice variants)
    'spliceai_ag': 0.12,
    'spliceai_al': 0.12,
    'spliceai_dg': 0.12,
    'spliceai_dl': 0.12,
    'spliceai_max': 0.25,  # Maximum SpliceAI score
    # Enhanced splice predictors
    'mmsplice': 0.10,        # Modular modeling of splicing
    'ada_score': 0.08,
    'rf_score': 0.08,
    'dbscsnv_ada': 0.06,
    'dbscsnv_rf': 0.06,
    # Conservation predictors
    'phastcons_vert': 0.02,
    'phastcons_mamm': 0.02,
    'phastcons_prim': 0.015,
    'phastcons_vert_ranked': 0.02,
    'phastcons_mamm_ranked': 0.02,
    'phastcons_prim_ranked': 0.015,
    'siphy': 0.015,
    # Additional missense predictors
    'metasvm': 0.09,
    'metasvm_ranked': 0.09,
    'metalr': 0.09,
    'metalr_ranked': 0.09,
    'mutationassessor': 0.06,
    'mutpred': 0.04,
    'mutpred_ranked': 0.04,
    'lrt': 0.03,
    'lrt_ranked': 0.03,
    'bayesdel_addaf': 0.06,
    'bayesdel_noaf': 0.06,
    # Functional predictors
    'fitcons': 0.05
}

# In silico predictor thresholds
INSILICO_THRESHOLDS = {
    'revel': {'pathogenic': 0.75, 'benign': 0.25},
    'cadd_phred': {'pathogenic': 25, 'benign': 10},
    'alphamissense': {'pathogenic': 0.564, 'benign': 0.34},
    'metarnn': {'pathogenic': 0.5, 'benign': 0.5},
    'metarnn_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'clinpred': {'pathogenic': 0.5, 'benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},
    'sift': {'pathogenic': 0.05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)
    'polyphen2': {'pathogenic': 0.85, 'benign': 0.15},
    'metalr': {'pathogenic': 0.5, 'benign': 0.5},
    'metalr_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'mutationassessor': {'pathogenic': 3.5, 'benign': 1.5},
    'mutpred': {'pathogenic': 0.7, 'benign': 0.3},
    'mutpred_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'lrt': {'pathogenic': 0.5, 'benign': 0.5},
    'lrt_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},
    # Functional predictors
    'fitcons': {'pathogenic': 0.7, 'benign': 0.3}
}

# In silico predictor thresholds': 0.5},
INSILICO_THRESHOLDS = {': 0.5},
    'revel': {'pathogenic': 0.75, 'benign': 0.25},M is inverted (lower = more pathogenic)
    'cadd_phred': {'pathogenic': 25, 'benign': 10}, 0.5},
    'alphamissense': {'pathogenic': 0.564, 'benign': 0.34},
    'metarnn': {'pathogenic': 0.5, 'benign': 0.5},},
    'metarnn_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'clinpred': {'pathogenic': 0.5, 'benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},VEAN is inverted (lower = more pathogenic)
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},': 5.0, 'benign': 2.0},
    'sift': {'pathogenic': 0.05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)0.5},
    'polyphen2': {'pathogenic': 0.85, 'benign': 0.15},   'phylop_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'mutationtaster': {'pathogenic': 0.5, 'benign': 0.5},    'phylop_mamm': {'pathogenic': 2.0, 'benign': 0.5},
    'mutationtaster_ranked': {'pathogenic': 0.5, 'benign': 0.5},ogenic': 0.5, 'benign': 0.5},
    'fathmm': {'pathogenic': -1.5, 'benign': 1.5},  # FATHMM is inverted (lower = more pathogenic)thogenic': 2.0, 'benign': 0.5},
    'fathmm_ranked': {'pathogenic': 0.5, 'benign': 0.5},ign': 0.5},
    # High priority predictor thresholds)
    'vest4': {'pathogenic': 0.7, 'benign': 0.3},
    'primateai': {'pathogenic': 0.8, 'benign': 0.2},.05},
    'esm1b': {'pathogenic': 0.5, 'benign': 0.5},
    'provean': {'pathogenic': -2.5, 'benign': -2.5},  # PROVEAN is inverted (lower = more pathogenic)05},
    'gerp_pp': {'pathogenic': 5.0, 'benign': 2.0},
    'phylop_vert': {'pathogenic': 2.0, 'benign': 0.5},
    'phylop_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phylop_mamm': {'pathogenic': 2.0, 'benign': 0.5},
    'phylop_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phylop_prim': {'pathogenic': 2.0, 'benign': 0.5},
    'phylop_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    # SpliceAI thresholds (splice-altering variants)': 0.5},
    'spliceai_ag': {'pathogenic': 0.2, 'benign': 0.05},c': 0.5, 'benign': 0.5},
    'spliceai_al': {'pathogenic': 0.2, 'benign': 0.05}, 'benign': 0.5},
    'spliceai_dg': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_dl': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_max': {'pathogenic': 0.5, 'benign': 0.1},
    # Enhanced splice predictors: 0.5},
    'mmsplice': {'pathogenic': 2.0, 'benign': 0.5},
    # Conservation predictors
    'phastcons_vert': {'pathogenic': 0.9, 'benign': 0.3},1.5},
    'phastcons_mamm': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_prim': {'pathogenic': 0.9, 'benign': 0.3},5},
    'phastcons_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phastcons_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},,
    'phastcons_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},.18},
    'siphy': {'pathogenic': 12.0, 'benign': 3.0},18},
    # Additional missense predictors
    'metasvm': {'pathogenic': 0.83, 'benign': 0.17},
    'metasvm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'metalr': {'pathogenic': 0.5, 'benign': 0.5},
    'metalr_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'mutationassessor': {'pathogenic': 3.5, 'benign': 1.5},
    'mutpred': {'pathogenic': 0.7, 'benign': 0.3},
    'mutpred_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'lrt': {'pathogenic': 0.5, 'benign': 0.5},},
    'lrt_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},
    # Functional predictorsn': -0.18},
    'fitcons': {'pathogenic': 0.7, 'benign': 0.3}0.16, 'benign': -0.18},
} Note: SIFT is inverted (lower = more pathogenic)

# In silico predictor thresholds': 0.5},
INSILICO_THRESHOLDS = {': 0.5},
    'revel': {'pathogenic': 0.75, 'benign': 0.25},M is inverted (lower = more pathogenic)
    'cadd_phred': {'pathogenic': 25, 'benign': 10}, 0.5},
    'alphamissense': {'pathogenic': 0.564, 'benign': 0.34},
    'metarnn': {'pathogenic': 0.5, 'benign': 0.5},},
    'metarnn_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'clinpred': {'pathogenic': 0.5, 'benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},VEAN is inverted (lower = more pathogenic)
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},': 5.0, 'benign': 2.0},
    'sift': {'pathogenic': 0.05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)0.5},
    'polyphen2': {'pathogenic': 0.85, 'benign': 0.15},   'phylop_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'mutationtaster': {'pathogenic': 0.5, 'benign': 0.5},    'phylop_mamm': {'pathogenic': 2.0, 'benign': 0.5},
    'mutationtaster_ranked': {'pathogenic': 0.5, 'benign': 0.5},ogenic': 0.5, 'benign': 0.5},
    'fathmm': {'pathogenic': -1.5, 'benign': 1.5},  # FATHMM is inverted (lower = more pathogenic)thogenic': 2.0, 'benign': 0.5},
    'fathmm_ranked': {'pathogenic': 0.5, 'benign': 0.5},ign': 0.5},
    # High priority predictor thresholds)
    'vest4': {'pathogenic': 0.7, 'benign': 0.3},
    'primateai': {'pathogenic': 0.8, 'benign': 0.2},.05},
    'esm1b': {'pathogenic': 0.5, 'benign': 0.5},
    'provean': {'pathogenic': -2.5, 'benign': -2.5},  # PROVEAN is inverted (lower = more pathogenic)05},
    'gerp_pp': {'pathogenic': 5.0, 'benign': 2.0},
    'phylop_vert': {'pathogenic': 2.0, 'benign': 0.5},
    'phylop_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phylop_mamm': {'pathogenic': 2.0, 'benign': 0.5},
    'phylop_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phylop_prim': {'pathogenic': 2.0, 'benign': 0.5},
    'phylop_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    # SpliceAI thresholds (splice-altering variants)
    'spliceai_ag': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_al': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_dg': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_dl': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_max': {'pathogenic': 0.5, 'benign': 0.1},
    # Enhanced splice predictors
    'mmsplice': {'pathogenic': 2.0, 'benign': 0.5},
    # Conservation predictors
    'phastcons_vert': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_mamm': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_prim': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phastcons_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phastcons_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'siphy': {'pathogenic': 12.0, 'benign': 3.0},
    # Additional missense predictors
    'metasvm': {'pathogenic': 0.83, 'benign': 0.17},
    'metasvm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'metalr': {'pathogenic': 0.5, 'benign': 0.5},
    'metalr_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'mutationassessor': {'pathogenic': 3.5, 'benign': 1.5},
    'mutpred': {'pathogenic': 0.7, 'benign': 0.3},
    'mutpred_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'lrt': {'pathogenic': 0.5, 'benign': 0.5},
    'lrt_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},
    # Functional predictorsn': 0.34},
    'fitcons': {'pathogenic': 0.7, 'benign': 0.3}benign': 0.5},
}0.5},

# In silico predictor thresholdsn': -0.18},
INSILICO_THRESHOLDS = {8},
    'revel': {'pathogenic': 0.75, 'benign': 0.25},SIFT is inverted (lower = more pathogenic)
    'cadd_phred': {'pathogenic': 25, 'benign': 10},15},
    'alphamissense': {'pathogenic': 0.564, 'benign': 0.34},
    'metarnn': {'pathogenic': 0.5, 'benign': 0.5},5, 'benign': 0.5},
    'metarnn_ranked': {'pathogenic': 0.5, 'benign': 0.5}, FATHMM is inverted (lower = more pathogenic)
    'clinpred': {'pathogenic': 0.5, 'benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18}, 0.7, 'benign': 0.3},
    'sift': {'pathogenic': 0.05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)2},
    'polyphen2': {'pathogenic': 0.85, 'benign': 0.15},   'esm1b': {'pathogenic': 0.5, 'benign': 0.5},
    'mutationtaster': {'pathogenic': 0.5, 'benign': 0.5},    'provean': {'pathogenic': -2.5, 'benign': -2.5},  # PROVEAN is inverted (lower = more pathogenic)
    'mutationtaster_ranked': {'pathogenic': 0.5, 'benign': 0.5},0, 'benign': 2.0},
    'fathmm': {'pathogenic': -1.5, 'benign': 1.5},  # FATHMM is inverted (lower = more pathogenic)thogenic': 2.0, 'benign': 0.5},
    'fathmm_ranked': {'pathogenic': 0.5, 'benign': 0.5},ign': 0.5},
    # High priority predictor thresholds5},
    'vest4': {'pathogenic': 0.7, 'benign': 0.3},},
    'primateai': {'pathogenic': 0.8, 'benign': 0.2},.5},
    'esm1b': {'pathogenic': 0.5, 'benign': 0.5},.5},
    'provean': {'pathogenic': -2.5, 'benign': -2.5},  # PROVEAN is inverted (lower = more pathogenic))
    'gerp_pp': {'pathogenic': 5.0, 'benign': 2.0},
    'phylop_vert': {'pathogenic': 2.0, 'benign': 0.5},
    'phylop_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phylop_mamm': {'pathogenic': 2.0, 'benign': 0.5},,
    'phylop_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phylop_prim': {'pathogenic': 2.0, 'benign': 0.5},
    'phylop_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    # SpliceAI thresholds (splice-altering variants)
    'spliceai_ag': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_al': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_dg': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_dl': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_max': {'pathogenic': 0.5, 'benign': 0.1},
    # Enhanced splice predictors
    'mmsplice': {'pathogenic': 2.0, 'benign': 0.5},
    # Conservation predictors
    'phastcons_vert': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_mamm': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_prim': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phastcons_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phastcons_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'siphy': {'pathogenic': 12.0, 'benign': 3.0},
    # Additional missense predictors
    'metasvm': {'pathogenic': 0.83, 'benign': 0.17},
    'metasvm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'metalr': {'pathogenic': 0.5, 'benign': 0.5},
    'metalr_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'mutationassessor': {'pathogenic': 3.5, 'benign': 1.5},
    'mutpred': {'pathogenic': 0.7, 'benign': 0.3},
    'mutpred_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'lrt': {'pathogenic': 0.5, 'benign': 0.5},
    'lrt_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},
    # Functional predictorsn': 0.34},
    'fitcons': {'pathogenic': 0.7, 'benign': 0.3}benign': 0.5},
}0.5},

# In silico predictor thresholdsn': -0.18},
INSILICO_THRESHOLDS = {8},
    'revel': {'pathogenic': 0.75, 'benign': 0.25},SIFT is inverted (lower = more pathogenic)
    'cadd_phred': {'pathogenic': 25, 'benign': 10},15},
    'alphamissense': {'pathogenic': 0.564, 'benign': 0.34},
    'metarnn': {'pathogenic': 0.5, 'benign': 0.5},5, 'benign': 0.5},
    'metarnn_ranked': {'pathogenic': 0.5, 'benign': 0.5}, FATHMM is inverted (lower = more pathogenic)
    'clinpred': {'pathogenic': 0.5, 'benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18}, 0.7, 'benign': 0.3},
    'sift': {'pathogenic': 0.05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)2},
    'polyphen2': {'pathogenic': 0.85, 'benign': 0.15},   'esm1b': {'pathogenic': 0.5, 'benign': 0.5},
    'mutationtaster': {'pathogenic': 0.5, 'benign': 0.5},    'provean': {'pathogenic': -2.5, 'benign': -2.5},  # PROVEAN is inverted (lower = more pathogenic)
    'mutationtaster_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'fathmm': {'pathogenic': -1.5, 'benign': 1.5},  # FATHMM is inverted (lower = more pathogenic)nic': 2.0, 'benign': 0.5},
    'fathmm_ranked': {'pathogenic': 0.5, 'benign': 0.5},rt_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    # High priority predictor thresholds
    'vest4': {'pathogenic': 0.7, 'benign': 0.3},'pathogenic': 0.5, 'benign': 0.5},
    'primateai': {'pathogenic': 0.8, 'benign': 0.2},hylop_prim': {'pathogenic': 2.0, 'benign': 0.5},
    'esm1b': {'pathogenic': 0.5, 'benign': 0.5},im_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'provean': {'pathogenic': -2.5, 'benign': -2.5},  # PROVEAN is inverted (lower = more pathogenic)
    'gerp_pp': {'pathogenic': 5.0, 'benign': 2.0},enic': 0.2, 'benign': 0.05},
    'phylop_vert': {'pathogenic': 2.0, 'benign': 0.5},pliceai_al': {'pathogenic': 0.2, 'benign': 0.05},
    'phylop_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},_dg': {'pathogenic': 0.2, 'benign': 0.05},
    'phylop_mamm': {'pathogenic': 2.0, 'benign': 0.5},
    'phylop_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},genic': 0.5, 'benign': 0.1},
    'phylop_prim': {'pathogenic': 2.0, 'benign': 0.5},Enhanced splice predictors
    'phylop_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},': {'pathogenic': 2.0, 'benign': 0.5},
    # SpliceAI thresholds (splice-altering variants)
    'spliceai_ag': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_al': {'pathogenic': 0.2, 'benign': 0.05},hastcons_mamm': 0.02,
    'spliceai_dg': {'pathogenic': 0.2, 'benign': 0.05},ns_prim': 0.015,
    'spliceai_dl': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_max': {'pathogenic': 0.5, 'benign': 0.1},
    # Enhanced splice predictors 0.015,
    'mmsplice': {'pathogenic': 2.0, 'benign': 0.5},
    # Conservation predictorsAdditional missense predictors
    'phastcons_vert': {'pathogenic': 0.9, 'benign': 0.3},0.09,
    'phastcons_mamm': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_prim': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phastcons_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phastcons_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},mutpred': 0.04,
    'siphy': {'pathogenic': 12.0, 'benign': 3.0},   'mutpred_ranked': 0.04,
    # Additional missense predictors    'lrt': 0.03,
    'metasvm': {'pathogenic': 0.83, 'benign': 0.17},
    'metasvm_ranked': {'pathogenic': 0.5, 'benign': 0.5},6,
    'metalr': {'pathogenic': 0.5, 'benign': 0.5},
    'metalr_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'mutationassessor': {'pathogenic': 3.5, 'benign': 1.5},
    'mutpred': {'pathogenic': 0.7, 'benign': 0.3},
    'mutpred_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'lrt': {'pathogenic': 0.5, 'benign': 0.5},lds
    'lrt_ranked': {'pathogenic': 0.5, 'benign': 0.5},NSILICO_THRESHOLDS = {
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},    'revel': {'pathogenic': 0.75, 'benign': 0.25},
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},': 10},
    # Functional predictors: {'pathogenic': 0.564, 'benign': 0.34},
    'fitcons': {'pathogenic': 0.7, 'benign': 0.3}genic': 0.5, 'benign': 0.5},
}pathogenic': 0.5, 'benign': 0.5},
nic': 0.5, 'benign': 0.5},
# In silico predictor thresholds'pathogenic': 0.16, 'benign': -0.18},
INSILICO_THRESHOLDS = {hogenic': 0.16, 'benign': -0.18},
    'revel': {'pathogenic': 0.75, 'benign': 0.25},c': 0.05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)
    'cadd_phred': {'pathogenic': 25, 'benign': 10},genic': 0.85, 'benign': 0.15},
    'alphamissense': {'pathogenic': 0.564, 'benign': 0.34},{'pathogenic': 0.5, 'benign': 0.5},
    'metarnn': {'pathogenic': 0.5, 'benign': 0.5},: {'pathogenic': 0.5, 'benign': 0.5},
    'metarnn_ranked': {'pathogenic': 0.5, 'benign': 0.5},genic': -1.5, 'benign': 1.5},  # FATHMM is inverted (lower = more pathogenic)
    'clinpred': {'pathogenic': 0.5, 'benign': 0.5},, 'benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},   # High priority predictor thresholds
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},    'vest4': {'pathogenic': 0.7, 'benign': 0.3},
    'sift': {'pathogenic': 0.05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)0.2},
    'polyphen2': {'pathogenic': 0.85, 'benign': 0.15},ogenic': 0.5, 'benign': 0.5},
    'mutationtaster': {'pathogenic': 0.5, 'benign': 0.5},
    'mutationtaster_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'fathmm': {'pathogenic': -1.5, 'benign': 1.5},  # FATHMM is inverted (lower = more pathogenic)
    'fathmm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    # High priority predictor thresholds
    'vest4': {'pathogenic': 0.7, 'benign': 0.3},
    'primateai': {'pathogenic': 0.8, 'benign': 0.2},
    'esm1b': {'pathogenic': 0.5, 'benign': 0.5},
    'provean': {'pathogenic': -2.5, 'benign': -2.5},  # PROVEAN is inverted (lower = more pathogenic)
    'gerp_pp': {'pathogenic': 5.0, 'benign': 2.0},
    'phylop_vert': {'pathogenic': 2.0, 'benign': 0.5},   'spliceai_al': {'pathogenic': 0.2, 'benign': 0.05},
    'phylop_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},    'spliceai_dg': {'pathogenic': 0.2, 'benign': 0.05},
    'phylop_mamm': {'pathogenic': 2.0, 'benign': 0.5}, 'benign': 0.05},
    'phylop_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},genic': 0.5, 'benign': 0.1},
    'phylop_prim': {'pathogenic': 2.0, 'benign': 0.5},    # Enhanced splice predictors
    'phylop_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},enign': 0.5},
    # SpliceAI thresholds (splice-altering variants)
    'spliceai_ag': {'pathogenic': 0.2, 'benign': 0.05},astcons_vert': 0.02,
    'spliceai_al': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_dg': {'pathogenic': 0.2, 'benign': 0.05},'phastcons_prim': 0.015,
    'spliceai_dl': {'pathogenic': 0.2, 'benign': 0.05},tcons_vert_ranked': 0.02,
    'spliceai_max': {'pathogenic': 0.5, 'benign': 0.1},
    # Enhanced splice predictors
    'mmsplice': {'pathogenic': 2.0, 'benign': 0.5},
    # Conservation predictors# Additional missense predictors
    'phastcons_vert': {'pathogenic': 0.9, 'benign': 0.3},': 0.09,
    'phastcons_mamm': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_prim': {'pathogenic': 0.9, 'benign': 0.3},talr': 0.09,
    'phastcons_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phastcons_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},'mutationassessor': 0.06,
    'phastcons_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'siphy': {'pathogenic': 12.0, 'benign': 3.0},
    # Additional missense predictors'lrt': 0.03,
    'metasvm': {'pathogenic': 0.83, 'benign': 0.17},
    'metasvm_ranked': {'pathogenic': 0.5, 'benign': 0.5},    'bayesdel_addaf': 0.06,
    'metalr': {'pathogenic': 0.5, 'benign': 0.5},
    'metalr_ranked': {'pathogenic': 0.5, 'benign': 0.5},redictors
    'mutationassessor': {'pathogenic': 3.5, 'benign': 1.5},
    'mutpred': {'pathogenic': 0.7, 'benign': 0.3},
    'mutpred_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'lrt': {'pathogenic': 0.5, 'benign': 0.5},sholds
    'lrt_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},: 0.75, 'benign': 0.25},
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},enic': 25, 'benign': 10},
    # Functional predictors'benign': 0.34},
    'fitcons': {'pathogenic': 0.7, 'benign': 0.3} 'benign': 0.5},
} 'benign': 0.5},
linpred': {'pathogenic': 0.5, 'benign': 0.5},
# In silico predictor thresholdspathogenic': 0.16, 'benign': -0.18},
INSILICO_THRESHOLDS = {
    'revel': {'pathogenic': 0.75, 'benign': 0.25},05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)
    'cadd_phred': {'pathogenic': 25, 'benign': 10},': 0.85, 'benign': 0.15},
    'alphamissense': {'pathogenic': 0.564, 'benign': 0.34},genic': 0.5, 'benign': 0.5},
    'metarnn': {'pathogenic': 0.5, 'benign': 0.5},{'pathogenic': 0.5, 'benign': 0.5},
    'metarnn_ranked': {'pathogenic': 0.5, 'benign': 0.5},-1.5, 'benign': 1.5},  # FATHMM is inverted (lower = more pathogenic)
    'clinpred': {'pathogenic': 0.5, 'benign': 0.5},enic': 0.5, 'benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},r thresholds
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},est4': {'pathogenic': 0.7, 'benign': 0.3},
    'sift': {'pathogenic': 0.05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)hogenic': 0.8, 'benign': 0.2},
    'polyphen2': {'pathogenic': 0.85, 'benign': 0.15},
    'mutationtaster': {'pathogenic': 0.5, 'benign': 0.5},,  # PROVEAN is inverted (lower = more pathogenic)
    'mutationtaster_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'fathmm': {'pathogenic': -1.5, 'benign': 1.5},  # FATHMM is inverted (lower = more pathogenic)
    'fathmm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    # High priority predictor thresholds
    'vest4': {'pathogenic': 0.7, 'benign': 0.3},
    'primateai': {'pathogenic': 0.8, 'benign': 0.2},
    'esm1b': {'pathogenic': 0.5, 'benign': 0.5},
    'provean': {'pathogenic': -2.5, 'benign': -2.5},  # PROVEAN is inverted (lower = more pathogenic)
    'gerp_pp': {'pathogenic': 5.0, 'benign': 2.0},
    'phylop_vert': {'pathogenic': 2.0, 'benign': 0.5},
    'phylop_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phylop_mamm': {'pathogenic': 2.0, 'benign': 0.5},
    'phylop_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phylop_prim': {'pathogenic': 2.0, 'benign': 0.5},
    'phylop_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    # SpliceAI thresholds (splice-altering variants)
    'spliceai_ag': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_al': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_dg': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_dl': {'pathogenic': 0.2, 'benign': 0.05},hastcons_vert_ranked': 0.02,
    'spliceai_max': {'pathogenic': 0.5, 'benign': 0.1},anked': 0.02,
    # Enhanced splice predictors
    'mmsplice': {'pathogenic': 2.0, 'benign': 0.5},
    # Conservation predictors
    'phastcons_vert': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_mamm': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_prim': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},etalr_ranked': 0.09,
    'phastcons_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},0.06,
    'phastcons_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'siphy': {'pathogenic': 12.0, 'benign': 3.0},
    # Additional missense predictors
    'metasvm': {'pathogenic': 0.83, 'benign': 0.17},
    'metasvm_ranked': {'pathogenic': 0.5, 'benign': 0.5},bayesdel_addaf': 0.06,
    'metalr': {'pathogenic': 0.5, 'benign': 0.5},   'bayesdel_noaf': 0.06,
    'metalr_ranked': {'pathogenic': 0.5, 'benign': 0.5},    # Functional predictors
    'mutationassessor': {'pathogenic': 3.5, 'benign': 1.5},
    'mutpred': {'pathogenic': 0.7, 'benign': 0.3},
    'mutpred_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'lrt': {'pathogenic': 0.5, 'benign': 0.5},
    'lrt_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},10},
    # Functional predictorsnign': 0.34},
    'fitcons': {'pathogenic': 0.7, 'benign': 0.3}
}thogenic': 0.5, 'benign': 0.5},
hogenic': 0.5, 'benign': 0.5},
# In silico predictor thresholds{'pathogenic': 0.16, 'benign': -0.18},
INSILICO_THRESHOLDS = {   'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},
    'revel': {'pathogenic': 0.75, 'benign': 0.25},    'sift': {'pathogenic': 0.05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)
    'cadd_phred': {'pathogenic': 25, 'benign': 10},c': 0.85, 'benign': 0.15},
    'alphamissense': {'pathogenic': 0.564, 'benign': 0.34},': {'pathogenic': 0.5, 'benign': 0.5},
    'metarnn': {'pathogenic': 0.5, 'benign': 0.5},': {'pathogenic': 0.5, 'benign': 0.5},
    'metarnn_ranked': {'pathogenic': 0.5, 'benign': 0.5},genic': -1.5, 'benign': 1.5},  # FATHMM is inverted (lower = more pathogenic)
    'clinpred': {'pathogenic': 0.5, 'benign': 0.5},: 0.5, 'benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},lds
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},: 0.7, 'benign': 0.3},
    'sift': {'pathogenic': 0.05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)': 0.2},
    'polyphen2': {'pathogenic': 0.85, 'benign': 0.15},c': 0.5, 'benign': 0.5},
    'mutationtaster': {'pathogenic': 0.5, 'benign': 0.5}, -2.5, 'benign': -2.5},  # PROVEAN is inverted (lower = more pathogenic)
    'mutationtaster_ranked': {'pathogenic': 0.5, 'benign': 0.5}, 'benign': 2.0},
    'fathmm': {'pathogenic': -1.5, 'benign': 1.5},  # FATHMM is inverted (lower = more pathogenic)'benign': 0.5},
    'fathmm_ranked': {'pathogenic': 0.5, 'benign': 0.5},thogenic': 0.5, 'benign': 0.5},
    # High priority predictor thresholds
    'vest4': {'pathogenic': 0.7, 'benign': 0.3}, 0.5, 'benign': 0.5},
    'primateai': {'pathogenic': 0.8, 'benign': 0.2},   'phylop_prim': {'pathogenic': 2.0, 'benign': 0.5},
    'esm1b': {'pathogenic': 0.5, 'benign': 0.5},    'phylop_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'provean': {'pathogenic': -2.5, 'benign': -2.5},  # PROVEAN is inverted (lower = more pathogenic)splice-altering variants)
    'gerp_pp': {'pathogenic': 5.0, 'benign': 2.0},thogenic': 0.2, 'benign': 0.05},
    'phylop_vert': {'pathogenic': 2.0, 'benign': 0.5},benign': 0.05},
    'phylop_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phylop_mamm': {'pathogenic': 2.0, 'benign': 0.5},ogenic': 0.2, 'benign': 0.05},
    'phylop_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},genic': 0.5, 'benign': 0.1},
    'phylop_prim': {'pathogenic': 2.0, 'benign': 0.5},
    'phylop_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    # SpliceAI thresholds (splice-altering variants)
    'spliceai_ag': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_al': {'pathogenic': 0.2, 'benign': 0.05},   'phastcons_mamm': 0.02,
    'spliceai_dg': {'pathogenic': 0.2, 'benign': 0.05},    'phastcons_prim': 0.015,
    'spliceai_dl': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_max': {'pathogenic': 0.5, 'benign': 0.1},ons_mamm_ranked': 0.02,
    # Enhanced splice predictorsnked': 0.015,
    'mmsplice': {'pathogenic': 2.0, 'benign': 0.5},
    # Conservation predictors
    'phastcons_vert': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_mamm': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_prim': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phastcons_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phastcons_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'siphy': {'pathogenic': 12.0, 'benign': 3.0},
    # Additional missense predictorsrt': 0.03,
    'metasvm': {'pathogenic': 0.83, 'benign': 0.17},
    'metasvm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'metalr': {'pathogenic': 0.5, 'benign': 0.5},
    'metalr_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'mutationassessor': {'pathogenic': 3.5, 'benign': 1.5},
    'mutpred': {'pathogenic': 0.7, 'benign': 0.3},
    'mutpred_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'lrt': {'pathogenic': 0.5, 'benign': 0.5}, In silico predictor thresholds
    'lrt_ranked': {'pathogenic': 0.5, 'benign': 0.5},INSILICO_THRESHOLDS = {
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},ic': 25, 'benign': 10},
    # Functional predictors
    'fitcons': {'pathogenic': 0.7, 'benign': 0.3}
}

# In silico predictor thresholds
INSILICO_THRESHOLDS = {
    'revel': {'pathogenic': 0.75, 'benign': 0.25},   'sift': {'pathogenic': 0.05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)
    'cadd_phred': {'pathogenic': 25, 'benign': 10},    'polyphen2': {'pathogenic': 0.85, 'benign': 0.15},
    'alphamissense': {'pathogenic': 0.564, 'benign': 0.34},enign': 0.5},
    'metarnn': {'pathogenic': 0.5, 'benign': 0.5},, 'benign': 0.5},
    'metarnn_ranked': {'pathogenic': 0.5, 'benign': 0.5},    'fathmm': {'pathogenic': -1.5, 'benign': 1.5},  # FATHMM is inverted (lower = more pathogenic)
    'clinpred': {'pathogenic': 0.5, 'benign': 0.5},benign': 0.5},
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},y predictor thresholds
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},
    'sift': {'pathogenic': 0.05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)
    'polyphen2': {'pathogenic': 0.85, 'benign': 0.15},
    'mutationtaster': {'pathogenic': 0.5, 'benign': 0.5},  # PROVEAN is inverted (lower = more pathogenic)
    'mutationtaster_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'fathmm': {'pathogenic': -1.5, 'benign': 1.5},  # FATHMM is inverted (lower = more pathogenic)
    'fathmm_ranked': {'pathogenic': 0.5, 'benign': 0.5},gn': 0.5},
    # High priority predictor thresholds
    'vest4': {'pathogenic': 0.7, 'benign': 0.3},benign': 0.5},
    'primateai': {'pathogenic': 0.8, 'benign': 0.2},enign': 0.5},
    'esm1b': {'pathogenic': 0.5, 'benign': 0.5},.5, 'benign': 0.5},
    'provean': {'pathogenic': -2.5, 'benign': -2.5},  # PROVEAN is inverted (lower = more pathogenic))
    'gerp_pp': {'pathogenic': 5.0, 'benign': 2.0},5},
    'phylop_vert': {'pathogenic': 2.0, 'benign': 0.5},   'spliceai_al': {'pathogenic': 0.2, 'benign': 0.05},
    'phylop_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},    'spliceai_dg': {'pathogenic': 0.2, 'benign': 0.05},
    'phylop_mamm': {'pathogenic': 2.0, 'benign': 0.5}, {'pathogenic': 0.2, 'benign': 0.05},
    'phylop_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},{'pathogenic': 0.5, 'benign': 0.1},
    'phylop_prim': {'pathogenic': 2.0, 'benign': 0.5},
    'phylop_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    # SpliceAI thresholds (splice-altering variants)
    'spliceai_ag': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_al': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_dg': {'pathogenic': 0.2, 'benign': 0.05},
    'spliceai_dl': {'pathogenic': 0.2, 'benign': 0.05},   'phastcons_vert_ranked': 0.02,
    'spliceai_max': {'pathogenic': 0.5, 'benign': 0.1},    'phastcons_mamm_ranked': 0.02,
    # Enhanced splice predictorsm_ranked': 0.015,
    'mmsplice': {'pathogenic': 2.0, 'benign': 0.5},
    # Conservation predictors
    'phastcons_vert': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_mamm': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_prim': {'pathogenic': 0.9, 'benign': 0.3},
    'phastcons_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phastcons_mamm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'phastcons_prim_ranked': {'pathogenic': 0.5, 'benign': 0.5},   'mutpred': 0.04,
    'siphy': {'pathogenic': 12.0, 'benign': 3.0},    'mutpred_ranked': 0.04,
    # Additional missense predictors,
    'metasvm': {'pathogenic': 0.83, 'benign': 0.17}, 0.03,
    'metasvm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'metalr': {'pathogenic': 0.5, 'benign': 0.5},
    'metalr_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'mutationassessor': {'pathogenic': 3.5, 'benign': 1.5},
    'mutpred': {'pathogenic': 0.7, 'benign': 0.3},
    'mutpred_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'lrt': {'pathogenic': 0.5, 'benign': 0.5}, In silico predictor thresholds
    'lrt_ranked': {'pathogenic': 0.5, 'benign': 0.5},INSILICO_THRESHOLDS = {
    'bayesdel_addaf': {'pathogenic': 0.16, 'benign': -0.18},0.75, 'benign': 0.25},
    'bayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},thogenic': 25, 'benign': 10},
    # Functional predictorsgenic': 0.564, 'benign': 0.34},
    'fitcons': {'pathogenic': 0.7, 'benign': 0.3}
}

# Gene-specific predictor preferences (original definition restored).18},
GENE_SPECIFIC_PREDICTORS = {ayesdel_noaf': {'pathogenic': 0.16, 'benign': -0.18},
    'BRCA1': {05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)
        'preferred_predictors': ['vest4', 'revel', 'clinpred'],  # Cancer-specific
        'weight_boost': 1.3
    },,
    'BRCA2': {FATHMM is inverted (lower = more pathogenic)
        'preferred_predictors': ['vest4', 'revel', 'clinpred'],athmm_ranked': {'pathogenic': 0.5, 'benign': 0.5},
        'weight_boost': 1.3y predictor thresholds
    },
    'TP53': {
        'preferred_predictors': ['vest4', 'alphamissense', 'esm1b'],
        'weight_boost': 1.2N is inverted (lower = more pathogenic)
    },gerp_pp': {'pathogenic': 5.0, 'benign': 2.0},
    'CFTR': {   'phylop_vert': {'pathogenic': 2.0, 'benign': 0.5},
        'preferred_predictors': ['cadd_phred', 'provean', 'sift'],    'phylop_vert_ranked': {'pathogenic': 0.5, 'benign': 0.5},
        'weight_boost': 1.1genic': 2.0, 'benign': 0.5},
    },d': {'pathogenic': 0.5, 'benign': 0.5},
    'TTN': {benign': 0.5},
        # Large gene - different predictor performance0.5, 'benign': 0.5},
        'preferred_predictors': ['revel', 'cadd_phred'],ring variants)
        'weight_boost': 1.0,enic': 0.2, 'benign': 0.05},
        'frequency_adjustment': 0.1  # Higher frequency toleranceic': 0.2, 'benign': 0.05},
    },   'spliceai_dg': {'pathogenic': 0.2, 'benign': 0.05},
    'CAMTA1': {    'spliceai_dl': {'pathogenic': 0.2, 'benign': 0.05},
        # Neurodevelopmental disorder gene
        'preferred_predictors': ['cadd_phred', 'revel', 'clinpred'],
        'weight_boost': 1.2,splice': {'pathogenic': 2.0, 'benign': 0.5},
        'frequency_adjustment': 0.0  # Standard frequency requirements
    }'phastcons_vert': 0.02,
}tcons_mamm': 0.02,

# Classification result colors
CLASSIFICATION_COLORS = {
    'Pathogenic': ('RED', 'BOLD'),
    'Likely Pathogenic': ('RED', None),'siphy': 0.015,
    'Likely Benign': ('GREEN', None),onal missense predictors
    'Benign': ('GREEN', 'BOLD'),
    'Variant of Uncertain Significance': ('YELLOW', 'BOLD'),tasvm_ranked': 0.09,
    'VUS': ('YELLOW', 'BOLD')'metalr': 0.09,
}

# Colorama color mappings for console output
COLORAMA_COLORS = {'mutpred_ranked': 0.04,
    'RED': '\033[91m',
    'GREEN': '\033[92m',,
    'YELLOW': '\033[93m',0.06,
    'BLUE': '\033[94m',
    'MAGENTA': '\033[95m',ors
    'CYAN': '\033[96m','fitcons': 0.05
    'WHITE': '\033[97m',
    'BOLD': '\033[1m',
    'UNDERLINE': '\033[4m',
    'END': '\033[0m',LDS = {
    'RESET': '\033[0m'  # Alias for END
}
se': {'pathogenic': 0.564, 'benign': 0.34},
# Section headers with formatted text and color0.5},
SECTION_HEADERS = { 0.5},
    'basic_info': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Variant Information{COLORAMA_COLORS['RESET']}",{'pathogenic': 0.5, 'benign': 0.5},
    'variant_info': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Variant Information{COLORAMA_COLORS['RESET']}",
    'population_data': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Population Frequency Data{COLORAMA_COLORS['RESET']}",
    'insilico_data': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}In Silico Prediction Scores{COLORAMA_COLORS['RESET']}",is inverted (lower = more pathogenic)
    'genetic_data': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Genetic Data{COLORAMA_COLORS['RESET']}",
    'functional_data': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Functional Data{COLORAMA_COLORS['RESET']}",
    'evidence_evaluation': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Evidence Evaluation{COLORAMA_COLORS['RESET']}",ter_ranked': {'pathogenic': 0.5, 'benign': 0.5},
    'classification': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Variant Classification{COLORAMA_COLORS['RESET']}",,  # FATHMM is inverted (lower = more pathogenic)
    'report_generation': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Report Generation{COLORAMA_COLORS['RESET']}",
    'analysis_complete': f"{COLORAMA_COLORS['GREEN']}{COLORAMA_COLORS['BOLD']}Analysis Complete{COLORAMA_COLORS['RESET']}"
}gn': 0.3},
 {'pathogenic': 0.8, 'benign': 0.2},
# Utility functions for colored output5},
from typing import Optional.5},  # PROVEAN is inverted (lower = more pathogenic)
ign': 2.0},
# Utility functions for colored output
def get_colored_message(message: str, color: str = 'WHITE', style: Optional[str] = None) -> str:c': 0.5, 'benign': 0.5},
    """0.5},
    Apply color and style to a message for console output.
    : {'pathogenic': 2.0, 'benign': 0.5},
    Args:nign': 0.5},
        message (str): The message to colorhresholds (splice-altering variants)
        color (str): Color name from COLORAMA_COLORS0.2, 'benign': 0.05},
        style (str): Style name from COLORAMA_COLORS (BOLD, UNDERLINE, etc.)
     0.05},
    Returns:
        str: Formatted message with color codes, 'benign': 0.1},
    """lice predictors
    colored_msg = COLORAMA_COLORS.get(color, '') + message'mmsplice': {'pathogenic': 2.0, 'benign': 0.5},
    
    if style and style in COLORAMA_COLORS:02,
        colored_msg = COLORAMA_COLORS[style] + colored_msg
    tcons_prim': 0.015,
    return colored_msg + COLORAMA_COLORS['END']d': 0.02,
'phastcons_mamm_ranked': 0.02,
# Test mode data for development and testing
TEST_MODE_DATA = {
    'basic_info': {
        'gene': 'CAMTA1',
        'variant_type': 'missense','metasvm_ranked': 0.09,
        'chromosome': '1',
        'position': '7249507',  # Real genomic position9,
        'ref_allele': 'A',r': 0.06,
        'alt_allele': 'G',,
        'amino_acid_change': 'p.Met107Val','mutpred_ranked': 0.04,
        'cdna_change': 'c.319A>G',
        'consequence': 'missense_variant'
    },
    'population_data': {
        'gnomad_af': 1.5e-5,  # Real gnomAD frequency 1.50e-05
        'gnomad_af_afr': 0.0,
        'gnomad_af_amr': 0.0,
        'gnomad_af_asj': 0.0,
        'gnomad_af_eas': 0.0,
        'gnomad_af_fin': 0.0,
        'gnomad_af_nfe': 0.0,'revel': {'pathogenic': 0.75, 'benign': 0.25},
        'gnomad_af_oth': 0.0
    },'alphamissense': {'pathogenic': 0.564, 'benign': 0.34},
    'insilico_data': {5, 'benign': 0.5},
        'revel_score': 0.395,     # Real REVEL score
        'cadd_phred': 24.3,       # Real CADD score
        'gerp_pp': 5.63,          # Real GERP++ scorehogenic': 0.16, 'benign': -0.18},
        'phylop_vert_ranked': 0.93579,  # Real PhyloP vertebrate ranked8},
        'phylop_mamm_ranked': 0.94714,  # Real PhyloP mammalian ranked05, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)
        'phylop_prim_ranked': 0.94297,  # Real PhyloP primate rankedphen2': {'pathogenic': 0.85, 'benign': 0.15},
        'alphamissense_score': 0.3896,  # Real AlphaMissense scorethogenic': 0.5, 'benign': 0.5},
        'clinpred_score': 0.349,         # Real ClinPred score'mutationtaster_ranked': {'pathogenic': 0.5, 'benign': 0.5},
        'bayesdel_noaf': 0.72752,       # Real BayesDel NoAF scoreinverted (lower = more pathogenic)
        'sift_score': 0.041,            # Real SIFT score
        'polyphen2_hdiv_score': 0.68658, # Real PolyPhen2 HDIV score# High priority predictor thresholds
        'mutationtaster_ranked': 0.51042, # Real MutationTaster ranked {'pathogenic': 0.7, 'benign': 0.3},
        'fathmm_ranked': 0.20382,       # Real FATHMM ranked, 'benign': 0.2},
        'vest4_score': 0.699,           # Real VEST4 score
        'esm1b': 0.62479,              # Real ESM1b scorebenign': -2.5},  # PROVEAN is inverted (lower = more pathogenic)
        'provean': -1.25,              # Real PROVEAN score,
        'metasvm': 0.03192,            # Real MetaSVM score,
        'metalr': 0.0946,              # Real MetaLR scoreenign': 0.5},
        'mutationassessor': 0.12,      # Real MutationAssessor scoreign': 0.5},
        'mutpred': 0.559,              # Real MutPred scorebenign': 0.5},
        'lrt': 0.000469                # Real LRT scoreathogenic': 2.0, 'benign': 0.5},
    },nign': 0.5},
    'genetic_data': {g variants)
        'inheritance': 'AD',  # CAMTA1 is autosomal dominant'benign': 0.05},
        'segregation': 'not_applicable',  # Only proband affected, no family data0.2, 'benign': 0.05},
        'de_novo': 'assumed',  # Real data: OMIM suggests de novo but not confirmed with parental testingeai_dg': {'pathogenic': 0.2, 'benign': 0.05},
        'maternity_confirmed': False,  # Real data: parental testing not performed 0.2, 'benign': 0.05},
        'paternity_confirmed': False,  # Real data: parental testing not performedbenign': 0.1},
        'compound_het': 'no'
    },
    'functional_data': {
        'functional_studies': 'inconclusive',  # Real data: inconclusive
        'phenotype_match': 'partial_match',    # Real data: partial matchphastcons_mamm': 0.02,
        'hpo_similarity_score': 0.75,         # Test HPO similarity score for PP4    'phastcons_prim': 0.015,
        'case_control': 'no'
    }
},
5,
# API configuration and endpoints
API_ENDPOINTS = {
    'ensembl_rest': 'https://rest.ensembl.org',
    'clinvar': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils',
    'dbsnp': 'https://api.ncbi.nlm.nih.gov/variation/v0',
    'gnomad': 'https://gnomad.broadinstitute.org/api',: 0.06,
    'lovd': 'https://databases.lovd.nl/shared',
    'pharmgkb': 'https://api.pharmgkb.org/v1',0.04,
    'varsome': 'https://varsome.com/variant/hg38',
    'timeout_seconds': 30,3,
    'max_retries': 3,6,
    'retry_delay': 1.06,
}tors

# Output formatting settings
OUTPUT_SETTINGS = {
    'max_line_length': 80,
    'indent_size': 4,
    'section_separator': '=' * 80,0.75, 'benign': 0.25},
    'subsection_separator': '-' * 40,ic': 25, 'benign': 10},
    'enable_colors': True,genic': 0.564, 'benign': 0.34},
    'timestamp_format': '%Y-%m-%d %H:%M:%S',.5, 'benign': 0.5},
    'decimal_places': 3,nic': 0.5, 'benign': 0.5},
    'score_format': '{:.3f}',0.5, 'benign': 0.5},
    'percentage_format': '{:.1%}',nic': 0.16, 'benign': -0.18},
    'cache_filename': 'api_cache.json',athogenic': 0.16, 'benign': -0.18},
    'max_cache_age_hours': 24,, 'benign': 0.05},  # Note: SIFT is inverted (lower = more pathogenic)
    'report_filename': 'variant_classification_report.txt',: 0.85, 'benign': 0.15},
    'log_filename': 'acmg_assistant.log'pathogenic': 0.5, 'benign': 0.5},
}pathogenic': 0.5, 'benign': 0.5},
athmm': {'pathogenic': -1.5, 'benign': 1.5},  # FATHMM is inverted (lower = more pathogenic)
# Input validation patternsd': {'pathogenic': 0.5, 'benign': 0.5},
VALIDATION_PATTERNS = {
    'gene_symbol': r'^[A-Z][A-Z0-9-]*$', 0.7, 'benign': 0.3},
    'chromosome': r'^(chr)?(1[0-9]|2[0-2]|[1-9]|X|Y|MT?)$',genic': 0.8, 'benign': 0.2},
    'position': r'^\d+$',0.5, 'benign': 0.5},
    'allele': r'^[ATCG]+$',: -2.5, 'benign': -2.5},  # PROVEAN is inverted (lower = more pathogenic)
    'amino_acid': r'^p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}$',: 5.0, 'benign': 2.0},
    'hgvs_cdna': r'^c\.\d+[ATCG]>[ATCG]$',  # Added cDNA validation pattern': 2.0, 'benign': 0.5},
    'frequency': r'^(0(\.\d+)?|1(\.0+)?)$',hogenic': 0.5, 'benign': 0.5},
    'score': r'^-?\d+(\.\d+)?$'': 2.0, 'benign': 0.5},
}

# Gene and variant aliases for user input flexibilityd': {'pathogenic': 0.5, 'benign': 0.5},
ALIASES = {ice-altering variants)
    'variant_types': {hogenic': 0.2, 'benign': 0.05},
        'missense': ['missense', 'missense_variant', 'amino_acid_substitution'],pliceai_al': {'pathogenic': 0.2, 'benign': 0.05},
        'nonsense': ['nonsense', 'stop_gained', 'premature_stop'],{'pathogenic': 0.2, 'benign': 0.05},
        'frameshift': ['frameshift', 'frameshift_variant', 'indel'],
        'splice_donor': ['splice_donor', 'splice_donor_variant', 'donor_splice'],genic': 0.5, 'benign': 0.1},
        'splice_acceptor': ['splice_acceptor', 'splice_acceptor_variant', 'acceptor_splice'],edictors
        'synonymous': ['synonymous', 'synonymous_variant', 'silent'],': 2.0, 'benign': 0.5},
        'intronic': ['intronic', 'intron_variant'],s
        'inframe_deletion': ['inframe_deletion', 'inframe_del'],
        'inframe_insertion': ['inframe_insertion', 'inframe_ins']
    },
    'inheritance_patterns': {02,
        'AD': ['AD', 'autosomal_dominant', 'dominant'],0.02,
        'AR': ['AR', 'autosomal_recessive', 'recessive'],nked': 0.015,
        'XLD': ['XLD', 'x_linked_dominant', 'x_dominant'],
        'XLR': ['XLR', 'x_linked_recessive', 'x_recessive'],e predictors
        'mitochondrial': ['mitochondrial', 'maternal', 'mtDNA']etasvm': 0.09,
    } 0.09,
}

# All possible variant consequences for validation,
ALL_VARIANT_CONSEQUENCES = [
    'missense_variant', 'nonsense_variant', 'stop_gained', 'frameshift_variant',4,
    'splice_donor_variant', 'splice_acceptor_variant', 'splice_region_variant',
    'synonymous_variant', 'intron_variant', 'regulatory_region_variant',
    'upstream_gene_variant', 'downstream_gene_variant', '5_prime_UTR_variant',
    '3_prime_UTR_variant', 'inframe_deletion', 'inframe_insertion',
    'start_lost', 'stop_lost', 'protein_altering_variant', 'coding_sequence_variant'
]

# Alias for compatibility with validators.py
VARIANT_CONSEQUENCES = ALL_VARIANT_CONSEQUENCESeferences (original definition restored)

# Error messages for user input validation
ERROR_MESSAGES = {: ['vest4', 'revel', 'clinpred'],  # Cancer-specific
    'invalid_numeric': "❌ Invalid input. Please enter a valid number.",  'weight_boost': 1.3
    'invalid_choice': "❌ Invalid choice. Please select from the available options.",
    'empty_input': "❌ Input cannot be empty. Please provide a value.",
    'invalid_gene': "❌ Invalid gene symbol format.",s': ['vest4', 'revel', 'clinpred'],
    'invalid_chromosome': "❌ Invalid chromosome format.",
    'invalid_position': "❌ Invalid genomic position.",
    'invalid_allele': "❌ Invalid allele sequence.",
    'invalid_frequency': "❌ Invalid frequency value (must be between 0 and 1).",ors': ['vest4', 'alphamissense', 'esm1b'],
    'invalid_score': "❌ Invalid score format.",.2
    'api_error': "❌ API request failed.",
    'file_not_found': "❌ File not found.",
    'parsing_error': "❌ Error parsing input data.",rs': ['cadd_phred', 'provean', 'sift'],
    'validation_failed': "❌ Data validation failed."1.1
}

# Success messagesrent predictor performance
SUCCESS_MESSAGES = { ['revel', 'cadd_phred'],
    'data_saved': "✅ Data saved successfully.",
    'analysis_complete': "✅ Variant analysis completed successfully.",: 0.1  # Higher frequency tolerance
    'report_generated': "✅ Report generated successfully.",
    'api_success': "✅ API request completed successfully.",
    'file_loaded': "✅ File loaded successfully.",
    'validation_passed': "✅ Data validation passed."': ['cadd_phred', 'revel', 'clinpred'],
}.2,
ent': 0.0  # Standard frequency requirements
# Warning messages
WARNING_MESSAGES = {
    'api_timeout': "⚠️ API request timed out. Using cached data if available.",
    'low_confidence': "⚠️ Low confidence in prediction results.",rs
    'missing_data': "⚠️ Some data points are missing.",
    'frequency_warning': "⚠️ Population frequency data may be incomplete.",OLD'),
    'predictor_warning': "⚠️ Limited in silico predictor data available.",, None),
    'gene_warning': "⚠️ Limited gene-specific data available."None),
},
D'),
# Info messages'BOLD')
INFO_MESSAGES = {
    'starting_analysis': "ℹ️ Starting variant analysis...",
    'loading_data': "ℹ️ Loading variant data...",rama color mappings for console output
    'calculating_scores': "ℹ️ Calculating predictor scores...", {
    'applying_criteria': "ℹ️ Applying ACMG/AMP criteria...",
    'generating_report': "ℹ️ Generating classification report...",
    'saving_results': "ℹ️ Saving analysis results..."
}
,
# AI Feature Status Messages
AI_STATUS_MESSAGES = {
    'literature_scanning': {
        'enabled': '🤖 AI Literature Scanning: ENABLED',m',
        'disabled': '📚 AI Literature Scanning: DISABLED',
        'analyzing': '🔍 AI Literature Analysis in progress...',ias for END
        'complete': '✅ AI Literature Analysis complete'
    },
    'phenotype_similarity': {d text and color
        'enabled': '🧬 AI Phenotype Matching: ENABLED',
        'disabled': '🔬 AI Phenotype Matching: DISABLED',basic_info': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Variant Information{COLORAMA_COLORS['RESET']}",
        'analyzing': '🔍 AI Phenotype Analysis in progress...',   'variant_info': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Variant Information{COLORAMA_COLORS['RESET']}",
        'complete': '✅ AI Phenotype Analysis complete'    'population_data': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Population Frequency Data{COLORAMA_COLORS['RESET']}",
    },AMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}In Silico Prediction Scores{COLORAMA_COLORS['RESET']}",
    'metascore': {LORS['RESET']}",
        'enabled': '⚡ Dynamic Metascore: ENABLED',nctional_data': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Functional Data{COLORAMA_COLORS['RESET']}",
        'disabled': '📊 Dynamic Metascore: DISABLED',BOLD']}Evidence Evaluation{COLORAMA_COLORS['RESET']}",
        'calculating': '🔍 Calculating Dynamic Metascore...','classification': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Variant Classification{COLORAMA_COLORS['RESET']}",
        'complete': '✅ Dynamic Metascore calculation complete'rt_generation': f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Report Generation{COLORAMA_COLORS['RESET']}",
    }A_COLORS['BOLD']}Analysis Complete{COLORAMA_COLORS['RESET']}"
}

# AI Feature Configurationility functions for colored output
AI_FEATURES_CONFIG = {import Optional
    'literature_scanning_enabled': True,
    'phenotype_similarity_enabled': True,ty functions for colored output
    'dynamic_metascore_enabled': True,d_message(message: str, color: str = 'WHITE', style: Optional[str] = None) -> str:
    'show_ai_status': True,
    'verbose_ai_output': Truetyle to a message for console output.
}

# Enhanced Metascore calculation with dynamic weighting   message (str): The message to color
def calculate_dynamic_vampp_score(insilico_scores, variant_type, allele_frequency=None, gene=None):    color (str): Color name from COLORAMA_COLORS
    """COLORAMA_COLORS (BOLD, UNDERLINE, etc.)
    Calculate Metascore with dynamic predictor weighting based on variant type.
    
    Args:ith color codes
        insilico_scores (dict): Available predictor scores
        variant_type (str): Type of variant (missense, nonsense, etc.)ge
        allele_frequency (float): Population allele frequency
        gene (str): Gene symbol for gene-specific adjustmentsMA_COLORS:
    
    Returns:
        dict: Enhanced Metascore results with detailed breakdown
    """
    
    # Get variant-specific weights, fallback to missense if not found_MODE_DATA = {
    variant_weights = VARIANT_TYPE_PREDICTOR_WEIGHTS.get(variant_type, 
                                                        VARIANT_TYPE_PREDICTOR_WEIGHTS['missense'])
    sense',
    # Initialize score calculation
    weighted_sum = 0.0    'position': '7249507',  # Real genomic position
    total_weight = 0.0
    predictor_details = {}
    used_predictors = []    'amino_acid_change': 'p.Met107Val',
    nge': 'c.319A>G',
    # Process each available predictor        'consequence': 'missense_variant'
    for predictor, score in insilico_scores.items():
        if predictor in variant_weights and score is not None:
            try: Real gnomAD frequency 1.50e-05
                # Normalize score to 0-1 range if needed.0,
                normalized_score = normalize_predictor_score(predictor, score)
                
                # Get variant-specific weight
                weight = variant_weights[predictor]  'gnomad_af_fin': 0.0,
                
                # Apply gene-specific adjustments if available.0
                if gene and gene in GENE_SPECIFIC_PREDICTORS:
                    gene_config = GENE_SPECIFIC_PREDICTORS[gene]
                    if predictor in gene_config.get('preferred_predictors', []):
                        weight *= gene_config.get('weight_boost', 1.0)  'cadd_phred': 24.3,       # Real CADD score
                         # Real GERP++ score
                # Calculate weighted contributiond': 0.93579,  # Real PhyloP vertebrate ranked
                contribution = normalized_score * weight.94714,  # Real PhyloP mammalian ranked
                weighted_sum += contributionimate ranked
                total_weight += weight
                  'clinpred_score': 0.349,         # Real ClinPred score
                # Store details for reporting752,       # Real BayesDel NoAF score
                predictor_details[predictor] = {41,            # Real SIFT score
                    'raw_score': score,core': 0.68658, # Real PolyPhen2 HDIV score
                    'normalized_score': normalized_score,
                    'weight': weight,
                    'contribution': contribution,   'vest4_score': 0.699,           # Real VEST4 score
                    'pathogenic_threshold': INSILICO_THRESHOLDS.get(predictor, {}).get('pathogenic', 0.5)       'esm1b': 0.62479,              # Real ESM1b score
                }        'provean': -1.25,              # Real PROVEAN score
                used_predictors.append(predictor)        # Real MetaSVM score
                ore
            except Exception as e:,      # Real MutationAssessor score
                # Log error but continue with other predictors: 0.559,              # Real MutPred score
                predictor_details[predictor] = {   # Real LRT score
                    'error': f"Processing error: {str(e)}",
                    'raw_score': score
                }nheritance': 'AD',  # CAMTA1 is autosomal dominant
    ble',  # Only proband affected, no family data
    # Calculate final Metascoremed',  # Real data: OMIM suggests de novo but not confirmed with parental testing
    if total_weight > 0:irmed': False,  # Real data: parental testing not performed
        vampp_score = weighted_sum / total_weightaternity_confirmed': False,  # Real data: parental testing not performed
    else:
        vampp_score = 0.0
    {
    # Apply frequency-based adjustmentsunctional_studies': 'inconclusive',  # Real data: inconclusive
    frequency_category = get_frequency_category(allele_frequency)ial_match',    # Real data: partial match
    frequency_adjustment = get_frequency_adjustment(variant_type, frequency_category)score': 0.75,         # Test HPO similarity score for PP4
    adjusted_vampp_score = max(0.0, min(1.0, vampp_score + frequency_adjustment)) 'no'
    
    # Calculate predictor consensus
    pathogenic_votes = 0
    benign_votes = 0 endpoints
    total_votes = 0
    tps://rest.ensembl.org',
    for predictor in used_predictors:ar': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils',
        if predictor in predictor_details:api.ncbi.nlm.nih.gov/variation/v0',
            details = predictor_details[predictor]nomad.broadinstitute.org/api',
            if 'normalized_score' in details:tabases.lovd.nl/shared',
                threshold = details['pathogenic_threshold']gkb': 'https://api.pharmgkb.org/v1',
                if details['normalized_score'] >= threshold:rsome.com/variant/hg38',
                    pathogenic_votes += 10,
                else:
                    benign_votes += 1_delay': 1.0
                total_votes += 1
    
    consensus = pathogenic_votes / total_votes if total_votes > 0 else 0.0tings
    TTINGS = {
    # Determine confidence levelax_line_length': 80,
    min_predictors = VARIANT_TYPE_MIN_PREDICTORS.get(variant_type, 5)4,
    if len(used_predictors) >= min_predictors and consensus >= 0.7:: '=' * 80,
        confidence = 'High'r': '-' * 40,
    elif len(used_predictors) >= 3 and consensus >= 0.6:ue,
        confidence = 'Medium'tamp_format': '%Y-%m-%d %H:%M:%S',
    else:3,
        confidence = 'Low'3f}',
    : '{:.1%}',
    # Get dynamic thresholds for this variant type and frequency_filename': 'api_cache.json',
    dynamic_thresholds = get_dynamic_thresholds(variant_type, frequency_category) 24,
    variant_classification_report.txt',
    return {mg_assistant.log'
        'vampp_score': vampp_score,
        'adjusted_vampp_score': adjusted_vampp_score,
        'variant_type': variant_type,rns
        'frequency_category': frequency_category,
        'frequency_adjustment': frequency_adjustment,_symbol': r'^[A-Z][A-Z0-9-]*$',
        'predictor_count': len(used_predictors),hromosome': r'^(chr)?(1[0-9]|2[0-2]|[1-9]|X|Y|MT?)$',
        'used_predictors': used_predictors,+$',
        'predictor_details': predictor_details,+$',
        'consensus': {[A-Z][a-z]{2}\d+[A-Z][a-z]{2}$',
            'pathogenic_votes': pathogenic_votes,\d+[ATCG]>[ATCG]$',  # Added cDNA validation pattern
            'benign_votes': benign_votes,ency': r'^(0(\.\d+)?|1(\.0+)?)$',
            'total_votes': total_votes,\.\d+)?$'
            'agreement': consensus
        },
        'confidence': confidence, variant aliases for user input flexibility
        'thresholds': dynamic_thresholds,
        'pp3_applicable': adjusted_vampp_score >= dynamic_thresholds['pp3'],
        'bp4_applicable': adjusted_vampp_score <= dynamic_thresholds['bp4'],issense', 'missense_variant', 'amino_acid_substitution'],
        'total_weight_used': total_weight,onsense': ['nonsense', 'stop_gained', 'premature_stop'],
        'gene_specific_adjustments': gene in GENE_SPECIFIC_PREDICTORS if gene else False': ['frameshift', 'frameshift_variant', 'indel'],
    }['splice_donor', 'splice_donor_variant', 'donor_splice'],
r': ['splice_acceptor', 'splice_acceptor_variant', 'acceptor_splice'],
# 1. VARIANT TYPE DYNAMIC WEIGHTINGsynonymous': ['synonymous', 'synonymous_variant', 'silent'],
# Based on ROC-AUC performance of predictors for different variant types  'intronic': ['intronic', 'intron_variant'],
VARIANT_TYPE_PREDICTOR_WEIGHTS = {n': ['inframe_deletion', 'inframe_del'],
    'missense': {ion': ['inframe_insertion', 'inframe_ins']
        # High-performance predictors for missense variants
        'revel': 0.25,ns': {
        'alphamissense': 0.20,D': ['AD', 'autosomal_dominant', 'dominant'],
        'cadd_phred': 0.15,autosomal_recessive', 'recessive'],
        'primateai': 0.12,x_linked_dominant', 'x_dominant'],
        'esm1b': 0.10,'x_linked_recessive', 'x_recessive'],
        'clinpred': 0.08,itochondrial': ['mitochondrial', 'maternal', 'mtDNA']
        'vest4': 0.05,
        'polyphen2': 0.03,
        'sift': 0.02,
        # Meta-predictorsible variant consequences for validation
        'metasvm': 0.04,ENCES = [
        'metalr': 0.04,'nonsense_variant', 'stop_gained', 'frameshift_variant',
        'mutationassessor': 0.03,nt', 'splice_acceptor_variant', 'splice_region_variant',
        'mutpred': 0.02,nymous_variant', 'intron_variant', 'regulatory_region_variant',
        'lrt': 0.02,pstream_gene_variant', 'downstream_gene_variant', '5_prime_UTR_variant',
        # Conservation scores (lower weight for missense)riant', 'inframe_deletion', 'inframe_insertion',
        'gerp_pp': 0.015,_lost', 'protein_altering_variant', 'coding_sequence_variant'
        'phylop_vert': 0.01,
        'phylop_mamm': 0.01,
        'phylop_prim': 0.01,r compatibility with validators.py
        'phastcons_vert': 0.01, ALL_VARIANT_CONSEQUENCES
        'phastcons_mamm': 0.01,
        'phastcons_prim': 0.01,er input validation
        # Functional predictionAGES = {
        'provean': 0.02,Invalid input. Please enter a valid number.",
        'bayesdel_addaf': 0.02, Invalid choice. Please select from the available options.",
        'bayesdel_noaf': 0.02,nput cannot be empty. Please provide a value.",
        'fathmm': 0.015,id_gene': "❌ Invalid gene symbol format.",
        'mutationtaster': 0.015some': "❌ Invalid chromosome format.",
    },"❌ Invalid genomic position.",
    'nonsense': {❌ Invalid allele sequence.",
        # For nonsense variants, conservation and splice predictors are more importantlid_frequency': "❌ Invalid frequency value (must be between 0 and 1).",
        'cadd_phred': 0.30,nvalid_score': "❌ Invalid score format.",
        'gerp_pp': 0.25,"❌ API request failed.",
        'phylop_vert': 0.15,❌ File not found.",
        'phylop_mamm': 0.10,Error parsing input data.",
        'phylop_prim': 0.08,: "❌ Data validation failed."
        'phastcons_vert': 0.05,
        'phastcons_mamm': 0.03,
        'phastcons_prim': 0.02,
        'spliceai_max': 0.02,  # Check for cryptic splice sites
        # Missense predictors less reliable for nonsensesaved': "✅ Data saved successfully.",
        'revel': 0.10,✅ Variant analysis completed successfully.",
        'alphamissense': 0.05,"✅ Report generated successfully.",
        'clinpred': 0.05PI request completed successfully.",
    },loaded': "✅ File loaded successfully.",
    'frameshift': {sed': "✅ Data validation passed."
        # Similar to nonsense but with some missense predictor utility
        'cadd_phred': 0.35,
        'gerp_pp': 0.20, messages
        'phylop_vert': 0.12,NG_MESSAGES = {
        'phylop_mamm': 0.08,   'api_timeout': "⚠️ API request timed out. Using cached data if available.",
        'phylop_prim': 0.06,    'low_confidence': "⚠️ Low confidence in prediction results.",
        'phastcons_vert': 0.04,e missing.",
        'phastcons_mamm': 0.03,a may be incomplete.",
        'phastcons_prim': 0.02,edictor_warning': "⚠️ Limited in silico predictor data available.",
        'spliceai_max': 0.05,ilable."
        'revel': 0.08,
        'alphamissense': 0.04,
        'clinpred': 0.03
    },_MESSAGES = {
    'splice_donor': {g_analysis': "ℹ️ Starting variant analysis...",
        # Splice predictors are primary for splice variantsiant data...",
        'spliceai_dg': 0.25,lculating_scores': "ℹ️ Calculating predictor scores...",
        'spliceai_dl': 0.25,criteria': "ℹ️ Applying ACMG/AMP criteria...",
        'spliceai_max': 0.15,t': "ℹ️ Generating classification report...",
        'mmsplice': 0.12, Saving analysis results..."
        'ada_score': 0.08,
        'rf_score': 0.08,
        'dbscsnv_ada': 0.05, Feature Status Messages
        'dbscsnv_rf': 0.05,
        'cadd_phred': 0.10,
        'gerp_pp': 0.08, AI Literature Scanning: ENABLED',
        'phylop_vert': 0.05,    'disabled': '📚 AI Literature Scanning: DISABLED',
        'phylop_mamm': 0.03,ure Analysis in progress...',
        'phylop_prim': 0.02,terature Analysis complete'
        'phastcons_vert': 0.02,
        'phastcons_mamm': 0.01,
        'phastcons_prim': 0.01atching: ENABLED',
    },henotype Matching: DISABLED',
    'splice_acceptor': {lysis in progress...',
        # Similar to splice_donor but with acceptor-specific weights
        'spliceai_ag': 0.25,
        'spliceai_al': 0.25,
        'spliceai_max': 0.15,
        'mmsplice': 0.12,    'disabled': '📊 Dynamic Metascore: DISABLED',
        'ada_score': 0.08,ing': '🔍 Calculating Dynamic Metascore...',
        'rf_score': 0.08,        'complete': '✅ Dynamic Metascore calculation complete'
        'dbscsnv_ada': 0.05,
        'dbscsnv_rf': 0.05,
        'cadd_phred': 0.10,
        'gerp_pp': 0.08,
        'phylop_vert': 0.05,
        'phylop_mamm': 0.03,
        'phylop_prim': 0.02,
        'phastcons_vert': 0.02,ynamic_metascore_enabled': True,
        'phastcons_mamm': 0.01,
        'phastcons_prim': 0.01
    },
    'synonymous': {
        # For synonymous variants, focus on conservation and splice impact
        'spliceai_max': 0.20,lculate_dynamic_vampp_score(insilico_scores, variant_type, allele_frequency=None, gene=None):
        'mmsplice': 0.15,
        'ada_score': 0.10,ic predictor weighting based on variant type.
        'rf_score': 0.10,
        'cadd_phred': 0.12,
        'gerp_pp': 0.10,
        'phylop_vert': 0.08,   variant_type (str): Type of variant (missense, nonsense, etc.)
        'phylop_mamm': 0.05,       allele_frequency (float): Population allele frequency
        'phylop_prim': 0.03,        gene (str): Gene symbol for gene-specific adjustments
        'phastcons_vert': 0.03,
        'phastcons_mamm': 0.02,
        'phastcons_prim': 0.02,nced Metascore results with detailed breakdown
        # Some missense predictors might still be useful
        'revel': 0.05,
        'alphamissense': 0.03,ific weights, fallback to missense if not found
        'clinpred': 0.02IANT_TYPE_PREDICTOR_WEIGHTS.get(variant_type, 
    },                                    VARIANT_TYPE_PREDICTOR_WEIGHTS['missense'])
    'intronic': {
        # For intronic variants, splice prediction is key   # Initialize score calculation
        'spliceai_max': 0.30,    weighted_sum = 0.0
        'mmsplice': 0.20, 0.0
        'ada_score': 0.15,
        'rf_score': 0.15,
        'dbscsnv_ada': 0.08,
        'dbscsnv_rf': 0.08,
        'cadd_phred': 0.08,
        'gerp_pp': 0.05,    if predictor in variant_weights and score is not None:
        'phylop_vert': 0.03,
        'phylop_mamm': 0.02,re to 0-1 range if needed
        'phylop_prim': 0.02,dictor_score(predictor, score)
        'phastcons_vert': 0.02,
        'phastcons_mamm': 0.01,fic weight
        'phastcons_prim': 0.01
    }
}specific adjustments if available
GENE_SPECIFIC_PREDICTORS:
# 2. DE NOVO ASSIGNMENT LOGIC
def assign_de_novo_criteria(de_novo_status, maternity_confirmed=None, paternity_confirmed=None):_predictors', []):
    """t *= gene_config.get('weight_boost', 1.0)
    Assign PS2 or PM6 based on de novo status and parental confirmation.            
    ed contribution
    Args:rmalized_score * weight
        de_novo_status (str): 'confirmed', 'assumed', or 'no'
        maternity_confirmed (bool): Whether maternity is confirmedt
        paternity_confirmed (bool): Whether paternity is confirmed
    
    Returns:
        dict: De novo criteria assignment
    """ore,
    result = {
        'PS2': False,ibution,
        'PM6': False,SHOLDS.get(predictor, {}).get('pathogenic', 0.5)
        'PS2_Very_Strong': False,
        'details': ''            used_predictors.append(predictor)
    }
    
    if de_novo_status == 'confirmed':                # Log error but continue with other predictors
        # PS2 requires both maternity and paternity confirmed= {
        if maternity_confirmed and paternity_confirmed:
            result['PS2'] = True
            result['details'] = 'De novo confirmed with both maternity and paternity testing'
        elif maternity_confirmed or paternity_confirmed:
            # Only one parent confirmed - still strong evidence but not full PS2core
            result['PS2'] = True
            result['details'] = 'De novo confirmed with single parent testing'hted_sum / total_weight
        else:
            # Confirmed but method unclear - assume PS2
            result['PS2'] = True
            result['details'] = 'De novo confirmed (parental testing method unclear)'ased adjustments
        frequency_category = get_frequency_category(allele_frequency)
    elif de_novo_status == 'assumed':e, frequency_category)
        # PM6 for assumed de novo (no parental testing)mpp_score + frequency_adjustment))
        result['PM6'] = True
        result['details'] = 'De novo assumed (no parental genetic testing performed)'
    
    elif de_novo_status == 'no' or de_novo_status == 'inherited':
        result['details'] = 'Variant is inherited or de novo status negative'
    
    return result
        if predictor in predictor_details:
# De novo assignment rules for different scenarios
DE_NOVO_ASSIGNMENT_RULES = {
    'confirmed_both_parents': {
        'criteria': 'PS2',
        'strength': 'Strong',
        'requirements': ['de_novo_status == "confirmed"', 'maternity_confirmed == True', 'paternity_confirmed == True'],                else:
        'description': 'De novo confirmed with both maternal and paternal testing'
    },
    'confirmed_single_parent': {
        'criteria': 'PS2',athogenic_votes / total_votes if total_votes > 0 else 0.0
        'strength': 'Strong',
        'requirements': ['de_novo_status == "confirmed"', 'maternity_confirmed == True OR paternity_confirmed == True'],
        'description': 'De novo confirmed with single parent testing'RIANT_TYPE_MIN_PREDICTORS.get(variant_type, 5)
    },ors and consensus >= 0.7:
    'assumed_clinical': {
        'criteria': 'PM6',en(used_predictors) >= 3 and consensus >= 0.6:
        'strength': 'Moderate','Medium'
        'requirements': ['de_novo_status == "assumed"'],
        'description': 'De novo assumed based on clinical assessment without parental testing'
    },
    'inherited_or_unknown': {iant type and frequency
        'criteria': None,c_thresholds = get_dynamic_thresholds(variant_type, frequency_category)
        'strength': None,
        'requirements': ['de_novo_status == "no" OR de_novo_status == "inherited"'],
        'description': 'Variant is inherited or de novo status not established'vampp_score,
    }score': adjusted_vampp_score,
}: variant_type,
requency_category': frequency_category,
# 3. DYNAMIC PP3/BP4 THRESHOLDSment': frequency_adjustment,
# Based on variant type and allele frequency categories
DYNAMIC_METASCORE_THRESHOLDS = {s': used_predictors,
    'missense': {ls': predictor_details,
        'ultra_rare': {  # AF ≤ 1e-5
            'pp3': 0.35,  'pathogenic_votes': pathogenic_votes,
            'bp4': 0.60benign_votes,
        },s': total_votes,
        'very_rare': {  # AF ≤ 1e-4: consensus
            'pp3': 0.45,
            'bp4': 0.50confidence,
        },hresholds': dynamic_thresholds,
        'moderate_rare': {  # 1e-4 < AF ≤ 1e-3': adjusted_vampp_score >= dynamic_thresholds['pp3'],
            'pp3': 0.55,],
            'bp4': 0.40sed': total_weight,
        },djustments': gene in GENE_SPECIFIC_PREDICTORS if gene else False
        'common': {  # AF > 1e-3
            'pp3': 0.70,
            'bp4': 0.30AMIC WEIGHTING
        }ant types
    },WEIGHTS = {
    'nonsense': {
        'ultra_rare': {ance predictors for missense variants
            'pp3': 0.25,  # Lower threshold for clear LoFrevel': 0.25,
            'bp4': 0.70   'alphamissense': 0.20,
        },
        'very_rare': {        'primateai': 0.12,
            'pp3': 0.30,
            'bp4': 0.65
        },
        'moderate_rare': {    'polyphen2': 0.03,
            'pp3': 0.40,
            'bp4': 0.55ctors
        },    'metasvm': 0.04,
        'common': {
            'pp3': 0.60,
            'bp4': 0.40
        }
    },wer weight for missense)
    'frameshift': {
        'ultra_rare': {
            'pp3': 0.30,
            'bp4': 0.65
        },stcons_vert': 0.01,
        'very_rare': {
            'pp3': 0.35,
            'bp4': 0.60    # Functional prediction
        },
        'moderate_rare': {
            'pp3': 0.45,
            'bp4': 0.50
        },
        'common': {
            'pp3': 0.65,
            'bp4': 0.35edictors are more important
        }d_phred': 0.30,
    },
    'splice_donor': {
        'ultra_rare': {    'phylop_mamm': 0.10,
            'pp3': 0.30,lop_prim': 0.08,
            'bp4': 0.65
        },
        'very_rare': {
            'pp3': 0.40,
            'bp4': 0.55   # Missense predictors less reliable for nonsense
        },        'revel': 0.10,
        'moderate_rare': {
            'pp3': 0.50,
            'bp4': 0.45
        },
        'common': { with some missense predictor utility
            'pp3': 0.65,
            'bp4': 0.35
        }
    },       'phylop_mamm': 0.08,
    'splice_acceptor': {        'phylop_prim': 0.06,
        'ultra_rare': {
            'pp3': 0.30,0.03,
            'bp4': 0.65 0.02,
        },.05,
        'very_rare': {
            'pp3': 0.40,0.04,
            'bp4': 0.55       'clinpred': 0.03
        },    },
        'moderate_rare': {
            'pp3': 0.50,e primary for splice variants
            'bp4': 0.45': 0.25,
        },5,
        'common': {
            'pp3': 0.65,2,
            'bp4': 0.35
        }ore': 0.08,
    },snv_ada': 0.05,
    'synonymous': {       'dbscsnv_rf': 0.05,
        'ultra_rare': {        'cadd_phred': 0.10,
            'pp3': 0.60,  # Higher threshold for synonymous
            'bp4': 0.40
        },        'phylop_mamm': 0.03,
        'very_rare': {
            'pp3': 0.70,rt': 0.02,
            'bp4': 0.30
        },
        'moderate_rare': {
            'pp3': 0.80, {
            'bp4': 0.20to splice_donor but with acceptor-specific weights
        },
        'common': {       'spliceai_al': 0.25,
            'pp3': 0.85,        'spliceai_max': 0.15,
            'bp4': 0.15
        }
    },
    'intronic': {
        'ultra_rare': {
            'pp3': 0.50,
            'bp4': 0.50
        },
        'very_rare': {  'phylop_mamm': 0.03,
            'pp3': 0.60,,
            'bp4': 0.40
        },
        'moderate_rare': {
            'pp3': 0.70,
            'bp4': 0.30
        },  # For synonymous variants, focus on conservation and splice impact
        'common': {': 0.20,
            'pp3': 0.80,
            'bp4': 0.20
        }
    }
}  'gerp_pp': 0.10,
': 0.08,
# 4. PP4 PHENOTYPE MATCH DYNAMIC ASSIGNMENT
def assign_pp4_phenotype_match(hpo_similarity_score):
    """
    Assign PP4 criteria based on HPO similarity score.
       'phastcons_prim': 0.02,
    Args:       # Some missense predictors might still be useful
        hpo_similarity_score (float): HPO similarity score (0-1)        'revel': 0.05,
    
    Returns:
        dict: PP4 assignment result
    """
    result = {nts, splice prediction is key
        'PP4': False,
        'strength': None,
        'details': ''
    }
    
    if hpo_similarity_score is None:  'dbscsnv_rf': 0.08,
        result['details'] = 'HPO similarity score not provided'08,
        return result
    ,
    if hpo_similarity_score >= 0.8:
        result['PP4'] = True.02,
        result['strength'] = 'Moderate'  # Upgrade to moderate for very high similarity02,
        result['details'] = f'High phenotype similarity (HPO score: {hpo_similarity_score:.3f})'.01,
    elif hpo_similarity_score >= 0.6: 0.01
        result['PP4'] = True
        result['strength'] = 'Supporting'
        result['details'] = f'Moderate phenotype similarity (HPO score: {hpo_similarity_score:.3f})'
    else:C
        result['PP4'] = Falsetatus, maternity_confirmed=None, paternity_confirmed=None):
        result['details'] = f'Low phenotype similarity (HPO score: {hpo_similarity_score:.3f} < 0.6)'
    de novo status and parental confirmation.
    return result

# PP4 phenotype assignment rules  de_novo_status (str): 'confirmed', 'assumed', or 'no'
PP4_PHENOTYPE_ASSIGNMENT = {irmed (bool): Whether maternity is confirmed
    'high_similarity': {d (bool): Whether paternity is confirmed
        'score_range': [0.8, 1.0],
        'criteria': 'PP4',
        'strength': 'Moderate',  # Upgrade for very high similarity
        'description': 'High phenotype similarity with gene-disease association'
    },esult = {
    'moderate_similarity': {       'PS2': False,
        'score_range': [0.6, 0.8],        'PM6': False,
        'criteria': 'PP4',
        'strength': 'Supporting', 'details': ''
        'description': 'Moderate phenotype similarity with gene-disease association'
    },
    'low_similarity': {_novo_status == 'confirmed':
        'score_range': [0.0, 0.6],ernity and paternity confirmed
        'criteria': None,
        'strength': None,
        'description': 'Insufficient phenotype similarity for PP4 application'        result['details'] = 'De novo confirmed with both maternity and paternity testing'
    } maternity_confirmed or paternity_confirmed:
}but not full PS2
     result['PS2'] = True
# Minimum predictor counts for confident variant type assessment'details'] = 'De novo confirmed with single parent testing'
VARIANT_TYPE_MIN_PREDICTORS = {
    'missense': 5,onfirmed but method unclear - assume PS2
    'nonsense': 3,
    'frameshift': 3,        result['details'] = 'De novo confirmed (parental testing method unclear)'
    'splice_donor': 4,
    'splice_acceptor': 4,ed':
    'synonymous': 3,de novo (no parental testing)
    'intronic': 3
} 'De novo assumed (no parental genetic testing performed)'

# Helper functions
def normalize_predictor_score(predictor, score):
    """Normalize predictor scores to 0-1 range for consistent weighting."""
    
    # Handle inverted predictors (lower = more pathogenic)
    inverted_predictors = ['sift', 'fathmm', 'provean']assignment rules for different scenarios
    
    if predictor in inverted_predictors:arents': {
        if predictor == 'sift':   'criteria': 'PS2',
            # SIFT: 0-1, lower = more pathogenic    'strength': 'Strong',
            return 1.0 - score'requirements': ['de_novo_status == "confirmed"', 'maternity_confirmed == True', 'paternity_confirmed == True'],
        elif predictor == 'fathmm':confirmed with both maternal and paternal testing'
            # FATHMM: typically -16 to +16, lower = more pathogenic
            normalized = max(0.0, min(1.0, (16 - score) / 32))
            return normalized
        elif predictor == 'provean':
            # PROVEAN: typically -14 to +14, lower = more pathogenic'requirements': ['de_novo_status == "confirmed"', 'maternity_confirmed == True OR paternity_confirmed == True'],
            normalized = max(0.0, min(1.0, (14 - score) / 28))ting'
            return normalized
    
    # Handle specific score rangesPM6',
    if predictor == 'cadd_phred':
        # CADD: 0-50+, higher = more pathogenic': ['de_novo_status == "assumed"'],
        return min(1.0, score / 50.0)inical assessment without parental testing'
    elif predictor == 'gerp_pp':
        # GERP++: -12 to +12, higher = more conservedknown': {
        return max(0.0, min(1.0, (score + 12) / 24))
    elif predictor in ['phylop_vert', 'phylop_mamm', 'phylop_prim']:
        # PhyloP: -20 to +20, higher = more conservedde_novo_status == "no" OR de_novo_status == "inherited"'],
        return max(0.0, min(1.0, (score + 20) / 40)) inherited or de novo status not established'
    elif predictor == 'mutationassessor':
        # MutationAssessor: 0-6, higher = more pathogenic
        return min(1.0, score / 6.0)
    PP3/BP4 THRESHOLDS
    # Default: assume 0-1 rangee frequency categories
    return max(0.0, min(1.0, score))METASCORE_THRESHOLDS = {

def get_frequency_category(allele_frequency):'ultra_rare': {  # AF ≤ 1e-5
    """Determine frequency category for threshold selection."""
    if allele_frequency is None or allele_frequency == 0:
        return 'ultra_rare'
    elif allele_frequency <= 1e-5:
        return 'ultra_rare''pp3': 0.45,
    elif allele_frequency <= 1e-4:
        return 'very_rare'
    elif allele_frequency <= 1e-3:_rare': {  # 1e-4 < AF ≤ 1e-3
        return 'moderate_rare'
    else:
        return 'common'

def get_frequency_adjustment(variant_type, frequency_category):
    """Get frequency-based score adjustment.""": 0.30
    adjustments = {
        'ultra_rare': -0.05,  # Slightly lower threshold
        'very_rare': -0.02,
        'moderate_rare': 0.0,'ultra_rare': {
        'common': 0.05       # Slightly higher threshold needed
    }
    return adjustments.get(frequency_category, 0.0)
y_rare': {
def get_dynamic_thresholds(variant_type, frequency_category):
    """Get dynamic PP3/BP4 thresholds based on variant type and frequency."""
    variant_thresholds = DYNAMIC_METASCORE_THRESHOLDS.get(variant_type, 
                                                     DYNAMIC_METASCORE_THRESHOLDS['missense'])e_rare': {
    return variant_thresholds.get(frequency_category, {'pp3': 0.5, 'bp4': 0.5})

# Enhanced ACMG evidence integration functions
def get_variant_type_evidence_modifiers(variant_type):
    """Get evidence strength modifiers based on variant type."""
    modifiers = {
        'missense': {
            'PP3': 1.0,  # Standard strength
            'BP4': 1.0,
            'PVS1': 0.0,  # Not applicableare': {
            'PP2': 1.2   # Boost for missense in constrained genes
        },
        'nonsense': {
            'PP3': 0.8,  # Reduced importance of in silico for clear LoF
            'BP4': 0.5,  # Very unlikely to be benign
            'PVS1': 1.0,  # Applicable
            'PP2': 0.0   # Not applicable
        },
        'frameshift': {    'pp3': 0.45,
            'PP3': 0.7,
            'BP4': 0.3,
            'PVS1': 1.0,
            'PP2': 0.0        'pp3': 0.65,
        }, 0.35
        'splice_donor': {        }
            'PP3': 1.1,  # Splice predictors are strong
            'BP4': 0.4,lice_donor': {
            'PVS1': 1.0,
            'PP2': 0.0        'pp3': 0.30,
        },   'bp4': 0.65
        'splice_acceptor': {
            'PP3': 1.1,    'very_rare': {
            'BP4': 0.4,'pp3': 0.40,
            'PVS1': 1.0,
            'PP2': 0.0 },
        },
        'synonymous': {.50,
            'PP3': 1.2,  # Important for synonymous to check splice impact
            'BP4': 1.0,
            'PVS1': 0.0,
            'PP2': 0.0       'pp3': 0.65,
        },        'bp4': 0.35
        'intronic': {
            'PP3': 1.3,  # Very important for intronic variants
            'BP4': 1.0,
            'PVS1': 0.0,
            'PP2': 0.0
        }
    }   },
    return modifiers.get(variant_type, modifiers['missense'])    'very_rare': {

def calculate_weighted_evidence_score(evidence_dict, variant_type):
    """Calculate weighted evidence score considering variant type."""
    modifiers = get_variant_type_evidence_modifiers(variant_type)
    
    pathogenic_score = 0
    benign_score = 0
    
    # Weight pathogenic evidence
    for criterion, details in evidence_dict.get('pathogenic_criteria', {}).items():
        if details.get('applies', False):
            base_weight = {,
                'Very Strong': 8,'synonymous': {
                'Strong': 4,
                'Moderate': 2,
                'Supporting': 1        'bp4': 0.40
            }.get(details.get('strength', 'Supporting'), 1)
            {
            modifier = modifiers.get(criterion, 1.0)
            pathogenic_score += base_weight * modifier
    
    # Weight benign evidence{
    for criterion, details in evidence_dict.get('benign_criteria', {}).items():
        if details.get('applies', False):        'bp4': 0.20
            base_weight = {
                'Stand-alone': 8,
                'Strong': 4,
                'Supporting': 1        'bp4': 0.15
            }.get(details.get('strength', 'Supporting'), 1)
            
            modifier = modifiers.get(criterion, 1.0)'intronic': {
            benign_score += base_weight * modifier
    
    return {        'bp4': 0.50
        'pathogenic_score': pathogenic_score,
        'benign_score': benign_score,
        'net_score': pathogenic_score - benign_score,        'pp3': 0.60,
        'confidence': min(pathogenic_score + benign_score, 10) / 10
    }
    'moderate_rare': {
# Statistical thresholds for various tests
STATISTICAL_THRESHOLDS = {
    'case_control_odds_ratio': 2.0,    },
    'case_control_p_value': 0.05,
    'segregation_lod_score': 3.0,
    'frequency_difference_threshold': 0.001,
    'minimum_case_count': 10,
    'minimum_control_count': 100
}

# Metascore thresholds (legacy compatibility)
VAMPP_SCORE_THRESHOLDS = {e):
    'pp3_threshold': 0.5,
    'bp4_threshold': 0.5,
    'high_confidence': 0.8,
    'low_confidence': 0.3
}po_similarity_score (float): HPO similarity score (0-1)

# Missing constants for backward compatibility
FUNCTIONAL_IMPACT_WEIGHTS = {}
PREDICTOR_TIERS = {}
ADVANCED_VAMPP_CONFIG = {}result = {
HIGH_IMPACT_VARIANT_CONFIG = {}
ML_ENSEMBLE_CONFIG = {}        'strength': None,
IMPACT_BASED_EVIDENCE_MODULATION = {
    'PP3': {},
    'BP4': {}    
}
 = 'HPO similarity score not provided'
# AI-POWERED LITERATURE SCANNING AND AUTOMATED EVIDENCE DETECTION
# =================================================================
 0.8:
# PubMed API configuration for literature scanning] = True
PUBMED_API_CONFIG = {th'] = 'Moderate'  # Upgrade to moderate for very high similarity
    'base_url': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils',tails'] = f'High phenotype similarity (HPO score: {hpo_similarity_score:.3f})'
    'email': 'cansevilmiss@gmail.com',  # Required for NCBI APIif hpo_similarity_score >= 0.6:
    'tool': 'ACMG_Assistant',
    'max_results': 50,trength'] = 'Supporting'
    'timeout': 30,       result['details'] = f'Moderate phenotype similarity (HPO score: {hpo_similarity_score:.3f})'
    'rate_limit': 3  # requests per second    else:
}
 result['details'] = f'Low phenotype similarity (HPO score: {hpo_similarity_score:.3f} < 0.6)'
# Literature scanning patterns for variant evidence
LITERATURE_SEARCH_PATTERNS = {return result
    'variant_position': [
        '{gene}[title] AND {position}[title]',
        '{gene}[title] AND {amino_acid_change}[title]',
        '{gene}[title] AND {cdna_change}[title]','high_similarity': {
        '{gene} AND c.{position} AND pathogenic',re_range': [0.8, 1.0],
        '{gene} AND p.{amino_acid_residue} AND variant'
    ], 'strength': 'Moderate',  # Upgrade for very high similarity
    'functional_studies': [ption': 'High phenotype similarity with gene-disease association'
        '{gene}[title] AND functional[title] AND (study OR assay OR experiment)',
        '{gene} AND functional AND (pathogenic OR benign OR neutral)',
        '{gene} AND "functional analysis" AND variant',
        '{gene} AND "in vitro" AND mutation',
        '{gene} AND "cell culture" AND variant'ng',
    ],henotype similarity with gene-disease association'
    'case_control': [
        '{gene} AND case AND control AND frequency', {
        '{gene} AND cohort AND variant AND pathogenic',   'score_range': [0.0, 0.6],
        '{gene} AND "case-control" AND association',    'criteria': None,
        '{gene} AND epidemiology AND variant''strength': None,
    ], application'
    'segregation': [
        '{gene} AND segregation AND family',
        '{gene} AND inheritance AND pedigree',
        '{gene} AND "family study" AND variant',m predictor counts for confident variant type assessment
        '{gene} AND cosegregation AND pathogenic'
    ]
}

# NLP patterns for extracting evidence from abstractsdonor': 4,
NLP_EVIDENCE_PATTERNS = {
    'pathogenic_evidence': [
        r'pathogenic.*variant',c': 3
        r'disease.*causing',
        r'deleterious.*effect',
        r'loss.*function',
        r'damaging.*mutation',
        r'causal.*variant',ighting."""
        r'functionally.*significant'
    ],re pathogenic)
    'benign_evidence': [
        r'benign.*variant',
        r'neutral.*effect',
        r'no.*functional.*impact',
        r'polymorphism',FT: 0-1, lower = more pathogenic
        r'normal.*function',
        r'not.*pathogenic',
        r'likely.*benign'
    ],lized = max(0.0, min(1.0, (16 - score) / 32))
    'functional_studies': [    return normalized
        r'functional.*stud(y|ies)','provean':
        r'in.*vitro.*assay',-14 to +14, lower = more pathogenic
        r'cell.*culture.*experiment',        normalized = max(0.0, min(1.0, (14 - score) / 28))
        r'protein.*function.*analysis', normalized
        r'enzymatic.*activity',    
        r'luciferase.*assay',
        r'rescue.*experiment'
    ],
    'same_position': [.0, score / 50.0)
        r'same.*position',== 'gerp_pp':
        r'identical.*residue',
        r'same.*amino.*acid',re + 12) / 24))
        r'previously.*reported.*position',_mamm', 'phylop_prim']:
        r'recurrent.*site'e conserved
    ]/ 40))
}
MutationAssessor: 0-6, higher = more pathogenic
def search_pubmed_literature_legacy(gene, variant_info, search_type='variant_position'):(1.0, score / 6.0)
    """
    Legacy search function for backward compatibility.
    
    Args:
        gene (str): Gene symbolequency_category(allele_frequency):
        variant_info (dict): Variant information including position, amino acid changerequency category for threshold selection."""
        search_type (str): Type of literature search to performquency == 0:
    
    Returns:
        dict: Literature search results with evidence extraction
    """allele_frequency <= 1e-4:
    import requests   return 'very_rare'
    import reelif allele_frequency <= 1e-3:
    import time
    from urllib.parse import quote    else:
    
    results = {
        'search_performed': True,
        'total_papers': 0,"""Get frequency-based score adjustment."""
        'relevant_papers': [],tments = {
        'evidence_found': { lower threshold
            'PS1_potential': False,  # Same amino acid change pathogenic
            'PM5_potential': False,  # Different change at same position    'moderate_rare': 0.0,
            'PS3_potential': False,  # Functional studies supporting pathogenicmon': 0.05       # Slightly higher threshold needed
            'BS3_potential': False,  # Functional studies supporting benign
            'PS4_potential': False,  # Case-control evidenceurn adjustments.get(frequency_category, 0.0)
            'PP1_potential': False   # Segregation evidence
        },get_dynamic_thresholds(variant_type, frequency_category):
        'search_queries': [],""
        'error': NoneMIC_METASCORE_THRESHOLDS.get(variant_type, 
    } DYNAMIC_METASCORE_THRESHOLDS['missense'])
    })
    try:
        # Prepare search queries
        patterns = LITERATURE_SEARCH_PATTERNS.get(search_type, [])type):
        position = variant_info.get('position', '')
        amino_acid_change = variant_info.get('amino_acid_change', '') = {
        cdna_change = variant_info.get('cdna_change', '')
                'PP3': 1.0,  # Standard strength
        # Extract amino acid residue number for broader searches
        amino_acid_residue = ''
        if amino_acid_change and 'p.' in amino_acid_change:        'PP2': 1.2   # Boost for missense in constrained genes
            import re
            match = re.search(r'p\.[A-Z][a-z]{2}(\d+)', amino_acid_change)
            if match:        'PP3': 0.8,  # Reduced importance of in silico for clear LoF
                amino_acid_residue = match.group(1)5,  # Very unlikely to be benign
        S1': 1.0,  # Applicable
        queries = []
        for pattern in patterns:
            query = pattern.format(
                gene=gene,
                position=position,
                amino_acid_change=amino_acid_change.replace('p.', ''),
                cdna_change=cdna_change.replace('c.', ''),
                amino_acid_residue=amino_acid_residue
            )   'splice_donor': {
            queries.append(query)        'PP3': 1.1,  # Splice predictors are strong
         0.4,
        results['search_queries'] = queries            'PVS1': 1.0,
        
        # Perform PubMed searches (mock implementation for now)
        # In real implementation, would use actual PubMed API calls
        for query in queries[:3]:  # Limit to first 3 queries for efficiency            'PP3': 1.1,
            time.sleep(1.0 / PUBMED_API_CONFIG['rate_limit'])  # Rate limiting
            
            # Mock response - in real implementation, would parse actual PubMed results
            mock_papers = [
                {
                    'pmid': '12345678',1.2,  # Important for synonymous to check splice impact
                    'title': f'{gene} variant functional analysis',
                    'abstract': f'This study demonstrates that the {amino_acid_change} variant in {gene} shows loss of function in cell culture assays.',
                    'year': '2023',
                    'relevance_score': 0.85  },
                }
            ]3,  # Very important for intronic variants
            
            results['relevant_papers'].extend(mock_papers) 0.0,
              'PP2': 0.0
        # Analyze abstracts for evidence (mock NLP analysis)
        if results['relevant_papers']:
            results['total_papers'] = len(results['relevant_papers'])
            
            # Mock evidence detection based on patterns
            for paper in results['relevant_papers']:
                abstract = paper.get('abstract', '').lower()
                
                # Check for functional study evidence   pathogenic_score = 0
                if any(re.search(pattern, abstract, re.IGNORECASE)     benign_score = 0
                       for pattern in NLP_EVIDENCE_PATTERNS['functional_studies']):
                    if any(re.search(pattern, abstract, re.IGNORECASE) eight pathogenic evidence
                           for pattern in NLP_EVIDENCE_PATTERNS['pathogenic_evidence']):teria', {}).items():
                        results['evidence_found']['PS3_potential'] = True    if details.get('applies', False):
                    elif any(re.search(pattern, abstract, re.IGNORECASE)    base_weight = {
                             for pattern in NLP_EVIDENCE_PATTERNS['benign_evidence']):
                        results['evidence_found']['BS3_potential'] = True
                
                # Check for same position evidence            'Supporting': 1
                if any(re.search(pattern, abstract, re.IGNORECASE) }.get(details.get('strength', 'Supporting'), 1)
                       for pattern in NLP_EVIDENCE_PATTERNS['same_position']):
                    if any(re.search(pattern, abstract, re.IGNORECASE)      modifier = modifiers.get(criterion, 1.0)
                           for pattern in NLP_EVIDENCE_PATTERNS['pathogenic_evidence']):hogenic_score += base_weight * modifier
                        results['evidence_found']['PS1_potential'] = True
                    else:
                        results['evidence_found']['PM5_potential'] = Truein evidence_dict.get('benign_criteria', {}).items():
        es', False):
    except Exception as e:
        results['error'] = str(e),
        results['search_performed'] = False
    rting': 1
    return resultsails.get('strength', 'Supporting'), 1)
       
def analyze_variant_literature(variant_data):        modifier = modifiers.get(criterion, 1.0)
    """    benign_score += base_weight * modifier
    Comprehensive literature analysis for automatic evidence detection.
    
    Args:
        variant_data: VariantData object with variant information
    'net_score': pathogenic_score - benign_score,
    Returns:_score, 10) / 10
        dict: Literature analysis results with evidence recommendations
    """
    gene = variant_data.basic_info.get('gene')tical thresholds for various tests
    variant_info = {
        'position': variant_data.basic_info.get('position'),l_odds_ratio': 2.0,
        'amino_acid_change': variant_data.basic_info.get('amino_acid_change'),
        'cdna_change': variant_data.basic_info.get('cdna_change')score': 3.0,
    }
    ': 10,
    # Perform different types of literature searches
    search_results = {
        'variant_position': search_pubmed_literature_legacy(gene, variant_info, 'variant_position'),
        'functional_studies': search_pubmed_literature_legacy(gene, variant_info, 'functional_studies'),ompatibility)
        'case_control': search_pubmed_literature_legacy(gene, variant_info, 'case_control'),
        'segregation': search_pubmed_literature_legacy(gene, variant_info, 'segregation')_threshold': 0.5,
    }
    
    # Consolidate evidence recommendations
    evidence_recommendations = {
        'PS1_recommended': False,
        'PM5_recommended': False,ward compatibility
        'PS3_recommended': False,
        'BS3_recommended': False,
        'PS4_recommended': False,
        'PP1_recommended': False,G = {}
        'literature_summary': '',
        'total_papers_found': 0,
        'confidence_level': 'Low'
    }
    
    total_papers = sum(result.get('total_papers', 0) for result in search_results.values())
    evidence_recommendations['total_papers_found'] = total_papers AUTOMATED EVIDENCE DETECTION
    
    # Aggregate evidence from all searches
    all_evidence = {}ature scanning
    for search_type, result in search_results.items():
        if result.get('search_performed') and not result.get('error'):e_url': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils',
            for evidence_type, found in result.get('evidence_found', {}).items():ired for NCBI API
                if found:l': 'ACMG_Assistant',
                    all_evidence[evidence_type] = True
    
    # Make recommendations based on aggregated evidence
    if all_evidence.get('PS1_potential'):
        evidence_recommendations['PS1_recommended'] = True
    
    if all_evidence.get('PM5_potential'):
        evidence_recommendations['PM5_recommended'] = Trueiant_position': [
    title]',
    if all_evidence.get('PS3_potential'):ND {amino_acid_change}[title]',
        evidence_recommendations['PS3_recommended'] = True
    AND pathogenic',
    if all_evidence.get('BS3_potential'):ue} AND variant'
        evidence_recommendations['BS3_recommended'] = True
    
    if all_evidence.get('PS4_potential'):ne}[title] AND functional[title] AND (study OR assay OR experiment)',
        evidence_recommendations['PS4_recommended'] = TrueR neutral)',
    ne} AND "functional analysis" AND variant',
    if all_evidence.get('PP1_potential'):
        evidence_recommendations['PP1_recommended'] = True
    
    # Set confidence based on number of papers and evidence foundntrol': [
    if total_papers >= 10 and len(all_evidence) >= 2:ncy',
        evidence_recommendations['confidence_level'] = 'High'genic',
    elif total_papers >= 5 and len(all_evidence) >= 1:
        evidence_recommendations['confidence_level'] = 'Medium'
    
    # Create literature summary
    evidence_types = [k.replace('_potential', '').replace('_recommended', '') 
                     for k in all_evidence.keys()]
    if evidence_types:
        evidence_recommendations['literature_summary'] = (ne} AND cosegregation AND pathogenic'
            f"Literature analysis of {total_papers} papers suggests potential evidence for: "
            f"{', '.join(evidence_types)}"
        )
    else:s for extracting evidence from abstracts
        evidence_recommendations['literature_summary'] = (EVIDENCE_PATTERNS = {
            f"Literature analysis of {total_papers} papers found limited specific evidence" [
        )
        r'disease.*causing',
    return evidence_recommendationsous.*effect',
        r'loss.*function',
# AUTOMATIC PHENOTYPE SIMILARITY SCORING
# ======================================

# HPO (Human Phenotype Ontology) similarity configuration
HPO_SIMILARITY_CONFIG = {'benign_evidence': [
    'api_endpoint': 'https://hpo.jax.org/api',  # HPO API endpoint
    'similarity_methods': ['jaccard', 'cosine', 'semantic'],al.*effect',
    'confidence_thresholds': {ional.*impact',
        'high': 0.8,sm',
        'medium': 0.6,*function',
        'low': 0.4
    },
    'gene_phenotype_database': 'https://hpo.jax.org/api/hpo/gene/',
    'timeout': 30
}

def calculate_hpo_similarity(patient_hpo_terms, gene_symbol):ulture.*experiment',
    """ein.*function.*analysis',
    Calculate HPO similarity between patient phenotype and gene-associated phenotypes.*activity',
    ase.*assay',
    Args:
        patient_hpo_terms (list): List of HPO terms for patient phenotype
        gene_symbol (str): Gene symbol to get associated phenotypes
    
    Returns:
        dict: HPO similarity analysis results
    """sly.*reported.*position',
    results = {rrent.*site'
        'similarity_score': 0.0,
        'method_used': 'mock',
        'patient_terms_count': len(patient_hpo_terms) if patient_hpo_terms else 0,
        'gene_terms_count': 0,
        'matching_terms': [],
        'confidence_level': 'Low',lity.
        'calculation_details': '',
        'error': None
    }): Gene symbol
    nt_info (dict): Variant information including position, amino acid change
    try:arch_type (str): Type of literature search to perform
        # Mock implementation - in real version would use HPO API
        # Get gene-associated HPO terms (mock data)
        gene_hpo_terms = get_gene_hpo_terms_mock(gene_symbol)erature search results with evidence extraction
        results['gene_terms_count'] = len(gene_hpo_terms)
        
        if patient_hpo_terms and gene_hpo_terms:
            # Calculate Jaccard similarity (mock calculation)
            patient_set = set(patient_hpo_terms)
            gene_set = set(gene_hpo_terms)
            
            intersection = patient_set.intersection(gene_set)ch_performed': True,
            union = patient_set.union(gene_set)total_papers': 0,
               'relevant_papers': [],
            results['matching_terms'] = list(intersection)    'evidence_found': {
            change pathogenic
            if union:            'PM5_potential': False,  # Different change at same position
                jaccard_similarity = len(intersection) / len(union)ogenic
                results['similarity_score'] = jaccard_similaritypporting benign
                rol evidence
                # Determine confidence level1_potential': False   # Segregation evidence
                if jaccard_similarity >= HPO_SIMILARITY_CONFIG['confidence_thresholds']['high']:    },
                    results['confidence_level'] = 'High'
                elif jaccard_similarity >= HPO_SIMILARITY_CONFIG['confidence_thresholds']['medium']:
                    results['confidence_level'] = 'Medium'}
                
                results['calculation_details'] = (
                    f"Jaccard similarity: {len(intersection)} matching terms / "
                    f"{len(union)} total unique terms = {jaccard_similarity:.3f}"    patterns = LITERATURE_SEARCH_PATTERNS.get(search_type, [])
                )riant_info.get('position', '')
        hange = variant_info.get('amino_acid_change', '')
    except Exception as e:riant_info.get('cdna_change', '')
        results['error'] = str(e)
        # Extract amino acid residue number for broader searches
    return results
d_change and 'p.' in amino_acid_change:
def get_gene_hpo_terms_mock(gene_symbol):
    """Mock function to get HPO terms associated with a gene."""re.search(r'p\.[A-Z][a-z]{2}(\d+)', amino_acid_change)
    # In real implementation, would query HPO database        if match:
    gene_hpo_map = {.group(1)
        'CAMTA1': [
            'HP:0001263',  # Global developmental delay
            'HP:0001250',  # Seizures
            'HP:0002119',  # Ventriculomegalyern.format(
            'HP:0001344',  # Absent speech
            'HP:0002540',  # Inability to walk=position,
            'HP:0000252'   # Microcephaly            amino_acid_change=amino_acid_change.replace('p.', ''),
        ],ange.replace('c.', ''),
        'BRCA1': [ino_acid_residue
            'HP:0003002',  # Breast carcinoma
            'HP:0100615',  # Ovarian neoplasmappend(query)
            'HP:0002664',  # Neoplasm    
            'HP:0000006'   # Autosomal dominant inheritanceeries'] = queries
        ],        
        'BRCA2': [
            'HP:0003002',  # Breast carcinomaubMed API calls
            'HP:0100615',  # Ovarian neoplasmn queries[:3]:  # Limit to first 3 queries for efficiency
            'HP:0002664',  # Neoplasmte limiting
            'HP:0000006'   # Autosomal dominant inheritance        
        ]actual PubMed results
    }
                {
    return gene_hpo_map.get(gene_symbol, [])
 f'{gene} variant functional analysis',
def auto_calculate_phenotype_similarity(variant_data, patient_hpo_terms=None): {amino_acid_change} variant in {gene} shows loss of function in cell culture assays.',
    """023',
    Automatically calculate phenotype similarity and update PP4 scoring.
                }
    Args:
        variant_data: VariantData object
        patient_hpo_terms (list): Optional HPO terms for patient
    
    Returns:
        dict: Updated phenotype similarity resultsers']:
    """rs'])
    gene = variant_data.basic_info.get('gene')
    tterns
    # If no HPO terms provided, try to extract from existing phenotype data        for paper in results['relevant_papers']:
    if not patient_hpo_terms:stract = paper.get('abstract', '').lower()
        # Mock extraction from phenotype_match field                
        phenotype_match = variant_data.functional_data.get('phenotype_match')
        if phenotype_match == 'specific_match':         if any(re.search(pattern, abstract, re.IGNORECASE) 
            patient_hpo_terms = ['HP:0001263', 'HP:0001250']  # Mock high matches']):
        elif phenotype_match == 'partial_match':                if any(re.search(pattern, abstract, re.IGNORECASE) 
            patient_hpo_terms = ['HP:0001263']  # Mock partial match                  for pattern in NLP_EVIDENCE_PATTERNS['pathogenic_evidence']):
        else:]['PS3_potential'] = True
            patient_hpo_terms = []e.IGNORECASE) 
                             for pattern in NLP_EVIDENCE_PATTERNS['benign_evidence']):
    # Calculate HPO similarity            results['evidence_found']['BS3_potential'] = True
    hpo_results = calculate_hpo_similarity(patient_hpo_terms, gene)
             # Check for same position evidence
    # Update PP4 assignment based on similarity scoreh(pattern, abstract, re.IGNORECASE) 
    pp4_assignment = assign_pp4_phenotype_match(hpo_results['similarity_score'])
                    if any(re.search(pattern, abstract, re.IGNORECASE) 
    # Combine results in NLP_EVIDENCE_PATTERNS['pathogenic_evidence']):
    results = {['evidence_found']['PS1_potential'] = True
        'hpo_similarity_score': hpo_results['similarity_score'],
        'hpo_confidence': hpo_results['confidence_level'],                    results['evidence_found']['PM5_potential'] = True
        'calculation_method': hpo_results['method_used'],
        'matching_terms_count': len(hpo_results['matching_terms']),
        'pp4_applicable': pp4_assignment['PP4'],results['error'] = str(e)
        'pp4_strength': pp4_assignment['strength'],
        'recommendation': hpo_results['calculation_details'],
        'auto_calculated': True
    }
    
    return results

# =============================================================================
# AI-POWERED LITERATURE SCANNING AND ANALYSIS
# =============================================================================nt information

# AI Literature Scanning Configuration
AI_LITERATURE_CONFIG = {h evidence recommendations
    'pubmed_api_base': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/',
    'max_results_per_query': 100,ene = variant_data.basic_info.get('gene')
    'search_timeout': 30,    variant_info = {
    'nlp_models': {
        'biobert': 'dmis-lab/biobert-base-cased-v1.1',nfo.get('amino_acid_change'),
        'scispacy': 'en_core_sci_lg',
        'clinicalbert': 'emilyalsentzer/Bio_ClinicalBERT'    }
    },
    'confidence_thresholds': { of literature searches
        'high': 0.85,
        'medium': 0.65,ch_pubmed_literature_legacy(gene, variant_info, 'variant_position'),
        'low': 0.45search_pubmed_literature_legacy(gene, variant_info, 'functional_studies'),
    },': search_pubmed_literature_legacy(gene, variant_info, 'case_control'),
    'evidence_keywords': {search_pubmed_literature_legacy(gene, variant_info, 'segregation')
        'PS1': ['same amino acid', 'identical residue', 'previously reported', 'pathogenic variant'],
        'PM5': ['different amino acid', 'same residue', 'missense change', 'pathogenic at position'],
        'PS3': ['functional study', 'in vitro', 'in vivo', 'damaging effect', 'protein function'],commendations
        'PS4': ['case-control', 'segregation', 'affected individuals', 'prevalence'],
        'PP1': ['cosegregation', 'family history', 'affected relatives', 'inheritance pattern'],
        'PP5': ['pathogenic', 'likely pathogenic', 'clinical significance', 'variant classification']
    }  'PS3_recommended': False,
}nded': False,
False,
def search_pubmed_literature(variant_info, evidence_type='PS1', max_results=50):False,
    """y': '',
    Search PubMed for literature evidence related to a variant.': 0,
       'confidence_level': 'Low'
    Args:   }
        variant_info (dict): Variant information including gene, position, change    
        evidence_type (str): Type of evidence to search for (PS1, PM5, PS3, PS4, etc.)
        max_results (int): Maximum number of results to returndence_recommendations['total_papers_found'] = total_papers
    
    Returns:# Aggregate evidence from all searches
        dict: Literature search results with relevance scoresvidence = {}
    """:
    results = {d') and not result.get('error'):
        'search_performed': True,'evidence_found', {}).items():
        'evidence_type': evidence_type,            if found:
        'search_terms': [],        all_evidence[evidence_type] = True
        'pubmed_results': [],
        'relevance_scores': [],ake recommendations based on aggregated evidence
        'evidence_found': False,ence.get('PS1_potential'):
        'confidence_level': 'None',ns['PS1_recommended'] = True
        'summary': '',
        'error': Noneial'):
    }PM5_recommended'] = True
    
    try:otential'):
        # Extract variant informationS3_recommended'] = True
        gene = variant_info.get('gene', '')
        aa_change = variant_info.get('aa_change', '')f all_evidence.get('BS3_potential'):
        cdna_change = variant_info.get('cdna_change', '')    evidence_recommendations['BS3_recommended'] = True
        
        # Build search terms based on evidence type):
        search_terms = []
        keywords = AI_LITERATURE_CONFIG['evidence_keywords'].get(evidence_type, [])
        
        # Base search terms
        if gene:
            search_terms.append(f'"{gene}"[Gene]')t confidence based on number of papers and evidence found
        if aa_change:ce) >= 2:
            search_terms.append(f'"{aa_change}"')idence_level'] = 'High'
        if cdna_change: total_papers >= 5 and len(all_evidence) >= 1:
            search_terms.append(f'"{cdna_change}"')el'] = 'Medium'
        
        # Add evidence-specific keywords
        for keyword in keywords:').replace('_recommended', '') 
            search_terms.append(f'"{keyword}"')         for k in all_evidence.keys()]
        :
        # Add variant-specific terms
        if evidence_type == 'PS1':gests potential evidence for: "
            if aa_change:
                # Search for same amino acid change
                position = ''.join(filter(str.isdigit, aa_change))
                if position:mary'] = (
                    search_terms.append(f'"{position}"')und limited specific evidence"
                    
        elif evidence_type == 'PM5':
            if aa_change:
                # Search for different changes at same position
                position = ''.join(filter(str.isdigit, aa_change))OTYPE SIMILARITY SCORING
                if position:======
                    search_terms.append(f'"{position}"')
                    search_terms.append('"missense"')
                    
        elif evidence_type == 'PS3':ndpoint
            search_terms.extend(['"functional analysis"', '"protein function"', '"in vitro"'])ilarity_methods': ['jaccard', 'cosine', 'semantic'],
            
        elif evidence_type == 'PS4':
            search_terms.extend(['"case control"', '"segregation"', '"affected"'])
        
        results['search_terms'] = search_terms
        ps://hpo.jax.org/api/hpo/gene/',
        # Mock PubMed API call (in real implementation, use requests)
        # This is a simulation of what would happen with real API
        search_query = ' AND '.join(search_terms)
        rms, gene_symbol):
        # Mock results based on evidence type and variant
        mock_results = _generate_mock_pubmed_results(variant_info, evidence_type)ulate HPO similarity between patient phenotype and gene-associated phenotypes.
        results['pubmed_results'] = mock_results
        
        # Analyze results for evidencePO terms for patient phenotype
        if mock_results:ciated phenotypes
            # Calculate relevance scores using NLP analysis
            relevance_scores = []
            for result in mock_results:dict: HPO similarity analysis results
                score = _calculate_literature_relevance(result, evidence_type, variant_info)
                relevance_scores.append(score)
            
            results['relevance_scores'] = relevance_scores
            if patient_hpo_terms else 0,
            # Determine if evidence was foundgene_terms_count': 0,
            max_score = max(relevance_scores) if relevance_scores else 0'matching_terms': [],
            thresholds = AI_LITERATURE_CONFIG['confidence_thresholds']: 'Low',
            ,
            if max_score >= thresholds['high']:    'error': None
                results['evidence_found'] = True
                results['confidence_level'] = 'High'    
            elif max_score >= thresholds['medium']:
                results['evidence_found'] = TrueAPI
                results['confidence_level'] = 'Medium'ock data)
            elif max_score >= thresholds['low']:
                results['evidence_found'] = True    results['gene_terms_count'] = len(gene_hpo_terms)
                results['confidence_level'] = 'Low'
            _hpo_terms and gene_hpo_terms:
            # Generate summary        # Calculate Jaccard similarity (mock calculation)
            results['summary'] = _generate_literature_summary(
                mock_results, relevance_scores, evidence_type
            )        
    nt_set.intersection(gene_set)
    except Exception as e:ient_set.union(gene_set)
        results['error'] = str(e)
    g_terms'] = list(intersection)
    return results

def _generate_mock_pubmed_results(variant_info, evidence_type):
    """Generate mock PubMed results for testing."""e'] = jaccard_similarity
    gene = variant_info.get('gene', 'UNKNOWN')            
    aa_change = variant_info.get('aa_change', '')evel
    thresholds']['high']:
    # Mock database of literature results
    mock_db = {            elif jaccard_similarity >= HPO_SIMILARITY_CONFIG['confidence_thresholds']['medium']:
        'CAMTA1': {ce_level'] = 'Medium'
            'PS1': [                
                {
                    'pmid': '12345678',             f"Jaccard similarity: {len(intersection)} matching terms / "
                    'title': f'Pathogenic {aa_change} variant in {gene} associated with developmental delay',similarity:.3f}"
                    'abstract': f'We report a pathogenic {aa_change} variant in {gene} that causes severe developmental delay...',            )
                    'authors': 'Smith J, et al.',
                    'journal': 'Am J Hum Genet',
                    'year': '2023'
                }
            ],esults
            'PM5': [
                {_gene_hpo_terms_mock(gene_symbol):
                    'pmid': '23456789',ed with a gene."""
                    'title': f'Novel missense variants in {gene} at position {aa_change[:4] if aa_change else ""}',# In real implementation, would query HPO database
                    'abstract': f'Different missense changes at the same position in {gene} show pathogenic effects...',
                    'authors': 'Johnson M, et al.',
                    'journal': 'Hum Mutat',
                    'year': '2022'
                }
            ],
            'PS3': [
                {HP:0000252'   # Microcephaly
                    'pmid': '34567890',
                    'title': f'Functional characterization of {gene} variants in cell culture',    'BRCA1': [
                    'abstract': f'In vitro studies demonstrate loss of function for {gene} variants...',
                    'authors': 'Davis R, et al.',
                    'journal': 'Mol Cell Biol',
                    'year': '2023'       'HP:0000006'   # Autosomal dominant inheritance
                }
            ]
        },       'HP:0003002',  # Breast carcinoma
        'BRCA1': {
            'PS1': [
                {       'HP:0000006'   # Autosomal dominant inheritance
                    'pmid': '45678901',    ]
                    'title': f'Recurrent {aa_change} mutation in BRCA1 in breast cancer families',
                    'abstract': 'The same amino acid change has been reported multiple times as pathogenic...',
                    'authors': 'Wilson K, et al.',
                    'journal': 'Cancer Res',
                    'year': '2023'ta, patient_hpo_terms=None):
                }""
            ]Automatically calculate phenotype similarity and update PP4 scoring.
        }
    }
    
    return mock_db.get(gene, {}).get(evidence_type, [])    patient_hpo_terms (list): Optional HPO terms for patient

def _calculate_literature_relevance(result, evidence_type, variant_info):
    """Calculate relevance score using mock NLP analysis."""    dict: Updated phenotype similarity results
    # Mock NLP scoring based on content analysis
    score = 0.0_info.get('gene')
    
    title = result.get('title', '').lower()g phenotype data
    abstract = result.get('abstract', '').lower()
    
    # Check for variant-specific termset('phenotype_match')
    gene = variant_info.get('gene', '').lower() 'specific_match':
    aa_change = variant_info.get('aa_change', '').lower():0001250']  # Mock high match
        elif phenotype_match == 'partial_match':
    if gene in title:patient_hpo_terms = ['HP:0001263']  # Mock partial match
        score += 0.3
    if gene in abstract:
        score += 0.2
    
    if aa_change and aa_change in title:ne)
        score += 0.4
    if aa_change and aa_change in abstract:
        score += 0.3'similarity_score'])
    
    # Check for evidence-specific keywords
    keywords = AI_LITERATURE_CONFIG['evidence_keywords'].get(evidence_type, [])esults = {
    for keyword in keywords:        'hpo_similarity_score': hpo_results['similarity_score'],
        if keyword.lower() in title:
            score += 0.2used'],
        if keyword.lower() in abstract:
            score += 0.1        'pp4_applicable': pp4_assignment['PP4'],
    rength'],
    # Recent publications get bonusion': hpo_results['calculation_details'],
    year = result.get('year', '2000')
    if year.isdigit() and int(year) >= 2020:
        score += 0.1
    
    return min(score, 1.0)
 =============================================================================
def _generate_literature_summary(results, scores, evidence_type):# AI-POWERED LITERATURE SCANNING AND ANALYSIS
    """Generate a summary of literature findings."""
    if not results:
        return f"No literature found for {evidence_type} evidence."
    ITERATURE_CONFIG = {
    high_score_count = sum(1 for score in scores if score >= 0.85)
    medium_score_count = sum(1 for score in scores if 0.65 <= score < 0.85)'max_results_per_query': 100,
    ch_timeout': 30,
    summary = f"Found {len(results)} publications for {evidence_type} evidence. "
    if high_score_count > 0:
        summary += f"{high_score_count} high-confidence matches. "
    if medium_score_count > 0:lsentzer/Bio_ClinicalBERT'
        summary += f"{medium_score_count} medium-confidence matches. "},
    nce_thresholds': {
    # Add specific findings
    if evidence_type == 'PS1': 'medium': 0.65,
        summary += "Same amino acid change previously reported as pathogenic."
    elif evidence_type == 'PM5':
        summary += "Different pathogenic changes at same amino acid position."
    elif evidence_type == 'PS3':        'PS1': ['same amino acid', 'identical residue', 'previously reported', 'pathogenic variant'],
        summary += "Functional studies support damaging effect.", 'missense change', 'pathogenic at position'],
    elif evidence_type == 'PS4': 'PS3': ['functional study', 'in vitro', 'in vivo', 'damaging effect', 'protein function'],
        summary += "Case-control data available."cted individuals', 'prevalence'],
        'PP1': ['cosegregation', 'family history', 'affected relatives', 'inheritance pattern'],
    return summaryPP5': ['pathogenic', 'likely pathogenic', 'clinical significance', 'variant classification']

def apply_ai_literature_evidence(variant_info, evidence_type='PS1'):
    """
    Apply AI-powered literature scanning to determine evidence applicability._type='PS1', max_results=50):
    
    Args:
        variant_info (dict): Variant information    





















































































































































































































































# METASCORE CALCULATION (RENAMED# =============================================================================    }        'auto_assigned': True        'calculation_method': 'Enhanced multi-method similarity',        'matching_terms': jaccard_results['matching_terms'],        'confidence_level': jaccard_results['confidence_level'],        'semantic_similarity': semantic_results['similarity_score'],        'cosine_similarity': cosine_results['similarity_score'],        'jaccard_similarity': jaccard_results['similarity_score'],        'combined_similarity_score': combined_score,        'pp4_strength': pp4_strength,        'pp4_applicable': pp4_applicable,    return {            pp4_strength = pp4_rules['low_confidence']        pp4_applicable = True    elif combined_score >= confidence_thresholds['low']:        pp4_strength = pp4_rules['medium_confidence']        pp4_applicable = True    elif combined_score >= confidence_thresholds['medium']:        pp4_strength = pp4_rules['high_confidence']        pp4_applicable = True    if combined_score >= confidence_thresholds['high']:        pp4_strength = 'None'    pp4_applicable = False        confidence_thresholds = ENHANCED_PHENOTYPE_CONFIG['confidence_thresholds']    pp4_rules = ENHANCED_PHENOTYPE_CONFIG['pp4_assignment_rules']    # Determine PP4 assignment        )        semantic_results['similarity_score'] * 0.3        cosine_results['similarity_score'] * 0.3 +        jaccard_results['similarity_score'] * 0.4 +    combined_score = (    # Combine results for final scoring        )        patient_hpo_terms, gene, method='semantic'    semantic_results = calculate_enhanced_phenotype_similarity(    )        patient_hpo_terms, gene, method='cosine'    cosine_results = calculate_enhanced_phenotype_similarity(    )        patient_hpo_terms, gene, method='jaccard'    jaccard_results = calculate_enhanced_phenotype_similarity(    # Calculate similarity using multiple methods                patient_hpo_terms = []        else:            patient_hpo_terms = ['HP:0001263', 'HP:0001250']        elif phenotype_match == 'partial_match':            patient_hpo_terms = ['HP:0001263', 'HP:0001250', 'HP:0002119']        if phenotype_match == 'specific_match':        phenotype_match = variant_data.functional_data.get('phenotype_match')    if not patient_hpo_terms:    # Extract HPO terms if not provided        gene = variant_data.basic_info.get('gene')    """        dict: Enhanced PP4 assignment results    Returns:            patient_hpo_terms (list): Patient's HPO terms        variant_data: VariantData object    Args:        Automatically assign PP4 using enhanced phenotype similarity.    """def auto_assign_pp4_with_enhanced_similarity(variant_data, patient_hpo_terms=None):    return min(total_similarity, 1.0)        total_similarity = direct_similarity + (related_score / 10)    direct_similarity = len(direct_matches) / max(len(terms1), len(terms2))    # Combine direct and related matches                        related_score += 0.1                if term1[:7] == term2[:7]:  # Same HP:0001xxx prefix                # Mock similarity based on term structure            if term1 != term2:        for term2 in terms2:    for term1 in terms1:    related_score = 0.0    # Mock related term matching        direct_matches = set(terms1).intersection(set(terms2))    # Simple mock: check for term overlaps and related terms            return 0.0    if not terms1 or not terms2:        # In real implementation, would use HPO ontology structure    # Mock semantic similarity calculation    """Calculate semantic similarity between HPO term sets."""def _calculate_semantic_similarity(terms1, terms2):    return results            results['error'] = str(e)    except Exception as e:                )            f"{len(gene_hpo_terms)} total gene terms)"            f"({len(results['matching_terms'])} matching terms, "            f"{method.title()} similarity: {score:.3f} "        results['calculation_details'] = (        # Generate detailed explanation                    results['confidence_level'] = 'Low'        elif score >= thresholds['low']:            results['confidence_level'] = 'Medium'        elif score >= thresholds['medium']:            results['confidence_level'] = 'High'        if score >= thresholds['high']:                thresholds = ENHANCED_PHENOTYPE_CONFIG['confidence_thresholds']        score = results['similarity_score']        # Determine confidence level                    results['weighted_score'] = min(weighted_score, 1.0)            ) / len(gene_hpo_terms)                for _ in results['matching_terms']                ENHANCED_PHENOTYPE_CONFIG['hpo_weights'].get('exact_match', 1.0)             weighted_score = sum(        if results['matching_terms']:        # Calculate weighted score                    results['semantic_similarity'] = semantic_score            results['similarity_score'] = semantic_score            semantic_score = _calculate_semantic_similarity(patient_hpo_terms, gene_hpo_terms)            # Mock semantic similarity using HPO hierarchy        elif method == 'semantic':                                results['matching_terms'] = list(intersection)                results['similarity_score'] = similarity_score                similarity_score = len(intersection) / (len(patient_set) * len(gene_set)) ** 0.5            if patient_set and gene_set:            intersection = patient_set.intersection(gene_set)            # Mock cosine similarity calculation        elif method == 'cosine':                                results['matching_terms'] = list(intersection)                results['similarity_score'] = similarity_score                similarity_score = len(intersection) / len(union)            if union:                        union = patient_set.union(gene_set)            intersection = patient_set.intersection(gene_set)        if method == 'jaccard':        # Calculate different similarity metrics                gene_set = set(gene_hpo_terms)        patient_set = set(patient_hpo_terms)                    return results            results['calculation_details'] = "Insufficient HPO terms for comparison"        if not gene_hpo_terms or not patient_hpo_terms:                gene_hpo_terms = get_gene_hpo_terms_mock(gene_symbol)        # Get gene-associated HPO terms    try:        }        'error': None        'semantic_similarity': 0.0,        'weighted_score': 0.0,        'matching_terms': [],        'calculation_details': '',        'confidence_level': 'None',        'similarity_score': 0.0,        'method_used': method,    results = {    """        dict: Enhanced similarity results    Returns:            method (str): Similarity calculation method        gene_symbol (str): Gene symbol        patient_hpo_terms (list): Patient's HPO terms    Args:        Calculate enhanced phenotype similarity using multiple methods.    """def calculate_enhanced_phenotype_similarity(patient_hpo_terms, gene_symbol, method='jaccard'):}    }        'sibling_term': 0.6        'child_term': 0.9,        'parent_term': 0.8,        'exact_match': 1.0,    'hpo_weights': {    },        'low_confidence': 'PP4_Supporting'        'medium_confidence': 'PP4',        'high_confidence': 'PP4_Strong',    'pp4_assignment_rules': {    },        'low': 0.25        'medium': 0.50,        'high': 0.75,    'confidence_thresholds': {    'default_method': 'jaccard',    'similarity_methods': ['jaccard', 'cosine', 'semantic'],ENHANCED_PHENOTYPE_CONFIG = {# Enhanced Phenotype Similarity Configuration# =============================================================================# ENHANCED AUTOMATIC PHENOTYPE SIMILARITY SCORING# =============================================================================    }        'ai_powered': True        'search_details': literature_results,        'literature_summary': literature_results['summary'],        'confidence_level': confidence_level,        'evidence_applicable': evidence_applicable,    return {                evidence_applicable = confidence_level in ['High', 'Medium', 'Low']            # Medium confidence acceptable for supporting evidence        elif evidence_type in ['PS3', 'PS4']:            evidence_applicable = confidence_level in ['High', 'Medium']            # Require high confidence for strong evidence        if evidence_type in ['PS1', 'PM5']:        # Different thresholds for different evidence types                confidence_level = literature_results['confidence_level']    if literature_results['evidence_found']:        confidence_level = 'None'    evidence_applicable = False    # Determine evidence applicability        literature_results = search_pubmed_literature(variant_info, evidence_type)    # Perform literature search    """        dict: Evidence evaluation results    Returns:            evidence_type (str): ACMG evidence type to evaluate    Args:
        variant_info (dict): Variant information including gene, position, change
        evidence_type (str): Type of evidence to search for (PS1, PM5, PS3, PS4, etc.)
        max_results (int): Maximum number of results to return
    
    Returns:
        dict: Literature search results with relevance scores
    """
    results = {
        'search_performed': True,
        'evidence_type': evidence_type,
        'search_terms': [],
        'pubmed_results': [],
        'relevance_scores': [],
        'evidence_found': False,
        'confidence_level': 'None',
        'summary': '',
        'error': None
    }
    
    try:
        # Extract variant information
        gene = variant_info.get('gene', '')
        aa_change = variant_info.get('aa_change', '')
        cdna_change = variant_info.get('cdna_change', '')
        
        # Build search terms based on evidence type
        search_terms = []
        keywords = AI_LITERATURE_CONFIG['evidence_keywords'].get(evidence_type, [])
        
        # Base search terms
        if gene:
            search_terms.append(f'"{gene}"[Gene]')
        if aa_change:
            search_terms.append(f'"{aa_change}"')
        if cdna_change:
            search_terms.append(f'"{cdna_change}"')
        
        # Add evidence-specific keywords
        for keyword in keywords:
            search_terms.append(f'"{keyword}"')
        
        # Add variant-specific terms
        if evidence_type == 'PS1':
            if aa_change:
                # Search for same amino acid change
                position = ''.join(filter(str.isdigit, aa_change))
                if position:
                    search_terms.append(f'"{position}"')
                    
        elif evidence_type == 'PM5':
            if aa_change:
                # Search for different changes at same position
                position = ''.join(filter(str.isdigit, aa_change))
                if position:
                    search_terms.append(f'"{position}"')
                    search_terms.append('"missense"')
                    
        elif evidence_type == 'PS3':
            search_terms.extend(['"functional analysis"', '"protein function"', '"in vitro"'])
            
        elif evidence_type == 'PS4':
            search_terms.extend(['"case control"', '"segregation"', '"affected"'])
        
        results['search_terms'] = search_terms
        
        # Mock PubMed API call (in real implementation, use requests)
        # This is a simulation of what would happen with real API
        search_query = ' AND '.join(search_terms)
        
        # Mock results based on evidence type and variant
        mock_results = _generate_mock_pubmed_results(variant_info, evidence_type)
        results['pubmed_results'] = mock_results
        
        # Analyze results for evidence
        if mock_results:
            # Calculate relevance scores using NLP analysis
            relevance_scores = []
            for result in mock_results:
                score = _calculate_literature_relevance(result, evidence_type, variant_info)
                relevance_scores.append(score)
            
            results['relevance_scores'] = relevance_scores
            
            # Determine if evidence was found
            max_score = max(relevance_scores) if relevance_scores else 0
            thresholds = AI_LITERATURE_CONFIG['confidence_thresholds']
            
            if max_score >= thresholds['high']:
                results['evidence_found'] = True
                results['confidence_level'] = 'High'
            elif max_score >= thresholds['medium']:
                results['evidence_found'] = True
                results['confidence_level'] = 'Medium'
            elif max_score >= thresholds['low']:
                results['evidence_found'] = True
                results['confidence_level'] = 'Low'
            
            # Generate summary
            results['summary'] = _generate_literature_summary(
                mock_results, relevance_scores, evidence_type
            )
    
    except Exception as e:
        results['error'] = str(e)
    
    return results

def _generate_mock_pubmed_results(variant_info, evidence_type):
    """Generate mock PubMed results for testing."""
    gene = variant_info.get('gene', 'UNKNOWN')
    aa_change = variant_info.get('aa_change', '')
    
    # Mock database of literature results
    mock_db = {
        'CAMTA1': {
            'PS1': [
                {
                    'pmid': '12345678',
                    'title': f'Pathogenic {aa_change} variant in {gene} associated with developmental delay',
                    'abstract': f'We report a pathogenic {aa_change} variant in {gene} that causes severe developmental delay...',
                    'authors': 'Smith J, et al.',
                    'journal': 'Am J Hum Genet',
                    'year': '2023'
                }
            ],
            'PM5': [
                {
                    'pmid': '23456789',
                    'title': f'Novel missense variants in {gene} at position {aa_change[:4] if aa_change else ""}',
                    'abstract': f'Different missense changes at the same position in {gene} show pathogenic effects...',
                    'authors': 'Johnson M, et al.',
                    'journal': 'Hum Mutat',
                    'year': '2022'
                }
            ],
            'PS3': [
                {
                    'pmid': '34567890',
                    'title': f'Functional characterization of {gene} variants in cell culture',
                    'abstract': f'In vitro studies demonstrate loss of function for {gene} variants...',
                    'authors': 'Davis R, et al.',
                    'journal': 'Mol Cell Biol',
                    'year': '2023'
                }
            ]
        },
        'BRCA1': {
            'PS1': [
                {
                    'pmid': '45678901',
                    'title': f'Recurrent {aa_change} mutation in BRCA1 in breast cancer families',
                    'abstract': 'The same amino acid change has been reported multiple times as pathogenic...',
                    'authors': 'Wilson K, et al.',
                    'journal': 'Cancer Res',
                    'year': '2023'
                }
            ]
        }
    }
    
    return mock_db.get(gene, {}).get(evidence_type, [])

def _calculate_literature_relevance(result, evidence_type, variant_info):
    """Calculate relevance score using mock NLP analysis."""
    # Mock NLP scoring based on content analysis
    score = 0.0
    
    title = result.get('title', '').lower()
    abstract = result.get('abstract', '').lower()
    
    # Check for variant-specific terms
    gene = variant_info.get('gene', '').lower()
    aa_change = variant_info.get('aa_change', '').lower()
    
    if gene in title:
        score += 0.3
    if gene in abstract:
        score += 0.2
    
    if aa_change and aa_change in title:
        score += 0.4
    if aa_change and aa_change in abstract:
        score += 0.3
    
    # Check for evidence-specific keywords
    keywords = AI_LITERATURE_CONFIG['evidence_keywords'].get(evidence_type, [])
    for keyword in keywords:
        if keyword.lower() in title:
            score += 0.2
        if keyword.lower() in abstract:
            score += 0.1
    
    # Recent publications get bonus
    year = result.get('year', '2000')
   