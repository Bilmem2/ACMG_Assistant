"""
Configuration constants for ACMG Variant Classification Assistant
===============================================================

This file contains all the constants, thresholds, and configuration
values used throughout the application.
"""

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

# In silico predictor weights for VAMPP-score-like metascore
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
    'dbscsnv_rf': 0.06
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
    # Other splice predictors
    'ada_score': {'pathogenic': 0.6, 'benign': 0.4},
    'rf_score': {'pathogenic': 0.6, 'benign': 0.4},
    'dbscsnv_ada': {'pathogenic': 0.6, 'benign': 0.4},
    'dbscsnv_rf': {'pathogenic': 0.6, 'benign': 0.4}
}

# Variant consequence types
# Variant consequence types (Sequence Ontology terms)
VARIANT_CONSEQUENCES = {
    'HIGH_IMPACT': [
        'chromosome_number_variation',
        'exon_loss_variant',
        'frameshift_variant',
        'rare_amino_acid_variant',
        'splice_acceptor_variant',
        'splice_donor_variant',
        'start_lost',
        'stop_gained',
        'stop_lost',
        'transcript_ablation'
    ],
    'MODERATE_IMPACT': [
        'missense_variant',
        'regulatory_region_ablation',
        'splice_region_variant',
        'protein_altering_variant',
        'inframe_deletion',
        'inframe_insertion',
        'conservative_inframe_deletion',
        'conservative_inframe_insertion',
        'disruptive_inframe_deletion',
        'disruptive_inframe_insertion'
    ],
    'LOW_IMPACT': [
        'synonymous_variant',
        'incomplete_terminal_codon_variant',
        'start_retained_variant',
        'stop_retained_variant',
        'NMD_transcript_variant'
    ],
    'MODIFIER': [
        'intron_variant',
        'intergenic_region',
        'upstream_gene_variant',
        'downstream_gene_variant',
        '3_prime_UTR_variant',
        '5_prime_UTR_variant',
        'non_coding_transcript_exon_variant',
        'non_coding_transcript_variant',
        'coding_sequence_variant',
        'mature_miRNA_variant',
        'regulatory_region_variant',
        'TF_binding_site_variant',
        'sequence_feature'
    ],
    'SPLICE_VARIANTS': [
        'splice_acceptor_variant',
        'splice_donor_variant', 
        'splice_region_variant'
    ],
    'REGULATORY': [
        '5_prime_UTR_variant',
        '3_prime_UTR_variant',
        'regulatory_region_variant',
        'TF_binding_site_variant',
        'upstream_gene_variant',
        'downstream_gene_variant'
    ]
}

# All possible variant consequences for user selection
ALL_VARIANT_CONSEQUENCES = [
    # High impact
    'frameshift_variant',
    'nonsense_variant',
    'stop_gained',
    'stop_lost',
    'start_lost',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'transcript_ablation',
    'exon_loss_variant',
    
    # Moderate impact  
    'missense_variant',
    'protein_altering_variant',
    'inframe_deletion',
    'inframe_insertion',
    'conservative_inframe_deletion',
    'conservative_inframe_insertion',
    'disruptive_inframe_deletion',
    'disruptive_inframe_insertion',
    'splice_region_variant',
    'regulatory_region_ablation',
    
    # Low impact
    'synonymous_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'incomplete_terminal_codon_variant',
    'NMD_transcript_variant',
    
    # Modifier/Intronic/Regulatory
    'intron_variant',
    'intergenic_region',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'upstream_gene_variant', 
    'downstream_gene_variant',
    'non_coding_transcript_exon_variant',
    'non_coding_transcript_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    'regulatory_region_variant',
    'TF_binding_site_variant',
    'sequence_feature',
    
    # Other/Unknown
    'other',
    'unknown'
]

# Population databases
POPULATION_DATABASES = {
    'gnomAD': {'url': 'https://gnomad.broadinstitute.org/', 'versions': ['v2.1.1', 'v3.1.2', 'v4.0']},
    'ExAC': {'url': 'http://exac.broadinstitute.org/', 'versions': ['v1.0']},
    '1000G': {'url': 'https://www.internationalgenome.org/', 'versions': ['phase3']},
    'ESP': {'url': 'https://evs.gs.washington.edu/EVS/', 'versions': ['v2']}
}

# Disease inheritance patterns
INHERITANCE_PATTERNS = {
    'AD': 'Autosomal Dominant',
    'AR': 'Autosomal Recessive',
    'XL': 'X-linked',
    'XLD': 'X-linked Dominant',
    'XLR': 'X-linked Recessive',
    'MT': 'Mitochondrial',
    'SMD': 'Somatic',
    'UNK': 'Unknown'
}

# Zygosity types
ZYGOSITY_TYPES = ['homozygous', 'heterozygous', 'hemizygous']

# Functional test types
FUNCTIONAL_TEST_TYPES = [
    'protein_function',
    'splicing',
    'promoter_activity',
    'dna_binding',
    'enzyme_activity',
    'cellular_localization',
    'protein_stability',
    'other'
]

# Statistical significance thresholds
STATISTICAL_THRESHOLDS = {
    'fisher_exact': 0.05,
    'chi_square': 0.05,
    'vampp_score': 0.5
}

# VAMPP-score thresholds for PP3/BP4 criteria
VAMPP_SCORE_THRESHOLDS = {
    'pp3': {
        'common_variants': 0.8,     # gnomAD AF > 1e-3
        'moderate_rare': 0.65,      # gnomAD AF > 1e-4
        'very_rare': 0.5,           # gnomAD AF <= 1e-4
    },
    'bp4': {
        'common_variants': 0.2,     # gnomAD AF > 1e-3
        'moderate_rare': 0.3,       # gnomAD AF > 1e-4
        'very_rare': 0.4,           # gnomAD AF <= 1e-4
    },
    'prevalence_adjustment': {
        'common_disease': 1e-3,     # Prevalence > 1e-3
        'rare_disease': 1e-5,       # Prevalence < 1e-5
        'adjustment_factor': 0.1    # Threshold adjustment
    }
}

# API endpoints
API_ENDPOINTS = {
    'ensembl': 'https://rest.ensembl.org',
    'clinvar': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils',
    'varsome': 'https://varsome.com/variant/hg38'
}

# File paths and output settings
OUTPUT_SETTINGS = {
    'report_filename': 'variant_classification_report.txt',
    'log_filename': 'variant_classification_log.txt',
    'cache_filename': 'api_cache.json',
    'max_cache_age_hours': 24
}

# Validation patterns
VALIDATION_PATTERNS = {
    'hgvs_cdna': r'^c\.\d+[ATCG]>[ATCG]$|^c\.\d+[_\-]\d+del|^c\.\d+[_\-]\d+ins|^c\.\d+[_\-]\d+dup',
    'hgvs_protein': r'^p\.[A-Z]{3}\d+[A-Z]{3}$|^p\.[A-Z]{3}\d+\*$|^p\.\*\d+[A-Z]{3}$|^p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}$|^p\.[A-Z][a-z]{2}\d+\*$',
    'chromosome': r'^(chr)?(1[0-9]|2[0-2]|[1-9]|X|Y|MT?)$',
    'gene_symbol': r'^[A-Z][A-Z0-9\-]*$',
    'allele_frequency': r'^0(\.\d+)?$|^1(\.0+)?$'
}

# Aliases for user input
ALIASES = {
    # Gene symbols
    'brca': 'BRCA1',
    'brca1': 'BRCA1',
    'brca2': 'BRCA2',
    'tp53': 'TP53',
    'apc': 'APC',
    'mlh1': 'MLH1',
    'msh2': 'MSH2',
    'msh6': 'MSH6',
    'pms2': 'PMS2',
    
    # Chromosomes
    'chr1': '1', 'chr2': '2', 'chr3': '3', 'chr4': '4', 'chr5': '5',
    'chr6': '6', 'chr7': '7', 'chr8': '8', 'chr9': '9', 'chr10': '10',
    'chr11': '11', 'chr12': '12', 'chr13': '13', 'chr14': '14', 'chr15': '15',
    'chr16': '16', 'chr17': '17', 'chr18': '18', 'chr19': '19', 'chr20': '20',
    'chr21': '21', 'chr22': '22', 'chrx': 'X', 'chry': 'Y',
    
    # Variant types
    'del': 'deletion',
    'ins': 'insertion',
    'dup': 'duplication',
    'snv': 'substitution',
    'indel': 'indel',
    'cnv': 'copy_number_variant',
    
    # Inheritance patterns
    'ad': 'AD',
    'ar': 'AR',
    'xl': 'XL',
    'xld': 'XLD',
    'xlr': 'XLR',
    'dom': 'AD',
    'rec': 'AR',
    'dominant': 'AD',
    'recessive': 'AR',
    'x-linked': 'XL',
    
    # Zygosity
    'het': 'heterozygous',
    'hom': 'homozygous',
    'hem': 'hemizygous',
    'hetero': 'heterozygous',
    'homo': 'homozygous',
    'hemi': 'hemizygous',
    
    # Yes/No responses
    'y': 'yes',
    'n': 'no',
    'true': 'yes',
    'false': 'no',
    '1': 'yes',
    '0': 'no',
    
    # Common abbreviations
    'na': 'not_available',
    'n/a': 'not_available',
    'unknown': 'not_available',
    'unk': 'not_available',
    'not specified': 'not_available',
    'ns': 'not_available'
}

# Test mode sample data
TEST_MODE_DATA = {
    'basic_info': {
        'gene': 'BRCA1',
        'chromosome': '17',
        'position': '41276045',
        'ref_allele': 'C',
        'alt_allele': 'T',
        'cdna_change': 'c.5266dupC',
        'protein_change': 'p.Gln1756Profs*74',
        'variant_type': 'frameshift',
        'consequence': 'frameshift_variant',
        'clinvar_status': 'pathogenic'
    },
    'population_data': {
        'gnomad_af': 0.0001,
        'gnomad_af_popmax': 0.0002,
        'gnomad_hom_count': 0,
        'gnomad_het_count': 3
    },
    'insilico_data': {
        'revel': 0.85,
        'cadd_phred': 24.1,
        'gerp_pp': 4.4,
        'phylop_vert_ranked': 0.8,
        'phylop_mamm_ranked': 0.75,
        'alphamissense': 0.7,
        'metarnn': 0.8,
        'clinpred': 0.9,
        'bayesdel_addaf': 0.3,
        'sift': 0.02,
        'polyphen2': 0.95,
        'mutationtaster': 0.88,
        'fathmm': -3.5,
        # High priority predictors
        'vest4': 0.78,
        'primateai': 0.82,
        'esm1b': 0.7,
        'provean': -4.2,
        'mmsplice': 1.5
    },
    'genetic_data': {
        'inheritance': 'AD',
        'zygosity': 'heterozygous',
        'phase': 'unknown'
    },
    'functional_data': {
        'segregation': 'cosegregates',
        'denovo': 'confirmed',
        'functional_studies': 'damaging'
    }
}

# Color codes for terminal output
COLORS = {
    'RED': '\033[91m',
    'GREEN': '\033[92m',
    'YELLOW': '\033[93m',
    'BLUE': '\033[94m',
    'MAGENTA': '\033[95m',
    'CYAN': '\033[96m',
    'WHITE': '\033[97m',
    'BOLD': '\033[1m',
    'UNDERLINE': '\033[4m',
    'END': '\033[0m'
}

# Colorama colors for cross-platform compatibility
try:
    from colorama import Fore, Back, Style
    COLORAMA_COLORS = {
        'RED': Fore.RED,
        'GREEN': Fore.GREEN,
        'YELLOW': Fore.YELLOW,
        'BLUE': Fore.BLUE,
        'MAGENTA': Fore.MAGENTA,
        'CYAN': Fore.CYAN,
        'WHITE': Fore.WHITE,
        'BOLD': Style.BRIGHT,
        'DIM': Style.DIM,
        'RESET': Style.RESET_ALL
    }
except ImportError:
    COLORAMA_COLORS = {key: '' for key in ['RED', 'GREEN', 'YELLOW', 'BLUE', 'MAGENTA', 'CYAN', 'WHITE', 'BOLD', 'DIM', 'RESET']}

# Color formatting functions
def get_colored_message(message_type, message):
    """Get a colored message based on type."""
    colors = {
        'error': f'{COLORAMA_COLORS["RED"]}❌ {message}{COLORAMA_COLORS["RESET"]}',
        'success': f'{COLORAMA_COLORS["GREEN"]}✅ {message}{COLORAMA_COLORS["RESET"]}',
        'warning': f'{COLORAMA_COLORS["YELLOW"]}⚠️ {message}{COLORAMA_COLORS["RESET"]}',
        'info': f'{COLORAMA_COLORS["CYAN"]}ℹ️ {message}{COLORAMA_COLORS["RESET"]}',
        'header': f'{COLORAMA_COLORS["CYAN"]}{COLORAMA_COLORS["BOLD"]}{message}{COLORAMA_COLORS["RESET"]}'
    }
    return colors.get(message_type, message)

# Error messages with color formatting
ERROR_MESSAGES = {
    'invalid_input': 'Invalid input provided. Please check your entry and try again.',
    'api_error': 'API request failed. Please check your internet connection.',
    'validation_error': 'Input validation failed. Please ensure all fields are correctly formatted.',
    'classification_error': 'Unable to classify variant with provided evidence.',
    'file_error': 'Error reading or writing file.',
    'network_error': 'Network connection error. Some features may be unavailable.'
}

# Success messages with color formatting
SUCCESS_MESSAGES = {
    'classification_complete': 'Variant classification completed successfully.',
    'report_generated': 'Classification report generated successfully.',
    'data_validated': 'All input data validated successfully.',
    'api_data_retrieved': 'External data retrieved successfully.'
}

# Warning messages with color formatting
WARNING_MESSAGES = {
    'missing_data': 'Some data fields are missing. Classification may be less accurate.',
    'low_confidence': 'Low confidence classification. Consider additional evidence.',
    'conflicting_evidence': 'Conflicting evidence detected. Manual review recommended.',
    'api_fallback': 'API unavailable, using cached or manual data.',
    'gene_specific_threshold': 'Using gene-specific thresholds for population frequency.'
}

# Info messages with color formatting
INFO_MESSAGES = {
    'vampp_score_calculated': 'VAMPP-score calculated from in silico predictors.',
    'fisher_test_performed': 'Fisher\'s exact test performed for population analysis.',
    'cache_used': 'Using cached API data.',
    'test_mode_active': 'Test mode active - using sample data.',
    'manual_override': 'Manual override applied.'
}

# Classification result colors
CLASSIFICATION_COLORS = {
    'Pathogenic': ('RED', 'BOLD'),
    'Likely Pathogenic': ('RED', None),
    'Likely Benign': ('GREEN', None),
    'Benign': ('GREEN', 'BOLD'),
    'Variant of Uncertain Significance': ('YELLOW', 'BOLD'),
    'VUS': ('YELLOW', 'BOLD')
}

# Section headers with color formatting
SECTION_HEADERS = {
    'welcome': '🧬 ACMG Variant Classification Assistant',
    'basic_info': '📋 Basic Variant Information',
    'population_data': '👥 Population Data',
    'insilico_data': '🔬 In Silico Prediction Scores',
    'genetic_data': '🧬 Genetic Information',
    'functional_data': '⚡ Functional Studies',
    'classification_result': '🎯 Classification Result',
    'evidence_summary': '📊 Evidence Summary',
    'report_generation': '📄 Report Generation'
}

# Version information
VERSION_INFO = {
    'version': '1.0.0',
    'release_date': '2025-07-06',
    'author': 'Can Sevilmiş',
    'license': 'MIT',
    'guidelines': ['ACMG/AMP 2015', 'ACMG/AMP 2023']
}
