"""
Test Scenarios Configuration
===========================

Contains test data scenarios for ACMG variant classification testing.

Author: Can Sevilmiş
License: MIT License
"""

# Test mode data - simplified version for now
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
    }
}

# Legacy TEST_MODE_DATA for backward compatibility (points to first scenario)
TEST_MODE_DATA = TEST_SCENARIOS['1_pvs1_nonsense']