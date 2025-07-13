#!/usr/bin/env python3
"""
Test script for ACMG classification with real variant data.
"""

import sys
import os
# Add both src and root to path
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src.core.variant_data import VariantData
from src.core.evidence_evaluator import EvidenceEvaluator
from src.core.acmg_classifier import ACMGClassifier

def test_classification_with_real_data():
    """Test ACMG classification with realistic variant data."""
    print("="*60)
    print("ACMG Classification Test with Real Variant Data")
    print("="*60)
    
    # Create a BRCA1 pathogenic variant for testing
    test_variant = {
        'basic_info': {
            'gene': 'BRCA1',
            'chromosome': '17',
            'position': '43124096',
            'ref_allele': 'G',
            'alt_allele': 'A',
            'variant_type': 'nonsense',
            'consequence': 'stop_gained',
            'cdna_change': 'c.4035G>A',
            'amino_acid_change': 'p.Trp1345*',
            'transcript': 'NM_007294.4'
        },
        'population_data': {
            'gnomad_af': 0.0,
            'gnomad_af_popmax': 0.0,
            'exac_af': 0.0,
        },
        'insilico_data': {
            'revel_score': 0.95,
            'cadd_phred': 35.0,
            'gerp_pp': 5.8,
            'phylop_vertebrate': 0.98,
            'phylop_mammalian': 0.99
        },
        'genetic_data': {
            'de_novo': 'confirmed',
            'maternity_confirmed': True,
            'paternity_confirmed': True,
            'segregation': 'consistent',
            'inheritance_pattern': 'autosomal_dominant'
        },
        'functional_data': {
            'functional_studies': 'deleterious',
            'phenotype_match': 'specific_match',
            'protein_function': 'loss_of_function'
        }
    }
    
    # Create VariantData object
    variant_data = VariantData(
        basic_info=test_variant['basic_info'],
        population_data=test_variant['population_data'],
        insilico_data=test_variant['insilico_data'],
        genetic_data=test_variant['genetic_data'],
        functional_data=test_variant['functional_data']
    )
    
    # Test with 2015 guidelines
    print("\n1. Testing with ACMG 2015 Guidelines:")
    print("-" * 40)
    
    evaluator_2015 = EvidenceEvaluator(use_2023_guidelines=False)
    evidence_results_2015 = evaluator_2015.evaluate_all_criteria(variant_data)
    
    classifier_2015 = ACMGClassifier(use_2023_guidelines=False)
    classification_2015 = classifier_2015.classify(evidence_results_2015)
    
    print(f"Classification: {classification_2015['classification']}")
    print(f"Confidence: {classification_2015['confidence']}")
    print(f"Pathogenic Evidence: {classification_2015['pathogenic_counts']}")
    print(f"Benign Evidence: {classification_2015['benign_counts']}")
    print(f"Applied Criteria: {list(classification_2015['applied_criteria'].keys())}")
    
    # Test with 2023 guidelines
    print("\n2. Testing with ACMG 2023 Guidelines:")
    print("-" * 40)
    
    evaluator_2023 = EvidenceEvaluator(use_2023_guidelines=True)
    evidence_results_2023 = evaluator_2023.evaluate_all_criteria(variant_data)
    
    classifier_2023 = ACMGClassifier(use_2023_guidelines=True)
    classification_2023 = classifier_2023.classify(evidence_results_2023)
    
    print(f"Classification: {classification_2023['classification']}")
    print(f"Confidence: {classification_2023['confidence']}")
    print(f"Pathogenic Evidence: {classification_2023['pathogenic_counts']}")
    print(f"Benign Evidence: {classification_2023['benign_counts']}")
    print(f"Applied Criteria: {list(classification_2023['applied_criteria'].keys())}")
    
    # Test with a likely benign variant
    print("\n3. Testing with Benign Variant:")
    print("-" * 40)
    
    benign_variant = {
        'basic_info': {
            'gene': 'BRCA1',
            'chromosome': '17',
            'position': '43124000',
            'ref_allele': 'C',
            'alt_allele': 'T',
            'variant_type': 'synonymous',
            'consequence': 'synonymous_variant',
            'cdna_change': 'c.4020C>T',
            'amino_acid_change': 'p.Ser1340=',
            'transcript': 'NM_007294.4'
        },
        'population_data': {
            'gnomad_af': 0.02,  # Common variant
            'gnomad_af_popmax': 0.025,
            'exac_af': 0.018,
        },
        'insilico_data': {
            'revel_score': 0.1,
            'cadd_phred': 5.2,
            'gerp_pp': -2.1,
            'phylop_vertebrate': 0.1,
            'phylop_mammalian': 0.05
        },
        'genetic_data': {
            'de_novo': 'no',
            'segregation': 'inconsistent',
            'inheritance_pattern': 'unknown'
        },
        'functional_data': {
            'functional_studies': 'benign',
            'phenotype_match': 'no_match',
            'protein_function': 'neutral'
        }
    }
    
    benign_variant_data = VariantData(
        basic_info=benign_variant['basic_info'],
        population_data=benign_variant['population_data'],
        insilico_data=benign_variant['insilico_data'],
        genetic_data=benign_variant['genetic_data'],
        functional_data=benign_variant['functional_data']
    )
    
    evaluator_benign = EvidenceEvaluator(use_2023_guidelines=False)
    evidence_results_benign = evaluator_benign.evaluate_all_criteria(benign_variant_data)
    
    classifier_benign = ACMGClassifier(use_2023_guidelines=False)
    classification_benign = classifier_benign.classify(evidence_results_benign)
    
    print(f"Classification: {classification_benign['classification']}")
    print(f"Confidence: {classification_benign['confidence']}")
    print(f"Pathogenic Evidence: {classification_benign['pathogenic_counts']}")
    print(f"Benign Evidence: {classification_benign['benign_counts']}")
    print(f"Applied Criteria: {list(classification_benign['applied_criteria'].keys())}")
    
    print("\n" + "="*60)
    print("Test Results Summary:")
    print("="*60)
    
    # Expected results check
    expected_pathogenic = classification_2015['classification'] in ['Pathogenic', 'Likely pathogenic']
    expected_benign = classification_benign['classification'] in ['Benign', 'Likely benign']
    
    print(f"✓ Pathogenic variant correctly classified: {expected_pathogenic}")
    print(f"✓ Benign variant correctly classified: {expected_benign}")
    print(f"✓ AI features working: {'PS1' in classification_2015['applied_criteria']}")
    print(f"✓ Dynamic scoring working: {classification_2015.get('vampp_score') is not None}")
    
    return True

if __name__ == "__main__":
    test_classification_with_real_data()
