#!/usr/bin/env python3
"""
Test script for CAMTA1 variant: c.319A>G (p.Met107Val)
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from core.variant_data import VariantData
from core.acmg_classifier import ACMGClassifier
from core.evidence_evaluator import EvidenceEvaluator
from utils.report_generator import ReportGenerator

def test_camta1_variant():
    """Test the CAMTA1 c.319A>G variant."""
    
    # Create variant data based on the report
    basic_info = {
        'gene': 'CAMTA1',
        'chromosome': '1',
        'position': 7249507,
        'ref_allele': 'A',
        'alt_allele': 'G',
        'cdna_change': 'c.319A>G',
        'protein_change': 'p.Met107Val',
        'variant_type': 'missense',
        'transcript': 'unknown',
        'vep_consequence': 'unknown'
    }
    
    population_data = {
        'gnomad_af': 1.505e-05,
        'gnomad_af_popmax': None,
        'gnomad_ac': None,
        'gnomad_an': None,
        'gnomad_hom': None,
        'gnomad_subpops': {},
        'disease_prevalence': 1e-06,
        'thousand_genomes_af': None,
        'exac_af': None,
        'esp_af': None
    }
    
    insilico_data = {
        'revel': 0.395,
        'cadd_phred': 24.3,
        'clinpred': 0.349,
        'bayesdel_addaf': 0.7275,
        'alphamissense': 0.3896,
        'mutationtaster': 0.9999,
        'polyphen2_hdiv': 0.6866,
        'sift': 0.041,
        'fathmm_xf': 0.904,
        'mutationassessor': 0.12,
        'provean': 0.3158,
        'mutpred': 0.559,
        'metolr': 0.0946,
        'esm1b': -8.208,
        'lrt': 0.44119,
        'phylop_vertebrates': 0.94,
        'phylop_mammals': 3.0,
        'phylop_primates': 3.0,
        'gerp_rs': 5.63,
        'combined_metascore': 0.3245
    }
    
    genetic_data = {
        'inheritance_pattern': 'ad',
        'zygosity': 'heterozygous',
        'allelic_state': 'unknown',
        'family_history': 'no',
        'segregation_analysis': 'no',
        'de_novo': 'yes',
        'parental_testing': 'unknown',
        'consanguinity': 'unknown',
        'siblings_affected': 'no',
        'notes': ''
    }
    
    functional_data = {
        'functional_test': 'unknown',
        'phenotype_match': 'moderate',
        'paternity_confirmed': None,
        'notes': ''
    }
    
    # Create variant data object
    variant_data = VariantData(
        basic_info=basic_info,
        population_data=population_data,
        insilico_data=insilico_data,
        genetic_data=genetic_data,
        functional_data=functional_data
    )
    
    print("🧬 Testing CAMTA1 c.319A>G (p.Met107Val) variant")
    print("="*60)
    
    # Evaluate evidence
    print("\n🔍 Evaluating evidence...")
    evidence_evaluator = EvidenceEvaluator(use_2023_guidelines=False)
    evidence_results = evidence_evaluator.evaluate_all_criteria(variant_data)
    
    # Classify variant
    print("\n⚖️  Classifying variant...")
    classifier = ACMGClassifier(use_2023_guidelines=False)
    classification_result = classifier.classify(evidence_results)
    
    # Generate report
    print("\n📄 Generating report...")
    report_generator = ReportGenerator()
    report_path = report_generator.generate_report(
        variant_data, evidence_results, classification_result
    )
    
    # Display results
    print("\n" + "="*60)
    print("🎯 CLASSIFICATION RESULTS")
    print("="*60)
    
    print(f"📊 FINAL CLASSIFICATION: {classification_result['classification']}")
    print(f"🔍 CONFIDENCE LEVEL: {classification_result['confidence']}")
    
    print(f"\n📋 APPLIED CRITERIA:")
    for category, criteria in classification_result['applied_criteria'].items():
        if criteria:
            print(f"  {category}: {', '.join(criteria)}")
    
    if classification_result.get('suggestions'):
        print(f"\n💡 SUGGESTIONS:")
        for suggestion in classification_result['suggestions']:
            print(f"  • {suggestion}")
    
    print(f"\n📄 DETAILED REPORT: {report_path}")
    print("="*60)
    
    return classification_result

if __name__ == "__main__":
    result = test_camta1_variant()
    print(f"\n✅ Test completed successfully!")
    print(f"Final classification: {result['classification']}")
