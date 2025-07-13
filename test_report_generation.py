#!/usr/bin/env python3
"""
Test Report Generation
=====================

Test script to debug report generation issues an    for predictor in readme_predictors:
        # Check if predictors are defined in constants
        variations = [
            predictor,
            f"{predictor}_score",
            f"{predictor}_ranked",
            f"{predictor}_ag", f"{predictor}_al", f"{predictor}_dg", f"{predictor}_dl",  # For SpliceAI
            f"{predictor}_vert", f"{predictor}_mamm", f"{predictor}_prim",  # For PhyloP/phastCons
            f"{predictor}_addaf", f"{predictor}_noaf"  # For BayesDel
        ]
        
        found = False
        for var in variations:
            if var in INSILICO_THRESHOLDS:
                defined_predictors.append(var)
                found = True
                break
                
        if not found:
            undefined_predictors.append(predictor)ctor functionality.
"""

import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from core.variant_data import VariantData
from core.acmg_classifier import ACMGClassifier
from core.evidence_evaluator import EvidenceEvaluator
from utils.report_generator import ReportGenerator
from config.constants import VERSION_INFO, INSILICO_THRESHOLDS, INSILICO_WEIGHTS

def test_basic_report_generation():
    """Test basic report generation with minimal data."""
    print("🧪 Testing Basic Report Generation...")
    
    # Create minimal variant data
    variant_data = VariantData(
        basic_info={
            'gene': 'TEST_GENE',
            'chromosome': '1',
            'position': '12345',
            'ref_allele': 'A',
            'alt_allele': 'T',
            'variant_type': 'missense',
            'consequence': 'missense_variant',
            'cdna_change': 'c.123A>T',
            'protein_change': 'p.Lys41Met'
        },
        population_data={
            'gnomad_af': 0.000001,
            'exac_af': 0.000002,
            'dbsnp_id': 'rs123456789'
        },
        insilico_data={
            'revel': 0.8,
            'cadd_phred': 28.5,
            'alphamissense': 0.75
        },
        genetic_data={
            'inheritance_mode': 'autosomal_dominant',
            'segregation': 'cosegregates'
        },
        functional_data={
            'functional_studies': 'none'
        }
    )
    
    # Initialize components
    evidence_evaluator = EvidenceEvaluator()
    classifier = ACMGClassifier()
    report_generator = ReportGenerator()
    
    try:
        # Evaluate evidence
        print("  - Evaluating evidence...")
        evidence_results = evidence_evaluator.evaluate_all_criteria(variant_data)
        
        # Classify variant
        print("  - Classifying variant...")
        classification_result = classifier.classify(evidence_results)
        
        # Generate report
        print("  - Generating report...")
        report_path = report_generator.generate_report(
            variant_data, evidence_results, classification_result
        )
        
        print(f"✅ Report generated successfully: {report_path}")
        
        # Check if file exists and has content
        if os.path.exists(report_path):
            with open(report_path, 'r', encoding='utf-8') as f:
                content = f.read()
                print(f"   Report size: {len(content)} characters")
                if len(content) > 0:
                    print("   Report contains content ✅")
                else:
                    print("   ❌ Report is empty!")
        else:
            print("   ❌ Report file was not created!")
            
        return True
        
    except Exception as e:
        print(f"❌ Error during report generation: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def test_predictor_coverage():
    """Test coverage of all predictors mentioned in README."""
    print("\n🧪 Testing Predictor Coverage...")
    
    # README'de bahsedilen tüm predictor'lar
    readme_predictors = [
        # Primary Metascores
        'revel', 'cadd_phred', 'alphamissense', 'metarnn', 'clinpred', 
        'bayesdel', 'metasvm', 'metalr',
        
        # Missense Predictors
        'sift', 'polyphen2', 'provean', 'vest4', 'esm1b', 'mutationtaster',
        'fathmm', 'mutationassessor', 'mutpred', 'lrt',
        
        # Conservation & Evolutionary
        'phylop', 'phastcons', 'gerp_pp', 'siphy',
        
        # Splice Site Prediction
        'spliceai', 'ada', 'rf',
        
        # Functional & Regulatory
        'fitcons'
    ]
    
    # Check if predictors are defined in constants
    defined_predictors = []
    undefined_predictors = []
    
    for predictor in readme_predictors:
        # Check in INSILICO_THRESHOLDS
        variations = [
            predictor,
            f"{predictor}_score",
            f"{predictor}_ranked",
            f"{predictor}_ag", f"{predictor}_al", f"{predictor}_dg", f"{predictor}_dl"  # For SpliceAI
        ]
        
        found = False
        for var in variations:
            if var in INSILICO_THRESHOLDS:
                defined_predictors.append(var)
                found = True
                break
                
        if not found:
            undefined_predictors.append(predictor)
    
    print(f"✅ Defined predictors ({len(defined_predictors)}):")
    for pred in sorted(defined_predictors):
        print(f"   - {pred}")
        
    if undefined_predictors:
        print(f"\n❌ Undefined predictors ({len(undefined_predictors)}):")
        for pred in sorted(undefined_predictors):
            print(f"   - {pred}")
    else:
        print(f"\n✅ All README predictors are defined!")
        
    return len(undefined_predictors) == 0

def test_predictor_weights():
    """Test predictor weights and VAMPP score calculation."""
    print("\n🧪 Testing Predictor Weights...")
    
    # Check if all predictors in thresholds have weights
    predictors_with_thresholds = set(INSILICO_THRESHOLDS.keys())
    predictors_with_weights = set(INSILICO_WEIGHTS.keys())
    
    missing_weights = predictors_with_thresholds - predictors_with_weights
    extra_weights = predictors_with_weights - predictors_with_thresholds
    
    if missing_weights:
        print(f"❌ Predictors with thresholds but no weights ({len(missing_weights)}):")
        for pred in sorted(missing_weights):
            print(f"   - {pred}")
    
    if extra_weights:
        print(f"⚠️  Predictors with weights but no thresholds ({len(extra_weights)}):")
        for pred in sorted(extra_weights):
            print(f"   - {pred}")
    
    if not missing_weights and not extra_weights:
        print("✅ All predictors have both thresholds and weights!")
    
    return len(missing_weights) == 0

def main():
    """Run all tests."""
    print("="*80)
    print("🧬 ACMG ASSISTANT - REPORT GENERATION & PREDICTOR TESTS")
    print("="*80)
    
    success = True
    
    # Test 1: Basic report generation
    success &= test_basic_report_generation()
    
    # Test 2: Predictor coverage
    success &= test_predictor_coverage()
    
    # Test 3: Predictor weights
    success &= test_predictor_weights()
    
    print("\n" + "="*80)
    if success:
        print("✅ ALL TESTS PASSED!")
    else:
        print("❌ SOME TESTS FAILED!")
    print("="*80)
    
    return success

if __name__ == "__main__":
    main()
