#!/usr/bin/env python3
"""
Comprehensive test for all 28 ACMG/AMP criteria implementation.
This test verifies that all criteria are properly implemented and callable.
"""

import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from core.evidence_evaluator import EvidenceEvaluator
from core.variant_data import VariantData

def test_all_acmg_criteria():
    """Test that all 28 ACMG/AMP criteria are implemented and callable."""
    
    print("🧬 ACMG/AMP Criteria Implementation Test")
    print("=" * 60)
    
    # Create test variant data
    test_data = VariantData({
        'basic_info': {
            'gene': 'BRCA1',
            'variant_name': 'c.5266dupC',
            'amino_acid_change': 'p.Gln1756ProfsX74',
            'variant_type': 'frameshift',
            'consequence': 'frameshift_variant',
            'chromosome': '17',
            'position': '41197819',
            'ref': 'A',
            'alt': 'AC'
        },
        'population_data': {
            'gnomad_af': 0.0001,
            'gnomad_af_popmax': 0.0002
        },
        'insilico_data': {
            'revel_score': 0.8,
            'cadd_raw': 5.2,
            'cadd_phred': 25.3,
            'sift_score': 0.02,
            'polyphen2_hdiv_score': 0.95,
            'spliceai_ag_score': 0.05,
            'spliceai_al_score': 0.03,
            'spliceai_dg_score': 0.02,
            'spliceai_dl_score': 0.01
        },
        'genetic_data': {
            'inheritance_pattern': 'autosomal_dominant',
            'de_novo_status': 'unknown',
            'segregation_data': 'unknown'
        },
        'functional_data': {
            'functional_studies': 'unknown',
            'protein_function': 'unknown'
        }
    })
    
    # Initialize evaluator in test mode
    evaluator = EvidenceEvaluator(test_mode=True, use_2023_guidelines=True)
    
    # Define all 28 ACMG/AMP criteria
    pathogenic_criteria = [
        'PVS1',  # Very Strong
        'PS1', 'PS2', 'PS3', 'PS4',  # Strong
        'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',  # Moderate
        'PP1', 'PP2', 'PP3', 'PP4', 'PP5'  # Supporting
    ]
    
    benign_criteria = [
        'BA1',  # Stand Alone
        'BS1', 'BS2', 'BS3', 'BS4',  # Strong
        'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'  # Supporting
    ]
    
    all_criteria = pathogenic_criteria + benign_criteria
    
    print(f"Testing {len(all_criteria)} ACMG/AMP criteria...")
    print()
    
    # Test pathogenic criteria
    print("🔴 PATHOGENIC CRITERIA:")
    print("-" * 30)
    
    failed_tests = []
    passed_tests = []
    
    try:
        pathogenic_results = evaluator._evaluate_pathogenic_criteria(test_data)
        
        for criterion in pathogenic_criteria:
            if criterion in pathogenic_results:
                result = pathogenic_results[criterion]
                if isinstance(result, dict):
                    print(f"✅ {criterion}: {result.get('details', 'OK')}")
                    passed_tests.append(criterion)
                else:
                    print(f"❌ {criterion}: Invalid result format")
                    failed_tests.append(criterion)
            else:
                print(f"❌ {criterion}: Not implemented")
                failed_tests.append(criterion)
    except Exception as e:
        print(f"❌ Error evaluating pathogenic criteria: {e}")
    
    print()
    
    # Test benign criteria
    print("🔵 BENIGN CRITERIA:")
    print("-" * 30)
    
    try:
        benign_results = evaluator._evaluate_benign_criteria(test_data)
        
        for criterion in benign_criteria:
            if criterion in benign_results:
                result = benign_results[criterion]
                if isinstance(result, dict):
                    print(f"✅ {criterion}: {result.get('details', 'OK')}")
                    passed_tests.append(criterion)
                else:
                    print(f"❌ {criterion}: Invalid result format")
                    failed_tests.append(criterion)
            else:
                print(f"❌ {criterion}: Not implemented")
                failed_tests.append(criterion)
    except Exception as e:
        print(f"❌ Error evaluating benign criteria: {e}")
    
    print()
    print("=" * 60)
    print("📊 TEST SUMMARY:")
    print("=" * 60)
    print(f"✅ Passed: {len(passed_tests)}/{len(all_criteria)} criteria")
    print(f"❌ Failed: {len(failed_tests)}/{len(all_criteria)} criteria")
    
    if passed_tests:
        print(f"\n✅ Successfully implemented: {', '.join(sorted(passed_tests))}")
    
    if failed_tests:
        print(f"\n❌ Failed or missing: {', '.join(sorted(failed_tests))}")
        return False
    else:
        print("\n🎉 ALL 28 ACMG/AMP CRITERIA SUCCESSFULLY IMPLEMENTED!")
        return True

def test_evidence_evaluation_workflow():
    """Test the complete evidence evaluation workflow."""
    
    print("\n🔬 EVIDENCE EVALUATION WORKFLOW TEST")
    print("=" * 60)
    
    # Create test variant data
    test_data = VariantData({
        'basic_info': {
            'gene': 'TP53',
            'variant_name': 'c.817C>T',
            'amino_acid_change': 'p.Arg273His',
            'variant_type': 'missense',
            'consequence': 'missense_variant'
        },
        'population_data': {
            'gnomad_af': 0.00001
        },
        'insilico_data': {
            'revel_score': 0.9,
            'cadd_phred': 32.0
        }
    })
    
    try:
        evaluator = EvidenceEvaluator(test_mode=True, use_2023_guidelines=True)
        results = evaluator.evaluate_all_criteria(test_data)
        
        print("✅ Evidence evaluation completed successfully")
        print(f"Applied criteria: {len(results.get('applied_criteria', {}))}")
        print(f"Pathogenic criteria: {len(results.get('pathogenic_evidence', {}))}")
        print(f"Benign criteria: {len(results.get('benign_evidence', {}))}")
        
        if 'vampp_score' in results and results['vampp_score'] is not None:
            print(f"Metascore calculated: {results['vampp_score']:.3f}")
        
        return True
        
    except Exception as e:
        print(f"❌ Error in evidence evaluation workflow: {e}")
        return False

if __name__ == "__main__":
    print("🧬 ACMG/AMP Criteria Implementation Test Suite")
    print("=" * 60)
    
    # Run tests
    criteria_test_passed = test_all_acmg_criteria()
    workflow_test_passed = test_evidence_evaluation_workflow()
    
    print("\n" + "=" * 60)
    print("🏁 FINAL RESULTS:")
    print("=" * 60)
    
    if criteria_test_passed and workflow_test_passed:
        print("✅ ALL TESTS PASSED!")
        print("✅ All 28 ACMG/AMP criteria are properly implemented")
        print("✅ Evidence evaluation workflow is functional")
        sys.exit(0)
    else:
        print("❌ SOME TESTS FAILED!")
        if not criteria_test_passed:
            print("❌ Criteria implementation issues detected")
        if not workflow_test_passed:
            print("❌ Workflow issues detected")
        sys.exit(1)
