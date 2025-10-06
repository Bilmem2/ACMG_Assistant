"""
ACMG 2023 Mode Test
===================

Quick test to verify ACMG 2023 mode works correctly with stricter thresholds.
"""

import sys
import os

# Set UTF-8 encoding for console output
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from core.variant_data import VariantData
from core.evidence_evaluator import EvidenceEvaluator
from core.acmg_classifier import ACMGClassifier


def test_ps4_2023_mode():
    """Test PS4 with ACMG 2023 stricter thresholds."""
    print("\n" + "="*80)
    print("TEST 1: PS4 with ACMG 2023 Mode (OR≥5.0, cases≥10, controls≥2000)")
    print("="*80)
    
    # Test case that meets ACMG 2023 criteria
    variant = VariantData()
    variant.basic_info = {'gene': 'APOE', 'variant_type': 'missense'}
    variant.population_data = {
        'case_control_data': {
            'cases': 250,
            'case_allele_count': 75,
            'controls': 5000,
            'control_allele_count': 25,
            'odds_ratio': 6.5,  # ACMG 2023: ≥5.0 required
            'p_value': 0.00001,
            'significant': True
        }
    }
    
    # Test with ACMG 2015 (should pass with OR=6.5)
    evaluator_2015 = EvidenceEvaluator(use_2023_guidelines=False, test_mode=True)
    result_2015 = evaluator_2015.evaluate_all_criteria(variant)
    ps4_2015 = result_2015['pathogenic_criteria'].get('PS4', {})
    
    print("\n📋 ACMG 2015 Mode:")
    print(f"  - PS4 Applied: {ps4_2015.get('applies', False)}")
    print(f"  - Details: {ps4_2015.get('details', 'N/A')}")
    
    # Test with ACMG 2023
    evaluator_2023 = EvidenceEvaluator(use_2023_guidelines=True, test_mode=True)
    result_2023 = evaluator_2023.evaluate_all_criteria(variant)
    ps4_2023 = result_2023['pathogenic_criteria'].get('PS4', {})
    
    print("\n📋 ACMG 2023 Mode:")
    print(f"  - PS4 Applied: {ps4_2023.get('applies', False)}")
    print(f"  - Details: {ps4_2023.get('details', 'N/A')}")
    
    # Verify
    success = ps4_2015.get('applies') and ps4_2023.get('applies')
    print(f"\n{'✅ PASS' if success else '❌ FAIL'}: Both modes should apply PS4 for OR=6.5")
    
    return success


def test_ps4_2023_fails_low_or():
    """Test PS4 with ACMG 2023 - should fail with OR=3.5."""
    print("\n" + "="*80)
    print("TEST 2: PS4 with ACMG 2023 Mode - Low OR (should fail)")
    print("="*80)
    
    variant = VariantData()
    variant.basic_info = {'gene': 'APOE', 'variant_type': 'missense'}
    variant.population_data = {
        'case_control_data': {
            'cases': 150,
            'controls': 2500,
            'odds_ratio': 3.5,  # ACMG 2023: needs ≥5.0 (FAIL)
            'p_value': 0.00001,
            'significant': True
        }
    }
    
    # Test with ACMG 2015 (should pass with OR=3.5)
    evaluator_2015 = EvidenceEvaluator(use_2023_guidelines=False, test_mode=True)
    result_2015 = evaluator_2015.evaluate_all_criteria(variant)
    ps4_2015 = result_2015['pathogenic_criteria'].get('PS4', {})
    
    print("\n📋 ACMG 2015 Mode:")
    print(f"  - PS4 Applied: {ps4_2015.get('applies', False)}")
    print(f"  - Expected: True (OR≥2.0)")
    
    # Test with ACMG 2023 (should FAIL with OR=3.5)
    evaluator_2023 = EvidenceEvaluator(use_2023_guidelines=True, test_mode=True)
    result_2023 = evaluator_2023.evaluate_all_criteria(variant)
    ps4_2023 = result_2023['pathogenic_criteria'].get('PS4', {})
    
    print("\n📋 ACMG 2023 Mode:")
    print(f"  - PS4 Applied: {ps4_2023.get('applies', False)}")
    print(f"  - Expected: False (OR≥5.0 required)")
    print(f"  - Details: {ps4_2023.get('details', 'N/A')}")
    
    # Verify
    success = ps4_2015.get('applies') and not ps4_2023.get('applies')
    print(f"\n{'✅ PASS' if success else '❌ FAIL'}: 2015 should pass, 2023 should fail for OR=3.5")
    
    return success


def test_pp1_2023_strength_modifiers():
    """Test PP1 with ACMG 2023 LOD-based strength modifiers."""
    print("\n" + "="*80)
    print("TEST 3: PP1 with ACMG 2023 LOD-based Strength Modifiers")
    print("="*80)
    
    # Test LOD=5.5 (should be PP1_Strong in 2023)
    variant = VariantData()
    variant.basic_info = {'gene': 'BRCA1', 'variant_type': 'missense'}
    variant.genetic_data = {
        'segregation_data': {
            'segregates': True,
            'lod_score': 5.5,
            'families': 5,
            'informative_meioses': 8
        }
    }
    
    # ACMG 2015
    evaluator_2015 = EvidenceEvaluator(use_2023_guidelines=False, test_mode=True)
    result_2015 = evaluator_2015.evaluate_all_criteria(variant)
    pp1_2015 = result_2015['pathogenic_criteria'].get('PP1', {})
    
    print("\n📋 ACMG 2015 Mode (LOD=5.5):")
    print(f"  - PP1 Applied: {pp1_2015.get('applies', False)}")
    print(f"  - Strength: {pp1_2015.get('strength', 'N/A')}")
    print(f"  - Expected: Moderate (simple threshold)")
    
    # ACMG 2023
    evaluator_2023 = EvidenceEvaluator(use_2023_guidelines=True, test_mode=True)
    result_2023 = evaluator_2023.evaluate_all_criteria(variant)
    pp1_2023 = result_2023['pathogenic_criteria'].get('PP1', {})
    
    print("\n📋 ACMG 2023 Mode (LOD=5.5):")
    print(f"  - PP1 Applied: {pp1_2023.get('applies', False)}")
    print(f"  - Strength: {pp1_2023.get('strength', 'N/A')}")
    print(f"  - Expected: Strong (LOD≥5.0)")
    print(f"  - Details: {pp1_2023.get('details', 'N/A')}")
    
    # Verify
    success = (
        pp1_2015.get('strength') == 'Moderate' and
        pp1_2023.get('strength') == 'Strong'
    )
    print(f"\n{'✅ PASS' if success else '❌ FAIL'}: 2023 mode should upgrade to Strong for LOD≥5.0")
    
    return success


def main():
    """Run ACMG 2023 mode tests."""
    print("""
================================================================================
                                                                            
                      ACMG 2023 MODE TEST SUITE                           
                                                                            
  Validates that ACMG 2023 guidelines work correctly alongside 2015 mode.  
  Tests stricter thresholds and new strength modifiers.                    
                                                                            
================================================================================
    """)
    
    results = []
    results.append(test_ps4_2023_mode())
    results.append(test_ps4_2023_fails_low_or())
    results.append(test_pp1_2023_strength_modifiers())
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    total = len(results)
    passed = sum(results)
    failed = total - passed
    
    print(f"\nTotal Tests: {total}")
    print(f"Passed: {passed} ({100*passed/total:.1f}%)")
    print(f"Failed: {failed} ({100*failed/total:.1f}%)")
    
    if failed == 0:
        print("\n✅ ALL ACMG 2023 MODE TESTS PASSED!")
    else:
        print(f"\n⚠️  {failed} test(s) failed")
    print("="*80)


if __name__ == '__main__':
    main()
