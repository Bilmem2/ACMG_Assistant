#!/usr/bin/env python3
"""
Test script to verify PS2/PM6 logic for de novo variants.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src.core.variant_data import VariantData
from src.core.evidence_evaluator import EvidenceEvaluator

def test_denovo_logic():
    """Test PS2/PM6 logic for different de novo scenarios."""
    
    evaluator = EvidenceEvaluator()
    
    # Test 1: Confirmed de novo with parental testing
    print("Test 1: Confirmed de novo with parental testing")
    print("-" * 50)
    confirmed_variant = VariantData(
        basic_info={'gene': 'TEST'},
        genetic_data={
            'de_novo': 'confirmed',
            'maternity_confirmed': True,
            'paternity_confirmed': True
        },
        population_data={},
        insilico_data={},
        functional_data={}
    )
    
    ps2_result = evaluator._evaluate_ps2(confirmed_variant)
    pm6_result = evaluator._evaluate_pm6(confirmed_variant)
    
    print(f"PS2 applies: {ps2_result['applies']}")
    print(f"PS2 details: {ps2_result['details']}")
    print(f"PM6 applies: {pm6_result['applies']}")
    print(f"PM6 details: {pm6_result['details']}")
    print()
    
    # Test 2: Assumed de novo (no parental testing)
    print("Test 2: Assumed de novo (no parental testing)")
    print("-" * 50)
    assumed_variant = VariantData(
        basic_info={'gene': 'TEST'},
        genetic_data={
            'de_novo': 'assumed',
            'maternity_confirmed': False,
            'paternity_confirmed': False
        },
        population_data={},
        insilico_data={},
        functional_data={}
    )
    
    ps2_result = evaluator._evaluate_ps2(assumed_variant)
    pm6_result = evaluator._evaluate_pm6(assumed_variant)
    
    print(f"PS2 applies: {ps2_result['applies']}")
    print(f"PS2 details: {ps2_result['details']}")
    print(f"PM6 applies: {pm6_result['applies']}")
    print(f"PM6 details: {pm6_result['details']}")
    print()
    
    # Test 3: Confirmed de novo but no parental testing
    print("Test 3: Confirmed de novo but no parental testing")
    print("-" * 50)
    confused_variant = VariantData(
        basic_info={'gene': 'TEST'},
        genetic_data={
            'de_novo': 'confirmed',
            'maternity_confirmed': False,
            'paternity_confirmed': False
        },
        population_data={},
        insilico_data={},
        functional_data={}
    )
    
    ps2_result = evaluator._evaluate_ps2(confused_variant)
    pm6_result = evaluator._evaluate_pm6(confused_variant)
    
    print(f"PS2 applies: {ps2_result['applies']}")
    print(f"PS2 details: {ps2_result['details']}")
    print(f"PM6 applies: {pm6_result['applies']}")
    print(f"PM6 details: {pm6_result['details']}")

if __name__ == "__main__":
    test_denovo_logic()
