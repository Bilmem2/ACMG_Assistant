#!/usr/bin/env python3
"""
Test interactive mode with simulated input
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from io import StringIO
from unittest.mock import patch
from src.core.evidence_evaluator import EvidenceEvaluator
from src.core.variant_data import VariantData

def test_interactive_ps1():
    """Test PS1 interactive evaluation."""
    
    # Create test data
    basic_info = {
        'gene': 'BRCA1',
        'amino_acid_change': 'p.Arg123Trp',
        'variant_type': 'missense'
    }
    
    variant_data = VariantData(
        basic_info=basic_info,
        population_data={},
        insilico_data={},
        genetic_data={},
        functional_data={}
    )
    
    # Test with mock input
    evaluator = EvidenceEvaluator(use_2023_guidelines=False, test_mode=False)
    
    # Mock user input "y" (yes)
    with patch('builtins.input', return_value='y'):
        result = evaluator._evaluate_ps1(variant_data)
    
    print("PS1 Interactive Test Results:")
    print(f"Applies: {result['applies']}")
    print(f"Details: {result['details']}")
    print(f"Strength: {result['strength']}")
    
    # Test with mock input "n" (no)
    with patch('builtins.input', return_value='n'):
        result = evaluator._evaluate_ps1(variant_data)
    
    print("\nPS1 Interactive Test Results (No):")
    print(f"Applies: {result['applies']}")
    print(f"Details: {result['details']}")

if __name__ == "__main__":
    test_interactive_ps1()
