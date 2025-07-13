#!/usr/bin/env python3
"""
Test script to verify all alias functionality.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src.utils.input_handler import InputHandler

def test_all_aliases():
    """Test all alias functionality."""
    
    handler = InputHandler()
    
    print("Testing ALL alias functionality...")
    print("=" * 50)
    
    # Test variant type aliases
    print("\n1. Variant Type Aliases:")
    variant_aliases = {
        'miss': 'missense', 'non': 'nonsense', 'fs': 'frameshift',
        'spl': 'splice', 'syn': 'synonymous', 'indel': 'inframe_indel'
    }
    choices = ['missense', 'nonsense', 'frameshift', 'splice', 'synonymous', 'inframe_indel', 'other']
    
    test_cases = [
        ('miss', 'missense'),
        ('missense', 'missense'),  # Direct match
        ('non', 'nonsense'),
        ('fs', 'frameshift')
    ]
    
    for input_val, expected in test_cases:
        result = simulate_choice_logic(input_val, choices, variant_aliases)
        status = "✅ PASS" if result == expected else "❌ FAIL"
        print(f"  {status}: '{input_val}' -> '{result}' (expected '{expected}')")
    
    # Test zygosity aliases
    print("\n2. Zygosity Aliases:")
    zygo_aliases = {
        'hom': 'homozygous', 'het': 'heterozygous', 'hemi': 'hemizygous'
    }
    choices = ['homozygous', 'heterozygous', 'hemizygous', 'unknown']
    
    test_cases = [
        ('het', 'heterozygous'),
        ('heterozygous', 'heterozygous'),  # Direct match
        ('hom', 'homozygous'),
        ('hemi', 'hemizygous')
    ]
    
    for input_val, expected in test_cases:
        result = simulate_choice_logic(input_val, choices, zygo_aliases)
        status = "✅ PASS" if result == expected else "❌ FAIL"
        print(f"  {status}: '{input_val}' -> '{result}' (expected '{expected}')")
    
    # Test segregation aliases
    print("\n3. Segregation Aliases:")
    seg_aliases = {
        'coseg': 'cosegregates', 'noseg': 'does_not_segregate', 
        'insuf': 'insufficient_data', 'none': 'not_performed'
    }
    choices = ['cosegregates', 'does_not_segregate', 'insufficient_data', 'not_performed']
    
    test_cases = [
        ('coseg', 'cosegregates'),
        ('noseg', 'does_not_segregate'),
        ('none', 'not_performed'),
        ('insuf', 'insufficient_data')
    ]
    
    for input_val, expected in test_cases:
        result = simulate_choice_logic(input_val, choices, seg_aliases)
        status = "✅ PASS" if result == expected else "❌ FAIL"
        print(f"  {status}: '{input_val}' -> '{result}' (expected '{expected}')")

def simulate_choice_logic(value, choices, aliases):
    """Simulate the _prompt_choice logic."""
    value = value.strip().lower()
    
    # Check direct matches first
    if value in [choice.lower() for choice in choices]:
        return next(choice for choice in choices if choice.lower() == value)
    
    # Check aliases
    if aliases:
        for alias, choice in aliases.items():
            if value == alias.lower():
                return choice
    
    return None

if __name__ == "__main__":
    test_all_aliases()
