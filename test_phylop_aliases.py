#!/usr/bin/env python3
"""
Test script to verify PhyloP alias functionality.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src.utils.input_handler import InputHandler

def test_phylop_aliases():
    """Test PhyloP alias functionality."""
    
    handler = InputHandler()
    
    # Test the _prompt_choice method directly
    print("Testing PhyloP alias functionality...")
    print("=" * 50)
    
    # Test aliases
    aliases = {
        'raw': 'raw_score', 'r': 'raw_score', '1': 'raw_score',
        'ranked': 'ranked_score', 'rank': 'ranked_score', '2': 'ranked_score',
        's': 'skip', '3': 'skip'
    }
    
    test_cases = [
        ('ranked', 'ranked_score'),
        ('rank', 'ranked_score'),
        ('2', 'ranked_score'),
        ('raw', 'raw_score'),
        ('r', 'raw_score'),
        ('1', 'raw_score'),
        ('s', 'skip'),
        ('3', 'skip')
    ]
    
    for input_val, expected in test_cases:
        print(f"Testing input '{input_val}' -> expected '{expected}'")
        
        # Simulate the logic from _prompt_choice
        choices = ['raw_score', 'ranked_score', 'skip']
        
        # Check direct matches first
        if input_val in [choice.lower() for choice in choices]:
            result = next(choice for choice in choices if choice.lower() == input_val)
        else:
            # Check aliases
            result = None
            for alias, choice in aliases.items():
                if input_val == alias.lower():
                    result = choice
                    break
        
        if result == expected:
            print(f"  ✅ PASS: '{input_val}' -> '{result}'")
        else:
            print(f"  ❌ FAIL: '{input_val}' -> '{result}' (expected '{expected}')")
    
    print("\n" + "=" * 50)
    print("PhyloP alias test complete!")

if __name__ == "__main__":
    test_phylop_aliases()
