#!/usr/bin/env python3
"""
Test script for protein validation function
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from utils.input_handler import InputHandler

def test_protein_validation():
    """Test protein validation function with various formats."""
    
    handler = InputHandler()
    
    # Test proteins from your example
    test_proteins = [
        'p.Met107Val',      # 3-letter to 3-letter
        'p.M107V',          # 1-letter to 1-letter
        'p.R123Q',          # 1-letter to 1-letter
        'p.Arg123Gln',      # 3-letter to 3-letter
        'p.Gln1756Profs*74',# Frameshift
        'p.Trp24*',         # Nonsense
        'p.Arg123fs',       # Frameshift (short)
        'p.Met1?',          # Start lost
        'p.Ter110GlnextTer17', # Stop extension
        'p.=',              # Synonymous
        'p.?',              # Unknown effect
        'p.0',              # No protein
        'invalid_format',   # Invalid
        'p.Met107',         # Incomplete
        'p.Met107VaI',      # Typo (capital i instead of l)
        'p.Xxx123Yyy',      # Invalid amino acids
    ]
    
    print("🧪 Testing Protein Validation Function")
    print("="*50)
    
    for protein in test_proteins:
        try:
            result = handler._validate_hgvs_protein(protein)
            status = "✅ Valid" if result else "❌ Invalid"
            print(f"{protein:<20} : {status}")
        except Exception as e:
            print(f"{protein:<20} : ❌ Error: {str(e)}")
    
    print("\n" + "="*50)
    print("✅ Protein validation test completed!")

if __name__ == "__main__":
    test_protein_validation()
