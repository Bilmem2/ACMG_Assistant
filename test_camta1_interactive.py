#!/usr/bin/env python3
"""
Interactive test script for CAMTA1 variant using real data
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from acmg_assistant import ACMGAssistant
from io import StringIO
from unittest.mock import patch

def test_camta1_interactive():
    """Test CAMTA1 variant with interactive input simulation."""
    
    # Simulate user inputs based on your provided data
    user_inputs = [
        # Basic info
        'CAMTA1',           # Gene name
        '1',                # Chromosome
        '7249507',          # Position  
        'A',                # Reference allele
        'G',                # Alternative allele
        'c.319A>G',         # cDNA change
        'p.Met107Val',      # Protein change
        'missense',         # Variant type
        '',                 # Transcript (empty)
        
        # Population data
        '1.505e-05',        # gnomAD AF
        '',                 # gnomAD AF popmax (empty)
        '',                 # gnomAD AC (empty)
        '',                 # gnomAD AN (empty)
        '',                 # gnomAD hom (empty)
        '1e-06',            # Disease prevalence
        '',                 # 1000G AF (empty)
        '',                 # ExAC AF (empty)
        '',                 # ESP AF (empty)
        
        # In silico scores
        'y',                # Want to enter scores
        '0.395',            # REVEL
        '24.3',             # CADD phred
        '0.349',            # ClinPred
        '0.7275',           # BayesDel addAF
        '0.3896',           # AlphaMissense
        '0.9999',           # MutationTaster
        '0.6866',           # PolyPhen2 HDiv
        '0.041',            # SIFT
        '0.904',            # FATHMM-XF
        '0.12',             # MutationAssessor
        '0.3158',           # PROVEAN
        '0.559',            # MutPred
        '0.0946',           # MetaLR
        '-8.208',           # ESM1b
        '0.44119',          # LRT
        'raw',              # PhyloP choice
        '0.94',             # PhyloP vertebrates
        '3.0',              # PhyloP mammals
        '3.0',              # PhyloP primates
        '5.63',             # GERP++
        
        # Genetic data
        'ad',               # Inheritance pattern
        'heterozygous',     # Zygosity
        '',                 # Allelic state (empty)
        'no',               # Family history
        'no',               # Segregation analysis
        'yes',              # De novo
        '',                 # Parental testing (empty)
        '',                 # Consanguinity (empty)
        'no',               # Siblings affected
        '',                 # Notes (empty)
        
        # Functional data
        '',                 # Functional test (empty)
        'moderate',         # Phenotype match
        '',                 # Paternity confirmed (empty)
        '',                 # Notes (empty)
        
        'n'                 # Don't analyze another variant
    ]
    
    print("🧬 Testing CAMTA1 c.319A>G (p.Met107Val) - Interactive Mode")
    print("="*70)
    print("📝 Simulating user inputs from your provided data...")
    print()
    
    # Create assistant in non-test mode
    assistant = ACMGAssistant(use_2023_guidelines=False, test_mode=False)
    
    # Mock user input
    with patch('builtins.input', side_effect=user_inputs):
        try:
            assistant.run()
        except SystemExit:
            pass  # Normal exit
    
    print("\n✅ Interactive test completed!")

if __name__ == "__main__":
    test_camta1_interactive()
