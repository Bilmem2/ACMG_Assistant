#!/usr/bin/env python3
"""
Test algorithm performance across different variant types
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from core.variant_data import VariantData
from core.acmg_classifier import ACMGClassifier
from core.evidence_evaluator import EvidenceEvaluator
from colorama import Fore, Style, init

# Initialize colorama
init()

def test_variant_types():
    """Test algorithm with different variant types"""
    
    print(f"{Fore.CYAN}🧬 Testing Algorithm Across Variant Types{Style.RESET_ALL}")
    print("=" * 60)
    
    # Test scenarios for different variant types
    test_cases = [
        {
            'name': 'Missense Variant (BRCA1)',
            'basic_info': {
                'gene': 'BRCA1',
                'variant_type': 'missense',
                'consequence': 'missense_variant',
                'chromosome': '17'
            },
            'insilico_data': {
                'revel_score': 0.8,
                'cadd_phred_score': 28,
                'alphamissense_score': 0.7,
                'sift_score': 0.02,
                'polyphen2_score': 0.95
            },
            'expected_pp3': True
        },
        {
            'name': 'Nonsense Variant (TP53)',
            'basic_info': {
                'gene': 'TP53',
                'variant_type': 'nonsense',
                'consequence': 'stop_gained',
                'chromosome': '17'
            },
            'insilico_data': {
                'cadd_phred_score': 35,
                'phylop_vert_score': 5.2
            },
            'expected_pvs1': True
        },
        {
            'name': 'Splice Site Variant (CFTR)',
            'basic_info': {
                'gene': 'CFTR',
                'variant_type': 'splice',
                'consequence': 'splice_donor_variant',
                'chromosome': '7'
            },
            'insilico_data': {
                'spliceai_dg_score': 0.95,
                'spliceai_dl_score': 0.02,
                'ada_score': 0.8,
                'rf_score': 0.9
            },
            'expected_pvs1': True
        },
        {
            'name': 'Intronic Variant with Splice Impact',
            'basic_info': {
                'gene': 'DMD',
                'variant_type': 'intronic',
                'consequence': 'intron_variant',
                'chromosome': 'X'
            },
            'insilico_data': {
                'spliceai_ag_score': 0.7,
                'spliceai_dg_score': 0.3,
                'phylop_vert_score': 4.2,
                'gerp_pp_score': 5.8
            },
            'expected_ps1_or_pm4': True
        },
        {
            'name': 'Synonymous Variant (Low Impact)',
            'basic_info': {
                'gene': 'LDLR',
                'variant_type': 'synonymous',
                'consequence': 'synonymous_variant',
                'chromosome': '19'
            },
            'insilico_data': {
                'spliceai_ag_score': 0.05,
                'spliceai_dg_score': 0.03,
                'phylop_vert_score': 1.2,
                'cadd_phred_score': 8
            },
            'expected_benign': True
        },
        {
            'name': 'Frameshift Variant (LOF)',
            'basic_info': {
                'gene': 'BRCA2',
                'variant_type': 'frameshift',
                'consequence': 'frameshift_variant',
                'chromosome': '13'
            },
            'insilico_data': {
                'cadd_phred_score': 32,
                'phylop_vert_score': 3.8
            },
            'expected_pvs1': True
        }
    ]
    
    # Initialize classifier
    classifier = ACMGClassifier()
    evidence_evaluator = EvidenceEvaluator()
    
    results = []
    
    for i, case in enumerate(test_cases, 1):
        print(f"\n{Fore.YELLOW}Test {i}: {case['name']}{Style.RESET_ALL}")
        print(f"Type: {case['basic_info']['variant_type']}")
        print(f"Consequence: {case['basic_info']['consequence']}")
        
        try:
            # Create variant data
            variant_data = VariantData()
            
            # Set the data properly
            for key, value in case['basic_info'].items():
                variant_data.basic_info[key] = value
            
            for key, value in case['insilico_data'].items():
                variant_data.insilico_data[key] = value
            
            # Initialize empty dictionaries for other required data
            if not variant_data.population_data:
                variant_data.population_data = {}
            if not variant_data.genetic_data:
                variant_data.genetic_data = {}
            if not variant_data.functional_data:
                variant_data.functional_data = {}
            
            # Evaluate evidence
            evidence = evidence_evaluator.evaluate_all_criteria(variant_data)
            
            # Check specific expectations
            expectations_met = []
            
            if case.get('expected_pvs1'):
                pvs1_applies = evidence.get('PVS1', {}).get('applies', False)
                expectations_met.append(f"PVS1: {'✅' if pvs1_applies else '❌'}")
            
            if case.get('expected_pp3'):
                pp3_applies = evidence.get('PP3', {}).get('applies', False)
                expectations_met.append(f"PP3: {'✅' if pp3_applies else '❌'}")
            
            if case.get('expected_ps1_or_pm4'):
                ps1_applies = evidence.get('PS1', {}).get('applies', False)
                pm4_applies = evidence.get('PM4', {}).get('applies', False)
                expectations_met.append(f"PS1/PM4: {'✅' if (ps1_applies or pm4_applies) else '❌'}")
            
            if case.get('expected_benign'):
                bp4_applies = evidence.get('BP4', {}).get('applies', False)
                bp7_applies = evidence.get('BP7', {}).get('applies', False)
                expectations_met.append(f"Benign evidence: {'✅' if (bp4_applies or bp7_applies) else '❌'}")
            
            # Print results
            applied_criteria = [code for code, ev in evidence.items() if ev.get('applies', False)]
            print(f"Applied criteria: {', '.join(applied_criteria) if applied_criteria else 'None'}")
            print(f"Expectations: {', '.join(expectations_met)}")
            
            # Try classification
            try:
                classification = classifier.classify(evidence)
                print(f"Classification: {Fore.GREEN}{classification['classification']}{Style.RESET_ALL}")
                results.append({
                    'name': case['name'],
                    'classification': classification['classification'],
                    'applied_criteria': applied_criteria,
                    'expectations_met': len([e for e in expectations_met if '✅' in e])
                })
            except Exception as e:
                print(f"Classification error: {e}")
                results.append({
                    'name': case['name'],
                    'classification': 'Error',
                    'applied_criteria': applied_criteria,
                    'expectations_met': 0
                })
                
        except Exception as e:
            print(f"{Fore.RED}❌ Error: {e}{Style.RESET_ALL}")
            results.append({
                'name': case['name'],
                'classification': 'Error',
                'applied_criteria': [],
                'expectations_met': 0
            })
    
    # Summary
    print(f"\n{Fore.CYAN}📊 Test Summary{Style.RESET_ALL}")
    print("=" * 60)
    
    total_tests = len(results)
    successful_tests = len([r for r in results if r['classification'] != 'Error'])
    
    print(f"Total tests: {total_tests}")
    print(f"Successful: {successful_tests}")
    print(f"Success rate: {(successful_tests/total_tests)*100:.1f}%")
    
    for result in results:
        color = Fore.GREEN if result['classification'] != 'Error' else Fore.RED
        print(f"  {color}{result['name']}: {result['classification']}{Style.RESET_ALL}")

if __name__ == "__main__":
    test_variant_types()
