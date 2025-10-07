#!/usr/bin/env python3
"""
Test BRCA1 PVS1 with multi-source validation (gnomAD + ClinGen)
"""

import sys
sys.path.insert(0, 'src')

from core.variant_data import VariantData
from core.evidence_evaluator import EvidenceEvaluator

def test_brca1_nonsense():
    """Test BRCA1 p.Gln356Ter nonsense variant"""
    
    print("="*80)
    print("🧬 Testing BRCA1 PVS1 with Multi-Source Validation")
    print("="*80)
    print("\n📋 Variant Details:")
    print("  Gene: BRCA1")
    print("  Type: Nonsense (stop_gained)")
    print("  HGVSp: p.Gln356Ter")
    print("  Expected: Pathogenic (PVS1 + PM2)")
    print("\n🔬 Testing with:")
    print("  ✓ gnomAD constraint API (population LOF tolerance)")
    print("  ✓ ClinGen gene-disease validity API (disease-specific LOF mechanism)")
    print("="*80)
    
    # Create variant data
    variant_data = VariantData()
    variant_data.basic_info = {
        'gene': 'BRCA1',
        'variant_type': 'nonsense',
        'consequence': 'stop_gained',
        'hgvs_p': 'p.Gln356Ter',
        'transcript': 'NM_007294.3'
    }
    
    variant_data.population_data = {
        'gnomad_af': 0.0,  # Absent in gnomAD
        'gnomad_ac': 0,
        'gnomad_an': 251456,
        'gnomad_hom': 0
    }
    
    variant_data.insilico_data = {
        'revel_score': None,
        'cadd_score': None,
        'sift_pred': None,
        'polyphen_pred': None
    }
    
    # Initialize evaluator with API enabled
    evaluator = EvidenceEvaluator(use_2023_guidelines=False, test_mode=False)
    
    # Evaluate PVS1
    print("\n🔍 Evaluating PVS1...")
    pvs1_result = evaluator._evaluate_pvs1(variant_data)
    
    print("\n📊 RESULTS:")
    print("="*80)
    print(f"  PVS1 Applies: {pvs1_result['applies']}")
    print(f"  Strength: {pvs1_result.get('strength', 'N/A')}")
    print(f"  Data Source: {pvs1_result.get('data_source', 'N/A')}")
    print(f"\n  Details:\n    {pvs1_result['details']}")
    
    if 'constraint_metrics' in pvs1_result:
        print("\n  📈 Constraint Metrics:")
        metrics = pvs1_result['constraint_metrics']
        print(f"    gnomAD pLI: {metrics.get('gnomad_pLI', 'N/A')}")
        print(f"    gnomAD LOEUF: {metrics.get('gnomad_LOEUF', 'N/A')}")
        print(f"    gnomAD Classification: {metrics.get('gnomad_classification', 'N/A')}")
        print(f"    ClinGen supports LOF: {metrics.get('clingen_supports_lof', 'N/A')}")
        print(f"    ClinGen mechanism: {metrics.get('clingen_mechanism', 'N/A')}")
        if 'clingen_diseases' in metrics:
            print(f"    ClinGen LOF diseases: {', '.join(metrics['clingen_diseases'][:3])}")
    
    print("="*80)
    
    # Expected vs Actual
    expected_applies = True  # PVS1 should apply for BRCA1 nonsense
    if pvs1_result['applies'] == expected_applies:
        print("\n✅ TEST PASSED: PVS1 correctly applied")
        print("   BRCA1 LOF variants are pathogenic for breast/ovarian cancer")
        print("   ClinGen should confirm LOF mechanism despite gnomAD showing population tolerance")
    else:
        print("\n❌ TEST FAILED: PVS1 should apply for BRCA1 nonsense")
        print("   BRCA1 is a tumor suppressor - LOF variants cause cancer")
    
    print("\n💡 Interpretation:")
    print("   BRCA1 shows LOF tolerance in general population (pLI=0.0)")
    print("   BUT has strong disease-specific LOF mechanism (haploinsufficiency)")
    print("   Multi-source validation resolves this by prioritizing ClinGen disease evidence")
    print("="*80)

if __name__ == '__main__':
    test_brca1_nonsense()
