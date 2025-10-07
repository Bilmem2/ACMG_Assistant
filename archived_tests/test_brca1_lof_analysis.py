#!/usr/bin/env python3
"""
Analyze BRCA1 interpretations for LOF mechanism evidence
"""

import requests
import json
from collections import Counter

def analyze_brca1_lof():
    """Analyze BRCA1 interpretations for LOF pathogenicity"""
    
    print("="*80)
    print("🧬 BRCA1 LOF Mechanism Analysis via ClinGen eRepo")
    print("="*80)
    
    base_url = "https://erepo.genome.network/evrepo/api"
    
    # Get BRCA1 interpretations
    print("\n📥 Fetching BRCA1 interpretations...")
    response = requests.get(
        f"{base_url}/interpretations",
        params={'gene': 'BRCA1', 'matchLimit': '50'},
        timeout=15
    )
    
    if response.status_code != 200:
        print(f"❌ Failed to fetch interpretations: {response.status_code}")
        return
    
    data = response.json()
    interpretations = data.get('variantInterpretations', [])
    
    print(f"   Found {len(interpretations)} interpretations")
    
    # Statistics
    pathogenic_count = 0
    benign_count = 0
    vus_count = 0
    pvs1_count = 0
    lof_variants = []
    evidence_codes = Counter()
    
    print("\n🔍 Analyzing detailed interpretations...")
    
    for i, interp_summary in enumerate(interpretations[:20], 1):  # First 20
        interp_id = interp_summary.get('@id', '')
        
        try:
            # Get full interpretation
            detail_response = requests.get(interp_id, timeout=10)
            
            if detail_response.status_code == 200:
                detail = detail_response.json()
                
                # Extract outcome
                outcome = detail.get('statementOutcome', {}).get('label', 'Unknown')
                
                # Count outcomes
                if 'Pathogenic' in outcome:
                    pathogenic_count += 1
                elif 'Benign' in outcome:
                    benign_count += 1
                elif 'Uncertain' in outcome or 'VUS' in outcome:
                    vus_count += 1
                
                # Extract variant info
                variant = detail.get('variant', {})
                variant_name = ''
                if 'relatedIdentifier' in variant:
                    identifiers = variant['relatedIdentifier']
                    if len(identifiers) > 0:
                        variant_name = identifiers[0].get('label', '')
                
                # Check evidence lines for PVS1
                evidence_lines = detail.get('evidenceLine', [])
                has_pvs1 = False
                
                for ev_line in evidence_lines:
                    if isinstance(ev_line, dict):
                        # Check for evidence codes
                        ev_items = ev_line.get('evidenceItem', [])
                        for item in ev_items:
                            if isinstance(item, dict):
                                ev_code_obj = item.get('evidenceCode', {})
                                if isinstance(ev_code_obj, dict):
                                    ev_label = ev_code_obj.get('label', '')
                                    if ev_label:
                                        evidence_codes[ev_label] += 1
                                        
                                        if 'PVS1' in ev_label:
                                            has_pvs1 = True
                                            pvs1_count += 1
                
                # If pathogenic with PVS1, likely LOF
                if 'Pathogenic' in outcome and has_pvs1:
                    lof_variants.append({
                        'variant': variant_name,
                        'outcome': outcome
                    })
                    print(f"   {i}. ✅ LOF Pathogenic: {variant_name[:60]} → {outcome}")
                
        except Exception as e:
            print(f"   {i}. ⚠️  Error processing {interp_id[-50:]}: {str(e)[:50]}")
    
    # Results
    print("\n" + "="*80)
    print("📊 ANALYSIS RESULTS:")
    print("="*80)
    print(f"\n📈 Classification Distribution:")
    print(f"   Pathogenic: {pathogenic_count}")
    print(f"   Benign: {benign_count}")
    print(f"   VUS: {vus_count}")
    
    print(f"\n🧬 LOF Evidence:")
    print(f"   Variants with PVS1: {pvs1_count}")
    print(f"   LOF Pathogenic variants: {len(lof_variants)}")
    
    print(f"\n📋 Top Evidence Codes:")
    for code, count in evidence_codes.most_common(10):
        print(f"   {code}: {count}")
    
    print(f"\n💡 Conclusion:")
    if pvs1_count >= 3:
        print(f"   ✅ BRCA1 shows strong LOF pathogenic mechanism")
        print(f"   {pvs1_count} variants with PVS1 (null variant) evidence")
        print(f"   LOF variants cause BRCA1-related cancer predisposition")
    else:
        print(f"   ⚠️  Limited PVS1 evidence in sample")
    
    print("="*80)

if __name__ == '__main__':
    analyze_brca1_lof()
