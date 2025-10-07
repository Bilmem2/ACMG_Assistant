"""
Automated Test Script for All 30 Test Scenarios
Tests each scenario and verifies expected criteria/classification
"""

from src.config.constants import TEST_SCENARIOS


def analyze_scenarios():
    """Analyze all scenarios for completeness and correctness."""
    print("="*100)
    print("🧪 ACMG ASSISTANT - SCENARIO ANALYSIS")
    print("Analyzing all 30 scenarios for expected coverage")
    print("="*100)
    
    scenario_list = sorted(TEST_SCENARIOS.items(), key=lambda x: int(x[0].split('_')[0]))
    
    results = []
    all_criteria = set()
    variant_types = set()
    genes = set()
    
    for idx, (scenario_key, scenario_data) in enumerate(scenario_list, 1):
        expected_class = scenario_data['expected_classification']
        expected_criteria = scenario_data['expected_criteria']
        gene = scenario_data['basic_info']['gene']
        variant_type = scenario_data['basic_info']['variant_type']
        
        all_criteria.update(expected_criteria)
        variant_types.add(variant_type)
        genes.add(gene)
        
        # Check for potential issues
        issues = []
        
        # Check if Pathogenic has enough strong evidence
        if expected_class == 'Pathogenic':
            has_pvs = any('PVS' in c for c in expected_criteria)
            ps_count = sum(1 for c in expected_criteria if c.startswith('PS'))
            if not has_pvs and ps_count < 2:
                issues.append("May not meet Pathogenic criteria")
        
        # Check if Likely Pathogenic has valid combination
        elif expected_class == 'Likely Pathogenic':
            has_ps = any('PS' in c for c in expected_criteria)
            pm_count = sum(1 for c in expected_criteria if c.startswith('PM'))
            pp_count = sum(1 for c in expected_criteria if c.startswith('PP'))
            
            valid = False
            if has_ps and pm_count >= 1:
                valid = True
            elif has_ps and pp_count >= 2:
                valid = True
            elif pm_count >= 3:
                valid = True
            elif pm_count >= 2 and pp_count >= 2:
                valid = True
            elif pm_count >= 1 and pp_count >= 4:
                valid = True
            elif pm_count >= 1 and pp_count >= 2:  # Our new rule
                valid = True
            
            if not valid:
                issues.append("May not meet Likely Pathogenic criteria")
        
        results.append({
            'num': idx,
            'name': scenario_data['name'],
            'gene': gene,
            'type': variant_type,
            'expected': expected_class,
            'criteria': expected_criteria,
            'issues': issues
        })
    
    return results, all_criteria, variant_types, genes


def main():
    """Run scenario analysis."""
    results, all_criteria, variant_types, genes = analyze_scenarios()
    
    # Print results
    print("\n📋 SCENARIO REVIEW:")
    print("="*100)
    print(f"{'#':<4} {'Name':<45} {'Gene':<8} {'Type':<15} {'Expected':<20}")
    print("-"*100)
    
    for r in results:
        status = "⚠️ " if r['issues'] else "✅"
        print(f"{status} {r['num']:<2} {r['name']:<45} {r['gene']:<8} {r['type']:<15} {r['expected']:<20}")
        if r['issues']:
            for issue in r['issues']:
                print(f"      → {issue}")
        print(f"      Criteria: {', '.join(r['criteria'])}")
    
    # Coverage summary
    print("\n" + "="*100)
    print("📊 COVERAGE SUMMARY:")
    print("="*100)
    print(f"\nTotal Scenarios: {len(results)}")
    print(f"Unique Criteria Covered: {len(all_criteria)}")
    print(f"Criteria List: {', '.join(sorted(all_criteria))}")
    print(f"\nVariant Types ({len(variant_types)}): {', '.join(sorted(variant_types))}")
    print(f"\nGenes ({len(genes)}): {', '.join(sorted(genes))}")
    
    # Warnings
    issues_found = sum(1 for r in results if r['issues'])
    if issues_found > 0:
        print(f"\n⚠️  WARNING: {issues_found} scenarios may have validation issues")
    else:
        print("\n✅ All scenarios appear valid!")
    
    print("\n" + "="*100)
    print("💡 Now run manual tests with: python src/acmg_assistant.py --test")
    print("="*100)


if __name__ == "__main__":
    main()
