#!/usr/bin/env python3
"""
Comprehensive test scenarios for ACMG Variant Classification Assistant
Testing different variant types to evaluate algorithm performance
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from core.variant_data import VariantData
from core.acmg_classifier import ACMGClassifier
from core.evidence_evaluator import EvidenceEvaluator
from utils.report_generator import ReportGenerator

def test_scenario_1_vus():
    """
    Scenario 1: VUS - Mixed evidence
    BRCA1 missense with conflicting in silico predictions
    """
    print("🧬 SCENARIO 1: VUS - Mixed Evidence")
    print("="*60)
    print("BRCA1 c.4327C>T (p.Arg1443Cys) - Mixed in silico predictions")
    
    basic_info = {
        'gene': 'BRCA1',
        'chromosome': '17',
        'position': 41244000,
        'ref_allele': 'C',
        'alt_allele': 'T',
        'cdna_change': 'c.4327C>T',
        'protein_change': 'p.Arg1443Cys',
        'variant_type': 'missense',
        'transcript': 'NM_007294.3',
        'vep_consequence': 'missense_variant'
    }
    
    population_data = {
        'gnomad_af': 3.2e-05,
        'gnomad_af_popmax': 4.1e-05,
        'disease_prevalence': 1e-04,  # BRCA1 breast cancer
        'gnomad_ac': 12,
        'gnomad_an': 375000,
        'gnomad_hom': 0
    }
    
    # Mixed in silico predictions (some pathogenic, some benign)
    insilico_data = {
        'revel': 0.42,           # Borderline
        'cadd_phred': 22.1,      # Moderate
        'clinpred': 0.35,        # Borderline
        'bayesdel_addaf': 0.02,  # Benign
        'alphamissense': 0.45,   # Borderline
        'mutationtaster': 0.6,   # Moderate pathogenic
        'polyphen2_hdiv': 0.3,   # Benign
        'sift': 0.12,            # Benign (SIFT inverted)
        'fathmm_xf': 0.35,       # Borderline
        'phylop_vertebrates': 1.2, # Moderately conserved
        'phylop_mammals': 2.1,
        'phylop_primates': 1.8,
        'gerp_rs': 4.2,
        'metarnn': 0.38,
        'provean': -1.2,         # Borderline
        'lrt': 0.45
    }
    
    genetic_data = {
        'inheritance_pattern': 'ad',
        'zygosity': 'heterozygous',
        'family_history': 'unknown',
        'segregation_analysis': 'unknown',
        'de_novo': 'unknown',
        'parental_testing': 'unknown'
    }
    
    functional_data = {
        'functional_test': 'unknown',
        'phenotype_match': 'partial',
    }
    
    return VariantData(basic_info, population_data, insilico_data, genetic_data, functional_data)

def test_scenario_2_likely_benign():
    """
    Scenario 2: Likely Benign - High frequency + benign predictions
    Common population variant with benign in silico predictions
    """
    print("\n🧬 SCENARIO 2: Likely Benign - Common Variant")
    print("="*60)
    print("CFTR c.1408A>G (p.Met470Val) - Common benign variant")
    
    basic_info = {
        'gene': 'CFTR',
        'chromosome': '7',
        'position': 117559590,
        'ref_allele': 'A',
        'alt_allele': 'G',
        'cdna_change': 'c.1408A>G',
        'protein_change': 'p.Met470Val',
        'variant_type': 'missense',
        'transcript': 'NM_000492.3',
        'vep_consequence': 'missense_variant'
    }
    
    population_data = {
        'gnomad_af': 0.28,           # Very common (28%)
        'gnomad_af_popmax': 0.32,
        'disease_prevalence': 1e-05, # CF prevalence
        'gnomad_ac': 150000,
        'gnomad_an': 535000,
        'gnomad_hom': 25000
    }
    
    # Predominantly benign in silico predictions
    insilico_data = {
        'revel': 0.12,           # Benign
        'cadd_phred': 8.2,       # Low
        'clinpred': 0.15,        # Benign
        'bayesdel_addaf': -0.35, # Benign
        'alphamissense': 0.25,   # Benign
        'mutationtaster': 0.2,   # Benign
        'polyphen2_hdiv': 0.05,  # Benign
        'sift': 0.82,            # Benign (SIFT inverted)
        'fathmm_xf': 0.15,       # Benign
        'phylop_vertebrates': 0.2, # Not conserved
        'phylop_mammals': 0.1,
        'phylop_primates': 0.3,
        'gerp_rs': 1.2,
        'metarnn': 0.18,
        'provean': 1.8,          # Benign
        'lrt': 0.22
    }
    
    genetic_data = {
        'inheritance_pattern': 'ar',
        'zygosity': 'heterozygous',
        'family_history': 'no',
        'segregation_analysis': 'no',
        'de_novo': 'no'
    }
    
    functional_data = {
        'functional_test': 'unknown',
        'phenotype_match': 'no_match'
    }
    
    return VariantData(basic_info, population_data, insilico_data, genetic_data, functional_data)

def test_scenario_3_pathogenic():
    """
    Scenario 3: Pathogenic - Strong evidence
    Known pathogenic variant with strong computational support
    """
    print("\n🧬 SCENARIO 3: Pathogenic - Strong Evidence")
    print("="*60)
    print("TP53 c.817C>T (p.Arg273Cys) - Known pathogenic hotspot")
    
    basic_info = {
        'gene': 'TP53',
        'chromosome': '17',
        'position': 7674220,
        'ref_allele': 'C',
        'alt_allele': 'T',
        'cdna_change': 'c.817C>T',
        'protein_change': 'p.Arg273Cys',
        'variant_type': 'missense',
        'transcript': 'NM_000546.5',
        'vep_consequence': 'missense_variant'
    }
    
    population_data = {
        'gnomad_af': 2.1e-06,        # Very rare
        'gnomad_af_popmax': 3.2e-06,
        'disease_prevalence': 1e-03, # Cancer predisposition
        'gnomad_ac': 3,
        'gnomad_an': 1425000,
        'gnomad_hom': 0
    }
    
    # Strong pathogenic in silico predictions
    insilico_data = {
        'revel': 0.89,           # High pathogenic
        'cadd_phred': 35.2,      # High
        'clinpred': 0.92,        # High pathogenic
        'bayesdel_addaf': 0.78,  # Pathogenic
        'alphamissense': 0.85,   # High pathogenic
        'mutationtaster': 0.98,  # High pathogenic
        'polyphen2_hdiv': 0.95,  # Probably damaging
        'sift': 0.01,            # Damaging (SIFT inverted)
        'fathmm_xf': 0.92,       # High pathogenic
        'phylop_vertebrates': 5.8, # Highly conserved
        'phylop_mammals': 6.2,
        'phylop_primates': 5.9,
        'gerp_rs': 5.98,
        'metarnn': 0.88,
        'provean': -6.8,         # Damaging
        'lrt': 0.92,
        'esm1b': -12.5           # Damaging
    }
    
    genetic_data = {
        'inheritance_pattern': 'ad',
        'zygosity': 'heterozygous',
        'family_history': 'yes',
        'segregation_analysis': 'yes',
        'de_novo': 'confirmed',
        'parental_testing': 'confirmed'
    }
    
    functional_data = {
        'functional_test': 'damaging',
        'phenotype_match': 'specific_match'
    }
    
    return VariantData(basic_info, population_data, insilico_data, genetic_data, functional_data)

def test_scenario_4_likely_pathogenic():
    """
    Scenario 4: Likely Pathogenic - Moderate evidence
    Rare missense with moderate pathogenic support
    """
    print("\n🧬 SCENARIO 4: Likely Pathogenic - Moderate Evidence")
    print("="*60)
    print("SCN5A c.5457G>A (p.Arg1819Gln) - Arrhythmia-associated")
    
    basic_info = {
        'gene': 'SCN5A',
        'chromosome': '3',
        'position': 38589532,
        'ref_allele': 'G',
        'alt_allele': 'A',
        'cdna_change': 'c.5457G>A',
        'protein_change': 'p.Arg1819Gln',
        'variant_type': 'missense',
        'transcript': 'NM_198056.2',
        'vep_consequence': 'missense_variant'
    }
    
    population_data = {
        'gnomad_af': 8.5e-06,
        'gnomad_af_popmax': 1.2e-05,
        'disease_prevalence': 5e-05, # Cardiac arrhythmia
        'gnomad_ac': 6,
        'gnomad_an': 705000,
        'gnomad_hom': 0
    }
    
    # Moderately pathogenic in silico predictions
    insilico_data = {
        'revel': 0.68,           # Moderate pathogenic
        'cadd_phred': 28.5,      # Moderate-high
        'clinpred': 0.72,        # Pathogenic
        'bayesdel_addaf': 0.45,  # Moderate pathogenic
        'alphamissense': 0.67,   # Pathogenic
        'mutationtaster': 0.85,  # Pathogenic
        'polyphen2_hdiv': 0.88,  # Probably damaging
        'sift': 0.02,            # Damaging
        'fathmm_xf': 0.75,       # Pathogenic
        'phylop_vertebrates': 4.2, # Conserved
        'phylop_mammals': 4.8,
        'phylop_primates': 4.5,
        'gerp_rs': 5.2,
        'metarnn': 0.71,
        'provean': -4.2,         # Damaging
        'lrt': 0.78
    }
    
    genetic_data = {
        'inheritance_pattern': 'ad',
        'zygosity': 'heterozygous',
        'family_history': 'yes',
        'segregation_analysis': 'partial',
        'de_novo': 'unknown',
        'parental_testing': 'unknown'
    }
    
    functional_data = {
        'functional_test': 'unknown',
        'phenotype_match': 'moderate'
    }
    
    return VariantData(basic_info, population_data, insilico_data, genetic_data, functional_data)

def test_scenario_5_benign():
    """
    Scenario 5: Benign - Strong benign evidence
    High frequency synonymous variant
    """
    print("\n🧬 SCENARIO 5: Benign - Synonymous High Frequency")
    print("="*60)
    print("LDLR c.1959C>T (p.Ile653=) - Common synonymous variant")
    
    basic_info = {
        'gene': 'LDLR',
        'chromosome': '19',
        'position': 11227255,
        'ref_allele': 'C',
        'alt_allele': 'T',
        'cdna_change': 'c.1959C>T',
        'protein_change': 'p.Ile653=',
        'variant_type': 'synonymous',
        'transcript': 'NM_000527.4',
        'vep_consequence': 'synonymous_variant'
    }
    
    population_data = {
        'gnomad_af': 0.15,           # Common (15%)
        'gnomad_af_popmax': 0.18,
        'disease_prevalence': 2e-04, # Familial hypercholesterolemia
        'gnomad_ac': 75000,
        'gnomad_an': 500000,
        'gnomad_hom': 8500
    }
    
    # Benign scores (not applicable for synonymous but testing algorithm)
    insilico_data = {
        'phylop_vertebrates': 0.8,   # Moderate conservation
        'phylop_mammals': 1.2,
        'phylop_primates': 0.9,
        'gerp_rs': 2.1
    }
    
    genetic_data = {
        'inheritance_pattern': 'ad',
        'zygosity': 'heterozygous',
        'family_history': 'no',
        'segregation_analysis': 'no',
        'de_novo': 'no'
    }
    
    functional_data = {
        'functional_test': 'unknown',
        'phenotype_match': 'no_match'
    }
    
    return VariantData(basic_info, population_data, insilico_data, genetic_data, functional_data)

def run_comprehensive_test():
    """Run all test scenarios and compare results."""
    
    scenarios = [
        ("VUS - Mixed Evidence", test_scenario_1_vus()),
        ("Likely Benign - Common Variant", test_scenario_2_likely_benign()),
        ("Pathogenic - Strong Evidence", test_scenario_3_pathogenic()),
        ("Likely Pathogenic - Moderate Evidence", test_scenario_4_likely_pathogenic()),
        ("Benign - Synonymous High Frequency", test_scenario_5_benign())
    ]
    
    results = []
    
    for scenario_name, variant_data in scenarios:
        print(f"\n{'='*80}")
        print(f"🧪 TESTING: {scenario_name}")
        print(f"{'='*80}")
        
        # Evaluate evidence
        evidence_evaluator = EvidenceEvaluator(use_2023_guidelines=False)
        evidence_results = evidence_evaluator.evaluate_all_criteria(variant_data)
        
        # Classify variant
        classifier = ACMGClassifier(use_2023_guidelines=False)
        classification_result = classifier.classify(evidence_results)
        
        # Store results
        vampp_score = evidence_results.get('vampp_score', 'N/A')
        results.append({
            'scenario': scenario_name,
            'variant': f"{variant_data.basic_info['gene']} {variant_data.basic_info['protein_change']}",
            'classification': classification_result['classification'],
            'confidence': classification_result['confidence'],
            'vampp_score': vampp_score,
            'applied_criteria': classification_result['applied_criteria']
        })
        
        # Display results
        print(f"📊 CLASSIFICATION: {classification_result['classification']}")
        print(f"🔍 CONFIDENCE: {classification_result['confidence']}")
        print(f"🧮 VAMPP-SCORE: {vampp_score}")
        print(f"📋 CRITERIA:")
        for category, criteria in classification_result['applied_criteria'].items():
            if criteria:
                print(f"  {category}: {', '.join(criteria)}")
    
    # Summary table
    print(f"\n{'='*100}")
    print("📈 COMPREHENSIVE TEST RESULTS SUMMARY")
    print(f"{'='*100}")
    print(f"{'Scenario':<35} {'Variant':<25} {'Classification':<18} {'VAMPP-Score':<12} {'Confidence'}")
    print(f"{'-'*100}")
    
    for result in results:
        vampp_display = f"{result['vampp_score']:.3f}" if isinstance(result['vampp_score'], float) else str(result['vampp_score'])
        print(f"{result['scenario']:<35} {result['variant']:<25} {result['classification']:<18} {vampp_display:<12} {result['confidence']}")
    
    print(f"{'-'*100}")
    print("✅ Comprehensive testing completed!")
    
    return results

if __name__ == "__main__":
    results = run_comprehensive_test()
    
    # Additional analysis
    print(f"\n📊 ALGORITHM PERFORMANCE ANALYSIS:")
    print(f"{'='*50}")
    
    vampp_scores = [r['vampp_score'] for r in results if isinstance(r['vampp_score'], float)]
    if vampp_scores:
        print(f"VAMPP-Score Range: {min(vampp_scores):.3f} - {max(vampp_scores):.3f}")
        print(f"Average VAMPP-Score: {sum(vampp_scores)/len(vampp_scores):.3f}")
    
    classifications = [r['classification'] for r in results]
    print(f"Classifications Distribution:")
    for classification in set(classifications):
        count = classifications.count(classification)
        print(f"  {classification}: {count}")
