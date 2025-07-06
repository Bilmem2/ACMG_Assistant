#!/usr/bin/env python3
"""
VUS (Variant of Uncertain Significance) Test Scenarios
====================================================

Comprehensive test scenarios for VUS classification with different in silico score combinations.
Tests the robustness of the metascore algorithm and ACMG classification.

Author: Can Sevilmiş
Last Updated: July 6, 2025
"""

import sys
import os
from colorama import init, Fore, Style
from datetime import datetime

# Initialize colorama
init(autoreset=True)

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.variant_data import VariantData
from core.evidence_evaluator import EvidenceEvaluator
from core.acmg_classifier import ACMGClassifier
from utils.report_generator import ReportGenerator

class VUSTestScenarios:
    """Test scenarios for VUS classification."""
    
    def __init__(self):
        self.evidence_evaluator = EvidenceEvaluator(use_2023_guidelines=False)
        self.classifier = ACMGClassifier(use_2023_guidelines=False)
        self.report_generator = ReportGenerator()
    
    def run_all_scenarios(self):
        """Run all VUS test scenarios."""
        print(f"\n{Fore.CYAN}{'='*80}")
        print(f"{Fore.CYAN}{Style.BRIGHT}🧬 VUS TEST SCENARIOS - METASCORE ALGORITHM VALIDATION")
        print(f"{Fore.CYAN}{'='*80}")
        print(f"{Fore.YELLOW}Test Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{Fore.CYAN}{'='*80}")
        
        scenarios = [
            self.scenario_1_high_conservation_low_pathogenicity,
            self.scenario_2_mixed_predictions,
            self.scenario_3_weak_evidence_combination,
            self.scenario_4_conservation_dominant,
            self.scenario_5_pathogenicity_dominant,
            self.scenario_6_conflicting_scores,
            self.scenario_7_extreme_conservation,
            self.scenario_8_minimal_evidence
        ]
        
        results = []
        
        for i, scenario_func in enumerate(scenarios, 1):
            print(f"\n{Fore.MAGENTA}{'='*60}")
            print(f"{Fore.MAGENTA}{Style.BRIGHT}📊 SCENARIO {i}")
            print(f"{Fore.MAGENTA}{'='*60}")
            
            result = scenario_func()
            results.append(result)
            
            # Brief summary
            print(f"{Fore.GREEN}✅ Classification: {result['classification']}")
            print(f"{Fore.CYAN}📈 Confidence: {result['confidence']}")
            print(f"{Fore.BLUE}🔍 Applied Criteria: {sum(len(criteria) for criteria in result['applied_criteria'].values())}")
        
        # Final summary
        self._print_summary(results)
        
        return results
    
    def scenario_1_high_conservation_low_pathogenicity(self):
        """
        Scenario 1: High conservation scores but low pathogenicity predictions
        Expected: Should lean towards VUS with slight pathogenic tendency due to conservation
        """
        print(f"{Fore.YELLOW}📝 Scenario 1: High Conservation, Low Pathogenicity")
        print(f"{Fore.WHITE}Gene: TP53 (tumor suppressor)")
        print(f"{Fore.WHITE}Variant: c.742C>T (p.Arg248Trp)")
        print(f"{Fore.WHITE}Known ClinVar: VUS")
        print(f"{Fore.WHITE}Characteristics: Strong conservation, weak pathogenicity predictors")
        
        variant_data = VariantData(
            basic_info={
                'gene': 'TP53',
                'chromosome': '17',
                'position': '7673776',
                'ref_allele': 'C',
                'alt_allele': 'T',
                'hgvs_c': 'c.742C>T',
                'hgvs_p': 'p.Arg248Trp',
                'variant_type': 'missense',
                'transcript': 'NM_000546.5',
                'exon': '7',
                'protein_change': 'p.Arg248Trp'
            },
            population_data={
                'gnomad_af': 0.00012,  # Rare but present
                'gnomad_af_popmax': 0.00089,
                'gnomad_an': 251456,
                'gnomad_ac': 30,
                'gnomad_hom': 0,
                'gnomad_hemi': 0,
                'exac_af': 0.00008,
                'thousand_genomes_af': 0.0001,
                'prevalence': 0.0001  # Common disease
            },
            insilico_data={
                # Strong conservation
                'phylop_100way': 8.2,  # Very high
                'phylop_30way': 7.8,   # Very high
                'phylop_17way': 6.9,   # High
                'phastcons_100way': 0.998,  # Very high
                'phastcons_30way': 0.995,   # Very high
                'phastcons_17way': 0.987,   # High
                'gerp_score': 5.8,     # Very high
                'siphy_score': 15.2,   # Very high
                
                # Weak pathogenicity predictors
                'sift_score': 0.12,        # Tolerated (>0.05)
                'sift_prediction': 'T',
                'polyphen2_hdiv_score': 0.654,  # Possibly damaging
                'polyphen2_hdiv_prediction': 'P',
                'polyphen2_hvar_score': 0.598,  # Possibly damaging
                'polyphen2_hvar_prediction': 'P',
                'lrt_score': 0.234,        # Neutral
                'lrt_prediction': 'N',
                'mutationtaster_score': 0.678,  # Disease causing
                'mutationtaster_prediction': 'D',
                'mutationassessor_score': 2.89,  # Medium
                'mutationassessor_prediction': 'M',
                'fathmm_score': -0.89,     # Tolerated
                'fathmm_prediction': 'T',
                'provean_score': -1.89,    # Neutral
                'provean_prediction': 'N',
                'vest4_score': 0.423,      # Below threshold
                'metasvm_score': 0.234,    # Tolerated
                'metasvm_prediction': 'T',
                'metalr_score': 0.198,     # Tolerated
                'metalr_prediction': 'T',
                'integrated_fitcons_score': 0.234,
                'gm12878_fitcons_score': 0.198,
                'h1_hesc_fitcons_score': 0.276,
                'huvec_fitcons_score': 0.198
            },
            genetic_data={
                'mode_of_inheritance': 'autosomal_dominant',
                'penetrance': 'high',
                'allelic_requirement': 'monoallelic',
                'disease_mechanism': 'loss_of_function',
                'variant_consequence': 'missense_variant',
                'variant_effect': 'uncertain'
            },
            functional_data={
                'protein_domain': 'DNA-binding domain',
                'protein_structure_affected': True,
                'catalytic_residue': False,
                'binding_site': True,
                'known_functional_domain': True,
                'previous_functional_studies': False,
                'protein_stability_affected': 'unknown',
                'enzymatic_activity_affected': 'unknown'
            }
        )
        
        return self._evaluate_and_classify(variant_data, "Scenario 1")
    
    def scenario_2_mixed_predictions(self):
        """
        Scenario 2: Mixed pathogenicity predictions with moderate conservation
        Expected: Should be classified as VUS with balanced evidence
        """
        print(f"{Fore.YELLOW}📝 Scenario 2: Mixed Pathogenicity Predictions")
        print(f"{Fore.WHITE}Gene: BRCA1 (breast cancer susceptibility)")
        print(f"{Fore.WHITE}Variant: c.2077G>A (p.Asp693Asn)")
        print(f"{Fore.WHITE}Known ClinVar: VUS")
        print(f"{Fore.WHITE}Characteristics: Mixed predictions, moderate conservation")
        
        variant_data = VariantData(
            basic_info={
                'gene': 'BRCA1',
                'chromosome': '17',
                'position': '41244618',
                'ref_allele': 'G',
                'alt_allele': 'A',
                'hgvs_c': 'c.2077G>A',
                'hgvs_p': 'p.Asp693Asn',
                'variant_type': 'missense',
                'transcript': 'NM_007294.3',
                'exon': '11',
                'protein_change': 'p.Asp693Asn'
            },
            population_data={
                'gnomad_af': 0.00023,
                'gnomad_af_popmax': 0.00156,
                'gnomad_an': 251456,
                'gnomad_ac': 58,
                'gnomad_hom': 0,
                'gnomad_hemi': 0,
                'exac_af': 0.00019,
                'thousand_genomes_af': 0.0002,
                'prevalence': 0.00125  # Breast cancer prevalence
            },
            insilico_data={
                # Moderate conservation
                'phylop_100way': 4.2,
                'phylop_30way': 3.8,
                'phylop_17way': 3.1,
                'phastcons_100way': 0.756,
                'phastcons_30way': 0.698,
                'phastcons_17way': 0.645,
                'gerp_score': 3.2,
                'siphy_score': 8.9,
                
                # Mixed pathogenicity predictors
                'sift_score': 0.02,        # Damaging
                'sift_prediction': 'D',
                'polyphen2_hdiv_score': 0.842,  # Probably damaging
                'polyphen2_hdiv_prediction': 'D',
                'polyphen2_hvar_score': 0.798,  # Probably damaging
                'polyphen2_hvar_prediction': 'D',
                'lrt_score': 0.789,        # Deleterious
                'lrt_prediction': 'D',
                'mutationtaster_score': 0.456,  # Borderline
                'mutationtaster_prediction': 'N',
                'mutationassessor_score': 2.34,  # Medium
                'mutationassessor_prediction': 'M',
                'fathmm_score': 0.45,      # Tolerated
                'fathmm_prediction': 'T',
                'provean_score': -3.67,    # Deleterious
                'provean_prediction': 'D',
                'vest4_score': 0.567,      # Moderate
                'metasvm_score': 0.678,    # Damaging
                'metasvm_prediction': 'D',
                'metalr_score': 0.543,     # Damaging
                'metalr_prediction': 'D',
                'integrated_fitcons_score': 0.456,
                'gm12878_fitcons_score': 0.398,
                'h1_hesc_fitcons_score': 0.512,
                'huvec_fitcons_score': 0.423
            },
            genetic_data={
                'mode_of_inheritance': 'autosomal_dominant',
                'penetrance': 'incomplete',
                'allelic_requirement': 'monoallelic',
                'disease_mechanism': 'loss_of_function',
                'variant_consequence': 'missense_variant',
                'variant_effect': 'uncertain'
            },
            functional_data={
                'protein_domain': 'BRCT domain',
                'protein_structure_affected': True,
                'catalytic_residue': False,
                'binding_site': False,
                'known_functional_domain': True,
                'previous_functional_studies': False,
                'protein_stability_affected': 'unknown',
                'enzymatic_activity_affected': 'unknown'
            }
        )
        
        return self._evaluate_and_classify(variant_data, "Scenario 2")
    
    def scenario_3_weak_evidence_combination(self):
        """
        Scenario 3: Weak evidence across all categories
        Expected: Should be classified as VUS with low confidence
        """
        print(f"{Fore.YELLOW}📝 Scenario 3: Weak Evidence Combination")
        print(f"{Fore.WHITE}Gene: LDLR (familial hypercholesterolemia)")
        print(f"{Fore.WHITE}Variant: c.1327A>G (p.Ile443Val)")
        print(f"{Fore.WHITE}Known ClinVar: VUS")
        print(f"{Fore.WHITE}Characteristics: Weak evidence across all categories")
        
        variant_data = VariantData(
            basic_info={
                'gene': 'LDLR',
                'chromosome': '19',
                'position': '11218090',
                'ref_allele': 'A',
                'alt_allele': 'G',
                'hgvs_c': 'c.1327A>G',
                'hgvs_p': 'p.Ile443Val',
                'variant_type': 'missense',
                'transcript': 'NM_000527.4',
                'exon': '9',
                'protein_change': 'p.Ile443Val'
            },
            population_data={
                'gnomad_af': 0.00089,  # Slightly higher frequency
                'gnomad_af_popmax': 0.00345,
                'gnomad_an': 251456,
                'gnomad_ac': 224,
                'gnomad_hom': 0,
                'gnomad_hemi': 0,
                'exac_af': 0.00076,
                'thousand_genomes_af': 0.0008,
                'prevalence': 0.002  # Familial hypercholesterolemia
            },
            insilico_data={
                # Weak conservation
                'phylop_100way': 1.8,
                'phylop_30way': 1.2,
                'phylop_17way': 0.9,
                'phastcons_100way': 0.432,
                'phastcons_30way': 0.345,
                'phastcons_17way': 0.298,
                'gerp_score': 1.1,
                'siphy_score': 3.2,
                
                # Weak pathogenicity predictors
                'sift_score': 0.089,       # Borderline
                'sift_prediction': 'T',
                'polyphen2_hdiv_score': 0.387,  # Benign
                'polyphen2_hdiv_prediction': 'B',
                'polyphen2_hvar_score': 0.234,  # Benign
                'polyphen2_hvar_prediction': 'B',
                'lrt_score': 0.123,        # Neutral
                'lrt_prediction': 'N',
                'mutationtaster_score': 0.234,  # Polymorphism
                'mutationtaster_prediction': 'P',
                'mutationassessor_score': 1.67,  # Low
                'mutationassessor_prediction': 'L',
                'fathmm_score': -1.23,     # Tolerated
                'fathmm_prediction': 'T',
                'provean_score': -0.89,    # Neutral
                'provean_prediction': 'N',
                'vest4_score': 0.234,      # Low
                'metasvm_score': 0.123,    # Tolerated
                'metasvm_prediction': 'T',
                'metalr_score': 0.089,     # Tolerated
                'metalr_prediction': 'T',
                'integrated_fitcons_score': 0.198,
                'gm12878_fitcons_score': 0.156,
                'h1_hesc_fitcons_score': 0.234,
                'huvec_fitcons_score': 0.178
            },
            genetic_data={
                'mode_of_inheritance': 'autosomal_dominant',
                'penetrance': 'incomplete',
                'allelic_requirement': 'monoallelic',
                'disease_mechanism': 'loss_of_function',
                'variant_consequence': 'missense_variant',
                'variant_effect': 'uncertain'
            },
            functional_data={
                'protein_domain': 'LDL-receptor class A domain',
                'protein_structure_affected': False,
                'catalytic_residue': False,
                'binding_site': False,
                'known_functional_domain': True,
                'previous_functional_studies': False,
                'protein_stability_affected': 'unknown',
                'enzymatic_activity_affected': 'unknown'
            }
        )
        
        return self._evaluate_and_classify(variant_data, "Scenario 3")
    
    def scenario_4_conservation_dominant(self):
        """
        Scenario 4: Extremely high conservation, mixed pathogenicity
        Expected: Should lean towards likely pathogenic due to conservation
        """
        print(f"{Fore.YELLOW}📝 Scenario 4: Conservation Dominant")
        print(f"{Fore.WHITE}Gene: MECP2 (Rett syndrome)")
        print(f"{Fore.WHITE}Variant: c.473C>T (p.Thr158Met)")
        print(f"{Fore.WHITE}Known ClinVar: VUS")
        print(f"{Fore.WHITE}Characteristics: Extremely high conservation")
        
        variant_data = VariantData(
            basic_info={
                'gene': 'MECP2',
                'chromosome': 'X',
                'position': '154030662',
                'ref_allele': 'C',
                'alt_allele': 'T',
                'hgvs_c': 'c.473C>T',
                'hgvs_p': 'p.Thr158Met',
                'variant_type': 'missense',
                'transcript': 'NM_004992.3',
                'exon': '4',
                'protein_change': 'p.Thr158Met'
            },
            population_data={
                'gnomad_af': 0.000032,  # Very rare
                'gnomad_af_popmax': 0.00023,
                'gnomad_an': 251456,
                'gnomad_ac': 8,
                'gnomad_hom': 0,
                'gnomad_hemi': 0,
                'exac_af': 0.000019,
                'thousand_genomes_af': 0.00002,
                'prevalence': 0.00001  # Rett syndrome prevalence
            },
            insilico_data={
                # Extremely high conservation
                'phylop_100way': 9.8,   # Maximum
                'phylop_30way': 9.5,    # Maximum
                'phylop_17way': 8.9,    # Very high
                'phastcons_100way': 1.0,    # Maximum
                'phastcons_30way': 1.0,     # Maximum
                'phastcons_17way': 0.998,   # Very high
                'gerp_score': 6.18,         # Maximum
                'siphy_score': 29.9,        # Very high
                
                # Mixed pathogenicity predictors
                'sift_score': 0.034,        # Damaging
                'sift_prediction': 'D',
                'polyphen2_hdiv_score': 0.567,  # Possibly damaging
                'polyphen2_hdiv_prediction': 'P',
                'polyphen2_hvar_score': 0.489,  # Possibly damaging
                'polyphen2_hvar_prediction': 'P',
                'lrt_score': 0.456,         # Neutral
                'lrt_prediction': 'N',
                'mutationtaster_score': 0.789,  # Disease causing
                'mutationtaster_prediction': 'D',
                'mutationassessor_score': 3.12,  # High
                'mutationassessor_prediction': 'H',
                'fathmm_score': -0.12,      # Tolerated
                'fathmm_prediction': 'T',
                'provean_score': -2.34,     # Deleterious
                'provean_prediction': 'D',
                'vest4_score': 0.678,       # Moderate
                'metasvm_score': 0.456,     # Borderline
                'metasvm_prediction': 'T',
                'metalr_score': 0.398,      # Tolerated
                'metalr_prediction': 'T',
                'integrated_fitcons_score': 0.567,
                'gm12878_fitcons_score': 0.489,
                'h1_hesc_fitcons_score': 0.623,
                'huvec_fitcons_score': 0.512
            },
            genetic_data={
                'mode_of_inheritance': 'x_linked_dominant',
                'penetrance': 'complete',
                'allelic_requirement': 'monoallelic',
                'disease_mechanism': 'loss_of_function',
                'variant_consequence': 'missense_variant',
                'variant_effect': 'uncertain'
            },
            functional_data={
                'protein_domain': 'Methyl-CpG binding domain',
                'protein_structure_affected': True,
                'catalytic_residue': False,
                'binding_site': True,
                'known_functional_domain': True,
                'previous_functional_studies': False,
                'protein_stability_affected': 'unknown',
                'enzymatic_activity_affected': 'unknown'
            }
        )
        
        return self._evaluate_and_classify(variant_data, "Scenario 4")
    
    def scenario_5_pathogenicity_dominant(self):
        """
        Scenario 5: Strong pathogenicity predictions, moderate conservation
        Expected: Should lean towards likely pathogenic due to pathogenicity scores
        """
        print(f"{Fore.YELLOW}📝 Scenario 5: Pathogenicity Dominant")
        print(f"{Fore.WHITE}Gene: CFTR (cystic fibrosis)")
        print(f"{Fore.WHITE}Variant: c.1822G>A (p.Gly608Ser)")
        print(f"{Fore.WHITE}Known ClinVar: VUS")
        print(f"{Fore.WHITE}Characteristics: Strong pathogenicity predictions")
        
        variant_data = VariantData(
            basic_info={
                'gene': 'CFTR',
                'chromosome': '7',
                'position': '117307161',
                'ref_allele': 'G',
                'alt_allele': 'A',
                'hgvs_c': 'c.1822G>A',
                'hgvs_p': 'p.Gly608Ser',
                'variant_type': 'missense',
                'transcript': 'NM_000492.3',
                'exon': '13',
                'protein_change': 'p.Gly608Ser'
            },
            population_data={
                'gnomad_af': 0.000067,  # Rare
                'gnomad_af_popmax': 0.00045,
                'gnomad_an': 251456,
                'gnomad_ac': 17,
                'gnomad_hom': 0,
                'gnomad_hemi': 0,
                'exac_af': 0.000052,
                'thousand_genomes_af': 0.00006,
                'prevalence': 0.0004  # Cystic fibrosis prevalence
            },
            insilico_data={
                # Moderate conservation
                'phylop_100way': 5.2,
                'phylop_30way': 4.8,
                'phylop_17way': 4.1,
                'phastcons_100way': 0.823,
                'phastcons_30way': 0.756,
                'phastcons_17way': 0.698,
                'gerp_score': 4.1,
                'siphy_score': 12.3,
                
                # Strong pathogenicity predictors
                'sift_score': 0.001,        # Damaging
                'sift_prediction': 'D',
                'polyphen2_hdiv_score': 0.998,  # Probably damaging
                'polyphen2_hdiv_prediction': 'D',
                'polyphen2_hvar_score': 0.987,  # Probably damaging
                'polyphen2_hvar_prediction': 'D',
                'lrt_score': 0.989,         # Deleterious
                'lrt_prediction': 'D',
                'mutationtaster_score': 0.999,  # Disease causing
                'mutationtaster_prediction': 'D',
                'mutationassessor_score': 4.67,  # High
                'mutationassessor_prediction': 'H',
                'fathmm_score': 3.45,       # Damaging
                'fathmm_prediction': 'D',
                'provean_score': -6.78,     # Deleterious
                'provean_prediction': 'D',
                'vest4_score': 0.892,       # High
                'metasvm_score': 0.956,     # Damaging
                'metasvm_prediction': 'D',
                'metalr_score': 0.934,      # Damaging
                'metalr_prediction': 'D',
                'integrated_fitcons_score': 0.789,
                'gm12878_fitcons_score': 0.823,
                'h1_hesc_fitcons_score': 0.856,
                'huvec_fitcons_score': 0.798
            },
            genetic_data={
                'mode_of_inheritance': 'autosomal_recessive',
                'penetrance': 'complete',
                'allelic_requirement': 'biallelic',
                'disease_mechanism': 'loss_of_function',
                'variant_consequence': 'missense_variant',
                'variant_effect': 'uncertain'
            },
            functional_data={
                'protein_domain': 'NBD1 domain',
                'protein_structure_affected': True,
                'catalytic_residue': False,
                'binding_site': True,
                'known_functional_domain': True,
                'previous_functional_studies': False,
                'protein_stability_affected': 'unknown',
                'enzymatic_activity_affected': 'unknown'
            }
        )
        
        return self._evaluate_and_classify(variant_data, "Scenario 5")
    
    def scenario_6_conflicting_scores(self):
        """
        Scenario 6: Conflicting in silico scores
        Expected: Should be classified as VUS due to conflicting evidence
        """
        print(f"{Fore.YELLOW}📝 Scenario 6: Conflicting Scores")
        print(f"{Fore.WHITE}Gene: SCN5A (Brugada syndrome)")
        print(f"{Fore.WHITE}Variant: c.2389G>A (p.Asp797Asn)")
        print(f"{Fore.WHITE}Known ClinVar: VUS")
        print(f"{Fore.WHITE}Characteristics: Highly conflicting in silico scores")
        
        variant_data = VariantData(
            basic_info={
                'gene': 'SCN5A',
                'chromosome': '3',
                'position': '38589446',
                'ref_allele': 'G',
                'alt_allele': 'A',
                'hgvs_c': 'c.2389G>A',
                'hgvs_p': 'p.Asp797Asn',
                'variant_type': 'missense',
                'transcript': 'NM_198056.2',
                'exon': '15',
                'protein_change': 'p.Asp797Asn'
            },
            population_data={
                'gnomad_af': 0.00015,
                'gnomad_af_popmax': 0.00098,
                'gnomad_an': 251456,
                'gnomad_ac': 38,
                'gnomad_hom': 0,
                'gnomad_hemi': 0,
                'exac_af': 0.00012,
                'thousand_genomes_af': 0.0001,
                'prevalence': 0.0002  # Brugada syndrome
            },
            insilico_data={
                # Conflicting conservation scores
                'phylop_100way': 7.8,   # High
                'phylop_30way': 2.1,    # Low
                'phylop_17way': 5.6,    # Moderate
                'phastcons_100way': 0.95,   # High
                'phastcons_30way': 0.23,    # Low
                'phastcons_17way': 0.67,    # Moderate
                'gerp_score': 5.2,          # High
                'siphy_score': 4.1,         # Low
                
                # Conflicting pathogenicity predictors
                'sift_score': 0.001,        # Damaging
                'sift_prediction': 'D',
                'polyphen2_hdiv_score': 0.156,  # Benign
                'polyphen2_hdiv_prediction': 'B',
                'polyphen2_hvar_score': 0.876,  # Probably damaging
                'polyphen2_hvar_prediction': 'D',
                'lrt_score': 0.123,         # Neutral
                'lrt_prediction': 'N',
                'mutationtaster_score': 0.987,  # Disease causing
                'mutationtaster_prediction': 'D',
                'mutationassessor_score': 1.23,  # Low
                'mutationassessor_prediction': 'L',
                'fathmm_score': 2.34,       # Damaging
                'fathmm_prediction': 'D',
                'provean_score': -0.56,     # Neutral
                'provean_prediction': 'N',
                'vest4_score': 0.789,       # High
                'metasvm_score': 0.234,     # Tolerated
                'metasvm_prediction': 'T',
                'metalr_score': 0.823,      # Damaging
                'metalr_prediction': 'D',
                'integrated_fitcons_score': 0.345,
                'gm12878_fitcons_score': 0.789,
                'h1_hesc_fitcons_score': 0.123,
                'huvec_fitcons_score': 0.656
            },
            genetic_data={
                'mode_of_inheritance': 'autosomal_dominant',
                'penetrance': 'incomplete',
                'allelic_requirement': 'monoallelic',
                'disease_mechanism': 'loss_of_function',
                'variant_consequence': 'missense_variant',
                'variant_effect': 'uncertain'
            },
            functional_data={
                'protein_domain': 'Voltage-gated sodium channel',
                'protein_structure_affected': True,
                'catalytic_residue': False,
                'binding_site': False,
                'known_functional_domain': True,
                'previous_functional_studies': False,
                'protein_stability_affected': 'unknown',
                'enzymatic_activity_affected': 'unknown'
            }
        )
        
        return self._evaluate_and_classify(variant_data, "Scenario 6")
    
    def scenario_7_extreme_conservation(self):
        """
        Scenario 7: Extreme conservation but benign predictions
        Expected: Should be interesting to see how algorithm handles this conflict
        """
        print(f"{Fore.YELLOW}📝 Scenario 7: Extreme Conservation vs Benign Predictions")
        print(f"{Fore.WHITE}Gene: ACTB (cytoskeletal beta-actin)")
        print(f"{Fore.WHITE}Variant: c.689A>G (p.Gln230Arg)")
        print(f"{Fore.WHITE}Known ClinVar: VUS")
        print(f"{Fore.WHITE}Characteristics: Maximum conservation, benign predictions")
        
        variant_data = VariantData(
            basic_info={
                'gene': 'ACTB',
                'chromosome': '7',
                'position': '5566778',
                'ref_allele': 'A',
                'alt_allele': 'G',
                'hgvs_c': 'c.689A>G',
                'hgvs_p': 'p.Gln230Arg',
                'variant_type': 'missense',
                'transcript': 'NM_001101.3',
                'exon': '4',
                'protein_change': 'p.Gln230Arg'
            },
            population_data={
                'gnomad_af': 0.00045,  # Moderate frequency
                'gnomad_af_popmax': 0.00234,
                'gnomad_an': 251456,
                'gnomad_ac': 113,
                'gnomad_hom': 0,
                'gnomad_hemi': 0,
                'exac_af': 0.00039,
                'thousand_genomes_af': 0.0004,
                'prevalence': 0.0001  # Baraitser-Winter syndrome
            },
            insilico_data={
                # Maximum conservation (actin is highly conserved)
                'phylop_100way': 9.99,  # Maximum
                'phylop_30way': 9.95,   # Maximum
                'phylop_17way': 9.89,   # Maximum
                'phastcons_100way': 1.0,    # Maximum
                'phastcons_30way': 1.0,     # Maximum
                'phastcons_17way': 1.0,     # Maximum
                'gerp_score': 6.18,         # Maximum
                'siphy_score': 29.9,        # Maximum
                
                # Benign pathogenicity predictors
                'sift_score': 0.234,        # Tolerated
                'sift_prediction': 'T',
                'polyphen2_hdiv_score': 0.087,  # Benign
                'polyphen2_hdiv_prediction': 'B',
                'polyphen2_hvar_score': 0.123,  # Benign
                'polyphen2_hvar_prediction': 'B',
                'lrt_score': 0.056,         # Neutral
                'lrt_prediction': 'N',
                'mutationtaster_score': 0.089,  # Polymorphism
                'mutationtaster_prediction': 'P',
                'mutationassessor_score': 0.89,  # Neutral
                'mutationassessor_prediction': 'N',
                'fathmm_score': -2.34,      # Tolerated
                'fathmm_prediction': 'T',
                'provean_score': -0.23,     # Neutral
                'provean_prediction': 'N',
                'vest4_score': 0.123,       # Low
                'metasvm_score': 0.067,     # Tolerated
                'metasvm_prediction': 'T',
                'metalr_score': 0.089,      # Tolerated
                'metalr_prediction': 'T',
                'integrated_fitcons_score': 0.134,
                'gm12878_fitcons_score': 0.098,
                'h1_hesc_fitcons_score': 0.156,
                'huvec_fitcons_score': 0.123
            },
            genetic_data={
                'mode_of_inheritance': 'autosomal_dominant',
                'penetrance': 'complete',
                'allelic_requirement': 'monoallelic',
                'disease_mechanism': 'loss_of_function',
                'variant_consequence': 'missense_variant',
                'variant_effect': 'uncertain'
            },
            functional_data={
                'protein_domain': 'Actin domain',
                'protein_structure_affected': False,
                'catalytic_residue': False,
                'binding_site': False,
                'known_functional_domain': True,
                'previous_functional_studies': False,
                'protein_stability_affected': 'unknown',
                'enzymatic_activity_affected': 'unknown'
            }
        )
        
        return self._evaluate_and_classify(variant_data, "Scenario 7")
    
    def scenario_8_minimal_evidence(self):
        """
        Scenario 8: Minimal evidence available
        Expected: Should be classified as VUS with very low confidence
        """
        print(f"{Fore.YELLOW}📝 Scenario 8: Minimal Evidence")
        print(f"{Fore.WHITE}Gene: ABCC8 (neonatal diabetes)")
        print(f"{Fore.WHITE}Variant: c.1234G>A (p.Gly412Ser)")
        print(f"{Fore.WHITE}Known ClinVar: VUS")
        print(f"{Fore.WHITE}Characteristics: Very limited evidence across all categories")
        
        variant_data = VariantData(
            basic_info={
                'gene': 'ABCC8',
                'chromosome': '11',
                'position': '17414431',
                'ref_allele': 'G',
                'alt_allele': 'A',
                'hgvs_c': 'c.1234G>A',
                'hgvs_p': 'p.Gly412Ser',
                'variant_type': 'missense',
                'transcript': 'NM_000352.3',
                'exon': '8',
                'protein_change': 'p.Gly412Ser'
            },
            population_data={
                'gnomad_af': 0.0012,  # Higher frequency
                'gnomad_af_popmax': 0.0078,
                'gnomad_an': 251456,
                'gnomad_ac': 302,
                'gnomad_hom': 1,
                'gnomad_hemi': 0,
                'exac_af': 0.0009,
                'thousand_genomes_af': 0.001,
                'prevalence': 0.00005  # Neonatal diabetes
            },
            insilico_data={
                # Very limited conservation data
                'phylop_100way': 0.5,
                'phylop_30way': 0.3,
                'phylop_17way': 0.1,
                'phastcons_100way': 0.156,
                'phastcons_30way': 0.089,
                'phastcons_17way': 0.067,
                'gerp_score': 0.2,
                'siphy_score': 1.1,
                
                # Very limited pathogenicity predictors
                'sift_score': 0.156,        # Tolerated
                'sift_prediction': 'T',
                'polyphen2_hdiv_score': 0.234,  # Benign
                'polyphen2_hdiv_prediction': 'B',
                'polyphen2_hvar_score': 0.189,  # Benign
                'polyphen2_hvar_prediction': 'B',
                'lrt_score': 0.089,         # Neutral
                'lrt_prediction': 'N',
                'mutationtaster_score': 0.123,  # Polymorphism
                'mutationtaster_prediction': 'P',
                'mutationassessor_score': 0.67,  # Neutral
                'mutationassessor_prediction': 'N',
                'fathmm_score': -1.89,      # Tolerated
                'fathmm_prediction': 'T',
                'provean_score': -0.45,     # Neutral
                'provean_prediction': 'N',
                'vest4_score': 0.089,       # Very low
                'metasvm_score': 0.034,     # Tolerated
                'metasvm_prediction': 'T',
                'metalr_score': 0.056,      # Tolerated
                'metalr_prediction': 'T',
                'integrated_fitcons_score': 0.067,
                'gm12878_fitcons_score': 0.045,
                'h1_hesc_fitcons_score': 0.089,
                'huvec_fitcons_score': 0.056
            },
            genetic_data={
                'mode_of_inheritance': 'autosomal_recessive',
                'penetrance': 'complete',
                'allelic_requirement': 'biallelic',
                'disease_mechanism': 'loss_of_function',
                'variant_consequence': 'missense_variant',
                'variant_effect': 'uncertain'
            },
            functional_data={
                'protein_domain': 'ABC transporter',
                'protein_structure_affected': False,
                'catalytic_residue': False,
                'binding_site': False,
                'known_functional_domain': True,
                'previous_functional_studies': False,
                'protein_stability_affected': 'unknown',
                'enzymatic_activity_affected': 'unknown'
            }
        )
        
        return self._evaluate_and_classify(variant_data, "Scenario 8")
    
    def _evaluate_and_classify(self, variant_data: VariantData, scenario_name: str) -> dict:
        """Evaluate evidence and classify variant for a scenario."""
        # Evaluate evidence
        evidence_results = self.evidence_evaluator.evaluate_all_criteria(variant_data)
        
        # Classify variant
        classification_result = self.classifier.classify(evidence_results)
        
        # Print detailed results
        print(f"\n{Fore.CYAN}🔍 Evidence Analysis:")
        vampp_score = evidence_results.get('vampp_score', 0)
        conservation_score = evidence_results.get('conservation_score', 0)
        pathogenicity_score = evidence_results.get('pathogenicity_score', 0)
        
        print(f"{Fore.WHITE}VAMPP Score: {vampp_score:.3f}")
        print(f"{Fore.WHITE}Conservation Score: {conservation_score:.3f}")
        print(f"{Fore.WHITE}Pathogenicity Score: {pathogenicity_score:.3f}")
        
        print(f"\n{Fore.GREEN}📊 Classification Results:")
        print(f"{Fore.WHITE}Classification: {classification_result['classification']}")
        print(f"{Fore.WHITE}Confidence: {classification_result['confidence']}")
        
        # Show applied criteria
        total_criteria = sum(len(criteria) for criteria in classification_result['applied_criteria'].values())
        print(f"{Fore.WHITE}Applied Criteria ({total_criteria} total):")
        for category, criteria in classification_result['applied_criteria'].items():
            if criteria:
                print(f"  {category}: {', '.join(criteria)}")
        
        # Generate report
        report_path = self.report_generator.generate_report(
            variant_data, evidence_results, classification_result
        )
        
        # Add scenario info to result
        classification_result['scenario_name'] = scenario_name
        classification_result['report_path'] = report_path
        classification_result['vampp_score'] = vampp_score
        classification_result['conservation_score'] = conservation_score
        classification_result['pathogenicity_score'] = pathogenicity_score
        
        return classification_result
    
    def _print_summary(self, results: list):
        """Print summary of all scenario results."""
        print(f"\n{Fore.CYAN}{'='*80}")
        print(f"{Fore.CYAN}{Style.BRIGHT}📊 VUS SCENARIOS SUMMARY")
        print(f"{Fore.CYAN}{'='*80}")
        
        # Classification distribution
        classifications = {}
        for result in results:
            classification = result['classification']
            if classification not in classifications:
                classifications[classification] = 0
            classifications[classification] += 1
        
        print(f"\n{Fore.YELLOW}{Style.BRIGHT}📈 Classification Distribution:")
        for classification, count in classifications.items():
            print(f"{Fore.WHITE}  {classification}: {count} scenario(s)")
        
        # Confidence analysis
        print(f"\n{Fore.YELLOW}{Style.BRIGHT}🔍 Confidence Analysis:")
        confidence_levels = {}
        for result in results:
            confidence = result['confidence']
            if confidence not in confidence_levels:
                confidence_levels[confidence] = 0
            confidence_levels[confidence] += 1
        
        for confidence, count in confidence_levels.items():
            print(f"{Fore.WHITE}  {confidence}: {count} scenario(s)")
        
        # Score analysis
        print(f"\n{Fore.YELLOW}{Style.BRIGHT}🎯 Score Analysis:")
        vampp_scores = [r['vampp_score'] for r in results]
        conservation_scores = [r['conservation_score'] for r in results]
        pathogenicity_scores = [r['pathogenicity_score'] for r in results]
        
        print(f"{Fore.WHITE}  VAMPP Score Range: {min(vampp_scores):.3f} - {max(vampp_scores):.3f}")
        print(f"{Fore.WHITE}  Conservation Score Range: {min(conservation_scores):.3f} - {max(conservation_scores):.3f}")
        print(f"{Fore.WHITE}  Pathogenicity Score Range: {min(pathogenicity_scores):.3f} - {max(pathogenicity_scores):.3f}")
        
        # Most interesting scenarios
        print(f"\n{Fore.YELLOW}{Style.BRIGHT}🎪 Most Interesting Scenarios:")
        
        # Highest VAMPP score
        highest_vampp = max(results, key=lambda x: x['vampp_score'])
        print(f"{Fore.WHITE}  Highest VAMPP Score: {highest_vampp['scenario_name']} ({highest_vampp['vampp_score']:.3f})")
        
        # Most criteria applied
        most_criteria = max(results, key=lambda x: sum(len(criteria) for criteria in x['applied_criteria'].values()))
        total_criteria = sum(len(criteria) for criteria in most_criteria['applied_criteria'].values())
        print(f"{Fore.WHITE}  Most Criteria Applied: {most_criteria['scenario_name']} ({total_criteria} criteria)")
        
        # Most conflicting (VUS with high scores)
        vus_results = [r for r in results if 'uncertain' in r['classification'].lower()]
        if vus_results:
            most_conflicting = max(vus_results, key=lambda x: x['vampp_score'])
            print(f"{Fore.WHITE}  Most Conflicting VUS: {most_conflicting['scenario_name']} ({most_conflicting['vampp_score']:.3f})")
        
        print(f"\n{Fore.GREEN}{Style.BRIGHT}✅ All {len(results)} scenarios completed successfully!")
        print(f"{Fore.CYAN}{'='*80}")


def main():
    """Main execution function."""
    test_scenarios = VUSTestScenarios()
    results = test_scenarios.run_all_scenarios()
    
    # Keep terminal open
    print(f"\n{Fore.CYAN}Press Enter to exit...")
    input()

if __name__ == "__main__":
    main()
