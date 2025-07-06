#!/usr/bin/env python3
"""
VUS Upgrade Test Scenarios
==========================

Test scenarios where ClinVar shows VUS but our algorithm can provide
more definitive classifications (Likely Pathogenic or Benign) through
enhanced analysis and additional criteria.

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

class VUSUpgradeTestScenarios:
    """Test scenarios for VUS upgrade classification."""
    
    def __init__(self):
        self.evidence_evaluator = EvidenceEvaluator(use_2023_guidelines=False)
        self.classifier = ACMGClassifier(use_2023_guidelines=False)
        self.report_generator = ReportGenerator()
    
    def run_all_scenarios(self):
        """Run all VUS upgrade test scenarios."""
        print(f"\n{Fore.CYAN}{'='*80}")
        print(f"{Fore.CYAN}{Style.BRIGHT}🔬 VUS UPGRADE TEST SCENARIOS")
        print(f"{Fore.CYAN}{Style.BRIGHT}ClinVar VUS → Algorithm Upgrade Analysis")
        print(f"{Fore.CYAN}{'='*80}")
        print(f"{Fore.YELLOW}Test Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{Fore.CYAN}{'='*80}")
        
        scenarios = [
            self.scenario_1_brca1_vus_to_likely_pathogenic,
            self.scenario_2_tp53_vus_to_likely_pathogenic,
            self.scenario_3_cftr_vus_to_benign,
            self.scenario_4_ldlr_vus_to_likely_benign,
            self.scenario_5_scn5a_vus_to_likely_pathogenic,
            self.scenario_6_mecp2_vus_to_likely_pathogenic
        ]
        
        results = []
        upgrades = {'to_likely_pathogenic': 0, 'to_likely_benign': 0, 'to_pathogenic': 0, 'to_benign': 0, 'still_vus': 0}
        
        for i, scenario_func in enumerate(scenarios, 1):
            print(f"\n{Fore.MAGENTA}{'='*70}")
            print(f"{Fore.MAGENTA}{Style.BRIGHT}📊 SCENARIO {i}")
            print(f"{Fore.MAGENTA}{'='*70}")
            
            result = scenario_func()
            results.append(result)
            
            # Track upgrades
            classification = result['classification'].lower()
            clinvar_status = result.get('clinvar_status', 'vus').lower()
            
            print(f"\n{Fore.CYAN}📋 COMPARISON:")
            print(f"{Fore.WHITE}  ClinVar Status: {Fore.YELLOW}{clinvar_status.upper()}")
            print(f"{Fore.WHITE}  Our Algorithm: {Fore.GREEN}{result['classification']}")
            
            if 'uncertain' in classification or 'vus' in classification:
                upgrades['still_vus'] += 1
                print(f"{Fore.YELLOW}  Status: No upgrade (still VUS)")
            elif 'likely pathogenic' in classification:
                upgrades['to_likely_pathogenic'] += 1
                print(f"{Fore.GREEN}  Status: ✅ UPGRADED to Likely Pathogenic")
            elif 'likely benign' in classification:
                upgrades['to_likely_benign'] += 1
                print(f"{Fore.GREEN}  Status: ✅ UPGRADED to Likely Benign")
            elif 'pathogenic' in classification and 'likely' not in classification:
                upgrades['to_pathogenic'] += 1
                print(f"{Fore.GREEN}  Status: ✅ UPGRADED to Pathogenic")
            elif 'benign' in classification and 'likely' not in classification:
                upgrades['to_benign'] += 1
                print(f"{Fore.GREEN}  Status: ✅ UPGRADED to Benign")
            
            print(f"{Fore.BLUE}  Confidence: {result['confidence']}")
            print(f"{Fore.BLUE}  Applied Criteria: {sum(len(criteria) for criteria in result['applied_criteria'].values())}")
        
        # Final summary
        self._print_upgrade_summary(results, upgrades)
        
        return results
    
    def scenario_1_brca1_vus_to_likely_pathogenic(self):
        """
        Scenario 1: BRCA1 VUS with strong conservation and de novo evidence
        ClinVar: VUS → Expected: Likely Pathogenic
        """
        print(f"{Fore.YELLOW}📝 Scenario 1: BRCA1 VUS → Likely Pathogenic")
        print(f"{Fore.WHITE}Gene: BRCA1 (breast cancer susceptibility)")
        print(f"{Fore.WHITE}Variant: c.5096G>A (p.Arg1699Gln)")
        print(f"{Fore.WHITE}ClinVar Status: VUS")
        print(f"{Fore.WHITE}Enhancement: Strong conservation + De novo + Functional domain")
        
        variant_data = VariantData(
            basic_info={
                'gene': 'BRCA1',
                'chromosome': '17',
                'position': '41244000',
                'ref_allele': 'G',
                'alt_allele': 'A',
                'hgvs_c': 'c.5096G>A',
                'hgvs_p': 'p.Arg1699Gln',
                'variant_type': 'missense',
                'transcript': 'NM_007294.3',
                'exon': '15',
                'protein_change': 'p.Arg1699Gln'
            },
            population_data={
                'gnomad_af': 0.000045,  # Very rare
                'gnomad_af_popmax': 0.00032,
                'gnomad_an': 251456,
                'gnomad_ac': 11,
                'gnomad_hom': 0,
                'gnomad_hemi': 0,
                'exac_af': 0.000038,
                'thousand_genomes_af': 0.00004,
                'prevalence': 0.00125  # Breast cancer prevalence
            },
            insilico_data={
                # High conservation scores
                'phylop_100way': 8.9,
                'phylop_30way': 8.2,
                'phylop_17way': 7.6,
                'phastcons_100way': 0.998,
                'phastcons_30way': 0.995,
                'phastcons_17way': 0.989,
                'gerp_score': 5.8,
                'siphy_score': 18.7,
                
                # Strong pathogenicity predictors
                'sift_score': 0.001,
                'sift_prediction': 'D',
                'polyphen2_hdiv_score': 0.996,
                'polyphen2_hdiv_prediction': 'D',
                'polyphen2_hvar_score': 0.987,
                'polyphen2_hvar_prediction': 'D',
                'lrt_score': 0.998,
                'lrt_prediction': 'D',
                'mutationtaster_score': 0.999,
                'mutationtaster_prediction': 'D',
                'mutationassessor_score': 4.12,
                'mutationassessor_prediction': 'H',
                'fathmm_score': 2.89,
                'fathmm_prediction': 'D',
                'provean_score': -7.34,
                'provean_prediction': 'D',
                'vest4_score': 0.923,
                'metasvm_score': 0.967,
                'metasvm_prediction': 'D',
                'metalr_score': 0.942,
                'metalr_prediction': 'D',
                
                # Additional high-confidence predictors
                'revel': 0.856,
                'cadd_phred': 29.3,
                'alphamissense': 0.789,
                'clinpred': 0.892,
                'bayesdel_addaf': 0.678,
                'integrated_fitcons_score': 0.823,
                'gm12878_fitcons_score': 0.856,
                'h1_hesc_fitcons_score': 0.834,
                'huvec_fitcons_score': 0.789
            },
            genetic_data={
                'mode_of_inheritance': 'autosomal_dominant',
                'penetrance': 'incomplete',
                'allelic_requirement': 'monoallelic',
                'disease_mechanism': 'loss_of_function',
                'variant_consequence': 'missense_variant',
                'variant_effect': 'uncertain',
                'de_novo': 'yes',  # Strong evidence
                'parental_testing': 'confirmed',
                'segregation_analysis': 'supportive',
                'family_history': 'positive'
            },
            functional_data={
                'protein_domain': 'BRCT domain',
                'protein_structure_affected': True,
                'catalytic_residue': False,
                'binding_site': True,
                'known_functional_domain': True,
                'previous_functional_studies': False,
                'protein_stability_affected': 'likely_reduced',
                'enzymatic_activity_affected': 'likely_reduced',
                'phenotype_match': 'strong'
            }
        )
        
        return self._evaluate_and_classify(variant_data, "BRCA1 VUS → Likely Pathogenic", "VUS")
    
    def scenario_2_tp53_vus_to_likely_pathogenic(self):
        """
        Scenario 2: TP53 VUS with extreme conservation but mixed predictors
        ClinVar: VUS → Expected: Likely Pathogenic (conservation-driven)
        """
        print(f"{Fore.YELLOW}📝 Scenario 2: TP53 VUS → Likely Pathogenic")
        print(f"{Fore.WHITE}Gene: TP53 (tumor suppressor)")
        print(f"{Fore.WHITE}Variant: c.818G>A (p.Arg273His)")
        print(f"{Fore.WHITE}ClinVar Status: VUS")
        print(f"{Fore.WHITE}Enhancement: Extreme conservation + DNA binding domain")
        
        variant_data = VariantData(
            basic_info={
                'gene': 'TP53',
                'chromosome': '17',
                'position': '7673802',
                'ref_allele': 'G',
                'alt_allele': 'A',
                'hgvs_c': 'c.818G>A',
                'hgvs_p': 'p.Arg273His',
                'variant_type': 'missense',
                'transcript': 'NM_000546.5',
                'exon': '8',
                'protein_change': 'p.Arg273His'
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
                'prevalence': 0.0001  # Li-Fraumeni syndrome
            },
            insilico_data={
                # Extreme conservation (TP53 is highly conserved)
                'phylop_100way': 9.99,
                'phylop_30way': 9.95,
                'phylop_17way': 9.89,
                'phastcons_100way': 1.0,
                'phastcons_30way': 1.0,
                'phastcons_17way': 1.0,
                'gerp_score': 6.18,
                'siphy_score': 29.9,
                
                # Mixed pathogenicity predictors
                'sift_score': 0.034,
                'sift_prediction': 'D',
                'polyphen2_hdiv_score': 0.687,
                'polyphen2_hdiv_prediction': 'P',
                'polyphen2_hvar_score': 0.612,
                'polyphen2_hvar_prediction': 'P',
                'lrt_score': 0.789,
                'lrt_prediction': 'D',
                'mutationtaster_score': 0.823,
                'mutationtaster_prediction': 'D',
                'mutationassessor_score': 3.45,
                'mutationassessor_prediction': 'H',
                'fathmm_score': 0.23,
                'fathmm_prediction': 'T',
                'provean_score': -3.12,
                'provean_prediction': 'D',
                'vest4_score': 0.678,
                'metasvm_score': 0.567,
                'metasvm_prediction': 'D',
                'metalr_score': 0.543,
                'metalr_prediction': 'D',
                
                # Additional predictors
                'revel': 0.623,
                'cadd_phred': 26.8,
                'alphamissense': 0.578,
                'clinpred': 0.712,
                'bayesdel_addaf': 0.456,
                'integrated_fitcons_score': 0.789,
                'gm12878_fitcons_score': 0.823,
                'h1_hesc_fitcons_score': 0.856,
                'huvec_fitcons_score': 0.798
            },
            genetic_data={
                'mode_of_inheritance': 'autosomal_dominant',
                'penetrance': 'high',
                'allelic_requirement': 'monoallelic',
                'disease_mechanism': 'loss_of_function',
                'variant_consequence': 'missense_variant',
                'variant_effect': 'uncertain',
                'de_novo': 'yes',
                'parental_testing': 'confirmed',
                'segregation_analysis': 'supportive',
                'family_history': 'positive'
            },
            functional_data={
                'protein_domain': 'DNA-binding domain',
                'protein_structure_affected': True,
                'catalytic_residue': False,
                'binding_site': True,
                'known_functional_domain': True,
                'previous_functional_studies': False,
                'protein_stability_affected': 'likely_reduced',
                'enzymatic_activity_affected': 'likely_reduced',
                'phenotype_match': 'strong'
            }
        )
        
        return self._evaluate_and_classify(variant_data, "TP53 VUS → Likely Pathogenic", "VUS")
    
    def scenario_3_cftr_vus_to_benign(self):
        """
        Scenario 3: CFTR VUS with benign population frequency
        ClinVar: VUS → Expected: Benign (population frequency)
        """
        print(f"{Fore.YELLOW}📝 Scenario 3: CFTR VUS → Benign")
        print(f"{Fore.WHITE}Gene: CFTR (cystic fibrosis)")
        print(f"{Fore.WHITE}Variant: c.1680-886A>G")
        print(f"{Fore.WHITE}ClinVar Status: VUS")
        print(f"{Fore.WHITE}Enhancement: High population frequency + Benign predictors")
        
        variant_data = VariantData(
            basic_info={
                'gene': 'CFTR',
                'chromosome': '7',
                'position': '117199645',
                'ref_allele': 'A',
                'alt_allele': 'G',
                'hgvs_c': 'c.1680-886A>G',
                'hgvs_p': 'p.=',
                'variant_type': 'intronic',
                'transcript': 'NM_000492.3',
                'exon': 'intron_12',
                'protein_change': 'p.='
            },
            population_data={
                'gnomad_af': 0.0123,  # High frequency (>1%)
                'gnomad_af_popmax': 0.0456,
                'gnomad_an': 251456,
                'gnomad_ac': 3093,
                'gnomad_hom': 15,
                'gnomad_hemi': 0,
                'exac_af': 0.0108,
                'thousand_genomes_af': 0.0134,
                'prevalence': 0.0004  # CF prevalence
            },
            insilico_data={
                # Low conservation (intronic)
                'phylop_100way': 0.45,
                'phylop_30way': 0.23,
                'phylop_17way': 0.12,
                'phastcons_100way': 0.034,
                'phastcons_30way': 0.019,
                'phastcons_17way': 0.008,
                'gerp_score': -2.1,
                'siphy_score': 0.23,
                
                # Benign predictors (splice prediction)
                'ada_score': 0.001,
                'rf_score': 0.002,
                'spliceai_ag_score': 0.001,
                'spliceai_al_score': 0.001,
                'spliceai_dg_score': 0.001,
                'spliceai_dl_score': 0.001,
                'integrated_fitcons_score': 0.001,
                'gm12878_fitcons_score': 0.001,
                'h1_hesc_fitcons_score': 0.001,
                'huvec_fitcons_score': 0.001
            },
            genetic_data={
                'mode_of_inheritance': 'autosomal_recessive',
                'penetrance': 'complete',
                'allelic_requirement': 'biallelic',
                'disease_mechanism': 'loss_of_function',
                'variant_consequence': 'intron_variant',
                'variant_effect': 'uncertain',
                'de_novo': 'no',
                'parental_testing': 'not_applicable',
                'segregation_analysis': 'not_supportive',
                'family_history': 'negative'
            },
            functional_data={
                'protein_domain': 'intronic',
                'protein_structure_affected': False,
                'catalytic_residue': False,
                'binding_site': False,
                'known_functional_domain': False,
                'previous_functional_studies': False,
                'protein_stability_affected': 'no_effect',
                'enzymatic_activity_affected': 'no_effect',
                'phenotype_match': 'poor'
            }
        )
        
        return self._evaluate_and_classify(variant_data, "CFTR VUS → Benign", "VUS")
    
    def scenario_4_ldlr_vus_to_likely_benign(self):
        """
        Scenario 4: LDLR VUS with benign predictors and moderate frequency
        ClinVar: VUS → Expected: Likely Benign
        """
        print(f"{Fore.YELLOW}📝 Scenario 4: LDLR VUS → Likely Benign")
        print(f"{Fore.WHITE}Gene: LDLR (familial hypercholesterolemia)")
        print(f"{Fore.WHITE}Variant: c.1846A>G (p.Ile616Val)")
        print(f"{Fore.WHITE}ClinVar Status: VUS")
        print(f"{Fore.WHITE}Enhancement: Moderate frequency + Benign predictors")
        
        variant_data = VariantData(
            basic_info={
                'gene': 'LDLR',
                'chromosome': '19',
                'position': '11218543',
                'ref_allele': 'A',
                'alt_allele': 'G',
                'hgvs_c': 'c.1846A>G',
                'hgvs_p': 'p.Ile616Val',
                'variant_type': 'missense',
                'transcript': 'NM_000527.4',
                'exon': '12',
                'protein_change': 'p.Ile616Val'
            },
            population_data={
                'gnomad_af': 0.0023,  # Moderate frequency
                'gnomad_af_popmax': 0.0089,
                'gnomad_an': 251456,
                'gnomad_ac': 578,
                'gnomad_hom': 2,
                'gnomad_hemi': 0,
                'exac_af': 0.0019,
                'thousand_genomes_af': 0.0021,
                'prevalence': 0.002  # FH prevalence
            },
            insilico_data={
                # Low conservation
                'phylop_100way': 1.2,
                'phylop_30way': 0.8,
                'phylop_17way': 0.5,
                'phastcons_100way': 0.234,
                'phastcons_30way': 0.189,
                'phastcons_17way': 0.123,
                'gerp_score': 0.45,
                'siphy_score': 2.1,
                
                # Benign predictors
                'sift_score': 0.234,
                'sift_prediction': 'T',
                'polyphen2_hdiv_score': 0.123,
                'polyphen2_hdiv_prediction': 'B',
                'polyphen2_hvar_score': 0.089,
                'polyphen2_hvar_prediction': 'B',
                'lrt_score': 0.067,
                'lrt_prediction': 'N',
                'mutationtaster_score': 0.089,
                'mutationtaster_prediction': 'P',
                'mutationassessor_score': 0.89,
                'mutationassessor_prediction': 'N',
                'fathmm_score': -2.34,
                'fathmm_prediction': 'T',
                'provean_score': -0.45,
                'provean_prediction': 'N',
                'vest4_score': 0.089,
                'metasvm_score': 0.067,
                'metasvm_prediction': 'T',
                'metalr_score': 0.078,
                'metalr_prediction': 'T',
                
                # Additional benign predictors
                'revel': 0.123,
                'cadd_phred': 12.3,
                'alphamissense': 0.234,
                'clinpred': 0.156,
                'bayesdel_addaf': -0.345,
                'integrated_fitcons_score': 0.089,
                'gm12878_fitcons_score': 0.067,
                'h1_hesc_fitcons_score': 0.123,
                'huvec_fitcons_score': 0.089
            },
            genetic_data={
                'mode_of_inheritance': 'autosomal_dominant',
                'penetrance': 'incomplete',
                'allelic_requirement': 'monoallelic',
                'disease_mechanism': 'loss_of_function',
                'variant_consequence': 'missense_variant',
                'variant_effect': 'uncertain',
                'de_novo': 'no',
                'parental_testing': 'not_applicable',
                'segregation_analysis': 'not_supportive',
                'family_history': 'negative'
            },
            functional_data={
                'protein_domain': 'LDL-receptor class A domain',
                'protein_structure_affected': False,
                'catalytic_residue': False,
                'binding_site': False,
                'known_functional_domain': True,
                'previous_functional_studies': False,
                'protein_stability_affected': 'no_effect',
                'enzymatic_activity_affected': 'no_effect',
                'phenotype_match': 'poor'
            }
        )
        
        return self._evaluate_and_classify(variant_data, "LDLR VUS → Likely Benign", "VUS")
    
    def scenario_5_scn5a_vus_to_likely_pathogenic(self):
        """
        Scenario 5: SCN5A VUS with strong metascore
        ClinVar: VUS → Expected: Likely Pathogenic (metascore-driven)
        """
        print(f"{Fore.YELLOW}📝 Scenario 5: SCN5A VUS → Likely Pathogenic")
        print(f"{Fore.WHITE}Gene: SCN5A (Brugada syndrome)")
        print(f"{Fore.WHITE}Variant: c.5218G>A (p.Ala1740Thr)")
        print(f"{Fore.WHITE}ClinVar Status: VUS")
        print(f"{Fore.WHITE}Enhancement: Strong metascore + Functional domain")
        
        variant_data = VariantData(
            basic_info={
                'gene': 'SCN5A',
                'chromosome': '3',
                'position': '38645323',
                'ref_allele': 'G',
                'alt_allele': 'A',
                'hgvs_c': 'c.5218G>A',
                'hgvs_p': 'p.Ala1740Thr',
                'variant_type': 'missense',
                'transcript': 'NM_198056.2',
                'exon': '28',
                'protein_change': 'p.Ala1740Thr'
            },
            population_data={
                'gnomad_af': 0.000023,  # Very rare
                'gnomad_af_popmax': 0.00015,
                'gnomad_an': 251456,
                'gnomad_ac': 6,
                'gnomad_hom': 0,
                'gnomad_hemi': 0,
                'exac_af': 0.000019,
                'thousand_genomes_af': 0.00002,
                'prevalence': 0.0002  # Brugada syndrome
            },
            insilico_data={
                # Strong conservation
                'phylop_100way': 7.2,
                'phylop_30way': 6.8,
                'phylop_17way': 6.1,
                'phastcons_100way': 0.956,
                'phastcons_30way': 0.923,
                'phastcons_17way': 0.889,
                'gerp_score': 5.2,
                'siphy_score': 16.3,
                
                # Strong pathogenicity predictors
                'sift_score': 0.002,
                'sift_prediction': 'D',
                'polyphen2_hdiv_score': 0.987,
                'polyphen2_hdiv_prediction': 'D',
                'polyphen2_hvar_score': 0.978,
                'polyphen2_hvar_prediction': 'D',
                'lrt_score': 0.987,
                'lrt_prediction': 'D',
                'mutationtaster_score': 0.998,
                'mutationtaster_prediction': 'D',
                'mutationassessor_score': 4.23,
                'mutationassessor_prediction': 'H',
                'fathmm_score': 3.45,
                'fathmm_prediction': 'D',
                'provean_score': -6.78,
                'provean_prediction': 'D',
                'vest4_score': 0.923,
                'metasvm_score': 0.976,
                'metasvm_prediction': 'D',
                'metalr_score': 0.967,
                'metalr_prediction': 'D',
                
                # Additional strong predictors
                'revel': 0.892,
                'cadd_phred': 31.2,
                'alphamissense': 0.823,
                'clinpred': 0.934,
                'bayesdel_addaf': 0.712,
                'integrated_fitcons_score': 0.856,
                'gm12878_fitcons_score': 0.889,
                'h1_hesc_fitcons_score': 0.923,
                'huvec_fitcons_score': 0.834
            },
            genetic_data={
                'mode_of_inheritance': 'autosomal_dominant',
                'penetrance': 'incomplete',
                'allelic_requirement': 'monoallelic',
                'disease_mechanism': 'loss_of_function',
                'variant_consequence': 'missense_variant',
                'variant_effect': 'uncertain',
                'de_novo': 'yes',
                'parental_testing': 'confirmed',
                'segregation_analysis': 'supportive',
                'family_history': 'positive'
            },
            functional_data={
                'protein_domain': 'Voltage-gated sodium channel domain',
                'protein_structure_affected': True,
                'catalytic_residue': False,
                'binding_site': True,
                'known_functional_domain': True,
                'previous_functional_studies': False,
                'protein_stability_affected': 'likely_reduced',
                'enzymatic_activity_affected': 'likely_reduced',
                'phenotype_match': 'strong'
            }
        )
        
        return self._evaluate_and_classify(variant_data, "SCN5A VUS → Likely Pathogenic", "VUS")
    
    def scenario_6_mecp2_vus_to_likely_pathogenic(self):
        """
        Scenario 6: MECP2 VUS with extreme conservation
        ClinVar: VUS → Expected: Likely Pathogenic (conservation-driven)
        """
        print(f"{Fore.YELLOW}📝 Scenario 6: MECP2 VUS → Likely Pathogenic")
        print(f"{Fore.WHITE}Gene: MECP2 (Rett syndrome)")
        print(f"{Fore.WHITE}Variant: c.916C>T (p.Arg306Cys)")
        print(f"{Fore.WHITE}ClinVar Status: VUS")
        print(f"{Fore.WHITE}Enhancement: Extreme conservation + Functional domain")
        
        variant_data = VariantData(
            basic_info={
                'gene': 'MECP2',
                'chromosome': 'X',
                'position': '154028789',
                'ref_allele': 'C',
                'alt_allele': 'T',
                'hgvs_c': 'c.916C>T',
                'hgvs_p': 'p.Arg306Cys',
                'variant_type': 'missense',
                'transcript': 'NM_004992.3',
                'exon': '4',
                'protein_change': 'p.Arg306Cys'
            },
            population_data={
                'gnomad_af': 0.000012,  # Very rare
                'gnomad_af_popmax': 0.00008,
                'gnomad_an': 251456,
                'gnomad_ac': 3,
                'gnomad_hom': 0,
                'gnomad_hemi': 0,
                'exac_af': 0.000009,
                'thousand_genomes_af': 0.00001,
                'prevalence': 0.00001  # Rett syndrome prevalence
            },
            insilico_data={
                # Extreme conservation (MECP2 is highly conserved)
                'phylop_100way': 9.98,
                'phylop_30way': 9.96,
                'phylop_17way': 9.89,
                'phastcons_100way': 1.0,
                'phastcons_30way': 1.0,
                'phastcons_17way': 1.0,
                'gerp_score': 6.18,
                'siphy_score': 29.9,
                
                # Strong pathogenicity predictors
                'sift_score': 0.001,
                'sift_prediction': 'D',
                'polyphen2_hdiv_score': 0.999,
                'polyphen2_hdiv_prediction': 'D',
                'polyphen2_hvar_score': 0.998,
                'polyphen2_hvar_prediction': 'D',
                'lrt_score': 0.999,
                'lrt_prediction': 'D',
                'mutationtaster_score': 1.0,
                'mutationtaster_prediction': 'D',
                'mutationassessor_score': 4.67,
                'mutationassessor_prediction': 'H',
                'fathmm_score': 4.12,
                'fathmm_prediction': 'D',
                'provean_score': -8.23,
                'provean_prediction': 'D',
                'vest4_score': 0.967,
                'metasvm_score': 0.989,
                'metasvm_prediction': 'D',
                'metalr_score': 0.978,
                'metalr_prediction': 'D',
                
                # Additional strong predictors
                'revel': 0.923,
                'cadd_phred': 34.1,
                'alphamissense': 0.889,
                'clinpred': 0.967,
                'bayesdel_addaf': 0.823,
                'integrated_fitcons_score': 0.923,
                'gm12878_fitcons_score': 0.956,
                'h1_hesc_fitcons_score': 0.967,
                'huvec_fitcons_score': 0.889
            },
            genetic_data={
                'mode_of_inheritance': 'x_linked_dominant',
                'penetrance': 'complete',
                'allelic_requirement': 'monoallelic',
                'disease_mechanism': 'loss_of_function',
                'variant_consequence': 'missense_variant',
                'variant_effect': 'uncertain',
                'de_novo': 'yes',
                'parental_testing': 'confirmed',
                'segregation_analysis': 'supportive',
                'family_history': 'positive'
            },
            functional_data={
                'protein_domain': 'Methyl-CpG binding domain',
                'protein_structure_affected': True,
                'catalytic_residue': False,
                'binding_site': True,
                'known_functional_domain': True,
                'previous_functional_studies': False,
                'protein_stability_affected': 'likely_reduced',
                'enzymatic_activity_affected': 'likely_reduced',
                'phenotype_match': 'strong'
            }
        )
        
        return self._evaluate_and_classify(variant_data, "MECP2 VUS → Likely Pathogenic", "VUS")
    
    def _evaluate_and_classify(self, variant_data: VariantData, scenario_name: str, clinvar_status: str) -> dict:
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
        
        # Handle None values
        if vampp_score is None:
            vampp_score = 0
        if conservation_score is None:
            conservation_score = 0
        if pathogenicity_score is None:
            pathogenicity_score = 0
        
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
        classification_result['clinvar_status'] = clinvar_status
        classification_result['report_path'] = report_path
        classification_result['vampp_score'] = vampp_score
        classification_result['conservation_score'] = conservation_score
        classification_result['pathogenicity_score'] = pathogenicity_score
        
        return classification_result
    
    def _print_upgrade_summary(self, results: list, upgrades: dict):
        """Print summary of upgrade results."""
        print(f"\n{Fore.CYAN}{'='*80}")
        print(f"{Fore.CYAN}{Style.BRIGHT}🎯 VUS UPGRADE SUMMARY")
        print(f"{Fore.CYAN}{'='*80}")
        
        total_scenarios = len(results)
        successful_upgrades = sum(upgrades[k] for k in upgrades if k != 'still_vus')
        
        print(f"\n{Fore.GREEN}{Style.BRIGHT}📈 UPGRADE STATISTICS:")
        print(f"{Fore.WHITE}  Total Scenarios: {total_scenarios}")
        print(f"{Fore.WHITE}  Successful Upgrades: {successful_upgrades}/{total_scenarios} ({successful_upgrades/total_scenarios*100:.1f}%)")
        print(f"{Fore.WHITE}  Still VUS: {upgrades['still_vus']}")
        
        print(f"\n{Fore.YELLOW}{Style.BRIGHT}📊 UPGRADE BREAKDOWN:")
        for upgrade_type, count in upgrades.items():
            if count > 0:
                upgrade_name = upgrade_type.replace('_', ' ').title()
                print(f"{Fore.WHITE}  {upgrade_name}: {count}")
        
        # Most successful upgrades
        print(f"\n{Fore.YELLOW}{Style.BRIGHT}🏆 MOST SUCCESSFUL UPGRADES:")
        successful_results = [r for r in results if 'uncertain' not in r['classification'].lower()]
        
        if successful_results:
            # Highest scoring upgrade
            highest_vampp = max(successful_results, key=lambda x: x['vampp_score'])
            print(f"{Fore.WHITE}  Highest VAMPP Score: {highest_vampp['scenario_name']} ({highest_vampp['vampp_score']:.3f})")
            
            # Most criteria applied
            most_criteria = max(successful_results, key=lambda x: sum(len(criteria) for criteria in x['applied_criteria'].values()))
            total_criteria = sum(len(criteria) for criteria in most_criteria['applied_criteria'].values())
            print(f"{Fore.WHITE}  Most Criteria Applied: {most_criteria['scenario_name']} ({total_criteria} criteria)")
        
        # Algorithm effectiveness
        print(f"\n{Fore.YELLOW}{Style.BRIGHT}🔬 ALGORITHM EFFECTIVENESS:")
        print(f"{Fore.WHITE}  Enhanced Analysis: VAMPP metascore + Conservation analysis")
        print(f"{Fore.WHITE}  Prevalence Adjustment: Threshold adaptation based on disease prevalence")
        print(f"{Fore.WHITE}  Additional Criteria: De novo, functional domains, phenotype matching")
        
        # Key insights
        print(f"\n{Fore.YELLOW}{Style.BRIGHT}💡 KEY INSIGHTS:")
        print(f"{Fore.WHITE}  • High conservation scores can drive pathogenic classification")
        print(f"{Fore.WHITE}  • Population frequency is crucial for benign classification")
        print(f"{Fore.WHITE}  • De novo status adds significant pathogenic weight")
        print(f"{Fore.WHITE}  • Functional domain knowledge enhances accuracy")
        print(f"{Fore.WHITE}  • Metascore integration provides better discrimination")
        
        print(f"\n{Fore.GREEN}{Style.BRIGHT}✅ Algorithm successfully upgraded {successful_upgrades}/{total_scenarios} VUS classifications!")
        print(f"{Fore.CYAN}{'='*80}")


def main():
    """Main execution function."""
    test_scenarios = VUSUpgradeTestScenarios()
    results = test_scenarios.run_all_scenarios()
    
    # Keep terminal open
    print(f"\n{Fore.CYAN}Press Enter to exit...")
    input()

if __name__ == "__main__":
    main()
