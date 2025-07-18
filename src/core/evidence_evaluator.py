class VariantData:
    def __init__(self, basic_info=None, population_data=None, insilico_data=None, genetic_data=None, functional_data=None, patient_phenotypes=None, clinvar_data=None):
        self.basic_info = basic_info or {}
        self.population_data = population_data or {}
        self.insilico_data = insilico_data or {}
        self.genetic_data = genetic_data or {}
        self.functional_data = functional_data or {}
        self.patient_phenotypes = patient_phenotypes
        self.clinvar_data = clinvar_data

    @property
    def gene(self):
        return self.basic_info.get('gene', None)

    @property
    def hgvs_c(self):
        return self.basic_info.get('hgvs_c', None)
class InframeAnalyzer:
    def __init__(self):
        import logging
        self.logger = logging.getLogger("InframeAnalyzer")
        self.critical_regions = self._load_critical_regions()
        self.repeat_regions = self._load_repeat_regions()
        self.domain_boundaries = self._load_domain_boundaries()

    def evaluate_inframe_deletion(self, variant_data):
        """Enhanced inframe deletion evaluation"""
        if self._affects_critical_region(variant_data):
            return 'PM4'
        if self._affects_functional_domain(variant_data):
            return 'PM4'
        if self._affects_structural_integrity(variant_data):
            return 'PM4'
        if self._in_repeat_region(variant_data):
            return 'BP3'
        return None

    def _affects_critical_region(self, variant_data):
        gene = variant_data.basic_info.get('gene', None)
        hgvs_c = variant_data.basic_info.get('hgvs_c', None)
        position = self._extract_position(hgvs_c)
        if position is None:
            return False
        if gene in self.critical_regions:
            for region in self.critical_regions[gene]:
                if region['start'] <= position <= region['end']:
                    return True
        return False

    def _affects_functional_domain(self, variant_data):
        return False

    def _affects_structural_integrity(self, variant_data):
        return False

    def _in_repeat_region(self, variant_data):
        return False

    def _extract_position(self, hgvs_c):
        import re
        if not hgvs_c:
            if hasattr(self, 'logger'):
                self.logger.warning("hgvs_c is None or empty")
            return None
        match = re.search(r'c\.(\d+)', hgvs_c)
        if match:
            return int(match.group(1))
        if hasattr(self, 'logger'):
            self.logger.warning(f"No position found in hgvs_c: {hgvs_c}")
        return None

    def _load_critical_regions(self):
        return {'TP53': [{'start': 730, 'end': 750}]}

    def _load_repeat_regions(self):
        return {}

    def _load_domain_boundaries(self):
        return {}
"""
Evidence Evaluator Module
=========================

This module evaluates ACMG evidence criteria based on variant data.
It implements both ACMG/AMP 2015 and 2023 guidelines.
"""

import numpy as np
import math
from typing import Dict, List, Optional, Any, Tuple
from scipy import stats
from config.constants import (
    EVIDENCE_WEIGHTS, GENE_SPECIFIC_THRESHOLDS, INSILICO_WEIGHTS,
    INSILICO_THRESHOLDS, VAMPP_SCORE_THRESHOLDS, STATISTICAL_THRESHOLDS,
    LOF_INTOLERANT_GENES, LOF_TOLERANT_GENES
)


class EvidenceEvaluator:
    """
    Evaluates ACMG evidence criteria for variant classification.
    
    This class implements the logic for evaluating all ACMG evidence criteria
    based on variant data, including population frequencies, in silico predictions,
    functional studies, and inheritance patterns.
    """
    
    def __init__(self, use_2023_guidelines: bool = False, test_mode: bool = False):
        """
        Initialize the evidence evaluator.
        
        Args:
            use_2023_guidelines (bool): Whether to use ACMG 2023 guidelines
            test_mode (bool): Whether to run in test mode (skip interactive prompts)
        """
        self.use_2023_guidelines = use_2023_guidelines
        self.test_mode = test_mode
        self.applied_criteria = {}
        self.evidence_details = {}
        # Yeni modüller entegre ediliyor
        from core.functional_studies_evaluator import FunctionalStudiesEvaluator
        from core.phenotype_matcher import PhenotypeMatcher
        self.functional_studies_evaluator = FunctionalStudiesEvaluator()
        self.phenotype_matcher = PhenotypeMatcher()
    
    def evaluate_all_criteria(self, variant_data) -> Dict[str, Any]:
        """
        Evaluate all ACMG criteria for a variant.
        
        Args:
            variant_data: VariantData object containing all variant information
            
        Returns:
            Dict[str, Any]: Evaluation results
        """
        results = {
            'pathogenic_criteria': {},
            'benign_criteria': {},
            'applied_criteria': {},
            'evidence_details': {},
            'vampp_score': None,
            'statistical_tests': {}
        }
        
        # Evaluate pathogenic criteria
        results['pathogenic_criteria'] = self._evaluate_pathogenic_criteria(variant_data)
        
        # Evaluate benign criteria
        results['benign_criteria'] = self._evaluate_benign_criteria(variant_data)
        
        # Calculate VAMPP-score-like metascore
        if variant_data.basic_info.get('variant_type') == 'missense':
            results['vampp_score'] = self._calculate_vampp_score(variant_data)
        
        # Perform statistical tests
        results['statistical_tests'] = self._perform_statistical_tests(variant_data)
        
        # Consolidate applied criteria
        results['applied_criteria'] = self._consolidate_applied_criteria(
            results['pathogenic_criteria'], 
            results['benign_criteria']
        )
        
        # Store evidence details
        results['evidence_details'] = self.evidence_details
        
        return results
    
    def _evaluate_pathogenic_criteria(self, variant_data) -> Dict[str, Any]:
        """Evaluate pathogenic evidence criteria."""
        pathogenic = {}
        
        # PVS1 - Null variant in LOF gene
        pathogenic['PVS1'] = self._evaluate_pvs1(variant_data)
        
        # PS1 - Same amino acid change as pathogenic variant
        pathogenic['PS1'] = self._evaluate_ps1(variant_data)
        
        # PS2 - De novo variant
        pathogenic['PS2'] = self._evaluate_ps2(variant_data)
        
        # PS3 - Functional studies
        pathogenic['PS3'] = self._evaluate_ps3(variant_data)
        
        # PS4 - Case-control data
        pathogenic['PS4'] = self._evaluate_ps4(variant_data)
        
        # PM1 - Mutational hotspot
        pathogenic['PM1'] = self._evaluate_pm1(variant_data)
        
        # PM2 - Absent from controls
        pathogenic['PM2'] = self._evaluate_pm2(variant_data)
        
        # PM3 - Recessive disorder, in trans
        pathogenic['PM3'] = self._evaluate_pm3(variant_data)
        
        # PM4 - Protein length changes
        pathogenic['PM4'] = self._evaluate_pm4(variant_data)
        
        # PM5 - Novel missense at pathogenic residue
        pathogenic['PM5'] = self._evaluate_pm5(variant_data)
        
        # PM6 - Assumed de novo
        pathogenic['PM6'] = self._evaluate_pm6(variant_data)
        
        # PP1 - Segregation
        pathogenic['PP1'] = self._evaluate_pp1(variant_data)
        
        # PP2 - Missense in gene with low benign rate
        pathogenic['PP2'] = self._evaluate_pp2(variant_data)
        
        # PP3 - In silico evidence
        pathogenic['PP3'] = self._evaluate_pp3(variant_data)
        
        # PP4 - Phenotype match
        pathogenic['PP4'] = self._evaluate_pp4(variant_data)
        
        # PP5 - Reputable source (if enabled)
        if self.use_2023_guidelines:
            pathogenic['PP5'] = self._evaluate_pp5(variant_data)
        
        return pathogenic
    
    def _evaluate_benign_criteria(self, variant_data) -> Dict[str, Any]:
        """Evaluate benign evidence criteria."""
        benign = {}
        
        # BA1 - High frequency
        benign['BA1'] = self._evaluate_ba1(variant_data)
        
        # BS1 - Allele frequency too high
        benign['BS1'] = self._evaluate_bs1(variant_data)
        
        # BS2 - Healthy individual
        benign['BS2'] = self._evaluate_bs2(variant_data)
        
        # BS3 - Functional studies
        benign['BS3'] = self._evaluate_bs3(variant_data)
        
        # BS4 - Lack of segregation
        benign['BS4'] = self._evaluate_bs4(variant_data)
        
        # BP1 - Missense in truncating gene
        benign['BP1'] = self._evaluate_bp1(variant_data)
        
        # BP2 - Observed in trans/cis
        benign['BP2'] = self._evaluate_bp2(variant_data)
        
        # BP3 - In-frame indel in repeat
        benign['BP3'] = self._evaluate_bp3(variant_data)
        
        # BP4 - Computational evidence suggests no impact
        benign['BP4'] = self._evaluate_bp4(variant_data)
        
        # BP5 - Variant found in case with alternate cause
        benign['BP5'] = self._evaluate_bp5(variant_data)
        
        # BP6 - Reputable source reports benign
        benign['BP6'] = self._evaluate_bp6(variant_data)
        
        # BP7 - Synonymous variant with no impact
        benign['BP7'] = self._evaluate_bp7(variant_data)
        
        return benign
    
    def _evaluate_pvs1(self, variant_data) -> Dict[str, Any]:
        """Evaluate PVS1 - Null variant in LOF gene."""
        result = {'applies': False, 'strength': 'Very Strong', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        consequence = variant_data.basic_info.get('consequence', '').lower()
        gene = variant_data.basic_info.get('gene', '').upper()
        
        # Check if it's a loss-of-function variant
        lof_types = ['nonsense', 'frameshift', 'splice_donor', 'splice_acceptor', 'start_lost']
        lof_consequences = ['stop_gained', 'frameshift_variant', 'splice_donor_variant', 
                           'splice_acceptor_variant', 'start_lost']
        
        if variant_type in lof_types or consequence in lof_consequences:
            # Check if gene is LOF intolerant (required for PVS1)
            if gene in LOF_INTOLERANT_GENES:
                result['applies'] = True
                result['details'] = f"Loss-of-function variant ({variant_type or consequence}) in LOF intolerant gene {gene}"
            elif gene in LOF_TOLERANT_GENES:
                result['applies'] = False
                result['details'] = f"Loss-of-function variant ({variant_type or consequence}) in LOF tolerant gene {gene} - PVS1 not applicable"
            else:
                # Unknown gene - conservative approach, don't apply PVS1
                result['applies'] = False
                result['details'] = f"Loss-of-function variant ({variant_type or consequence}) in gene {gene} with unknown LOF tolerance - manual review required"
        
        # Check for intronic variants with high SpliceAI scores (splice-altering)
        elif variant_type == 'intronic':
            # Get SpliceAI scores
            spliceai_scores = []
            for score_key in ['spliceai_ag_score', 'spliceai_al_score', 'spliceai_dg_score', 'spliceai_dl_score']:
                if score_key in variant_data.insilico_data and variant_data.insilico_data[score_key] is not None:
                    score = variant_data.insilico_data[score_key]
                    spliceai_scores.append(score)
            
            # If any SpliceAI score is very high (>0.5), consider it splice-altering
            if spliceai_scores and max(spliceai_scores) > 0.5:
                # Check if gene is LOF intolerant
                if gene in LOF_INTOLERANT_GENES:
                    result['applies'] = True
                    result['strength'] = 'Strong'  # Use PS1 strength for splice-altering variants
                    result['details'] = f"Intronic variant with high SpliceAI score (max: {max(spliceai_scores):.3f}) in LOF intolerant gene {gene}"
                else:
                    result['applies'] = False
                    result['details'] = f"Intronic variant with high SpliceAI score (max: {max(spliceai_scores):.3f}) in gene {gene} - LOF tolerance unknown"
            elif spliceai_scores and max(spliceai_scores) > 0.2:
                result['applies'] = False
                result['details'] = f"Intronic variant with moderate SpliceAI score (max: {max(spliceai_scores):.3f}) - consider PS1 or PM4"
            else:
                result['details'] = "Intronic variant with low splice impact prediction"
        
        return result
    
    def _evaluate_ps1(self, variant_data) -> Dict[str, Any]:
        """Evaluate PS1 - Same amino acid change as pathogenic variant."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        gene = variant_data.basic_info.get('gene')
        aa_change = variant_data.basic_info.get('amino_acid_change')
        
        if gene and aa_change:
            if self.test_mode:
                # Test mode: Return default result without interaction
                result['details'] = f"Test mode: Manual literature review required for PS1 evaluation"
                result['manual_review'] = True
            else:
                # Interactive mode: Ask user for literature review results
                result = self._evaluate_ps1_interactive(variant_data, gene, aa_change)
        else:
            result['details'] = "Insufficient variant information for PS1 evaluation"
        
        return result
    
    def _evaluate_ps1_interactive(self, variant_data, gene, aa_change) -> Dict[str, Any]:
        """Interactive PS1 evaluation with user input."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        print(f"\n🔍 PS1 Evaluation: {gene} {aa_change}")
        print("─" * 50)
        print("QUESTION: Have you found the exact same amino acid change reported as pathogenic?")
        print()
        print("📋 Search recommendations:")
        print(f"   • PubMed: '{gene} {aa_change} pathogenic'")
        print(f"   • ClinVar: '{gene} {aa_change}'")
        print(f"   • HGMD: Same amino acid change")
        print(f"   • Literature: Functional studies on this change")
        print()
        print("💡 PS1 applies if the EXACT same amino acid change has been reported as pathogenic")
        print("   with sufficient evidence in the literature.")
        print()
        
        while True:
            choice = input("PS1 Evidence Found? (y/n/u for unknown): ").strip().lower()
            if choice in ['y', 'yes']:
                result['applies'] = True
                result['details'] = f"Same amino acid change {aa_change} in {gene} reported as pathogenic (PS1)"
                print(f"✅ PS1 applies: {result['details']}")
                break
            elif choice in ['n', 'no']:
                result['details'] = f"No evidence found for same amino acid change {aa_change} in {gene}"
                print(f"❌ PS1 does not apply: {result['details']}")
                break
            elif choice in ['u', 'unknown']:
                result['details'] = f"Evidence unclear for {aa_change} in {gene} - requires further investigation"
                print(f"⚠️  PS1 unclear: {result['details']}")
                break
            else:
                print("❌ Please enter 'y' for yes, 'n' for no, or 'u' for unknown")
        
        return result
    
    def _evaluate_ps2(self, variant_data) -> Dict[str, Any]:
        """Evaluate PS2 - De novo variant with enhanced logic."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        # Import the new de novo assignment logic
        # from config.constants import assign_de_novo_criteria
        
        # Try both genetic_data and functional_data for de novo status
        denovo_status = (variant_data.genetic_data.get('de_novo') or 
                        variant_data.genetic_data.get('denovo') or
                        variant_data.functional_data.get('de_novo') or
                        variant_data.functional_data.get('denovo'))
        
        # Get parental confirmation status
        maternity_confirmed = variant_data.genetic_data.get('maternity_confirmed')
        paternity_confirmed = variant_data.genetic_data.get('paternity_confirmed')
        
        # Mock de novo assignment logic
        # de_novo_assignment = assign_de_novo_criteria(
        #     denovo_status, 
        #     maternity_confirmed, 
        #     paternity_confirmed
        # )
        
        # Apply PS2 only for confirmed de novo WITH parental testing
        if denovo_status == 'confirmed' and maternity_confirmed and paternity_confirmed:
            result['applies'] = True
            if self.use_2023_guidelines:
                result['strength'] = 'Very Strong'
            result['details'] = 'De novo variant confirmed with parental testing'
        elif denovo_status == 'assumed':
            # For assumed de novo, use PM6 instead of PS2
            result['applies'] = False
            result['details'] = 'De novo assumed but not confirmed - see PM6'
        elif denovo_status == 'confirmed' and not (maternity_confirmed and paternity_confirmed):
            # Marked as confirmed but no parental testing - treat as assumed
            result['applies'] = False
            result['details'] = 'De novo confirmed but no parental testing - see PM6'
        else:
            result['details'] = 'No evidence of de novo variant'
        
        return result
    
    def _evaluate_ps3(self, variant_data) -> Dict[str, Any]:
        # Automated assignment for PS3/BS3 using FunctionalStudiesEvaluator
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        ps3_bs3 = self.functional_studies_evaluator.evaluate_functional_evidence(variant_data)
        if ps3_bs3 == 'PS3':
            result['applies'] = True
            result['details'] = "Functional studies support a damaging effect (PS3)"
        elif ps3_bs3 == 'BS3':
            result['details'] = "Functional studies support a benign effect (BS3)"
        else:
            result['details'] = "No or inconclusive functional studies data"
        return result
    
    def _evaluate_ps4(self, variant_data) -> Dict[str, Any]:
        """Evaluate PS4 - The prevalence of the variant in affected individuals is significantly increased compared to the prevalence in controls."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        gene = variant_data.basic_info.get('gene')
        variant_name = variant_data.basic_info.get('variant_name', 'variant')
        
        if gene:
            if self.test_mode:
                result['details'] = f"Test mode: Case-control study data not available for {gene} {variant_name}"
                result['manual_review'] = True
            else:
                # Interactive evaluation for case-control data
                result = self._evaluate_ps4_interactive(variant_data, gene, variant_name)
        else:
            result['details'] = "Insufficient variant information for PS4 evaluation"
        
        return result
    
    def _evaluate_ps4_interactive(self, variant_data, gene, variant_name) -> Dict[str, Any]:
        """Interactive PS4 evaluation with user input."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        print(f"\n🔍 PS4 Evaluation: {gene} {variant_name}")
        print("─" * 50)
        print("QUESTION: Is there case-control data showing significantly higher frequency in affected individuals?")
        print()
        print("📋 Search recommendations:")
        print(f"   • PubMed: Case-control studies for {gene}")
        print(f"   • ClinVar: Review submitted interpretations")
        print(f"   • GWAS studies for this condition")
        print(f"   • Population genetics studies")
        print()
        print("💡 PS4 applies if variant frequency is significantly higher in affected")
        print("   individuals compared to controls with appropriate statistical significance.")
        print()
        
        while True:
            choice = input("PS4 Evidence Found? (y/n/u for unknown): ").strip().lower()
            if choice in ['y', 'yes']:
                result['applies'] = True
                result['details'] = f"Case-control data shows significant enrichment in affected individuals for {gene} (PS4)"
                print(f"✅ PS4 applies: {result['details']}")
                break
            elif choice in ['n', 'no']:
                result['details'] = f"No significant case-control evidence found for {gene}"
                print(f"❌ PS4 does not apply: {result['details']}")
                break
            elif choice in ['u', 'unknown']:
                result['details'] = f"Case-control data unclear or insufficient for {gene}"
                print(f"⚠️  PS4 unclear: {result['details']}")
                break
            else:
                print("❌ Please enter 'y' for yes, 'n' for no, or 'u' for unknown")
        
        return result
    
    def _evaluate_pm1(self, variant_data) -> Dict[str, Any]:
        """Evaluate PM1 - Located in a mutational hot spot and/or critical functional domain."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        gene = variant_data.basic_info.get('gene')
        aa_change = variant_data.basic_info.get('amino_acid_change')
        
        if gene and aa_change:
            # Define known hotspot regions for common genes
            hotspot_genes = {
                'TP53': ['DNA_binding_domain', 'tetramerization_domain'],
                'KRAS': ['GTPase_domain'],
                'PIK3CA': ['helical_domain', 'kinase_domain'],
                'BRAF': ['kinase_domain'],
                'EGFR': ['kinase_domain', 'tyrosine_kinase_domain']
            }
            
            if gene.upper() in hotspot_genes:
                result['details'] = f"Gene {gene} has known hotspot regions - manual domain analysis required"
                result['manual_review'] = True
                result['search_recommendations'] = [
                    f"Check if variant is in known functional domain for {gene}",
                    f"Review protein structure and critical domains",
                    f"Check ClinVar for clustering of pathogenic variants",
                    f"Review literature for hotspot regions in {gene}"
                ]
                result['guidance'] = "PM1 applies if variant is in a mutational hotspot or critical functional domain"
            else:
                result['details'] = f"Manual domain analysis required for {gene}"
                result['manual_review'] = True
                result['search_recommendations'] = [
                    f"Identify functional domains for {gene}",
                    f"Check for clustering of pathogenic variants",
                    f"Review protein structure databases",
                    f"Check UniProt for domain annotations"
                ]
                result['guidance'] = "PM1 applies if variant is in a mutational hotspot or critical functional domain"
        else:
            result['details'] = "Insufficient variant information for PM1 evaluation"
        
        return result
    
    def _evaluate_pm2(self, variant_data) -> Dict[str, Any]:
        """Evaluate PM2 - Absent from controls."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        # Check population frequencies
        pop_data = variant_data.population_data
        gnomad_af = pop_data.get('gnomad_af')
        
        # PM2 applies if variant is absent or very rare in population databases
        if gnomad_af is None or gnomad_af == 0:
            result['applies'] = True
            result['details'] = "Variant absent from population databases"
        elif gnomad_af < 0.0001:  # Very rare
            result['applies'] = True
            result['details'] = f"Variant very rare in population (gnomAD AF: {gnomad_af})"
        else:
            result['details'] = f"Variant present in population databases (gnomAD AF: {gnomad_af})"
        
        return result
    
    def _evaluate_pm3(self, variant_data) -> Dict[str, Any]:
        """Evaluate PM3 - For recessive disorders, detected in trans with a pathogenic variant."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        gene = variant_data.basic_info.get('gene')
        variant_name = variant_data.basic_info.get('variant_name', 'variant')
        
        if gene:
            if self.test_mode:
                result['details'] = f"Test mode: Trans analysis not available for {gene} {variant_name}"
                result['manual_review'] = True
            else:
                # Interactive evaluation for trans analysis
                result = self._evaluate_pm3_interactive(variant_data, gene, variant_name)
        else:
            result['details'] = "Insufficient variant information for PM3 evaluation"
        
        return result
    
    def _evaluate_pm3_interactive(self, variant_data, gene, variant_name) -> Dict[str, Any]:
        """Interactive PM3 evaluation with user input."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        print(f"\n🔍 PM3 Evaluation: {gene} {variant_name}")
        print("─" * 50)
        print("QUESTION: For recessive disorders, is this variant detected in trans with a pathogenic variant?")
        print()
        print("📋 Search recommendations:")
        print(f"   • Check inheritance pattern for {gene}")
        print(f"   • Review family segregation data")
        print(f"   • Check for compound heterozygous variants")
        print(f"   • Verify parental origin if available")
        print()
        print("💡 PM3 applies only for recessive disorders when variant is in trans")
        print("   with a known pathogenic variant.")
        print()
        
        while True:
            choice = input("PM3 Evidence Found? (y/n/u for unknown): ").strip().lower()
            if choice in ['y', 'yes']:
                result['applies'] = True
                result['details'] = f"Variant in trans with pathogenic variant in recessive disorder for {gene} (PM3)"
                print(f"✅ PM3 applies: {result['details']}")
                break
            elif choice in ['n', 'no']:
                result['details'] = f"No trans pathogenic variant found or not recessive disorder for {gene}"
                print(f"❌ PM3 does not apply: {result['details']}")
                break
            elif choice in ['u', 'unknown']:
                result['details'] = f"Trans analysis unclear or insufficient data for {gene}"
                print(f"⚠️  PM3 unclear: {result['details']}")
                break
            else:
                print("❌ Please enter 'y' for yes, 'n' for no, or 'u' for unknown")
        
        return result
    
    def _evaluate_pm4(self, variant_data) -> Dict[str, Any]:
        """Evaluate PM4 - Protein length changes due to inframe indels."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        consequence = variant_data.basic_info.get('consequence', '').lower()
        
        # Check for inframe indels
        if variant_type == 'inframe_indel' or 'inframe' in consequence:
            result['applies'] = True
            result['details'] = "Inframe indel causing protein length change"
        elif variant_type in ['insertion', 'deletion']:
            # Could be inframe - need to check if length is multiple of 3
            result['details'] = "Indel variant - manual review required to determine if inframe"
        else:
            result['details'] = "No protein length changes detected"
        
        return result
    
    def _evaluate_pm5(self, variant_data) -> Dict[str, Any]:
        """Evaluate PM5 - Novel missense change at same residue as pathogenic variant."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        gene = variant_data.basic_info.get('gene')
        aa_change = variant_data.basic_info.get('amino_acid_change')
        
        if gene and aa_change:
            if self.test_mode:
                # Test mode: Return default result without interaction
                result['details'] = f"Test mode: Manual literature review required for PM5 evaluation"
                result['manual_review'] = True
            else:
                # Interactive mode: Ask user for literature review results
                result = self._evaluate_pm5_interactive(variant_data, gene, aa_change)
        else:
            result['details'] = "Insufficient variant information for PM5 evaluation"
        
        return result
    
    def _evaluate_pm5_interactive(self, variant_data, gene, aa_change) -> Dict[str, Any]:
        """Interactive PM5 evaluation with user input."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        print(f"\n🔍 PM5 Evaluation: {gene} {aa_change}")
        print("─" * 50)
        print("QUESTION: Have you found a DIFFERENT pathogenic missense change at the same amino acid position?")
        print()
        print("📋 Search recommendations:")
        print(f"   • ClinVar: Same position, different amino acid")
        print(f"   • HGMD: Same residue, different change")
        print(f"   • Literature: Functional studies on this position")
        print(f"   • Verify the current change is novel (not previously reported)")
        print()
        print("💡 PM5 applies if a DIFFERENT pathogenic missense change has been reported")
        print("   at the same amino acid position.")
        print()
        
        while True:
            choice = input("PM5 Evidence Found? (y/n/u for unknown): ").strip().lower()
            if choice in ['y', 'yes']:
                result['applies'] = True
                result['details'] = f"Different pathogenic missense change found at same position in {gene} (PM5)"
                print(f"✅ PM5 applies: {result['details']}")
                break
            elif choice in ['n', 'no']:
                result['details'] = f"No different pathogenic changes found at same position in {gene}"
                print(f"❌ PM5 does not apply: {result['details']}")
                break
            elif choice in ['u', 'unknown']:
                result['details'] = f"Evidence unclear for same position variants in {gene}"
                print(f"⚠️  PM5 unclear: {result['details']}")
                break
            else:
                print("❌ Please enter 'y' for yes, 'n' for no, or 'u' for unknown")
        
        return result
    
    def _evaluate_pm6(self, variant_data) -> Dict[str, Any]:
        """Evaluate PM6 - Assumed de novo."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        # Try both genetic_data and functional_data for de novo status
        denovo_status = (variant_data.genetic_data.get('de_novo') or 
                        variant_data.genetic_data.get('denovo') or
                        variant_data.functional_data.get('de_novo') or
                        variant_data.functional_data.get('denovo'))
        
        # Get parental confirmation status
        maternity_confirmed = variant_data.genetic_data.get('maternity_confirmed', False)
        paternity_confirmed = variant_data.genetic_data.get('paternity_confirmed', False)
        
        # Apply PM6 for assumed de novo (without parental confirmation)
        if denovo_status == 'assumed':
            result['applies'] = True
            result['details'] = 'De novo variant assumed but not confirmed with parental testing'
        elif denovo_status == 'confirmed' and not (maternity_confirmed and paternity_confirmed):
            # Even if marked as confirmed, if no parental testing, treat as assumed
            result['applies'] = True
            result['details'] = 'De novo variant without parental confirmation (assumed)'
        else:
            result['details'] = 'No evidence of de novo variant or confirmed with parental testing'
        
        return result
    
    def _evaluate_pp1(self, variant_data) -> Dict[str, Any]:
        """Evaluate PP1 - Cosegregation with disease in multiple affected family members."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        # Check segregation data
        segregation = variant_data.genetic_data.get('segregation')
        
        if segregation == 'cosegregates':
            result['applies'] = True
            result['details'] = "Variant cosegregates with disease in family (PP1)"
        elif segregation == 'does_not_segregate':
            result['details'] = "Variant does not cosegregate with disease - see BS4"
        elif segregation == 'insufficient_data':
            result['details'] = "Insufficient segregation data for PP1 evaluation"
        elif segregation == 'not_performed':
            result['details'] = "Segregation analysis not performed"
        else:
            result['details'] = "Segregation analysis required for PP1 evaluation"
            result['manual_review'] = True
            result['search_recommendations'] = [
                "Analyze affected family members for variant presence",
                "Check unaffected family members for variant absence", 
                "Consider statistical significance of segregation",
                "Review inheritance pattern consistency"
            ]
            result['guidance'] = "PP1 applies if variant cosegregates with disease in multiple affected family members"
        
        return result
    
    def _evaluate_pp2(self, variant_data) -> Dict[str, Any]:
        """Evaluate PP2 - Missense variant in a gene that has a low rate of benign missense variation."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        gene = variant_data.basic_info.get('gene')
        
        if variant_type == 'missense':
            # Define genes known to have low benign missense variation
            low_benign_genes = [
                'BRCA1', 'BRCA2', 'TP53', 'ATM', 'CHEK2', 'PALB2',
                'MLH1', 'MSH2', 'MSH6', 'PMS2', 'EPCAM',
                'APC', 'MUTYH', 'NTHL1', 'POLD1', 'POLE',
                'STK11', 'PTEN', 'CDH1', 'CTNNA1',
                'RB1', 'NF1', 'NF2', 'TSC1', 'TSC2',
                'VHL', 'MEN1', 'RET', 'SDHB', 'SDHC', 'SDHD'
            ]
            
            if gene and gene.upper() in low_benign_genes:
                result['applies'] = True
                result['details'] = f"Missense variant in {gene} (gene with low rate of benign missense variation) (PP2)"
            else:
                result['details'] = f"Gene {gene} not recognized as having low benign missense variation rate"
                result['manual_review'] = True
                result['search_recommendations'] = [
                    f"Check gnomAD constraint metrics for {gene}",
                    f"Review missense Z-score and constraint for {gene}",
                    f"Compare benign vs pathogenic missense rates in ClinVar",
                    f"Check gene-specific guidelines for {gene}"
                ]
                result['guidance'] = "PP2 applies if the gene has a low rate of benign missense variation and high constraint"
        else:
            result['details'] = "Not a missense variant"
        
        return result
    
    def _evaluate_pp3(self, variant_data) -> Dict[str, Any]:
        """Evaluate PP3 - Computational evidence suggests pathogenic impact."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        # Check in silico predictors
        insilico_data = variant_data.insilico_data
        
        # Count pathogenic predictions
        pathogenic_count = 0
        total_count = 0
        
        predictor_thresholds = {
            'cadd_phred': 20,
            'revel': 0.5,
            'sift': 0.05,  # Lower is more damaging
            'polyphen2': 0.5,
            'mutation_taster': 0.5,
            'dann': 0.5,
            'fathmm': 0.5,  # Lower is more damaging
        }
        
        for predictor, threshold in predictor_thresholds.items():
            score = insilico_data.get(predictor)
            if score is not None:
                total_count += 1
                if predictor in ['sift', 'fathmm']:
                    # For SIFT and FATHMM, lower values are more damaging
                    if score < threshold:
                        pathogenic_count += 1
                else:
                    # For other predictors, higher values are more damaging
                    if score > threshold:
                        pathogenic_count += 1
        
        # Apply PP3 if majority of predictors suggest pathogenic
        if total_count > 0:
            pathogenic_ratio = pathogenic_count / total_count
            if pathogenic_ratio >= 0.7:  # 70% or more predict pathogenic
                result['applies'] = True
                result['details'] = f"Computational evidence suggests pathogenic impact ({pathogenic_count}/{total_count} predictors)"
            else:
                result['details'] = f"Computational evidence inconclusive ({pathogenic_count}/{total_count} predictors suggest pathogenic)"
        else:
            result['details'] = "No computational prediction scores available"
        
        return result
    
    def _evaluate_pp4(self, variant_data) -> Dict[str, Any]:
        # Automated assignment for PP4/BP5 using PhenotypeMatcher
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        patient_phenotypes = getattr(variant_data, 'patient_phenotypes', None)
        gene = getattr(variant_data, 'gene', None)
        pp4_bp5 = self.phenotype_matcher.evaluate_phenotype_match(variant_data, patient_phenotypes)
        if pp4_bp5 == 'PP4':
            result['applies'] = True
            result['details'] = f"Patient phenotype highly specific for {gene}-related disease (PP4)"
        elif pp4_bp5 == 'BP5':
            result['details'] = f"Patient phenotype not consistent with {gene}-related disease (BP5)"
        else:
            result['details'] = "No or inconclusive phenotype data"
        return result
    
    def _evaluate_pp5(self, variant_data) -> Dict[str, Any]:
        """Evaluate PP5 - Reputable source recently reports variant as pathogenic."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        # Check if ClinVar data is available
        clinvar_data = getattr(variant_data, 'clinvar_data', None)
        
        if clinvar_data and clinvar_data.get('clinical_significance'):
            significance = clinvar_data['clinical_significance'].lower()
            if 'pathogenic' in significance and 'benign' not in significance:
                result['applies'] = True
                result['details'] = f"Reputable source (ClinVar) reports variant as pathogenic (PP5)"
            elif 'benign' in significance:
                result['details'] = f"Reputable source reports variant as benign - see BP6"
            else:
                result['details'] = f"ClinVar classification: {significance}"
        else:
            result['details'] = "Reputable source analysis required for PP5 evaluation"
            result['manual_review'] = True
            result['search_recommendations'] = [
                "Check ClinVar for clinical significance",
                "Review recent literature reports",
                "Check laboratory databases (HGMD, LOVD)",
                "Verify source credibility and recency"
            ]
            result['guidance'] = "PP5 applies if a reputable source recently reports the variant as pathogenic"
        
        return result
    
    def _evaluate_ba1(self, variant_data) -> Dict[str, Any]:
        """Evaluate BA1 - High allele frequency."""
        result = {'applies': False, 'strength': 'Stand-alone', 'details': ''}
        
        # Check population frequencies
        pop_data = variant_data.population_data
        gnomad_af = pop_data.get('gnomad_af')
        gene = variant_data.basic_info.get('gene', '').upper()
        
        # Get gene-specific threshold or use default
        gene_thresholds = GENE_SPECIFIC_THRESHOLDS.get(gene, GENE_SPECIFIC_THRESHOLDS['default'])
        ba1_threshold = gene_thresholds.get('BA1', 0.05)
        
        # BA1 applies if variant is common in population
        if gnomad_af is not None and gnomad_af > ba1_threshold:
            result['applies'] = True
            result['details'] = f"Variant common in population (gnomAD AF: {gnomad_af}, threshold: {ba1_threshold})"
        else:
            result['details'] = f"Variant not common in population (gnomAD AF: {gnomad_af or 'N/A'}, threshold: {ba1_threshold})"
        
        return result
    
    def _evaluate_bs1(self, variant_data) -> Dict[str, Any]:
        """Evaluate BS1 - Allele frequency higher than expected for disorder."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        # Check population frequencies
        pop_data = variant_data.population_data
        gnomad_af = pop_data.get('gnomad_af')
        gene = variant_data.basic_info.get('gene', '').upper()
        
        # Get gene-specific threshold or use default
        gene_thresholds = GENE_SPECIFIC_THRESHOLDS.get(gene, GENE_SPECIFIC_THRESHOLDS['default'])
        bs1_threshold = gene_thresholds.get('BS1', 0.01)
        
        # BS1 applies if variant frequency is higher than expected for disorder
        if gnomad_af is not None and gnomad_af > bs1_threshold:
            result['applies'] = True
            result['details'] = f"Variant frequency higher than expected for disorder (gnomAD AF: {gnomad_af}, threshold: {bs1_threshold})"
        else:
            result['details'] = f"Variant frequency not higher than expected (gnomAD AF: {gnomad_af or 'N/A'}, threshold: {bs1_threshold})"
        
        return result
    
    def _evaluate_bs2(self, variant_data) -> Dict[str, Any]:
        """Evaluate BS2 - Observed in a healthy adult individual for a recessive disorder."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        # Check case-control data
        case_control = variant_data.functional_data.get('case_control', 'not_available')
        inheritance = variant_data.genetic_data.get('inheritance', 'unknown')
        
        if case_control == 'observed_in_healthy' and inheritance in ['recessive', 'autosomal_recessive']:
            result['applies'] = True
            result['details'] = "Variant observed in healthy individual for recessive disorder (BS2)"
        elif case_control == 'observed_in_healthy' and inheritance not in ['recessive', 'autosomal_recessive']:
            result['details'] = "Variant observed in healthy individual but inheritance not recessive"
        elif case_control == 'not_observed_in_healthy':
            result['details'] = "Variant not observed in healthy individuals"
        else:
            result['details'] = "Case-control analysis required for BS2 evaluation"
            result['manual_review'] = True
            result['search_recommendations'] = [
                "Check if variant is observed in healthy individuals",
                "Confirm inheritance pattern is recessive",
                "Review population databases for healthy carriers",
                "Consider age-related penetrance"
            ]
            result['guidance'] = "BS2 applies if variant is observed in healthy adults for a recessive disorder"
        
        return result
    
    def _evaluate_bs3(self, variant_data) -> Dict[str, Any]:
        # Automated assignment for BS3 using FunctionalStudiesEvaluator
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        ps3_bs3 = self.functional_studies_evaluator.evaluate_functional_evidence(variant_data)
        if ps3_bs3 == 'BS3':
            result['applies'] = True
            result['details'] = "Functional studies support a benign effect (BS3)"
        else:
            result['details'] = "No or inconclusive functional studies data"
        return result
    
    def _evaluate_bs4(self, variant_data) -> Dict[str, Any]:
        """Evaluate BS4 - Lack of segregation in affected members of a family."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        # Check segregation data
        segregation = variant_data.genetic_data.get('segregation')
        
        if segregation == 'does_not_segregate':
            result['applies'] = True
            result['details'] = "Variant does not segregate with disease in family (BS4)"
        elif segregation == 'cosegregates':
            result['details'] = "Variant cosegregates with disease - see PP1"
        elif segregation == 'insufficient_data':
            result['details'] = "Insufficient segregation data for BS4 evaluation"
        elif segregation == 'not_performed':
            result['details'] = "Segregation analysis not performed"
        else:
            result['details'] = "Segregation analysis required for BS4 evaluation"
            result['manual_review'] = True
            result['search_recommendations'] = [
                "Check for variant presence in unaffected family members",
                "Verify variant absence in affected family members",
                "Consider statistical significance of non-segregation",
                "Review family structure and inheritance pattern"
            ]
            result['guidance'] = "BS4 applies if variant does not segregate with disease in affected family members"
        
        return result
    
    def _evaluate_bp1(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP1 - Missense variant in a gene for which primarily truncating variants cause disease."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        gene = variant_data.basic_info.get('gene')
        
        if variant_type == 'missense':
            # Define genes where primarily truncating variants cause disease
            truncating_mechanism_genes = [
                'APC', 'BRCA1', 'BRCA2', 'NF1', 'NF2', 'TSC1', 'TSC2',
                'MLH1', 'MSH2', 'MSH6', 'PMS2', 'EPCAM',
                'RB1', 'WT1', 'VHL', 'SMAD4', 'BMPR1A',
                'STK11', 'CDH1', 'PTEN', 'CDKN2A'
            ]
            
            if gene and gene.upper() in truncating_mechanism_genes:
                result['applies'] = True
                result['details'] = f"Missense variant in {gene} (gene where truncating variants primarily cause disease) (BP1)"
            else:
                result['details'] = f"Gene {gene} not recognized as primarily truncating mechanism"
                result['manual_review'] = True
                result['search_recommendations'] = [
                    f"Review ClinVar pathogenic variants for {gene}",
                    f"Check if truncating variants are the primary mechanism",
                    f"Review gene-specific guidelines for {gene}",
                    f"Check OMIM and literature for mechanism of disease"
                ]
                result['guidance'] = "BP1 applies if truncating variants are the primary disease mechanism for this gene"
        else:
            result['details'] = "Not a missense variant"
        
        return result
    
    def _evaluate_bp2(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP2 - Observed in trans with pathogenic variant for fully penetrant dominant gene/disorder."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        gene = variant_data.basic_info.get('gene')
        variant_name = variant_data.basic_info.get('variant_name', 'variant')
        
        if gene:
            if self.test_mode:
                result['details'] = f"Test mode: Trans analysis not available for {gene} {variant_name}"
                result['manual_review'] = True
            else:
                # Interactive evaluation for trans analysis
                result = self._evaluate_bp2_interactive(variant_data, gene, variant_name)
        else:
            result['details'] = "Insufficient variant information for BP2 evaluation"
        
        return result
    
    def _evaluate_bp2_interactive(self, variant_data, gene, variant_name) -> Dict[str, Any]:
        """Interactive BP2 evaluation with user input."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        print(f"\n🔍 BP2 Evaluation: {gene} {variant_name}")
        print("─" * 50)
        print("QUESTION: For dominant disorders, is this variant observed in trans with a pathogenic variant?")
        print()
        print("📋 Search recommendations:")
        print(f"   • Check inheritance pattern for {gene} (must be dominant)")
        print(f"   • Review family segregation data")
        print(f"   • Check for compound heterozygous variants")
        print(f"   • Verify parental origin if available")
        print()
        print("💡 BP2 applies only for fully penetrant dominant disorders when variant")
        print("   is observed in trans with a known pathogenic variant.")
        print()
        
        while True:
            choice = input("BP2 Evidence Found? (y/n/u for unknown): ").strip().lower()
            if choice in ['y', 'yes']:
                result['applies'] = True
                result['details'] = f"Variant in trans with pathogenic variant in dominant disorder for {gene} (BP2)"
                print(f"✅ BP2 applies: {result['details']}")
                break
            elif choice in ['n', 'no']:
                result['details'] = f"No trans pathogenic variant found or not dominant disorder for {gene}"
                print(f"❌ BP2 does not apply: {result['details']}")
                break
            elif choice in ['u', 'unknown']:
                result['details'] = f"Trans analysis unclear or insufficient data for {gene}"
                print(f"⚠️  BP2 unclear: {result['details']}")
                break
            else:
                print("❌ Please enter 'y' for yes, 'n' for no, or 'u' for unknown")
        
        return result
    
    def _evaluate_bp3(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP3 - In-frame deletions/insertions in repetitive region."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        consequence = variant_data.basic_info.get('consequence', '').lower()
        
        # Check if it's an in-frame indel
        if variant_type == 'inframe_indel' or 'inframe' in consequence:
            # This would ideally check if the variant is in a repetitive region
            # For now, we'll apply it to in-frame indels and ask for user confirmation
            if self.test_mode:
                result['details'] = "Test mode: In-frame indel - manual review required for repetitive region confirmation"
            else:
                # Interactive evaluation for repetitive region
                result = self._evaluate_bp3_interactive(variant_data)
        else:
            result['details'] = "Not an in-frame indel variant"
        
        return result
    
    def _evaluate_bp3_interactive(self, variant_data) -> Dict[str, Any]:
        """Interactive BP3 evaluation for repetitive regions."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        gene = variant_data.basic_info.get('gene')
        variant_type = variant_data.basic_info.get('variant_type')
        
        print(f"\n🔍 BP3 Evaluation: {gene} {variant_type}")
        print("─" * 50)
        print("QUESTION: Is this in-frame indel located in a repetitive region without known function?")
        print()
        print("📋 Consider:")
        print("   • Repetitive sequences (tandem repeats, homopolymers)")
        print("   • Regions with no known functional importance")
        print("   • Protein domains that tolerate length variation")
        print("   • Areas outside critical functional domains")
        print()
        print("💡 BP3 applies if the indel is in a repetitive region without known function.")
        print()
        
        while True:
            choice = input("Is this in a repetitive region? (y/n/u for unknown): ").strip().lower()
            if choice in ['y', 'yes']:
                result['applies'] = True
                result['details'] = f"In-frame indel in repetitive region without known function (BP3)"
                print(f"✅ BP3 applies: {result['details']}")
                break
            elif choice in ['n', 'no']:
                result['details'] = f"In-frame indel not in repetitive region"
                print(f"❌ BP3 does not apply: {result['details']}")
                break
            elif choice in ['u', 'unknown']:
                result['details'] = f"Repetitive region status unclear - manual review required"
                print(f"⚠️  BP3 unclear: {result['details']}")
                break
            else:
                print("❌ Please enter 'y' for yes, 'n' for no, or 'u' for unknown")
        
        return result
    
    def _evaluate_bp4(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP4 - Computational evidence suggests no impact."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        # Check in silico predictors
        insilico_data = variant_data.insilico_data
        
        # Count benign predictions
        benign_count = 0
        total_count = 0
        
        predictor_thresholds = {
            'cadd_phred': 20,
            'revel': 0.5,
            'sift': 0.05,  # Lower is more damaging
            'polyphen2': 0.5,
            'mutation_taster': 0.5,
            'dann': 0.5,
            'fathmm': 0.5,  # Lower is more damaging
        }
        
        for predictor, threshold in predictor_thresholds.items():
            score = insilico_data.get(predictor)
            if score is not None:
                total_count += 1
                if predictor in ['sift', 'fathmm']:
                    # For SIFT and FATHMM, higher values are more benign
                    if score > threshold:
                        benign_count += 1
                else:
                    # For other predictors, lower values are more benign
                    if score < threshold:
                        benign_count += 1
        
        # Apply BP4 if majority of predictors suggest benign
        if total_count > 0:
            benign_ratio = benign_count / total_count
            if benign_ratio >= 0.7:  # 70% or more predict benign
                result['applies'] = True
                result['details'] = f"Computational evidence suggests no impact ({benign_count}/{total_count} predictors)"
            else:
                result['details'] = f"Computational evidence inconclusive ({benign_count}/{total_count} predictors suggest benign)"
        else:
            result['details'] = "No computational prediction scores available"
        
        return result
    
    def _evaluate_bp5(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP5 - Variant found in a case with an alternate molecular basis for disease."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        gene = variant_data.basic_info.get('gene')
        variant_name = variant_data.basic_info.get('variant_name', 'variant')
        
        if gene:
            if self.test_mode:
                result['details'] = f"Test mode: Alternate molecular basis analysis not available for {gene} {variant_name}"
                result['manual_review'] = True
            else:
                # Interactive evaluation for alternate cause
                result = self._evaluate_bp5_interactive(variant_data, gene, variant_name)
        else:
            result['details'] = "Insufficient variant information for BP5 evaluation"
        
        return result
    
    def _evaluate_bp5_interactive(self, variant_data, gene, variant_name) -> Dict[str, Any]:
        """Interactive BP5 evaluation with user input."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        print(f"\n🔍 BP5 Evaluation: {gene} {variant_name}")
        print("─" * 50)
        print("QUESTION: Has an alternate molecular basis for disease been identified in this case?")
        print()
        print("📋 Search recommendations:")
        print(f"   • Check for other pathogenic variants in the case")
        print(f"   • Review full genetic testing results")
        print(f"   • Check for copy number variants (CNVs)")
        print(f"   • Consider other genes causing similar phenotype")
        print()
        print("💡 BP5 applies if another variant or molecular cause explains")
        print("   the patient's phenotype, making this variant less likely causal.")
        print()
        
        while True:
            choice = input("BP5 Evidence Found? (y/n/u for unknown): ").strip().lower()
            if choice in ['y', 'yes']:
                result['applies'] = True
                result['details'] = f"Alternate molecular basis for disease identified - {gene} variant less likely causal (BP5)"
                print(f"✅ BP5 applies: {result['details']}")
                break
            elif choice in ['n', 'no']:
                result['details'] = f"No alternate molecular basis identified for {gene}"
                print(f"❌ BP5 does not apply: {result['details']}")
                break
            elif choice in ['u', 'unknown']:
                result['details'] = f"Alternate cause analysis unclear for {gene}"
                print(f"⚠️  BP5 unclear: {result['details']}")
                break
            else:
                print("❌ Please enter 'y' for yes, 'n' for no, or 'u' for unknown")
        
        return result
    
    def _evaluate_bp6(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP6 - Reputable source recently reports variant as benign but evidence is not available."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        gene = variant_data.basic_info.get('gene')
        variant_name = variant_data.basic_info.get('variant_name', 'variant')
        
        if gene:
            if self.test_mode:
                result['details'] = f"Test mode: Reputable source analysis not available for {gene} {variant_name}"
                result['manual_review'] = True
            else:
                # Interactive evaluation for reputable source
                result = self._evaluate_bp6_interactive(variant_data, gene, variant_name)
        else:
            result['details'] = "Insufficient variant information for BP6 evaluation"
        
        return result
    
    def _evaluate_bp6_interactive(self, variant_data, gene, variant_name) -> Dict[str, Any]:
        """Interactive BP6 evaluation with user input."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        print(f"\n🔍 BP6 Evaluation: {gene} {variant_name}")
        print("─" * 50)
        print("QUESTION: Has a reputable source recently reported this variant as benign?")
        print()
        print("📋 Check reputable sources:")
        print(f"   • ClinVar: Expert panel classifications")
        print(f"   • Professional lab reports (Quest, LabCorp, etc.)")
        print(f"   • Published clinical guidelines")
        print(f"   • Disease-specific databases")
        print()
        print("💡 BP6 applies when a reputable source classifies variant as benign")
        print("   but detailed evidence is not available to review.")
        print()
        
        while True:
            choice = input("BP6 Evidence Found? (y/n/u for unknown): ").strip().lower()
            if choice in ['y', 'yes']:
                result['applies'] = True
                result['details'] = f"Reputable source reports {gene} variant as benign (BP6)"
                print(f"✅ BP6 applies: {result['details']}")
                break
            elif choice in ['n', 'no']:
                result['details'] = f"No reputable source benign classification found for {gene}"
                print(f"❌ BP6 does not apply: {result['details']}")
                break
            elif choice in ['u', 'unknown']:
                result['details'] = f"Reputable source classification unclear for {gene}"
                print(f"⚠️  BP6 unclear: {result['details']}")
                break
            else:
                print("❌ Please enter 'y' for yes, 'n' for no, or 'u' for unknown")
        
        return result
    
    def _evaluate_bp7(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP7 - Synonymous variant with no predicted splice impact."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        consequence = variant_data.basic_info.get('consequence', '').lower()
        
        # Check if it's a synonymous variant
        if variant_type == 'synonymous' or 'synonymous' in consequence:
            # Check SpliceAI scores to ensure no splice impact
            insilico_data = variant_data.insilico_data
            spliceai_scores = []
            
            for score_key in ['spliceai_ag_score', 'spliceai_al_score', 'spliceai_dg_score', 'spliceai_dl_score']:
                if score_key in insilico_data and insilico_data[score_key] is not None:
                    spliceai_scores.append(insilico_data[score_key])
            
            if spliceai_scores:
                max_spliceai = max(spliceai_scores)
                if max_spliceai < 0.1:  # Very low splice impact
                    result['applies'] = True
                    result['details'] = f"Synonymous variant with no predicted splice impact (SpliceAI max: {max_spliceai:.3f}) (BP7)"
                else:
                    result['details'] = f"Synonymous variant but potential splice impact (SpliceAI max: {max_spliceai:.3f}) - BP7 does not apply"
            else:
                # No SpliceAI data - apply BP7 for synonymous variants
                result['applies'] = True
                result['details'] = "Synonymous variant with no splice prediction data available (BP7)"
        else:
            result['details'] = "Not a synonymous variant"
        
        return result
    
    def _calculate_vampp_score(self, variant_data) -> float:
        """Calculate VAMPP-like score (now called Computational Metascore)."""
        print("    ⚡ Computational Metascore: Calculating based on in silico predictors...")
        
        # Get all available scores
        insilico_data = variant_data.insilico_data
        
        # Calculate metascore based on available predictors
        scores = []
        weights = []
        
        # Define predictor weights and thresholds
        predictor_config = {
            'cadd_phred': {'weight': 0.2, 'threshold': 20},
            'revel': {'weight': 0.2, 'threshold': 0.5},
            'sift': {'weight': 0.15, 'threshold': 0.05, 'inverse': True},  # Lower is more damaging
            'polyphen2': {'weight': 0.15, 'threshold': 0.5},
            'mutation_taster': {'weight': 0.1, 'threshold': 0.5},
            'dann': {'weight': 0.1, 'threshold': 0.5},
            'fathmm': {'weight': 0.1, 'threshold': 0.5, 'inverse': True}
        }
        
        for predictor, config in predictor_config.items():
            score = insilico_data.get(predictor)
            if score is not None:
                # Normalize score to 0-1 scale
                if config.get('inverse', False):
                    # For SIFT and FATHMM, lower values are more damaging
                    normalized = max(0, min(1, (config['threshold'] - score) / config['threshold']))
                else:
                    # For most predictors, higher values are more damaging
                    normalized = max(0, min(1, score / config['threshold']))
                
                scores.append(normalized)
                weights.append(config['weight'])
        
        # Calculate weighted average
        if scores:
            metascore = sum(s * w for s, w in zip(scores, weights)) / sum(weights)
        else:
            metascore = 0.0
            print("    ⚠️  No in silico scores available for metascore calculation")
        
        print(f"    ✅ Computational Metascore: {metascore:.3f} (based on {len(scores)} predictors)")
        
        return metascore
    
    def _perform_statistical_tests(self, variant_data) -> Dict[str, Any]:
        """Perform statistical tests."""
        return {'tests_performed': [], 'results': {}}
    
    def _consolidate_applied_criteria(self, pathogenic, benign) -> Dict[str, Any]:
        """Consolidate applied criteria."""
        applied = {}
        for category, criteria in [('pathogenic', pathogenic), ('benign', benign)]:
            for criterion, result in criteria.items():
                # Check both 'applicable' and 'applies' for compatibility
                if result.get('applicable', False) or result.get('applies', False):
                    applied[criterion] = result
        return applied
