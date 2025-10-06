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
        from utils.statistical_utils import StatisticalAnalyzer
        self.functional_studies_evaluator = FunctionalStudiesEvaluator()
        self.phenotype_matcher = PhenotypeMatcher()
        self.statistical_analyzer = StatisticalAnalyzer()
    
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
        hgvs_p = variant_data.basic_info.get('hgvs_p')
        
        # CRITICAL FIX: Check ClinVar data first
        clinvar = variant_data.clinvar_data or {}
        
        # Check if same AA change is marked as pathogenic in ClinVar
        if clinvar.get('same_aa_pathogenic'):
            result['applies'] = True
            result['details'] = f"Same amino acid change {aa_change or hgvs_p} reported as pathogenic in ClinVar (PS1)"
            return result
        
        # Check pathogenic_variants list for matching amino acid change
        pathogenic_variants = clinvar.get('pathogenic_variants', [])
        if pathogenic_variants and hgvs_p:
            for pv in pathogenic_variants:
                if pv.get('hgvs_p') == hgvs_p:
                    result['applies'] = True
                    result['details'] = f"Same amino acid change {hgvs_p} in {gene} reported as pathogenic (PS1)"
                    return result
        
        # If no automated data, check in test mode or interactive mode
        if gene and (aa_change or hgvs_p):
            if self.test_mode:
                # Test mode: Return default result without interaction
                result['details'] = f"Test mode: No ClinVar data for same AA change"
                result['manual_review'] = True
            else:
                # Interactive mode: Ask user for literature review results
                result = self._evaluate_ps1_interactive(variant_data, gene, aa_change or hgvs_p)
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
        """Evaluate PS2 - De novo variant with enhanced logic and strict 2023 upgrade rules."""
        result = {
            'applies': False, 
            'strength': 'Strong', 
            'details': '',
            'confidence': 'very_low',
            'data_source': 'user_input'
        }
        
        genetic_data = variant_data.genetic_data or {}
        
        # CRITICAL FIX: Check inheritance_pattern for de_novo
        inheritance_pattern = genetic_data.get('inheritance_pattern', '').lower()
        
        # Get parental confirmation status
        maternity_confirmed = genetic_data.get('maternity_confirmed', False)
        paternity_confirmed = genetic_data.get('paternity_confirmed', False)
        family_history = genetic_data.get('family_history', True)
        
        # Try both genetic_data and functional_data for de novo status (legacy support)
        denovo_status = (genetic_data.get('de_novo') or 
                        genetic_data.get('denovo') or
                        variant_data.functional_data.get('de_novo') or
                        variant_data.functional_data.get('denovo'))
        
        # Check if inheritance pattern indicates de novo
        is_de_novo = (inheritance_pattern == 'de_novo' or 
                     denovo_status == 'confirmed' or
                     (inheritance_pattern == 'de_novo' and not family_history))
        
        if not is_de_novo:
            result['details'] = 'No evidence of de novo variant'
            return result
        
        # STRICT 2023 RULES: PS2_Very_Strong only with BOTH parents confirmed
        if maternity_confirmed and paternity_confirmed:
            result['applies'] = True
            result['confidence'] = 'high'
            result['data_source'] = 'parental_testing'
            
            # Upgrade to Very Strong in 2023 guidelines
            if self.use_2023_guidelines:
                result['strength'] = 'Very Strong'
                result['details'] = 'De novo variant confirmed with maternity AND paternity testing (PS2_Very_Strong per ACMG 2023)'
            else:
                result['details'] = 'De novo variant confirmed with parental testing (PS2)'
                
        elif maternity_confirmed or paternity_confirmed:
            # Only one parent confirmed - downgrade to PM6
            result['applies'] = False
            result['confidence'] = 'medium'
            result['data_source'] = 'partial_parental_testing'
            result['details'] = 'De novo confirmed but only one parent tested - see PM6 instead of PS2'
            
        else:
            # De novo indicated but not confirmed with parental testing
            result['applies'] = False
            result['confidence'] = 'low'
            result['data_source'] = 'clinical_assumption'
            result['details'] = 'De novo assumed but not confirmed with parental testing - see PM6'
        
        return result
    
    def _evaluate_ps3(self, variant_data) -> Dict[str, Any]:
        """Evaluate PS3 - Functional studies supportive of damaging effect."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        functional_data = variant_data.functional_data or {}
        
        # CRITICAL FIX: Check has_functional_studies flag first
        has_studies = functional_data.get('has_functional_studies', False)
        functional_impact = functional_data.get('functional_impact', '').lower()
        
        if has_studies:
            # Check functional_impact field
            if functional_impact in ['damaging', 'pathogenic', 'loss_of_function', 'lof']:
                result['applies'] = True
                result['details'] = "Functional studies support a damaging effect (PS3)"
                return result
            elif functional_impact in ['benign', 'neutral', 'normal']:
                result['details'] = "Functional studies support benign effect - see BS3"
                return result
            
            # Check individual functional_studies array
            functional_studies = functional_data.get('functional_studies', [])
            if functional_studies:
                for study in functional_studies:
                    study_result = study.get('result', '').lower()
                    if 'loss' in study_result or 'damaging' in study_result or 'pathogenic' in study_result:
                        result['applies'] = True
                        result['details'] = f"Functional studies show {study_result} (PS3)"
                        return result
        
        # Fall back to FunctionalStudiesEvaluator
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
        """Evaluate PS4 - Case-control data with automated statistical analysis."""
        result = {
            'applies': False, 
            'strength': 'Strong', 
            'details': '',
            'confidence': 'very_low',
            'data_source': 'none'
        }
        
        gene = variant_data.basic_info.get('gene')
        variant_name = variant_data.basic_info.get('variant_name', 'variant')
        
        # CRITICAL FIX: Check population_data for case_control_data
        population_data = variant_data.population_data or {}
        case_control = population_data.get('case_control_data', {})
        
        # Check if simplified case-control data is available
        if case_control.get('significant'):
            odds_ratio = case_control.get('odds_ratio', 0)
            p_value = case_control.get('p_value', 1.0)
            cases = case_control.get('cases', 0)
            controls = case_control.get('controls', 0)
            
            # Select thresholds based on guidelines version
            if self.use_2023_guidelines:
                from config.constants import STATISTICAL_THRESHOLDS_2023
                thresholds = STATISTICAL_THRESHOLDS_2023
                guidelines_label = 'ACMG 2023'
            else:
                from config.constants import STATISTICAL_THRESHOLDS_2015
                thresholds = STATISTICAL_THRESHOLDS_2015
                guidelines_label = 'ACMG 2015'
            
            min_or = thresholds.get('case_control_odds_ratio')
            min_cases = thresholds.get('case_control_min_cases')
            min_controls = thresholds.get('case_control_min_controls')
            
            # Validate criteria
            meets_criteria = (
                odds_ratio >= min_or and
                p_value < 0.05 and
                cases >= min_cases and
                controls >= min_controls
            )
            
            if meets_criteria:
                result['applies'] = True
                result['data_source'] = 'case_control_study'
                result['confidence'] = 'high'
                if self.use_2023_guidelines:
                    result['details'] = f"Case-control study shows significant enrichment (OR={odds_ratio}, p={p_value}, n_cases={cases}, n_controls={controls}) - Meets {guidelines_label} criteria (PS4)"
                else:
                    result['details'] = f"Case-control study shows significant enrichment (OR={odds_ratio}, p={p_value}) (PS4)"
            else:
                result['applies'] = False
                result['confidence'] = 'low'
                result['details'] = f"Case-control data insufficient for {guidelines_label} (OR={odds_ratio} [need ≥{min_or}], cases={cases} [need ≥{min_cases}], controls={controls} [need ≥{min_controls}])"
            
            return result
        
        # Check if case-control data is available in functional_data (legacy format)
        functional_data = variant_data.functional_data or {}
        
        cases_with = functional_data.get('cases_with_variant')
        cases_total = functional_data.get('total_cases')
        controls_with = functional_data.get('controls_with_variant')
        controls_total = functional_data.get('total_controls')
        
        # If numerical data is available, perform automated Fisher's exact test
        if all(v is not None for v in [cases_with, cases_total, controls_with, controls_total]):
            fisher_result = self.statistical_analyzer.calculate_fishers_exact(
                cases_with, cases_total, controls_with, controls_total
            )
            
            if fisher_result['valid']:
                result['data_source'] = 'case_control_study'
                result['confidence'] = fisher_result['confidence']
                result['statistical_test'] = fisher_result
                
                if fisher_result['significant']:
                    result['applies'] = True
                    result['details'] = (
                        f"Case-control analysis shows significant enrichment in cases: "
                        f"{fisher_result['interpretation']}"
                    )
                else:
                    result['details'] = (
                        f"Case-control data available but not significant: "
                        f"{fisher_result['interpretation']}"
                    )
                
                # Store in evidence_details for report
                self.evidence_details['PS4_statistical_analysis'] = fisher_result
                
                return result
            else:
                result['details'] = f"Case-control data invalid: {fisher_result.get('error', 'Unknown error')}"
                return result
        
        # Fall back to interactive evaluation if no numerical data
        if gene:
            if self.test_mode:
                result['details'] = f"Test mode: Case-control study data not available for {gene} {variant_name}"
                result['data_source'] = 'test_mode'
                result['confidence'] = 'very_low'
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
        functional_data = variant_data.functional_data or {}
        
        # CRITICAL FIX: Check functional_data for hotspot/domain flags
        in_hotspot = functional_data.get('in_hotspot', False)
        in_functional_domain = functional_data.get('in_functional_domain', False)
        hotspot_name = functional_data.get('hotspot_name', '')
        domain_name = functional_data.get('domain_name', '')
        benign_variation = functional_data.get('benign_variation_in_domain', True)
        
        # Extract position from HGVS if available for API check
        position = None
        hgvs_p = variant_data.basic_info.get('hgvs_p', '')
        if hgvs_p:
            import re
            match = re.search(r'p\.\w+(\d+)', hgvs_p)
            if match:
                position = int(match.group(1))
        
        # Check if variant is in hotspot
        if in_hotspot:
            result['applies'] = True
            result['details'] = f"Variant in mutational hotspot"
            if hotspot_name:
                result['details'] += f" ({hotspot_name})"
            result['details'] += " (PM1)"
            return result
        
        # Check if variant is in critical functional domain WITHOUT benign variation
        if in_functional_domain and not benign_variation:
            result['applies'] = True
            result['details'] = f"Variant in critical functional domain"
            if domain_name:
                result['details'] += f" ({domain_name})"
            result['details'] += " without benign variation (PM1)"
            return result
        
        # If no manual flags set, try API-based hotspot detection
        if gene and not self.test_mode:
            # Import here to handle both interactive and test modes
            try:
                from config.constants import COLORAMA_COLORS
            except ImportError:
                # Fallback if running from different context
                COLORAMA_COLORS = {
                    'CYAN': '\033[96m', 'GREEN': '\033[92m', 
                    'YELLOW': '\033[93m', 'RESET': '\033[0m'
                }
            
            try:
                from utils.domain_api_client import DomainAPIClient
                
                # USER FEEDBACK: Inform about API check
                print(f"\n{COLORAMA_COLORS['CYAN']}🔍 Checking hotspot databases for {gene}", end='')
                if position:
                    print(f" position {position}", end='')
                print(f"...{COLORAMA_COLORS['RESET']}")
                
                domain_client = DomainAPIClient(cache_enabled=True)
                hotspot_info = domain_client.get_hotspot_info(gene, position)
                
                # USER FEEDBACK: Show what was found
                source = hotspot_info.get('source', 'none')
                if source != 'none':
                    print(f"{COLORAMA_COLORS['GREEN']}✓ Data retrieved from: {source}{COLORAMA_COLORS['RESET']}")
                    
                    # Show hotspot details if available
                    if 'hotspot_details' in hotspot_info:
                        details = hotspot_info['hotspot_details']
                        print(f"  📊 Tumor count: {details.get('tumor_count', 'N/A')}")
                        print(f"  📊 Mutation count: {details.get('mutation_count', 'N/A')}")
                    
                    # Show domain information
                    domains = hotspot_info.get('domains', [])
                    if domains:
                        print(f"  🧬 Domains found: {', '.join(domains[:3])}")
                        if len(domains) > 3:
                            print(f"     ... and {len(domains)-3} more")
                else:
                    print(f"{COLORAMA_COLORS['YELLOW']}⚠ No external database information available{COLORAMA_COLORS['RESET']}")
                
                # Check if position is a known hotspot
                if position and hotspot_info.get('position_is_hotspot'):
                    result['applies'] = True
                    result['details'] = f"Variant at known hotspot residue {position} in {gene}"
                    if hotspot_info.get('source'):
                        result['details'] += f" (Source: {hotspot_info['source']})"
                    result['details'] += " (PM1)"
                    result['api_source'] = source
                    result['external_data'] = hotspot_info
                    
                    # USER FEEDBACK: Clear success message
                    print(f"{COLORAMA_COLORS['GREEN']}✅ PM1 APPLIES: Position {position} is a documented hotspot{COLORAMA_COLORS['RESET']}")
                    return result
                
                # Check if gene has known hotspots (suggest manual review)
                if hotspot_info.get('is_hotspot_gene'):
                    domains = hotspot_info.get('domains', [])
                    hotspot_residues = hotspot_info.get('hotspot_residues', [])
                    
                    result['manual_review'] = True
                    result['details'] = f"Gene {gene} has known hotspot regions"
                    if domains:
                        result['details'] += f" (Domains: {', '.join(domains[:3])})"
                    if hotspot_residues and position:
                        closest = min(hotspot_residues, key=lambda x: abs(x - position))
                        if abs(closest - position) <= 5:
                            result['details'] += f" - Near hotspot residue {closest}"
                    result['details'] += f" (Source: {hotspot_info['source']})"
                    result['api_source'] = source
                    result['external_data'] = hotspot_info
                    
                    # USER FEEDBACK: Manual review suggestion
                    print(f"{COLORAMA_COLORS['YELLOW']}⚡ Manual review suggested: {gene} has documented hotspot regions{COLORAMA_COLORS['RESET']}")
                    if hotspot_residues:
                        print(f"   Known hotspot residues: {', '.join(map(str, hotspot_residues[:5]))}")
                        if len(hotspot_residues) > 5:
                            print(f"   ... and {len(hotspot_residues)-5} more")
                    return result
            
            except ImportError:
                print(f"{COLORAMA_COLORS['YELLOW']}⚠ Domain API module not available - using fallback{COLORAMA_COLORS['RESET']}")
            except Exception as e:
                # API error - don't fail, just fall back
                print(f"{COLORAMA_COLORS['YELLOW']}⚠ API check failed ({str(e)[:50]}) - using fallback data{COLORAMA_COLORS['RESET']}")
                pass
        
        # Legacy fallback: If flags not set, check for known hotspot genes
        if gene:
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
            else:
                result['details'] = f"No hotspot or functional domain data available"
        else:
            result['details'] = "Insufficient variant information for PM1 evaluation"
        
        return result
    
    def _evaluate_pm2(self, variant_data) -> Dict[str, Any]:
        """Evaluate PM2 - Absent from controls."""
        result = {
            'applies': False, 
            'strength': 'Moderate', 
            'details': '',
            'confidence': 'high',
            'data_source': 'population_database'
        }
        
        # Check population frequencies
        pop_data = variant_data.population_data
        gnomad_af = pop_data.get('gnomad_af')
        
        # CRITICAL FIX: PM2 should NOT apply if frequency is high enough for BA1/BS1
        # BA1 threshold: >5%, BS1 threshold: >1% (gene-specific)
        gene = variant_data.basic_info.get('gene', '').upper()
        gene_thresholds = GENE_SPECIFIC_THRESHOLDS.get(gene, GENE_SPECIFIC_THRESHOLDS['default'])
        bs1_threshold = gene_thresholds.get('BS1', 0.01)
        
        # Check if frequency is too high for PM2
        if gnomad_af is not None and gnomad_af >= bs1_threshold:
            # Frequency too high - variant is NOT absent from controls
            result['applies'] = False
            result['details'] = f"Variant present in population at significant frequency (gnomAD AF: {gnomad_af})"
            result['confidence'] = 'high'
            return result
        
        # Check if absent from controls flag is explicitly set
        if pop_data.get('absent_from_controls') is True:
            result['applies'] = True
            result['details'] = "Variant absent from population databases (PM2)"
            return result
        
        # PM2 applies if variant is absent or extremely rare
        if gnomad_af is None:
            # No frequency data - could be absent or just not in database
            result['applies'] = False
            result['details'] = "No population frequency data available"
            result['confidence'] = 'low'
        elif gnomad_af == 0.0:
            result['applies'] = True
            result['details'] = "Variant absent from population databases (gnomAD AF: 0.0) (PM2)"
        elif gnomad_af < 0.0001:  # Extremely rare (< 0.01%)
            result['applies'] = True
            result['details'] = f"Variant extremely rare in population (gnomAD AF: {gnomad_af}) (PM2)"
        else:
            result['confidence'] = 'medium'
            result['details'] = f"Variant present in population databases (gnomAD AF: {gnomad_af})"
        
        return result
    
    def _evaluate_pm3(self, variant_data) -> Dict[str, Any]:
        """Evaluate PM3 - For recessive disorders, detected in trans with a pathogenic variant."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        gene = variant_data.basic_info.get('gene')
        variant_name = variant_data.basic_info.get('variant_name', 'variant')
        
        # CRITICAL FIX: Check genetic_data for trans phasing information
        genetic_data = variant_data.genetic_data or {}
        phase = genetic_data.get('phase', '').lower()
        inheritance = genetic_data.get('inheritance_pattern', '').lower()
        other_variant = genetic_data.get('other_variant')
        
        # Check if variant is in trans with pathogenic variant in recessive disorder
        if phase == 'trans' and inheritance == 'recessive' and other_variant:
            other_classification = other_variant.get('classification', '').lower()
            if 'pathogenic' in other_classification:
                result['applies'] = True
                result['details'] = f"Variant in trans with pathogenic variant ({other_variant.get('hgvs_c', 'unknown')}) in recessive disorder (PM3)"
                return result
        
        if gene:
            if self.test_mode:
                result['details'] = f"Test mode: No trans pathogenic variant detected for {gene} {variant_name}"
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
        hgvs_p = variant_data.basic_info.get('hgvs_p')
        
        # CRITICAL FIX: Check ClinVar data for same residue different AA
        clinvar = variant_data.clinvar_data or {}
        
        # Check if same residue has different pathogenic AA change
        if clinvar.get('same_residue_different_aa'):
            result['applies'] = True
            result['details'] = f"Different pathogenic missense change at same residue (PM5)"
            return result
        
        # Check pathogenic_at_residue list
        pathogenic_at_residue = clinvar.get('pathogenic_at_residue', [])
        if pathogenic_at_residue and hgvs_p:
            # Extract residue position from hgvs_p (e.g., p.Arg273Cys -> 273)
            import re
            match = re.search(r'p\.\w+(\d+)\w+', hgvs_p)
            if match:
                current_position = match.group(1)
                
                # Check if any pathogenic variant at same position but different AA
                for pv in pathogenic_at_residue:
                    pv_hgvs = pv.get('hgvs_p', '')
                    if current_position in pv_hgvs and pv_hgvs != hgvs_p:
                        result['applies'] = True
                        result['details'] = f"Different pathogenic change at same residue: {pv_hgvs} (PM5)"
                        return result
        
        # If no automated data, use interactive mode
        if gene and (aa_change or hgvs_p):
            if self.test_mode:
                result['details'] = f"Test mode: No ClinVar data for same residue different AA"
                result['manual_review'] = True
            else:
                result = self._evaluate_pm5_interactive(variant_data, gene, aa_change or hgvs_p)
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
        
        genetic_data = variant_data.genetic_data or {}
        
        # CRITICAL FIX: Check multiple fields for assumed de novo
        inheritance_pattern = genetic_data.get('inheritance_pattern', '').lower()
        assumed_de_novo = genetic_data.get('assumed_de_novo', False)
        
        # Get parental confirmation status
        maternity_confirmed = genetic_data.get('maternity_confirmed', False)
        paternity_confirmed = genetic_data.get('paternity_confirmed', False)
        family_history = genetic_data.get('family_history', True)
        
        # Try legacy de novo status fields
        denovo_status = (genetic_data.get('de_novo') or 
                        genetic_data.get('denovo') or
                        variant_data.functional_data.get('de_novo') or
                        variant_data.functional_data.get('denovo'))
        
        # Check if de novo but not confirmed with BOTH parents
        is_de_novo = (inheritance_pattern == 'de_novo' or 
                     assumed_de_novo or
                     denovo_status in ['assumed', 'confirmed'])
        
        if not is_de_novo:
            result['details'] = 'No evidence of de novo variant'
            return result
        
        # PM6 applies if de novo but NOT confirmed with both parents
        both_parents_confirmed = maternity_confirmed and paternity_confirmed
        
        if assumed_de_novo or denovo_status == 'assumed':
            # Explicitly marked as assumed
            result['applies'] = True
            result['details'] = 'De novo variant assumed but not confirmed with parental testing (PM6)'
        elif is_de_novo and not both_parents_confirmed:
            # De novo indicated but parental testing incomplete
            result['applies'] = True
            if maternity_confirmed or paternity_confirmed:
                result['details'] = 'De novo with only one parent confirmed (PM6)'
            else:
                result['details'] = 'De novo assumed without parental testing (PM6)'
        else:
            # Both parents confirmed - should use PS2 instead
            result['details'] = 'De novo confirmed with both parents - see PS2 instead of PM6'
        
        return result
    
    def _evaluate_pp1(self, variant_data) -> Dict[str, Any]:
        """Evaluate PP1 - Cosegregation with disease with automated LOD scoring."""
        result = {
            'applies': False, 
            'strength': 'Supporting', 
            'details': '',
            'confidence': 'very_low',
            'data_source': 'none'
        }
        
        # CRITICAL FIX: Check segregation_data dictionary
        genetic_data = variant_data.genetic_data or {}
        segregation_data = genetic_data.get('segregation_data')
        
        # First check if segregation_data dictionary exists
        if segregation_data and isinstance(segregation_data, dict):
            segregates = segregation_data.get('segregates')
            lod_score = segregation_data.get('lod_score')
            
            # If explicit segregation status is provided
            if segregates is True and lod_score:
                result['applies'] = True
                result['data_source'] = 'segregation_analysis'
                
                # ACMG 2023: LOD-based strength modifiers
                if self.use_2023_guidelines:
                    if lod_score >= 5.0:
                        result['strength'] = 'Strong'
                        result['confidence'] = 'high'
                        result['details'] = f"Variant cosegregates with disease (LOD={lod_score:.2f}) - Strong segregation evidence (PP1_Strong per ACMG 2023)"
                    elif lod_score >= 3.0:
                        result['strength'] = 'Moderate'
                        result['confidence'] = 'high'
                        result['details'] = f"Variant cosegregates with disease (LOD={lod_score:.2f}) - Moderate segregation evidence (PP1_Moderate per ACMG 2023)"
                    else:  # 1.5 <= LOD < 3.0
                        result['strength'] = 'Supporting'
                        result['confidence'] = 'medium'
                        result['details'] = f"Variant cosegregates with disease (LOD={lod_score:.2f}) (PP1)"
                else:
                    # ACMG 2015: Simple strength assignment
                    result['confidence'] = 'high' if lod_score >= 3.0 else 'medium'
                    result['details'] = f"Variant cosegregates with disease (LOD={lod_score}) (PP1)"
                    if lod_score >= 3.0:
                        result['strength'] = 'Moderate'
                        result['details'] += " (strong segregation evidence)"
                
                return result
            elif segregates is False:
                result['details'] = "Variant does not cosegregate with disease - see BS4"
                return result
        
        # Check if structured segregation families data is available
        segregation_families = genetic_data.get('segregation_families')
        
        if segregation_families and isinstance(segregation_families, list):
            # Automated LOD score calculation
            lod_result = self.statistical_analyzer.calculate_lod_score(segregation_families)
            
            if lod_result['valid']:
                result['data_source'] = 'segregation_analysis'
                result['confidence'] = lod_result['confidence']
                result['lod_analysis'] = lod_result
                
                if lod_result['applies'] and lod_result['strength'] in ['supporting', 'strong']:
                    result['applies'] = True
                    result['details'] = lod_result['interpretation']
                    
                    # Potentially upgrade to PP1_Strong if LOD is very high
                    if lod_result['lod_score'] >= 3.0:
                        result['strength'] = 'Moderate'  # Could consider Strong with more families
                        result['details'] += " (strong segregation evidence)"
                else:
                    result['details'] = f"Segregation data insufficient: {lod_result['interpretation']}"
                
                # Store for report
                self.evidence_details['PP1_lod_analysis'] = lod_result
                
                return result
            else:
                result['details'] = f"Segregation analysis failed: {lod_result.get('error', 'Unknown error')}"
                return result
        
        # Fall back to simple segregation status check
        segregation = genetic_data.get('segregation')
        
        if segregation == 'cosegregates':
            result['applies'] = True
            result['confidence'] = 'low'
            result['data_source'] = 'qualitative_report'
            result['details'] = "Variant cosegregates with disease in family (PP1) - recommend quantitative LOD analysis"
        elif segregation == 'does_not_segregate':
            result['confidence'] = 'low'
            result['data_source'] = 'qualitative_report'
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
                "Calculate LOD score with structured family data",
                "Minimum 3 families recommended for statistical power"
            ]
            result['guidance'] = "PP1 applies if variant cosegregates with disease (LOD ≥1.5 for Supporting, ≥3.0 for Moderate/Strong)"
        
        return result
    
    def _evaluate_pp2(self, variant_data) -> Dict[str, Any]:
        """Evaluate PP2 - Missense variant in a gene that has a low rate of benign missense variation."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        gene = variant_data.basic_info.get('gene')
        variant_name = variant_data.basic_info.get('variant_name', 'variant')
        
        if variant_type == 'missense':
            # CRITICAL: Do NOT apply PP2 if population frequency suggests benign
            population_data = variant_data.population_data or {}
            gnomad_af = population_data.get('gnomad_af')
            
            # If allele frequency is high (>0.001 = 0.1%), this is likely benign missense
            if gnomad_af is not None and gnomad_af > 0.001:  # 0.1% threshold
                result['details'] = f"Not applying PP2: variant frequency too high ({gnomad_af:.4f}) for pathogenic missense"
                return result
            
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
                # In test mode, don't auto-apply PP2 (too many false positives)
                if self.test_mode:
                    result['details'] = f"Test mode: PP2 requires manual review for {gene} {variant_name}"
                    result['manual_review'] = True
                    result['guidance'] = "PP2 applies if the gene has low benign missense rate. Use cautiously."
                else:
                    # Interactive mode: Ask user
                    result = self._evaluate_pp2_interactive(variant_data, gene, variant_name)
            else:
                result['details'] = f"Gene {gene} not in low-benign-missense gene list"
                result['manual_review'] = True
        else:
            result['details'] = "Not a missense variant"
        
        return result
    
    def _evaluate_pp2_interactive(self, variant_data, gene, variant_name) -> Dict[str, Any]:
        """Interactive PP2 evaluation with user input."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        print(f"\n🔍 PP2 Evaluation: {gene} {variant_name}")
        print("─" * 50)
        print(f"QUESTION: Does {gene} have a LOW rate of benign missense variation?")
        print()
        print("📋 Check these resources:")
        print(f"   • gnomAD constraint metrics (missense Z-score)")
        print(f"   • ClinVar pathogenic vs benign missense ratio")
        print(f"   • Gene-specific guidelines")
        print(f"   • OMIM and literature")
        print()
        print("⚠️  WARNING: PP2 is commonly over-applied. Use only if:")
        print("   1. Gene has documented low benign missense variation")
        print("   2. Missense variants are a COMMON disease mechanism")
        print("   3. Gene has high missense constraint (Z-score > 2)")
        print()
        
        while True:
            choice = input("Apply PP2? (y/n/u for unknown): ").strip().lower()
            if choice in ['y', 'yes']:
                result['applies'] = True
                result['details'] = f"Missense in {gene} (low benign missense rate) (PP2)"
                print(f"✅ PP2 applies: {result['details']}")
                break
            elif choice in ['n', 'no']:
                result['details'] = f"PP2 not applicable for {gene}"
                print(f"❌ PP2 does not apply: {result['details']}")
                break
            elif choice in ['u', 'unknown']:
                result['details'] = f"PP2 unclear - requires gene-specific research for {gene}"
                print(f"⚠️  PP2 unclear: {result['details']}")
                break
            else:
                print("❌ Please enter 'y' for yes, 'n' for no, or 'u' for unknown")
        
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
        # Automated assignment for PP4 using PhenotypeMatcher or genetic_data
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        gene = variant_data.basic_info.get('gene')
        patient_phenotypes = getattr(variant_data, 'patient_phenotypes', None)
        
        # CRITICAL FIX: Check genetic_data for phenotype_specificity
        genetic_data = variant_data.genetic_data or {}
        phenotype_specificity = genetic_data.get('phenotype_specificity')
        
        # If explicit phenotype specificity is provided
        if phenotype_specificity == 'high' and patient_phenotypes:
            result['applies'] = True
            result['details'] = f"Patient phenotype highly specific for {gene}-related disease (PP4)"
            return result
        
        # Use PhenotypeMatcher for automated evaluation
        if patient_phenotypes and gene:
            pp4_bp5 = self.phenotype_matcher.evaluate_phenotype_match(variant_data, patient_phenotypes)
            if pp4_bp5 == 'PP4':
                result['applies'] = True
                result['details'] = f"Patient phenotype highly specific for {gene}-related disease (PP4)"
            elif pp4_bp5 == 'BP5':
                result['details'] = f"Patient phenotype not consistent with {gene}-related disease (BP5)"
            else:
                result['details'] = "No or inconclusive phenotype data"
        else:
            result['details'] = "No phenotype data or gene information available"
        
        return result
    
    def _evaluate_pp5(self, variant_data) -> Dict[str, Any]:
        """Evaluate PP5 - Reputable source with stricter validation requirements."""
        from config.constants import REPUTABLE_SOURCE_REQUIREMENTS
        from datetime import datetime
        
        result = {
            'applies': False, 
            'strength': 'Supporting', 
            'details': '',
            'confidence': 'very_low',
            'data_source': 'none'
        }
        
        # Check if ClinVar data is available
        clinvar_data = getattr(variant_data, 'clinvar_data', None)
        
        if clinvar_data and clinvar_data.get('clinical_significance'):
            significance = clinvar_data['clinical_significance'].lower()
            review_status = clinvar_data.get('review_status', '').lower()
            submitter = clinvar_data.get('submitter', '')
            last_evaluated = clinvar_data.get('last_evaluated')
            
            # Validate reputable source requirements
            is_expert_panel = any(panel.lower() in submitter.lower() 
                                 for panel in REPUTABLE_SOURCE_REQUIREMENTS['expert_panels'])
            
            # Check review status (stars)
            star_mapping = {
                'practice guideline': 4,
                'reviewed by expert panel': 3,
                'criteria provided, multiple submitters': 2,
                'criteria provided, single submitter': 1,
                'no assertion criteria provided': 0
            }
            
            stars = 0
            for status, star_count in star_mapping.items():
                if status in review_status:
                    stars = star_count
                    break
            
            # Check age of classification
            classification_recent = True
            if last_evaluated:
                try:
                    eval_date = datetime.fromisoformat(last_evaluated.replace('Z', '+00:00'))
                    age_years = (datetime.now() - eval_date).days / 365.25
                    classification_recent = age_years <= REPUTABLE_SOURCE_REQUIREMENTS['max_age_years']
                except:
                    classification_recent = False  # Assume old if can't parse
            
            # Apply strict criteria
            meets_requirements = (
                (is_expert_panel or stars >= REPUTABLE_SOURCE_REQUIREMENTS['min_stars']) and
                classification_recent
            )
            
            if 'pathogenic' in significance and 'benign' not in significance:
                if meets_requirements:
                    result['applies'] = True
                    result['confidence'] = 'high' if is_expert_panel else 'medium'
                    result['data_source'] = 'expert_panel' if is_expert_panel else 'clinvar_multi_submitter'
                    result['details'] = (
                        f"Reputable source ({submitter}, {stars}⭐) reports variant as pathogenic "
                        f"(PP5, reviewed: {last_evaluated or 'unknown'})"
                    )
                else:
                    result['confidence'] = 'low'
                    result['data_source'] = 'clinvar_insufficient_review'
                    result['details'] = (
                        f"ClinVar reports pathogenic but does not meet PP5 criteria "
                        f"(stars: {stars}, expert panel: {is_expert_panel}, recent: {classification_recent})"
                    )
            elif 'benign' in significance:
                result['details'] = f"Reputable source reports variant as benign - see BP6"
            else:
                result['details'] = f"ClinVar classification: {significance}"
        else:
            result['details'] = "Reputable source analysis required for PP5 evaluation"
            result['manual_review'] = True
            result['search_recommendations'] = [
                "Check ClinVar for expert panel classifications",
                "Review recent literature reports (within 5 years)",
                "Check laboratory databases (HGMD, LOVD)",
                "Verify source credibility (expert panels preferred)",
                "Confirm classification recency"
            ]
            result['guidance'] = (
                "PP5 applies if reputable source (expert panel or ≥2⭐ ClinVar with recent review) "
                "reports variant as pathogenic"
            )
        
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
        
        # CRITICAL FIX: Check if expected_max_af is provided (from test data or disease prevalence)
        expected_max_af = pop_data.get('expected_max_af')
        disease_prevalence = pop_data.get('disease_prevalence')
        
        # Calculate expected max AF from disease prevalence if available
        if disease_prevalence is not None and expected_max_af is None:
            # For dominant: expected_max_af ≈ disease prevalence
            # For recessive: expected_max_af ≈ sqrt(disease prevalence) * 2
            # Use conservative dominant calculation as default
            expected_max_af = disease_prevalence
        
        # Use provided expected_max_af or fall back to gene-specific threshold
        if expected_max_af is not None:
            bs1_threshold = expected_max_af
        else:
            # Get gene-specific threshold or use default
            gene_thresholds = GENE_SPECIFIC_THRESHOLDS.get(gene, GENE_SPECIFIC_THRESHOLDS['default'])
            bs1_threshold = gene_thresholds.get('BS1', 0.01)
        
        # BS1 applies if variant frequency is higher than expected for disorder
        # Must be significantly higher (not just slightly above)
        if gnomad_af is not None:
            # Check if BA1 threshold is already exceeded (>5%)
            ba1_threshold = gene_thresholds.get('BA1', 0.05) if expected_max_af is None else 0.05
            
            if gnomad_af > ba1_threshold:
                # BA1 takes precedence - don't apply BS1
                result['applies'] = False
                result['details'] = f"Variant frequency exceeds BA1 threshold (gnomAD AF: {gnomad_af})"
            elif gnomad_af > bs1_threshold:
                result['applies'] = True
                result['details'] = f"Variant frequency higher than expected for disorder (gnomAD AF: {gnomad_af}, expected max: {bs1_threshold})"
            else:
                result['details'] = f"Variant frequency not higher than expected (gnomAD AF: {gnomad_af}, expected max: {bs1_threshold})"
        else:
            result['details'] = f"No population frequency data available"
        
        return result
    
    def _evaluate_bs2(self, variant_data) -> Dict[str, Any]:
        """Evaluate BS2 - Observed in a healthy adult individual for a recessive disorder."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        # CRITICAL FIX: Check genetic_data for observed_in_healthy
        genetic_data = variant_data.genetic_data or {}
        observed_in_healthy = genetic_data.get('observed_in_healthy', False)
        zygosity = genetic_data.get('zygosity', '')
        age_observed = genetic_data.get('age_observed')
        disease_onset_age = genetic_data.get('disease_onset_age')
        penetrance = genetic_data.get('penetrance', 1.0)
        
        # BS2 requires: healthy individual + appropriate zygosity + age considerations
        if not observed_in_healthy:
            result['details'] = "Variant not observed in healthy individuals"
            return result
        
        # Check if observed at appropriate age (must be past expected disease onset)
        age_appropriate = True
        if age_observed and disease_onset_age:
            age_appropriate = age_observed > disease_onset_age
        
        # Check zygosity matches disorder type
        # Recessive: homozygous, Dominant: heterozygous, X-linked: hemizygous
        valid_zygosity = zygosity.lower() in ['homozygous', 'heterozygous', 'hemizygous']
        
        if observed_in_healthy and age_appropriate and valid_zygosity:
            result['applies'] = True
            detail_parts = [f"Variant observed in healthy individual ({zygosity})"]
            if age_observed:
                detail_parts.append(f"at age {age_observed}")
            if disease_onset_age:
                detail_parts.append(f"(expected onset: {disease_onset_age})")
            result['details'] = " ".join(detail_parts) + " (BS2)"
        elif observed_in_healthy and not age_appropriate:
            result['details'] = f"Observed in healthy individual but age ({age_observed}) not past expected onset ({disease_onset_age})"
        elif observed_in_healthy and not valid_zygosity:
            result['details'] = f"Observed in healthy individual but zygosity ({zygosity}) unclear"
        else:
            result['details'] = "Case-control analysis required for BS2 evaluation"
            result['manual_review'] = True
        
        return result
    
    def _evaluate_bs3(self, variant_data) -> Dict[str, Any]:
        """Evaluate BS3 - Well-established functional studies show no damaging effect."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        # CRITICAL FIX: Check functional_data first
        functional_data = variant_data.functional_data or {}
        has_studies = functional_data.get('has_functional_studies', False)
        functional_impact = functional_data.get('functional_impact', '').lower()
        functional_studies = functional_data.get('functional_studies', [])
        
        # Check if functional_impact explicitly indicates benign
        if has_studies and functional_impact == 'benign':
            result['applies'] = True
            details = "Functional studies show no damaging effect"
            
            # Add study details if available
            if functional_studies:
                study_count = len(functional_studies)
                details += f" ({study_count} {'study' if study_count == 1 else 'studies'})"
                # Add first study description
                if functional_studies[0].get('description'):
                    details += f": {functional_studies[0]['description']}"
            
            result['details'] = details + " (BS3)"
            return result
        
        # Fall back to FunctionalStudiesEvaluator for complex cases
        ps3_bs3 = self.functional_studies_evaluator.evaluate_functional_evidence(variant_data)
        if ps3_bs3 == 'BS3':
            result['applies'] = True
            result['details'] = "Functional studies support a benign effect (BS3)"
        elif has_studies:
            result['details'] = f"Functional studies available but impact unclear (reported: {functional_impact})"
        else:
            result['details'] = "No functional studies data available"
        
        return result
    
    def _evaluate_bs4(self, variant_data) -> Dict[str, Any]:
        """Evaluate BS4 - Lack of segregation with automated LOD analysis."""
        result = {
            'applies': False, 
            'strength': 'Strong', 
            'details': '',
            'confidence': 'very_low',
            'data_source': 'none'
        }
        
        # CRITICAL FIX: Check segregation_data dictionary
        genetic_data = variant_data.genetic_data or {}
        segregation_data = genetic_data.get('segregation_data')
        
        # First check if segregation_data dictionary exists
        if segregation_data and isinstance(segregation_data, dict):
            segregates = segregation_data.get('segregates')
            lod_score = segregation_data.get('lod_score')
            
            # If explicit non-segregation is provided
            if segregates is False and lod_score:
                result['applies'] = True
                result['confidence'] = 'high' if lod_score <= -2.0 else 'medium'
                result['data_source'] = 'segregation_analysis'
                result['details'] = f"Variant does not segregate with disease (LOD={lod_score}) (BS4)"
                return result
            elif segregates is True:
                result['details'] = "Variant cosegregates with disease - see PP1"
                return result
        
        # Check if structured segregation families data is available
        segregation_families = genetic_data.get('segregation_families')
        
        if segregation_families and isinstance(segregation_families, list):
            # Automated LOD score calculation
            lod_result = self.statistical_analyzer.calculate_lod_score(segregation_families)
            
            if lod_result['valid']:
                result['data_source'] = 'segregation_analysis'
                result['confidence'] = lod_result['confidence']
                result['lod_analysis'] = lod_result
                
                # BS4 applies for negative LOD (non-segregation)
                if lod_result['applies'] and lod_result['strength'] == 'benign_strong':
                    result['applies'] = True
                    result['details'] = lod_result['interpretation'] + " (BS4)"
                else:
                    result['details'] = f"Segregation analysis: {lod_result['interpretation']}"
                
                # Store for report
                self.evidence_details['BS4_lod_analysis'] = lod_result
                
                return result
            else:
                result['details'] = f"Segregation analysis failed: {lod_result.get('error', 'Unknown error')}"
                return result
        
        # Fall back to simple segregation status check
        segregation = genetic_data.get('segregation')
        
        if segregation == 'does_not_segregate':
            result['applies'] = True
            result['confidence'] = 'low'
            result['data_source'] = 'qualitative_report'
            result['details'] = "Variant does not segregate with disease in family (BS4) - recommend quantitative LOD analysis"
        elif segregation == 'cosegregates':
            result['confidence'] = 'low'
            result['data_source'] = 'qualitative_report'
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
                "Calculate LOD score with structured family data",
                "Review family structure and inheritance pattern"
            ]
            result['guidance'] = "BS4 applies if variant does not segregate with disease (negative LOD score)"
        
        return result
    
    def _evaluate_bp1(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP1 - Missense variant in a gene for which primarily truncating variants cause disease."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        gene = variant_data.basic_info.get('gene')
        
        if variant_type == 'missense':
            # Define genes where primarily truncating variants cause disease
            # NOTE: This list is conservative - only includes genes where truncating is PRIMARY mechanism
            truncating_mechanism_genes = [
                'DMD',  # Duchenne muscular dystrophy
                'NF1', 'NF2',  # Neurofibromatosis
                'TSC1', 'TSC2',  # Tuberous sclerosis
                'APC',  # Familial adenomatous polyposis
                'MLH1', 'MSH2', 'MSH6', 'PMS2', 'EPCAM',  # Lynch syndrome
                'RB1',  # Retinoblastoma
                'WT1',  # Wilms tumor
                'VHL',  # Von Hippel-Lindau
                'SMAD4', 'BMPR1A',  # Juvenile polyposis
                'STK11',  # Peutz-Jeghers
                'CDKN2A'  # Melanoma
                # NOTE: BRCA1, BRCA2, CDH1, PTEN removed - they have significant missense pathogenic variants
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
        
        # CRITICAL FIX: Check genetic_data for cis/trans phasing information
        genetic_data = variant_data.genetic_data or {}
        phase = genetic_data.get('phase', '').lower()
        inheritance = genetic_data.get('inheritance_pattern', '').lower()
        other_variant = genetic_data.get('other_variant')
        patient_unaffected = genetic_data.get('patient_unaffected', False)
        
        # BP2 applies if: (1) variant in cis with pathogenic OR (2) in trans with pathogenic in dominant disorder
        # NOTE: Test data shows cis case (same chromosome, patient unaffected)
        if other_variant:
            other_classification = other_variant.get('classification', '').lower()
            if 'pathogenic' in other_classification:
                # Case 1: Cis with pathogenic (patient should be unaffected or only mildly affected)
                if phase == 'cis' and patient_unaffected:
                    result['applies'] = True
                    result['details'] = f"Variant in cis with pathogenic variant ({other_variant.get('hgvs_c', 'unknown')}) and patient unaffected (BP2)"
                    return result
                # Case 2: Trans with pathogenic in fully penetrant dominant disorder
                elif phase == 'trans' and inheritance == 'dominant' and patient_unaffected:
                    result['applies'] = True
                    result['details'] = f"Variant in trans with pathogenic variant in dominant disorder and patient unaffected (BP2)"
                    return result
        
        if gene:
            if self.test_mode:
                result['details'] = f"Test mode: No cis/trans pathogenic variant configuration detected for {gene} {variant_name}"
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
        if variant_type == 'inframe_indel' or 'inframe' in consequence or variant_type in ['inframe_insertion', 'inframe_deletion']:
            # CRITICAL FIX: Check functional_data for repeat region information
            functional_data = variant_data.functional_data or {}
            in_repeat = functional_data.get('in_repeat_region', False)
            repeat_functional = functional_data.get('repeat_functional', True)  # Default: assume functional
            
            # BP3 applies if in repeat region WITHOUT known function
            if in_repeat and not repeat_functional:
                result['applies'] = True
                repeat_type = functional_data.get('repeat_type', 'unknown')
                result['details'] = f"In-frame indel in repetitive region ({repeat_type}) without known function (BP3)"
                return result
            elif in_repeat and repeat_functional:
                result['details'] = "In-frame indel in repetitive region but has known function"
                return result
            
            # This would ideally check if the variant is in a repetitive region
            # For now, we'll apply it to in-frame indels and ask for user confirmation
            if self.test_mode:
                result['details'] = "Test mode: In-frame indel - no repeat region data available"
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
        
        # CRITICAL FIX: Check both prediction labels (_pred) and scores
        # Check prediction labels first (more explicit)
        pred_mapping = {
            'sift_pred': {'benign': ['T', 'TOLERATED'], 'damaging': ['D', 'DELETERIOUS']},
            'polyphen2_hvar_pred': {'benign': ['B', 'BENIGN'], 'damaging': ['D', 'PROBABLY_DAMAGING', 'P']},
            'polyphen2_hdiv_pred': {'benign': ['B', 'BENIGN'], 'damaging': ['D', 'PROBABLY_DAMAGING', 'P']},
            'mutation_taster_pred': {'benign': ['N', 'P', 'POLYMORPHISM'], 'damaging': ['D', 'A', 'DISEASE_CAUSING']},
        }
        
        for pred_key, values in pred_mapping.items():
            pred_value = insilico_data.get(pred_key)
            if pred_value:
                total_count += 1
                if pred_value.upper() in [v.upper() for v in values['benign']]:
                    benign_count += 1
        
        # Check scores
        predictor_thresholds = {
            'cadd_phred': 20,
            'revel_score': 0.5,
            'revel': 0.5,
            'sift_score': 0.05,  # Higher is more benign (>0.05 = tolerated)
            'sift': 0.05,
            'polyphen2_score': 0.5,  # Lower is more benign (<0.5 = benign)
            'polyphen2': 0.5,
            'mutation_taster': 0.5,
            'dann_score': 0.5,
            'dann': 0.5,
            'fathmm_score': 0.5,  # Higher is more benign
            'fathmm': 0.5,
        }
        
        for predictor, threshold in predictor_thresholds.items():
            score = insilico_data.get(predictor)
            if score is not None:
                total_count += 1
                if predictor in ['sift', 'sift_score', 'fathmm', 'fathmm_score']:
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
        
        # CRITICAL FIX: Check genetic_data for alternate cause information
        genetic_data = variant_data.genetic_data or {}
        alternate_found = genetic_data.get('alternate_cause_found', False)
        alternate_variant = genetic_data.get('alternate_variant')
        
        # If alternate cause is explicitly documented
        if alternate_found and alternate_variant:
            alternate_gene = alternate_variant.get('gene', 'unknown')
            alternate_classification = alternate_variant.get('classification', 'unknown')
            explains = alternate_variant.get('explains_phenotype', False)
            
            if 'pathogenic' in alternate_classification.lower() and explains:
                result['applies'] = True
                result['details'] = f"Alternate molecular basis identified ({alternate_gene} {alternate_classification}) explaining phenotype (BP5)"
                return result
        
        if gene:
            if self.test_mode:
                result['details'] = f"Test mode: No alternate molecular basis documented for {gene} {variant_name}"
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
        """Evaluate BP6 - Reputable source reports variant as benign."""
        from config.constants import REPUTABLE_SOURCE_REQUIREMENTS
        from datetime import datetime
        
        result = {
            'applies': False, 
            'strength': 'Supporting', 
            'details': '',
            'confidence': 'very_low',
            'data_source': 'none'
        }
        
        # CRITICAL FIX: Check clinvar_data for benign classification
        clinvar_data = variant_data.clinvar_data or {}
        classification = clinvar_data.get('classification', '').lower()
        
        # Quick check: if explicitly marked as Benign/Likely benign
        if classification in ['benign', 'likely benign', 'likely_benign']:
            result['applies'] = True
            result['data_source'] = 'clinvar'
            result['confidence'] = 'high'
            result['details'] = f"Reputable source reports variant as {classification.title()} (BP6)"
            return result
        
        # Fall back to detailed ClinVar validation
        clinvar_data = getattr(variant_data, 'clinvar_data', None)
        
        if clinvar_data and clinvar_data.get('clinical_significance'):
            significance = clinvar_data['clinical_significance'].lower()
            review_status = clinvar_data.get('review_status', '').lower()
            submitter = clinvar_data.get('submitter', '')
            last_evaluated = clinvar_data.get('last_evaluated')
            
            # Validate reputable source requirements
            is_expert_panel = any(panel.lower() in submitter.lower() 
                                 for panel in REPUTABLE_SOURCE_REQUIREMENTS['expert_panels'])
            
            # Check review status (stars)
            star_mapping = {
                'practice guideline': 4,
                'reviewed by expert panel': 3,
                'criteria provided, multiple submitters': 2,
                'criteria provided, single submitter': 1,
                'no assertion criteria provided': 0
            }
            
            stars = 0
            for status, star_count in star_mapping.items():
                if status in review_status:
                    stars = star_count
                    break
            
            # Check age of classification
            classification_recent = True
            if last_evaluated:
                try:
                    eval_date = datetime.fromisoformat(last_evaluated.replace('Z', '+00:00'))
                    age_years = (datetime.now() - eval_date).days / 365.25
                    classification_recent = age_years <= REPUTABLE_SOURCE_REQUIREMENTS['max_age_years']
                except:
                    classification_recent = False
            
            # Apply strict criteria
            meets_requirements = (
                (is_expert_panel or stars >= REPUTABLE_SOURCE_REQUIREMENTS['min_stars']) and
                classification_recent
            )
            
            if 'benign' in significance and 'pathogenic' not in significance:
                if meets_requirements:
                    result['applies'] = True
                    result['confidence'] = 'high' if is_expert_panel else 'medium'
                    result['data_source'] = 'expert_panel' if is_expert_panel else 'clinvar_multi_submitter'
                    result['details'] = (
                        f"Reputable source ({submitter}, {stars}⭐) reports variant as benign "
                        f"(BP6, reviewed: {last_evaluated or 'unknown'})"
                    )
                else:
                    result['confidence'] = 'low'
                    result['data_source'] = 'clinvar_insufficient_review'
                    result['details'] = (
                        f"ClinVar reports benign but does not meet BP6 criteria "
                        f"(stars: {stars}, expert panel: {is_expert_panel}, recent: {classification_recent})"
                    )
            elif 'pathogenic' in significance:
                result['details'] = f"Reputable source reports variant as pathogenic - see PP5"
            else:
                result['details'] = f"ClinVar classification: {significance}"
        else:
            # Fall back to interactive if no ClinVar data
            gene = variant_data.basic_info.get('gene')
            variant_name = variant_data.basic_info.get('variant_name', 'variant')
            
            if gene and not self.test_mode:
                result = self._evaluate_bp6_interactive(variant_data, gene, variant_name)
            else:
                result['details'] = "Reputable source analysis required for BP6 evaluation"
                result['manual_review'] = True
                result['search_recommendations'] = [
                    "Check ClinVar for expert panel benign classifications",
                    "Review certified lab reports (≥2⭐)",
                    "Check published clinical guidelines (within 5 years)",
                    "Verify source credibility (expert panels preferred)"
                ]
                result['guidance'] = (
                    "BP6 applies if reputable source (expert panel or ≥2⭐ ClinVar with recent review) "
                    "reports variant as benign"
                )
        
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
