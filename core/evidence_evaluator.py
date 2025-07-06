"""
Evidence Evaluator Module
=========================

This module evaluates ACMG evidence criteria based on variant data.
It implements both ACMG/AMP 2015 and 2023 guidelines.
"""

import numpy as np
from typing import Dict, List, Optional, Any, Tuple
from scipy import stats
from config.constants import (
    EVIDENCE_WEIGHTS, GENE_SPECIFIC_THRESHOLDS, INSILICO_WEIGHTS,
    INSILICO_THRESHOLDS, VARIANT_CONSEQUENCES, STATISTICAL_THRESHOLDS,
    VAMPP_SCORE_THRESHOLDS
)


class EvidenceEvaluator:
    """
    Evaluates ACMG evidence criteria for variant classification.
    
    This class implements the logic for evaluating all ACMG evidence criteria
    based on variant data, including population frequencies, in silico predictions,
    functional studies, and inheritance patterns.
    """
    
    def __init__(self, use_2023_guidelines: bool = False):
        """
        Initialize the evidence evaluator.
        
        Args:
            use_2023_guidelines (bool): Whether to use ACMG 2023 guidelines
        """
        self.use_2023_guidelines = use_2023_guidelines
        self.applied_criteria = {}
        self.evidence_details = {}
    
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
        
        # BP4 - Computational evidence
        benign['BP4'] = self._evaluate_bp4(variant_data)
        
        # BP5 - Alternate molecular basis
        benign['BP5'] = self._evaluate_bp5(variant_data)
        
        # BP6 - Reputable source (if enabled)
        if self.use_2023_guidelines:
            benign['BP6'] = self._evaluate_bp6(variant_data)
        
        # BP7 - Synonymous variant
        benign['BP7'] = self._evaluate_bp7(variant_data)
        
        return benign
    
    def _evaluate_pvs1(self, variant_data) -> Dict[str, Any]:
        """Evaluate PVS1 - Null variant in LOF gene."""
        result = {'applies': False, 'strength': 'Very Strong', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        consequence = variant_data.basic_info.get('consequence', '').lower()
        
        # Check if it's a loss-of-function variant
        lof_types = ['nonsense', 'frameshift', 'splice_donor', 'splice_acceptor', 'start_lost']
        lof_consequences = ['stop_gained', 'frameshift_variant', 'splice_donor_variant', 
                           'splice_acceptor_variant', 'start_lost']
        
        if variant_type in lof_types or consequence in lof_consequences:
            # TODO: Add gene-specific LOF mechanism check
            # For now, assume LOF is a known mechanism for the gene
            result['applies'] = True
            result['details'] = f"Loss-of-function variant ({variant_type or consequence})"
        
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
                result['applies'] = True
                result['strength'] = 'Strong'  # Use PS1 strength for splice-altering variants
                result['details'] = f"Intronic variant with high SpliceAI score (max: {max(spliceai_scores):.3f})"
            elif spliceai_scores and max(spliceai_scores) > 0.2:
                result['applies'] = False
                result['details'] = f"Intronic variant with moderate SpliceAI score (max: {max(spliceai_scores):.3f}) - consider PS1 or PM4"
            else:
                result['details'] = "Intronic variant with low splice impact prediction"
        
        return result
    
    def _evaluate_ps1(self, variant_data) -> Dict[str, Any]:
        """Evaluate PS1 - Same amino acid change as pathogenic variant."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        # This would require database lookup for known pathogenic variants
        # at the same amino acid position
        result['details'] = "Requires manual verification of pathogenic variants at same position"
        
        return result
    
    def _evaluate_ps2(self, variant_data) -> Dict[str, Any]:
        """Evaluate PS2 - De novo variant."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        # Try both genetic_data and functional_data for de novo status
        denovo_status = (variant_data.genetic_data.get('de_novo') or 
                        variant_data.genetic_data.get('denovo') or
                        variant_data.functional_data.get('de_novo') or
                        variant_data.functional_data.get('denovo'))
        
        if denovo_status == 'confirmed' or denovo_status == 'yes':
            if self.use_2023_guidelines:
                result['applies'] = True
                result['strength'] = 'Very Strong'
                result['details'] = "Confirmed de novo variant (2023 guidelines - Very Strong)"
            else:
                result['applies'] = True
                result['details'] = "Confirmed de novo variant"
        elif denovo_status == 'assumed':
            # This would be PM6 instead
            result['details'] = "Assumed de novo - see PM6"
        
        return result
    
    def _evaluate_ps3(self, variant_data) -> Dict[str, Any]:
        """Evaluate PS3 - Functional studies supportive of damaging effect."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        functional_studies = variant_data.functional_data.get('functional_studies')
        
        if functional_studies == 'damaging':
            result['applies'] = True
            result['details'] = "Functional studies show damaging effect"
        elif functional_studies == 'benign':
            result['details'] = "Functional studies show benign effect - see BS3"
        elif functional_studies == 'inconclusive':
            result['details'] = "Functional studies inconclusive"
        
        return result
    
    def _evaluate_ps4(self, variant_data) -> Dict[str, Any]:
        """Evaluate PS4 - Case-control data."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        if variant_data.functional_data.get('case_control') == 'yes':
            # Perform Fisher's exact test
            cases_with = variant_data.functional_data.get('cases_with_variant', 0)
            total_cases = variant_data.functional_data.get('total_cases', 0)
            controls_with = variant_data.functional_data.get('controls_with_variant', 0)
            total_controls = variant_data.functional_data.get('total_controls', 0)
            
            if all(x is not None for x in [cases_with, total_cases, controls_with, total_controls]):
                p_value = self._fisher_exact_test(cases_with, total_cases, controls_with, total_controls)
                
                if p_value < STATISTICAL_THRESHOLDS['fisher_exact']:
                    result['applies'] = True
                    result['details'] = f"Significant case-control association (p={p_value:.2e})"
                else:
                    result['details'] = f"Non-significant case-control data (p={p_value:.3f})"
        
        return result
    
    def _evaluate_pm1(self, variant_data) -> Dict[str, Any]:
        """Evaluate PM1 - Mutational hotspot."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        # This would require domain/hotspot annotation
        result['details'] = "Requires manual assessment of mutational hotspot"
        
        return result
    
    def _evaluate_pm2(self, variant_data) -> Dict[str, Any]:
        """Evaluate PM2 - Absent from controls."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        # Get gene-specific thresholds
        gene = variant_data.basic_info.get('gene', 'default')
        thresholds = GENE_SPECIFIC_THRESHOLDS.get(gene, GENE_SPECIFIC_THRESHOLDS['default'])
        
        gnomad_af = variant_data.population_data.get('gnomad_af')
        
        if gnomad_af is not None:
            if gnomad_af <= thresholds['PM2']:
                result['applies'] = True
                result['details'] = f"Absent or very rare in gnomAD (AF={gnomad_af})"
            else:
                result['details'] = f"Present in gnomAD (AF={gnomad_af})"
        else:
            result['details'] = "No population frequency data available"
        
        return result
    
    def _evaluate_pm3(self, variant_data) -> Dict[str, Any]:
        """Evaluate PM3 - In trans with pathogenic variant."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        inheritance = variant_data.genetic_data.get('inheritance')
        compound_het = variant_data.genetic_data.get('compound_het')
        
        if inheritance == 'AR' and compound_het == 'yes':
            # This would require assessment of the second variant
            result['details'] = "Requires manual assessment of second variant pathogenicity"
        
        return result
    
    def _evaluate_pm4(self, variant_data) -> Dict[str, Any]:
        """Evaluate PM4 - Protein length changes."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        
        if variant_type == 'inframe_indel':
            # Check if it's in a non-repeat region
            result['details'] = "In-frame indel - requires assessment of repeat region"
        elif 'stop_loss' in variant_type or 'nonstop' in variant_type:
            result['applies'] = True
            result['details'] = "Stop-loss variant causing protein extension"
        
        return result
    
    def _evaluate_pm5(self, variant_data) -> Dict[str, Any]:
        """Evaluate PM5 - Novel missense at pathogenic residue."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        # This would require database lookup for known pathogenic variants
        # at the same position
        result['details'] = "Requires manual verification of pathogenic variants at same position"
        
        return result
    
    def _evaluate_pm6(self, variant_data) -> Dict[str, Any]:
        """Evaluate PM6 - Assumed de novo."""
        result = {'applies': False, 'strength': 'Moderate', 'details': ''}
        
        # Try both genetic_data and functional_data for de novo status
        denovo_status = (variant_data.genetic_data.get('de_novo') or 
                        variant_data.genetic_data.get('denovo') or
                        variant_data.functional_data.get('de_novo') or
                        variant_data.functional_data.get('denovo'))
        
        if denovo_status == 'assumed':
            result['applies'] = True
            result['details'] = "Assumed de novo (paternity/maternity not confirmed)"
        
        return result
    
    def _evaluate_pp1(self, variant_data) -> Dict[str, Any]:
        """Evaluate PP1 - Segregation with disease."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        segregation = variant_data.functional_data.get('segregation')
        
        if segregation == 'cosegregates':
            result['applies'] = True
            result['details'] = "Cosegregates with disease in family"
        elif segregation == 'does_not_segregate':
            result['details'] = "Does not segregate with disease - see BS4"
        
        return result
    
    def _evaluate_pp2(self, variant_data) -> Dict[str, Any]:
        """Evaluate PP2 - Missense in gene with low benign rate."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        
        if variant_type == 'missense':
            # This would require gene-specific missense tolerance data
            result['details'] = "Requires gene-specific missense tolerance assessment"
        
        return result
    
    def _evaluate_pp3(self, variant_data) -> Dict[str, Any]:
        """Evaluate PP3 - In silico evidence."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        # Only for missense variants
        if variant_data.basic_info.get('variant_type') != 'missense':
            result['details'] = "PP3 only applies to missense variants"
            return result
        
        # Calculate VAMPP-score-like metascore
        vampp_score = self._calculate_vampp_score(variant_data)
        
        if vampp_score is not None:
            # Use VAMPP-score approach with gene-specific and disease prevalence considerations
            gnomad_af = variant_data.population_data.get('gnomad_af', 0)
            disease_prevalence = variant_data.population_data.get('disease_prevalence', 1e-4)
            
            # VAMPP-score thresholds based on allele frequency
            if gnomad_af > 1e-3:  # Common variants
                vampp_threshold = VAMPP_SCORE_THRESHOLDS['pp3']['common_variants']
            elif gnomad_af > 1e-4:  # Moderately rare variants
                vampp_threshold = VAMPP_SCORE_THRESHOLDS['pp3']['moderate_rare']
            else:  # Very rare variants
                vampp_threshold = VAMPP_SCORE_THRESHOLDS['pp3']['very_rare']
            
            # Consider disease prevalence for threshold adjustment
            prevalence_thresholds = VAMPP_SCORE_THRESHOLDS['prevalence_adjustment']
            adjustment = prevalence_thresholds['adjustment_factor']
            
            if disease_prevalence > prevalence_thresholds['common_disease']:  # Common diseases
                vampp_threshold += adjustment
            elif disease_prevalence < prevalence_thresholds['rare_disease']:  # Very rare diseases
                vampp_threshold -= adjustment
            
            # Ensure threshold is within reasonable bounds
            vampp_threshold = max(0.5, min(0.9, vampp_threshold))
            
            if vampp_score >= vampp_threshold:
                result['applies'] = True
                result['details'] = f"VAMPP-like metascore: {vampp_score:.3f} (≥{vampp_threshold:.3f})"
            else:
                result['details'] = f"VAMPP-like metascore: {vampp_score:.3f} (<{vampp_threshold:.3f})"
        else:
            # Fallback to traditional approach if VAMPP-score cannot be calculated
            pathogenic_count = 0
            total_count = 0
            
            for predictor, score in variant_data.insilico_data.items():
                if score is not None and predictor in INSILICO_THRESHOLDS:
                    total_count += 1
                    threshold = INSILICO_THRESHOLDS[predictor]['pathogenic']
                    
                    if predictor in ['sift', 'fathmm']:
                        # SIFT and FATHMM are inverted (lower = more pathogenic)
                        if score <= threshold:
                            pathogenic_count += 1
                    else:
                        if score >= threshold:
                            pathogenic_count += 1
            
            if total_count > 0:
                pathogenic_ratio = pathogenic_count / total_count
                if pathogenic_ratio >= 0.7:  # At least 70% pathogenic predictions
                    result['applies'] = True
                    result['details'] = f"Multiple pathogenic predictions ({pathogenic_count}/{total_count})"
                else:
                    result['details'] = f"Mixed or benign predictions ({pathogenic_count}/{total_count})"
        
        return result
    
    def _evaluate_pp4(self, variant_data) -> Dict[str, Any]:
        """Evaluate PP4 - Phenotype match."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        phenotype_match = variant_data.functional_data.get('phenotype_match')
        
        if phenotype_match == 'specific_match':
            result['applies'] = True
            result['details'] = "Phenotype highly specific for gene"
        elif phenotype_match == 'partial_match':
            result['details'] = "Partial phenotype match"
        elif phenotype_match == 'no_match':
            result['details'] = "No phenotype match"
        
        return result
    
    def _evaluate_pp5(self, variant_data) -> Dict[str, Any]:
        """Evaluate PP5 - Reputable source (2023 guidelines)."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        # Check ClinVar status from basic info (user input or API)
        clinvar_status = variant_data.basic_info.get('clinvar_status', '').lower()
        
        # Also check ClinVar data if available
        if not clinvar_status and hasattr(variant_data, 'clinvar_data'):
            clinvar_significance = variant_data.clinvar_data.get('significance', '').lower()
            clinvar_status = clinvar_significance
        
        if 'pathogenic' in clinvar_status and 'conflicting' not in clinvar_status:
            result['applies'] = True
            result['details'] = f"ClinVar classification: {clinvar_status}"
        elif clinvar_status:
            result['details'] = f"ClinVar classification: {clinvar_status} (not supporting pathogenicity)"
        else:
            result['details'] = "No ClinVar data available"
        
        return result
    
    def _evaluate_ba1(self, variant_data) -> Dict[str, Any]:
        """Evaluate BA1 - High allele frequency."""
        result = {'applies': False, 'strength': 'Stand-alone', 'details': ''}
        
        # Get gene-specific thresholds
        gene = variant_data.basic_info.get('gene', 'default')
        thresholds = GENE_SPECIFIC_THRESHOLDS.get(gene, GENE_SPECIFIC_THRESHOLDS['default'])
        
        gnomad_af = variant_data.population_data.get('gnomad_af')
        
        if gnomad_af is not None:
            if gnomad_af >= thresholds['BA1']:
                result['applies'] = True
                result['details'] = f"High allele frequency in gnomAD (AF={gnomad_af})"
            else:
                result['details'] = f"Allele frequency below BA1 threshold (AF={gnomad_af})"
        
        return result
    
    def _evaluate_bs1(self, variant_data) -> Dict[str, Any]:
        """Evaluate BS1 - Allele frequency too high for disorder."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        # Get gene-specific thresholds
        gene = variant_data.basic_info.get('gene', 'default')
        thresholds = GENE_SPECIFIC_THRESHOLDS.get(gene, GENE_SPECIFIC_THRESHOLDS['default'])
        
        gnomad_af = variant_data.population_data.get('gnomad_af')
        
        if gnomad_af is not None:
            if gnomad_af >= thresholds['BS1']:
                result['applies'] = True
                result['details'] = f"Allele frequency too high for disorder (AF={gnomad_af})"
            else:
                result['details'] = f"Allele frequency below BS1 threshold (AF={gnomad_af})"
        
        return result
    
    def _evaluate_bs2(self, variant_data) -> Dict[str, Any]:
        """Evaluate BS2 - Observed in healthy individual."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        # This would require specific information about healthy individuals
        result['details'] = "Requires manual assessment of healthy individual data"
        
        return result
    
    def _evaluate_bs3(self, variant_data) -> Dict[str, Any]:
        """Evaluate BS3 - Functional studies show no damaging effect."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        functional_studies = variant_data.functional_data.get('functional_studies')
        
        if functional_studies == 'benign':
            result['applies'] = True
            result['details'] = "Functional studies show no damaging effect"
        elif functional_studies == 'damaging':
            result['details'] = "Functional studies show damaging effect - see PS3"
        
        return result
    
    def _evaluate_bs4(self, variant_data) -> Dict[str, Any]:
        """Evaluate BS4 - Lack of segregation."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        segregation = variant_data.functional_data.get('segregation')
        
        if segregation == 'does_not_segregate':
            result['applies'] = True
            result['details'] = "Lack of segregation in affected family members"
        elif segregation == 'cosegregates':
            result['details'] = "Cosegregates with disease - see PP1"
        
        return result
    
    def _evaluate_bp1(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP1 - Missense in gene where truncating variants cause disease."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        
        if variant_type == 'missense':
            # This would require gene-specific mechanism information
            result['details'] = "Requires gene-specific mechanism assessment"
        
        return result
    
    def _evaluate_bp2(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP2 - Observed in trans/cis with pathogenic variant."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        # This would require phasing information
        result['details'] = "Requires phasing information and pathogenic variant assessment"
        
        return result
    
    def _evaluate_bp3(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP3 - In-frame indel in repetitive region."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        
        if variant_type == 'inframe_indel':
            result['details'] = "In-frame indel - requires repeat region assessment"
        
        return result
    
    def _evaluate_bp4(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP4 - Computational evidence suggests no impact."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        # Only for missense variants
        if variant_data.basic_info.get('variant_type') != 'missense':
            result['details'] = "BP4 only applies to missense variants"
            return result
        
        # Calculate VAMPP-score-like metascore
        vampp_score = self._calculate_vampp_score(variant_data)
        
        if vampp_score is not None:
            # Use VAMPP-score approach with prevalence considerations
            gnomad_af = variant_data.population_data.get('gnomad_af', 0)
            disease_prevalence = variant_data.population_data.get('disease_prevalence', 1e-4)
            
            # VAMPP-score thresholds for benign classification
            if gnomad_af > 1e-3:  # Common variants
                vampp_threshold = VAMPP_SCORE_THRESHOLDS['bp4']['common_variants']
            elif gnomad_af > 1e-4:  # Moderately rare variants
                vampp_threshold = VAMPP_SCORE_THRESHOLDS['bp4']['moderate_rare']
            else:  # Very rare variants
                vampp_threshold = VAMPP_SCORE_THRESHOLDS['bp4']['very_rare']
            
            # Consider disease prevalence for threshold adjustment
            prevalence_thresholds = VAMPP_SCORE_THRESHOLDS['prevalence_adjustment']
            adjustment = prevalence_thresholds['adjustment_factor']
            
            if disease_prevalence > prevalence_thresholds['common_disease']:  # Common diseases
                vampp_threshold -= adjustment
            elif disease_prevalence < prevalence_thresholds['rare_disease']:  # Very rare diseases
                vampp_threshold += adjustment
            
            # Ensure threshold is within reasonable bounds
            vampp_threshold = max(0.1, min(0.5, vampp_threshold))
            
            if vampp_score <= vampp_threshold:
                result['applies'] = True
                result['details'] = f"VAMPP-like metascore: {vampp_score:.3f} (≤{vampp_threshold:.3f})"
            else:
                result['details'] = f"VAMPP-like metascore: {vampp_score:.3f} (>{vampp_threshold:.3f})"
        else:
            # Fallback to traditional approach if VAMPP-score cannot be calculated
            benign_count = 0
            total_count = 0
            
            for predictor, score in variant_data.insilico_data.items():
                if score is not None and predictor in INSILICO_THRESHOLDS:
                    total_count += 1
                    threshold = INSILICO_THRESHOLDS[predictor]['benign']
                    
                    if predictor in ['sift', 'fathmm']:
                        # SIFT and FATHMM are inverted (higher = more benign)
                        if score > threshold:
                            benign_count += 1
                    else:
                        if score <= threshold:
                            benign_count += 1
            
            if total_count > 0:
                benign_ratio = benign_count / total_count
                if benign_ratio >= 0.7:  # At least 70% benign predictions
                    result['applies'] = True
                    result['details'] = f"Multiple benign predictions ({benign_count}/{total_count})"
                else:
                    result['details'] = f"Mixed or pathogenic predictions ({benign_count}/{total_count})"
        
        return result
    
    def _evaluate_bp5(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP5 - Alternate molecular basis for disease."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        # This would require information about other causal variants
        result['details'] = "Requires assessment of alternate molecular basis"
        
        return result
    
    def _evaluate_bp6(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP6 - Reputable source (2023 guidelines)."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        # Check ClinVar status from basic info (user input or API)
        clinvar_status = variant_data.basic_info.get('clinvar_status', '').lower()
        
        # Also check ClinVar data if available
        if not clinvar_status and hasattr(variant_data, 'clinvar_data'):
            clinvar_significance = variant_data.clinvar_data.get('significance', '').lower()
            clinvar_status = clinvar_significance
        
        if 'benign' in clinvar_status and 'conflicting' not in clinvar_status:
            result['applies'] = True
            result['details'] = f"ClinVar classification: {clinvar_status}"
        elif clinvar_status:
            result['details'] = f"ClinVar classification: {clinvar_status} (not supporting benign classification)"
        else:
            result['details'] = "No ClinVar data available"
        
        return result
    
    def _evaluate_bp7(self, variant_data) -> Dict[str, Any]:
        """Evaluate BP7 - Synonymous variant with no splice impact."""
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        
        if variant_type == 'synonymous':
            # This would require splice prediction assessment
            result['details'] = "Synonymous variant - requires splice impact assessment"
        
        return result
    
    def _calculate_vampp_score(self, variant_data) -> Optional[float]:
        """
        Calculate VAMPP-score-like metascore for variants.
        
        Based on the original VAMPP-score approach:
        - Uses weighted combination of in silico predictors
        - Gene-specific performance weighting (simplified here)
        - Min-max normalization to 0-1 range
        - For non-missense variants, uses available conservation/functional scores
        """
        # For intronic and other non-missense variants, use conservation scores if available
        if variant_data.basic_info.get('variant_type') not in ['missense', 'intronic', 'synonymous']:
            return None
        
        variant_type = variant_data.basic_info.get('variant_type')
        
        # For intronic variants, use conservation and splice prediction scores
        if variant_type == 'intronic':
            conservation_scores = []
            conservation_weights = []
            
            # Use conservation scores
            for score_key in ['phylop_100way', 'phylop_30way', 'phylop_17way', 'gerp_score', 'siphy_score']:
                if score_key in variant_data.insilico_data and variant_data.insilico_data[score_key] is not None:
                    score = variant_data.insilico_data[score_key]
                    
                    # Normalize conservation scores
                    if score_key.startswith('phylop'):
                        if score <= 1.0:  # Likely ranked score (0-1)
                            normalized_score = score
                        else:  # Raw score, normalize from -20,10 to 0,1
                            normalized_score = (score + 20) / 30
                    elif score_key == 'gerp_score':
                        # GERP++: normalize from -12.3,6.17 to 0,1
                        normalized_score = (score + 12.3) / 18.47
                    elif score_key == 'siphy_score':
                        # SiPhy: normalize from 0,29.9 to 0,1
                        normalized_score = min(score / 29.9, 1.0)
                    else:
                        normalized_score = score
                    
                    # Ensure within 0-1 range
                    normalized_score = max(0, min(1, normalized_score))
                    conservation_scores.append(normalized_score)
                    conservation_weights.append(0.2)  # Equal weights for conservation
            
            # Use splice prediction scores if available
            splice_scores = []
            splice_weights = []
            
            # SpliceAI scores (most important for intronic variants)
            spliceai_scores = []
            for score_key in ['spliceai_ag_score', 'spliceai_al_score', 'spliceai_dg_score', 'spliceai_dl_score']:
                if score_key in variant_data.insilico_data and variant_data.insilico_data[score_key] is not None:
                    score = variant_data.insilico_data[score_key]
                    normalized_score = max(0, min(1, score))
                    spliceai_scores.append(normalized_score)
            
            # If SpliceAI scores are available, use the maximum as it indicates strongest splice effect
            if spliceai_scores:
                max_spliceai = max(spliceai_scores)
                splice_scores.append(max_spliceai)
                splice_weights.append(0.4)  # High weight for SpliceAI
            
            # Other splice prediction scores
            for score_key in ['ada_score', 'rf_score', 'dbscsnv_ada_score', 'dbscsnv_rf_score']:
                if score_key in variant_data.insilico_data and variant_data.insilico_data[score_key] is not None:
                    score = variant_data.insilico_data[score_key]
                    # These are typically 0-1 range
                    normalized_score = max(0, min(1, score))
                    splice_scores.append(normalized_score)
                    splice_weights.append(0.15)  # Lower weights for other splice predictors
            
            # Combine all scores
            all_scores = conservation_scores + splice_scores
            all_weights = conservation_weights + splice_weights
            
            if not all_scores:
                return 0.0  # Return 0 instead of None for intronic variants
            
            # Calculate weighted average
            weighted_sum = sum(score * weight for score, weight in zip(all_scores, all_weights))
            total_weight = sum(all_weights)
            
            if total_weight > 0:
                return weighted_sum / total_weight
            else:
                return 0.0
        
        # Mapping from data keys to weight keys
        key_mapping = {
            'sift_score': 'sift',
            'polyphen2_hdiv_score': 'polyphen2',
            'polyphen2_hvar_score': 'polyphen2', 
            'mutationtaster_score': 'mutationtaster',
            'fathmm_score': 'fathmm',
            'gerp_score': 'gerp_pp',
            'phylop_100way': 'phylop_vert',
            'phylop_30way': 'phylop_mamm',
            'phylop_17way': 'phylop_prim',
            'lrt_score': 'lrt',
            'provean_score': 'provean',
            'vest4_score': 'vest4',
            'metasvm_score': 'metasvm',
            'metalr_score': 'metalr',
            'mutationassessor_score': 'mutationassessor',
            'cadd_phred': 'cadd_phred',
            'revel_score': 'revel',
            'alphamissense_score': 'alphamissense',
            'clinpred_score': 'clinpred',
            'bayesdel_addaf_score': 'bayesdel_addaf',
            'bayesdel_noaf_score': 'bayesdel_noaf',
            'metarnn_score': 'metarnn'
        }
        
        # Priority predictors with higher weights (based on general performance)
        predictor_scores = []
        predictor_weights = []
        
        for data_key, score in variant_data.insilico_data.items():
            if score is not None and data_key in key_mapping:
                weight_key = key_mapping[data_key]
                if weight_key in INSILICO_WEIGHTS:
                    # Get base weight
                    weight = INSILICO_WEIGHTS[weight_key]
                    
                    # Normalize scores to 0-1 range (pathogenic direction)
                    if data_key == 'sift_score':
                        # SIFT is inverted (lower = more pathogenic)
                        normalized_score = 1 - score
                    elif data_key == 'cadd_phred':
                        # CADD phred: normalize using sigmoid-like transformation
                        normalized_score = min(score / 40, 1.0)  # More realistic upper bound
                    elif data_key in ['bayesdel_addaf_score', 'bayesdel_noaf_score']:
                        # BayesDel: normalize from -1,1 to 0,1
                        normalized_score = (score + 1) / 2
                    elif data_key == 'gerp_score':
                        # GERP++: normalize from -12.3,6.17 to 0,1
                        normalized_score = (score + 12.3) / 18.47
                    elif data_key in ['phylop_100way', 'phylop_30way', 'phylop_17way']:
                        # PhyloP scores: check if they are raw or ranked
                        if score <= 1.0:  # Likely ranked score (0-1)
                            normalized_score = score
                        else:  # Raw score, normalize from -20,10 to 0,1
                            normalized_score = (score + 20) / 30
                    elif data_key == 'fathmm_score':
                        # FATHMM: negative values are pathogenic, normalize from -16,16 to 0,1
                        normalized_score = 1 - ((score + 16) / 32)
                    elif data_key == 'mutationtaster_score':
                        # MutationTaster: depends on type (score or ranked)
                        if score <= 1.0:
                            normalized_score = score
                        else:
                            # If raw score, normalize
                            normalized_score = min(score / 10, 1.0)
                    elif data_key in ['revel_score', 'alphamissense_score', 'clinpred_score', 
                                     'polyphen2_hdiv_score', 'polyphen2_hvar_score',
                                     'metarnn_score', 'mutationassessor_score', 'vest4_score',
                                     'metasvm_score', 'metalr_score', 'lrt_score']:
                        # These are already 0-1 range or normalized
                        if data_key == 'provean_score':
                            # PROVEAN: negative values are pathogenic, normalize from -14,14 to 0,1
                            normalized_score = 1 - ((score + 14) / 28)
                        elif data_key == 'mutationassessor_score':
                            # MutationAssessor: normalize from -5,5 to 0,1
                            normalized_score = (score + 5) / 10
                        elif data_key == 'lrt_score':
                            # LRT: already 0-1 range
                            normalized_score = score
                        else:
                            # Already 0-1 range
                            normalized_score = score
                    elif data_key == 'provean_score':
                        # PROVEAN: negative values are pathogenic, normalize from -14,14 to 0,1
                        normalized_score = 1 - ((score + 14) / 28)
                    else:
                        # Default: assume 0-1 range
                        normalized_score = score
                    
                    # Ensure normalized score is within 0-1 range
                    normalized_score = max(0, min(1, normalized_score))
                    
                    predictor_scores.append(normalized_score)
                    predictor_weights.append(weight)
        
        if not predictor_scores:
            return None
        
        # Calculate weighted average (VAMPP-like approach)
        weighted_sum = sum(score * weight for score, weight in zip(predictor_scores, predictor_weights))
        total_weight = sum(predictor_weights)
        
        if total_weight > 0:
            raw_vampp_score = weighted_sum / total_weight
            
            # Apply sigmoid transformation for better distribution
            # This mimics the ranked score transformation in original VAMPP
            import math
            vampp_score = 1 / (1 + math.exp(-5 * (raw_vampp_score - 0.5)))
            
            return vampp_score
        else:
            return None
    
    def _perform_statistical_tests(self, variant_data) -> Dict[str, Any]:
        """Perform statistical tests on variant data."""
        results = {}
        
        # Fisher's exact test for case-control data
        if variant_data.functional_data.get('case_control') == 'yes':
            cases_with = variant_data.functional_data.get('cases_with_variant')
            total_cases = variant_data.functional_data.get('total_cases')
            controls_with = variant_data.functional_data.get('controls_with_variant')
            total_controls = variant_data.functional_data.get('total_controls')
            
            if all(x is not None for x in [cases_with, total_cases, controls_with, total_controls]):
                p_value = self._fisher_exact_test(cases_with, total_cases, controls_with, total_controls)
                results['fisher_exact'] = {
                    'p_value': p_value,
                    'significant': p_value < STATISTICAL_THRESHOLDS['fisher_exact']
                }
        
        return results
    
    def _fisher_exact_test(self, cases_with: int, total_cases: int, 
                          controls_with: int, total_controls: int) -> float:
        """Perform Fisher's exact test."""
        try:
            # Create contingency table
            table = np.array([
                [cases_with, total_cases - cases_with],
                [controls_with, total_controls - controls_with]
            ])
            
            # Perform Fisher's exact test
            test_result = stats.fisher_exact(table)
            
            # Extract p-value (second element of the result)
            if isinstance(test_result, tuple) and len(test_result) >= 2:
                p_val = test_result[1]  # type: ignore
                return float(p_val)  # type: ignore
            else:
                return 1.0  # Default p-value if test fails
        except Exception:
            return 1.0  # Default p-value if test fails
    
    def _consolidate_applied_criteria(self, pathogenic: Dict[str, Any], 
                                    benign: Dict[str, Any]) -> Dict[str, List[str]]:
        """Consolidate applied criteria from pathogenic and benign evaluations."""
        applied = {
            'pathogenic': [],
            'benign': []
        }
        
        # Collect applied pathogenic criteria
        for criterion, result in pathogenic.items():
            if result.get('applies', False):
                applied['pathogenic'].append(criterion)
        
        # Collect applied benign criteria
        for criterion, result in benign.items():
            if result.get('applies', False):
                applied['benign'].append(criterion)
        
        return applied
