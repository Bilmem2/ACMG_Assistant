# =============================================================================
# LEGACY VariantData Class
# =============================================================================
# NOTE: This is a LEGACY implementation kept for backward compatibility.
# CANONICAL IMPLEMENTATION: Use core/variant_data.py (dataclass version) instead.
# TODO: Migrate all usages to core/variant_data.VariantData and remove this class.
# =============================================================================
from typing import Optional


class VariantData:
    """
    Legacy VariantData class for backward compatibility.
    
    DEPRECATED: Use core.variant_data.VariantData (dataclass) instead.
    This class is maintained only for existing code that depends on it.
    
    Args:
        basic_info: Basic variant information (gene, position, alleles, etc.)
        population_data: Population frequency data from gnomAD, ExAC, etc.
        insilico_data: In silico prediction scores (REVEL, CADD, etc.)
        genetic_data: Genetic context (inheritance, zygosity, etc.)
        functional_data: Functional study results
        patient_phenotypes: Patient phenotype information
        clinvar_data: ClinVar database information
    """
    def __init__(self, basic_info=None, population_data=None, insilico_data=None, genetic_data=None, functional_data=None, patient_phenotypes=None, clinvar_data=None):
        self.basic_info = basic_info or {}
        self.population_data = population_data or {}
        self.insilico_data = insilico_data or {}
        self.genetic_data = genetic_data or {}
        self.functional_data = functional_data or {}
        self.patient_phenotypes = patient_phenotypes
        self.clinvar_data = clinvar_data

    @property
    def gene(self) -> Optional[str]:
        """Get gene symbol from basic_info."""
        return self.basic_info.get('gene', None)

    @property
    def hgvs_c(self) -> Optional[str]:
        """Get HGVS cDNA notation from basic_info."""
        return self.basic_info.get('hgvs_c', None)
# =============================================================================
# InframeAnalyzer Class
# =============================================================================
# NOTE: This class provides specialized analysis for in-frame deletions/insertions.
# Some methods are placeholder implementations pending full integration.
# =============================================================================
class InframeAnalyzer:
    """
    Analyzer for in-frame insertion/deletion variants.
    
    Evaluates whether in-frame indels affect critical functional regions,
    which influences PM4 (protein length changes) and BP3 (repeat region) criteria.
    
    Attributes:
        critical_regions: Dict mapping genes to critical region coordinates
        repeat_regions: Dict mapping genes to repeat region coordinates
        domain_boundaries: Dict mapping genes to functional domain boundaries
    """
    def __init__(self):
        import logging
        self.logger = logging.getLogger("InframeAnalyzer")
        self.critical_regions = self._load_critical_regions()
        self.repeat_regions = self._load_repeat_regions()
        self.domain_boundaries = self._load_domain_boundaries()

    def evaluate_inframe_deletion(self, variant_data) -> Optional[str]:
        """
        Evaluate in-frame deletion for ACMG criteria applicability.
        
        Args:
            variant_data: VariantData object containing variant information
            
        Returns:
            str or None: 'PM4' if deletion affects critical region/domain,
                        'BP3' if in repeat region, None otherwise
        """
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
        
        # Initialize API client for gnomAD/ClinVar integration
        from utils.api_client import APIClient
        from config.constants import API_SETTINGS
        
        try:
            self.api_client = APIClient(cache_enabled=True)
            self.api_enabled = API_SETTINGS.get('enabled', True)
        except Exception as e:
            print(f"⚠️  Warning: Could not initialize API client: {str(e)}")
            self.api_client = None
            self.api_enabled = False
        
        # Initialize multi-source API clients for predictors and population data
        self._init_multi_source_clients()
        
        # Yeni modüller entegre ediliyor
        from core.functional_studies_evaluator import FunctionalStudiesEvaluator
        from core.phenotype_matcher import PhenotypeMatcher
        from core.missense_evaluator import MissenseEvaluator
        from utils.statistical_utils import StatisticalAnalyzer
        self.functional_studies_evaluator = FunctionalStudiesEvaluator()
        self.phenotype_matcher = PhenotypeMatcher()
        self.statistical_analyzer = StatisticalAnalyzer()
        self.missense_evaluator = MissenseEvaluator()
        
        # Interactive evidence collector for literature-based criteria
        # Initialized lazily when needed (to avoid prompts in automated mode)
        self._interactive_collector = None
        self._manual_evidence = None
    
    def _init_multi_source_clients(self) -> None:
        """
        Initialize multi-source API clients for predictor and population data.
        
        These clients implement the "fetch once, interpret many" pattern where
        all external data is fetched at the start of evaluation, then evaluators
        work as pure interpreters of pre-fetched data.
        
        Uses a shared ResultCache instance for validated caching across all
        API clients.
        """
        from config.constants import API_SETTINGS
        
        try:
            from utils.predictor_api_client import PredictorAPIClient, PopulationAPIClient
            
            # Try to initialize shared ResultCache for validated caching
            result_cache = None
            try:
                from utils.cache import ResultCache
                result_cache = ResultCache(enabled=not self.test_mode)
            except ImportError:
                pass  # Cache module not available - continue without caching
            
            self.predictor_client = PredictorAPIClient(
                api_enabled=API_SETTINGS.get('enabled', True),
                timeout=API_SETTINGS.get('timeout', 15),
                test_mode=self.test_mode,
                result_cache=result_cache
            )
            self.population_client = PopulationAPIClient(
                api_enabled=API_SETTINGS.get('enabled', True),
                timeout=API_SETTINGS.get('timeout', 15),
                test_mode=self.test_mode,
                result_cache=result_cache
            )
            
            # Store cache reference for potential direct access
            self._result_cache = result_cache
            
        except ImportError as e:
            print(f"⚠️  Warning: Multi-source API clients not available: {str(e)}")
            self.predictor_client = None
            self.population_client = None
            self._result_cache = None
        except Exception as e:
            print(f"⚠️  Warning: Could not initialize multi-source clients: {str(e)}")
            self.predictor_client = None
            self.population_client = None
            self._result_cache = None
    
    def _fetch_external_data(self, variant_data) -> None:
        """
        Pre-fetch all external data for a variant before evaluation.
        
        This method implements the "fetch once, interpret many" pattern.
        All external predictor and population data is fetched here and
        stored in the variant_data object for use by evaluators.
        
        This approach:
        - Minimizes API calls by fetching once per variant
        - Enables offline testing with mocked data
        - Makes evaluators pure functions that interpret pre-fetched data
        - Provides graceful degradation when sources are unavailable
        
        Args:
            variant_data: VariantData object to populate with external data
        """
        # Skip if data already fetched (avoid redundant API calls)
        if hasattr(variant_data, 'predictor_scores') and variant_data.predictor_scores:
            return
        
        # Skip if insilico_data already present (backward compatibility)
        # This allows tests to provide their own predictor data without API override
        insilico_data = getattr(variant_data, 'insilico_data', {}) or {}
        if insilico_data:
            # Don't fetch from API if test data is already provided
            return
        
        basic_info = variant_data.basic_info or {}
        chrom = basic_info.get('chromosome')
        pos = basic_info.get('position')
        ref = basic_info.get('ref_allele')
        alt = basic_info.get('alt_allele')
        gene = basic_info.get('gene')
        
        # Fetch predictor scores from multi-source API
        if self.predictor_client and chrom and pos and ref and alt:
            try:
                variant_data.predictor_scores = self.predictor_client.get_predictor_scores(
                    chrom=str(chrom),
                    pos=int(pos),
                    ref=str(ref),
                    alt=str(alt),
                    gene=gene
                )
                
                # Log available predictors
                if variant_data.predictor_scores:
                    available = sum(
                        1 for s in variant_data.predictor_scores.values() 
                        if hasattr(s, 'value') and s.value is not None
                    )
                    if available > 0:
                        print(f"📊 Fetched {available} predictor scores from external APIs")
                        
            except Exception as e:
                print(f"⚠️  Warning: Could not fetch predictor scores: {str(e)}")
                variant_data.predictor_scores = None
        
        # Fetch population statistics from multi-source API
        if self.population_client and chrom and pos and ref and alt:
            try:
                variant_data.population_stats = self.population_client.get_population_stats(
                    chrom=str(chrom),
                    pos=int(pos),
                    ref=str(ref),
                    alt=str(alt)
                )
                
                # Log population data availability
                if variant_data.population_stats:
                    sources = list(variant_data.population_stats.keys())
                    if sources:
                        print(f"📊 Fetched population data from: {', '.join(sources)}")
                        
            except Exception as e:
                print(f"⚠️  Warning: Could not fetch population stats: {str(e)}")
                variant_data.population_stats = None
    
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
        
        # Pre-fetch external data (predictors, population frequencies)
        # This implements "fetch once, interpret many" pattern
        self._fetch_external_data(variant_data)
        
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
        
        # PP5 - Reputable source (available in both 2015 and 2023)
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
        """
        Evaluate PVS1 - Null variant in LOF gene.
        
        Now enhanced with:
        1. gnomAD constraint API integration for dynamic LOF intolerance determination
        2. ClinGen Dosage Sensitivity Map for haploinsufficiency evidence
        3. PVS1 strength modulation based on HI score (Quick Win #4)
        """
        result = {'applies': False, 'strength': 'Very Strong', 'details': '', 'data_source': 'manual'}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        consequence = variant_data.basic_info.get('consequence', '').lower()
        gene = variant_data.basic_info.get('gene', '').upper()
        
        # Check if it's a loss-of-function variant
        lof_types = ['nonsense', 'frameshift', 'splice_donor', 'splice_acceptor', 'start_lost', 'stop_lost']
        lof_consequences = ['stop_gained', 'frameshift_variant', 'splice_donor_variant', 
                           'splice_acceptor_variant', 'start_lost', 'stop_lost']
        
        if variant_type in lof_types or consequence in lof_consequences:
            # ENHANCED: Multi-source validation (gnomAD + ClinGen + Dosage Sensitivity)
            lof_status = self._check_lof_intolerance(gene)
            clingen_status = self._check_clingen_validity(gene)
            dosage_status = self._check_clingen_dosage_sensitivity(gene)  # Quick Win #4: NEW
            
            # Combine all three sources for decision
            if lof_status['source'] == 'gnomAD' or clingen_status['source'] == 'ClinGen' or dosage_status['source'] == 'ClinGen Dosage':
                # API-based determination (MULTI-SOURCE)
                pLI_val = lof_status.get('pLI', 0)
                LOEUF_val = lof_status.get('LOEUF', 0)
                
                # Check for conflict with manual list
                in_manual_intolerant = gene in LOF_INTOLERANT_GENES
                in_manual_tolerant = gene in LOF_TOLERANT_GENES
                
                # Decision logic: ClinGen LOF mechanism OR gnomAD intolerant OR dosage HI
                apply_pvs1 = False
                data_sources = []
                details_parts = []
                
                # ClinGen takes priority for disease-specific mechanism
                if clingen_status.get('supports_lof_pathogenicity'):
                    apply_pvs1 = True
                    lof_diseases = clingen_status.get('lof_diseases', [])
                    disease_str = lof_diseases[0] if lof_diseases else "disease"
                    details_parts.append(f"LOF mechanism established for {disease_str} (ClinGen)")
                    data_sources.append(f"ClinGen ({clingen_status.get('confidence', 'unknown')} confidence)")
                
                # gnomAD constraint as secondary evidence
                if lof_status['is_lof_intolerant']:
                    apply_pvs1 = True
                    details_parts.append(f"pLI={pLI_val:.3f}, LOEUF={LOEUF_val:.3f} (gnomAD)")
                    data_sources.append(f"gnomAD v4 ({lof_status['confidence']} confidence)")
                elif lof_status['classification'] == 'LOF_tolerant' and not apply_pvs1:
                    # Only if ClinGen didn't support LOF
                    details_parts.append(f"pLI={pLI_val:.3f}, LOEUF={LOEUF_val:.3f} → LOF tolerant (gnomAD)")
                    data_sources.append(f"gnomAD v4 ({lof_status['confidence']} confidence)")
                
                # QUICK WIN #4: ClinGen Dosage Sensitivity - Modulate PVS1 strength
                hi_score = dosage_status.get('haploinsufficiency_score')
                pvs1_recommendation = dosage_status.get('pvs1_recommendation', 'insufficient_data')
                
                if hi_score is not None:
                    if hi_score == 3:
                        # Sufficient evidence: Keep PVS1 Very Strong
                        result['strength'] = 'Very Strong'
                        details_parts.append(f"HI Score=3 (Sufficient) - PVS1 Very Strong maintained")
                    elif hi_score == 2:
                        # Some evidence: Downgrade to PS1 Strong
                        result['strength'] = 'Strong'
                        details_parts.append(f"HI Score=2 (Some) - PVS1 downgraded to PS1 Strong")
                    elif hi_score == 1:
                        # Little evidence: Downgrade to PM2 Moderate
                        result['strength'] = 'Moderate'
                        details_parts.append(f"HI Score=1 (Little) - PVS1 downgraded to PM2 Moderate")
                    elif hi_score == 0 or hi_score == 40:
                        # No/unlikely HI: Consider not applying
                        apply_pvs1 = False
                        details_parts.append(f"HI Score={hi_score} (No/Unlikely evidence) - PVS1 not applied")
                    
                    data_sources.append(f"ClinGen Dosage (HI={hi_score})")
                else:
                    # No dosage data: Keep original PVS1 strength
                    details_parts.append("HI data unavailable - default PVS1 strength")
                
                # Conflict warnings
                conflict_warning = ""
                if in_manual_intolerant and not apply_pvs1:
                    conflict_warning = f" ⚠️ NOTE: Manual list classifies {gene} as LOF intolerant, but API data suggests tolerant"
                elif in_manual_tolerant and apply_pvs1:
                    conflict_warning = f" ⚠️ NOTE: Manual list classifies {gene} as LOF tolerant, but API data suggests intolerant"
                
                # Special case: ClinGen supports LOF but gnomAD says tolerant (e.g., BRCA1)
                if clingen_status.get('supports_lof_pathogenicity') and lof_status['classification'] == 'LOF_tolerant':
                    conflict_warning += f" 📘 CONTEXT: Gene is LOF tolerant in general population but LOF variants cause specific diseases (tumor suppressor, etc.)"
                
                # Build result
                if apply_pvs1:
                    result['applies'] = True
                    result['details'] = (
                        f"Loss-of-function variant ({variant_type or consequence}) in {gene}. "
                        f"{' | '.join(details_parts)}{conflict_warning}"
                    )
                    result['data_source'] = ' + '.join(data_sources)
                    result['constraint_metrics'] = {
                        'gnomad_pLI': pLI_val,
                        'gnomad_LOEUF': LOEUF_val,
                        'gnomad_classification': lof_status.get('classification'),
                        'clingen_supports_lof': clingen_status.get('supports_lof_pathogenicity'),
                        'clingen_mechanism': clingen_status.get('primary_mechanism'),
                        'clingen_diseases': clingen_status.get('lof_diseases', []),
                        'clingen_hi_score': hi_score,
                        'pvs1_recommendation': pvs1_recommendation
                    }
                else:
                    result['applies'] = False
                    result['details'] = (
                        f"Loss-of-function variant ({variant_type or consequence}) in {gene}. "
                        f"{' | '.join(details_parts)} - PVS1 not applicable{conflict_warning}"
                    )
                    result['data_source'] = ' + '.join(data_sources) if data_sources else 'API'
                    result['constraint_metrics'] = {
                        'gnomad_pLI': pLI_val,
                        'gnomad_LOEUF': LOEUF_val,
                        'gnomad_classification': lof_status.get('classification'),
                        'clingen_supports_lof': clingen_status.get('supports_lof_pathogenicity'),
                        'clingen_mechanism': clingen_status.get('primary_mechanism'),
                        'clingen_hi_score': hi_score,
                        'pvs1_recommendation': pvs1_recommendation
                    }
            else:
                # Fallback to manual hard-coded list (API unavailable)
                result = self._check_lof_manual(gene, variant_type, consequence, result)
                result['details'] += " [gnomAD/ClinGen APIs unavailable - using manual classification]"
        
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
                # Check LOF intolerance (API or manual)
                lof_status = self._check_lof_intolerance(gene)
                
                if lof_status['is_lof_intolerant']:
                    result['applies'] = True
                    result['strength'] = 'Strong'  # Use PS1 strength for splice-altering variants
                    result['details'] = f"Intronic variant with high SpliceAI score (max: {max(spliceai_scores):.3f}) in LOF intolerant gene {gene}"
                    result['data_source'] = lof_status['source']
                else:
                    result['applies'] = False
                    result['details'] = f"Intronic variant with high SpliceAI score (max: {max(spliceai_scores):.3f}) in gene {gene} - LOF tolerance unknown"
            elif spliceai_scores and max(spliceai_scores) > 0.2:
                result['applies'] = False
                result['details'] = f"Intronic variant with moderate SpliceAI score (max: {max(spliceai_scores):.3f}) - consider PS1 or PM4"
            else:
                result['details'] = "Intronic variant with low splice impact prediction"
        
        return result
    
    def _check_lof_intolerance(self, gene: str) -> Dict[str, Any]:
        """
        Check if gene is LOF intolerant using gnomAD API or fallback to manual list.
        
        Args:
            gene: Gene symbol (e.g., 'BRCA1')
            
        Returns:
            Dict with:
                - is_lof_intolerant (bool or None)
                - classification (str): 'LOF_intolerant', 'LOF_tolerant', 'uncertain', 'unknown'
                - source (str): 'gnomAD' or 'manual'
                - pLI, LOEUF, etc. (if from gnomAD)
                - confidence (str)
        """
        from config.constants import API_SETTINGS
        
        # Try API if enabled
        if API_SETTINGS.get('enabled', True) and hasattr(self, 'api_client'):
            try:
                constraint_data = self.api_client.get_gene_constraint(gene)
                
                if 'error' not in constraint_data:
                    return {
                        'is_lof_intolerant': constraint_data['is_lof_intolerant'],
                        'classification': constraint_data['classification'],
                        'source': 'gnomAD',
                        'confidence': constraint_data['confidence'],
                        'pLI': constraint_data.get('pLI'),
                        'LOEUF': constraint_data.get('LOEUF'),
                        'oe_lof': constraint_data.get('oe_lof')
                    }
            except Exception as e:
                print(f"⚠️  gnomAD API failed for {gene}: {str(e)}, falling back to manual list")
        
        # Fallback to manual list
        return {
            'is_lof_intolerant': None,
            'classification': 'unknown',
            'source': 'manual',
            'confidence': 'low'
        }
    
    def _check_lof_manual(self, gene: str, variant_type: str, consequence: str, result: Dict) -> Dict:
        """Check LOF intolerance using hard-coded gene lists (fallback method)."""
        if gene in LOF_INTOLERANT_GENES:
            result['applies'] = True
            result['details'] = f"Loss-of-function variant ({variant_type or consequence}) in LOF intolerant gene {gene} (manual list)"
            result['data_source'] = 'manual (hard-coded list)'
        elif gene in LOF_TOLERANT_GENES:
            result['applies'] = False
            result['details'] = f"Loss-of-function variant ({variant_type or consequence}) in LOF tolerant gene {gene} - PVS1 not applicable (manual list)"
            result['data_source'] = 'manual (hard-coded list)'
        else:
            # Unknown gene - conservative approach, don't apply PVS1
            result['applies'] = False
            result['details'] = f"Loss-of-function variant ({variant_type or consequence}) in gene {gene} with unknown LOF tolerance - manual review required"
            result['data_source'] = 'unknown (not in manual list, API unavailable)'
        
        return result
    
    def _check_clingen_validity(self, gene: str, disease: Optional[str] = None) -> Dict[str, Any]:
        """
        Check if gene has ClinGen gene-disease validity evidence for LOF mechanism.
        
        Args:
            gene: Gene symbol (e.g., 'BRCA1')
            disease: Optional disease name for filtering
        
        Returns:
            Dict with:
                - supports_lof_pathogenicity (bool or None)
                - source (str): 'ClinGen' or 'manual'
                - confidence (str)
                - lof_diseases (list): List of diseases with LOF mechanism
                - primary_mechanism (str): LOF mechanism type if found
        """
        from config.constants import API_SETTINGS
        
        # Try API if enabled
        if API_SETTINGS.get('enabled', True) and hasattr(self, 'api_client'):
            try:
                clingen_data = self.api_client.get_clingen_gene_validity(gene, disease)
                
                if 'error' not in clingen_data:
                    return {
                        'supports_lof_pathogenicity': clingen_data.get('supports_lof_pathogenicity'),
                        'source': 'ClinGen',
                        'confidence': clingen_data.get('confidence', 'unknown'),
                        'lof_diseases': clingen_data.get('lof_diseases', []),
                        'primary_mechanism': clingen_data.get('primary_mechanism'),
                        'curations': clingen_data.get('curations', [])
                    }
            except Exception as e:
                print(f"⚠️  ClinGen API failed for {gene}: {str(e)}")
        
        # Fallback: No ClinGen data available
        return {
            'supports_lof_pathogenicity': None,
            'source': 'manual',
            'confidence': 'low',
            'lof_diseases': [],
            'primary_mechanism': None
        }
    
    def _check_clingen_dosage_sensitivity(self, gene: str) -> Dict[str, Any]:
        """
        Check ClinGen Dosage Sensitivity Map for haploinsufficiency evidence.
        
        Quick Win #4 implementation: Query ClinGen for HI/TS scores to modulate PVS1 strength.
        
        Args:
            gene: Gene symbol (e.g., 'BRCA1', 'TP53')
        
        Returns:
            Dict with:
                - haploinsufficiency_score (int or None): 0-3 or 40
                - triplosensitivity_score (int or None): 0-3 or 40
                - pvs1_recommendation (str): Strength recommendation
                - source (str): 'ClinGen Dosage' or 'manual'
                - confidence (str): Confidence level
        """
        from config.constants import API_SETTINGS
        
        # Try API if enabled
        if API_SETTINGS.get('enabled', True) and hasattr(self, 'api_client'):
            try:
                dosage_data = self.api_client.get_clingen_dosage_sensitivity(gene)
                
                if 'error' not in dosage_data:
                    return {
                        'haploinsufficiency_score': dosage_data.get('haploinsufficiency_score'),
                        'triplosensitivity_score': dosage_data.get('triplosensitivity_score'),
                        'pvs1_recommendation': dosage_data.get('pvs1_recommendation'),
                        'source': 'ClinGen Dosage',
                        'confidence': dosage_data.get('confidence', 'unknown'),
                        'haploinsufficiency_description': dosage_data.get('haploinsufficiency_description'),
                        'clingen_gene_url': dosage_data.get('clingen_gene_url')
                    }
            except Exception as e:
                print(f"⚠️  ClinGen Dosage API failed for {gene}: {str(e)}")
        
        # Fallback: No dosage data available
        return {
            'haploinsufficiency_score': None,
            'triplosensitivity_score': None,
            'pvs1_recommendation': 'insufficient_data',
            'source': 'manual',
            'confidence': 'none'
        }

    
    def _evaluate_ps1(self, variant_data) -> Dict[str, Any]:
        """Evaluate PS1 - Same amino acid change as pathogenic variant."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        gene = variant_data.basic_info.get('gene')
        aa_change = variant_data.basic_info.get('amino_acid_change')
        hgvs_p = variant_data.basic_info.get('hgvs_p')
        
        # **QUICK WIN #3: Query ClinVar for variants at same position**
        if self.api_client and gene and hgvs_p:
            try:
                clinvar_search = self.api_client.search_clinvar_variants_at_position(gene, hgvs_p)
                
                if clinvar_search.get('same_aa_pathogenic'):
                    result['applies'] = True
                    result['details'] = f"✓ PS1 applies: Same amino acid change {hgvs_p} in {gene} reported as pathogenic in ClinVar"
                    result['source'] = 'ClinVar E-utilities'
                    return result
                else:
                    result['details'] = f"No pathogenic variants with same amino acid change found in ClinVar for {gene} {hgvs_p}"
                    return result
            except Exception as e:
                print(f"⚠️ ClinVar search failed: {e}")
        
        # FALLBACK: Check ClinVar data from variant_data
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
        """
        Enhanced Interactive PS1 evaluation with intelligent guidance system.
        
        Features:
        - Context-aware search recommendations
        - Evidence quality assessment
        - Decision tree guidance
        - Real-time validation feedback
        """
        result = {'applies': False, 'strength': 'Strong', 'details': '', 'guidance_used': True}
        
        print(f"\n🧬 ACMG PS1 Interactive Guidance System")
        print("=" * 60)
        print(f"🎯 Evaluating: {gene} {aa_change}")
        print(f"📖 Criterion: PS1 - Same amino acid change as pathogenic variant")
        print()
        
        # STEP 1: Context-aware search strategy
        print("🔍 STEP 1: INTELLIGENT SEARCH STRATEGY")
        print("─" * 40)
        
        # Gene-specific search recommendations
        gene_context = self._get_gene_context(gene)
        if gene_context['is_cancer_gene']:
            print(f"🎯 {gene} is a known cancer gene - Focus on:")
            print(f"   • Somatic mutation databases (COSMIC, TCGA)")
            print(f"   • Tumor suppressor/oncogene literature")
        elif gene_context['is_mendelian_gene']:
            print(f"🎯 {gene} is associated with Mendelian disease - Focus on:")
            print(f"   • Inherited variant databases (ClinVar, HGMD)")
            print(f"   • Family segregation studies")
        
        print(f"\n📚 Recommended search terms:")
        print(f"   🔸 Primary: '{gene} {aa_change} pathogenic'")
        print(f"   🔸 Alternative: '{gene} {aa_change} disease-causing'")
        print(f"   🔸 Functional: '{gene} {aa_change} functional analysis'")
        print(f"   🔸 Position: '{gene} codon {self._extract_codon_number(aa_change)} mutation'")
        
        print(f"\n🗄️  Database search priority:")
        print(f"   1️⃣  ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/")
        print(f"   2️⃣  HGMD Professional (if available)")
        print(f"   3️⃣  PubMed: https://pubmed.ncbi.nlm.nih.gov/")
        print(f"   4️⃣  LOVD: https://www.lovd.nl/")
        print(f"   5️⃣  Gene-specific databases")
        
        # STEP 2: Evidence quality assessment guide
        print(f"\n📊 STEP 2: EVIDENCE QUALITY ASSESSMENT")
        print("─" * 40)
        print("💡 PS1 requires HIGH-QUALITY evidence. Check for:")
        print()
        print("✅ STRONG EVIDENCE (Apply PS1):")
        print("   • Multiple independent reports of pathogenicity")
        print("   • Functional studies confirming damage")
        print("   • Clear disease association with segregation")
        print("   • Expert panel classification (ClinGen, etc.)")
        print()
        print("⚠️  MODERATE EVIDENCE (Consider PM5 instead):")
        print("   • Single case report without functional data")
        print("   • Conflicting interpretations in databases")
        print("   • Limited segregation data")
        print()
        print("❌ WEAK EVIDENCE (Do not apply PS1):")
        print("   • Only computational predictions")
        print("   • Uncertain significance classifications")
        print("   • No independent validation")
        
        # STEP 3: Interactive decision tree
        print(f"\n🌳 STEP 3: DECISION TREE GUIDANCE")
        print("─" * 40)
        
        # Question 1: Basic evidence check
        while True:
            print(f"\n❓ Q1: Did you find ANY reports of {aa_change} in {gene}?")
            q1 = input("   Enter: (y)es / (n)o / (h)elp: ").strip().lower()
            
            if q1 in ['h', 'help']:
                print(f"\n💡 HELP: Search tips for {gene} {aa_change}:")
                print(f"   • Try alternative amino acid notations (3-letter vs 1-letter)")
                print(f"   • Search for the codon position without specific change")
                print(f"   • Check if gene has alternative names/symbols")
                print(f"   • Look for functional domain information")
                continue
            elif q1 in ['n', 'no']:
                result['details'] = f"No reports found for {aa_change} in {gene} - PS1 does not apply"
                print(f"\n❌ PS1 Result: Does not apply")
                print(f"💡 Recommendation: Consider PM5 if different AA change at same position is pathogenic")
                return result
            elif q1 in ['y', 'yes']:
                break
            else:
                print("❌ Please enter 'y', 'n', or 'h'")
        
        # Question 2: Pathogenicity evidence
        while True:
            print(f"\n❓ Q2: Are these reports classified as PATHOGENIC/LIKELY PATHOGENIC?")
            print("   (Not VUS, not conflicting interpretations)")
            q2 = input("   Enter: (y)es / (n)o / (c)onflicting / (h)elp: ").strip().lower()
            
            if q2 in ['h', 'help']:
                print(f"\n💡 HELP: Pathogenicity assessment:")
                print(f"   ✅ PATHOGENIC: Clear disease-causing evidence")
                print(f"   ✅ LIKELY PATHOGENIC: Strong evidence, high confidence")
                print(f"   ❌ VUS: Uncertain significance, insufficient evidence")
                print(f"   ❌ CONFLICTING: Mixed pathogenic/benign interpretations")
                print(f"   💡 Look for: Expert review panels, functional validation")
                continue
            elif q2 in ['n', 'no']:
                result['details'] = f"Reports of {aa_change} in {gene} not classified as pathogenic - PS1 does not apply"
                print(f"\n❌ PS1 Result: Does not apply")
                return result
            elif q2 in ['c', 'conflicting']:
                result['details'] = f"Conflicting interpretations for {aa_change} in {gene} - PS1 requires manual expert review"
                print(f"\n⚠️  PS1 Result: Requires expert review")
                print(f"💡 Recommendation: Consult clinical geneticist or molecular pathologist")
                return result
            elif q2 in ['y', 'yes']:
                break
            else:
                print("❌ Please enter 'y', 'n', 'c', or 'h'")
        
        # Question 3: Evidence strength
        while True:
            print(f"\n❓ Q3: How strong is the pathogenic evidence?")
            print("   (s)trong - Multiple independent sources, functional data")
            print("   (m)oderate - Single well-documented case, some functional data")
            print("   (w)eak - Limited evidence, mostly computational")
            q3 = input("   Enter: (s)trong / (m)oderate / (w)eak / (h)elp: ").strip().lower()
            
            if q3 in ['h', 'help']:
                print(f"\n💡 HELP: Evidence strength criteria:")
                print(f"   🔥 STRONG: ≥2 independent studies + functional validation")
                print(f"   🔶 MODERATE: 1 good study + some supporting evidence")
                print(f"   🔸 WEAK: Case reports only, no functional data")
                print(f"   💡 Consider: Study quality, sample size, replication")
                continue
            elif q3 in ['s', 'strong']:
                result['applies'] = True
                result['strength'] = 'Strong'
                result['details'] = f"Strong evidence: Same amino acid change {aa_change} in {gene} reported as pathogenic with robust supporting data (PS1)"
                result['confidence'] = 'high'
                break
            elif q3 in ['m', 'moderate']:
                result['applies'] = True
                result['strength'] = 'Strong'  # PS1 is always Strong when applied
                result['details'] = f"Moderate evidence: Same amino acid change {aa_change} in {gene} reported as pathogenic (PS1)"
                result['confidence'] = 'medium'
                print(f"\n⚠️  Note: Consider additional supporting evidence (PM1, PM5, PP3)")
                break
            elif q3 in ['w', 'weak']:
                result['applies'] = False
                result['details'] = f"Weak evidence for {aa_change} in {gene} - insufficient for PS1"
                print(f"\n❌ PS1 Result: Does not apply (evidence too weak)")
                print(f"💡 Recommendation: Consider PM5 or PP3 instead")
                return result
            else:
                print("❌ Please enter 's', 'm', 'w', or 'h'")
        
        # STEP 4: Final validation and recommendations
        print(f"\n✅ STEP 4: FINAL VALIDATION")
        print("─" * 40)
        if result['applies']:
            print(f"🎉 PS1 APPLIES: {result['details']}")
            print(f"💪 Strength: {result['strength']} (ACMG Standard)")
            print(f"📊 Confidence: {result.get('confidence', 'medium')}")
            
            print(f"\n🔗 SYNERGISTIC CRITERIA TO CONSIDER:")
            print(f"   • PM1: If in mutational hotspot/critical domain")
            print(f"   • PP3: If computational predictors support pathogenicity")
            print(f"   • PP4: If patient phenotype matches gene-disease association")
            
        print(f"\n📝 DOCUMENTATION REMINDER:")
        print(f"   • Record specific publications/database entries used")
        print(f"   • Note evidence quality and any limitations")
        print(f"   • Consider periodic re-evaluation as new data emerges")
        
        return result
    
    def _get_gene_context(self, gene: str) -> Dict[str, Any]:
        """Get contextual information about gene for targeted search recommendations."""
        # This could be enhanced with API calls to gene databases
        cancer_genes = {'TP53', 'BRCA1', 'BRCA2', 'APC', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'KRAS', 'PIK3CA'}
        mendelian_genes = {'CFTR', 'DMD', 'FBN1', 'LDLR', 'PKU', 'HBB', 'F8', 'F9'}
        
        return {
            'is_cancer_gene': gene.upper() in cancer_genes,
            'is_mendelian_gene': gene.upper() in mendelian_genes,
            'gene_symbol': gene.upper()
        }
    
    def _extract_codon_number(self, aa_change: str) -> str:
        """Extract codon number from amino acid change notation."""
        import re
        if not aa_change:
            return "unknown"
        
        # Handle different formats: p.Arg273His, R273H, etc.
        match = re.search(r'(\d+)', aa_change)
        if match:
            return match.group(1)
        return "unknown"
    
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
        
        # Interactive guidance for PM1 if no automated data available
        if gene and not self.test_mode and not result.get('applies') and not result.get('manual_review'):
            result = self._evaluate_pm1_interactive(variant_data, gene, position)
        
        # Test mode or no gene: return negative result without API calls
        elif gene:
            # In test mode, use GeneSpecificRules which queries remote APIs only
            try:
                from core.gene_specific_rules import GeneSpecificRules
                rules = GeneSpecificRules()
                pm1_result = rules.evaluate_pm1(gene=gene, position=position)
                
                if pm1_result.applies:
                    result['applies'] = True
                    result['details'] = pm1_result.reason
                    result['api_source'] = pm1_result.source
                    result['confidence'] = pm1_result.confidence
                    if pm1_result.evidence_code == 'PM1_supporting':
                        result['strength'] = 'Supporting'
                else:
                    result['details'] = pm1_result.reason or "No hotspot/domain data from remote APIs"
            except ImportError:
                result['details'] = "No hotspot or functional domain data available"
            except Exception as e:
                result['details'] = f"PM1 evaluation error: {str(e)[:50]}"
        else:
            result['details'] = "Insufficient variant information for PM1 evaluation"
        
        return result
    
    def _evaluate_pm1_interactive(self, variant_data, gene: str, position: Optional[int]) -> Dict[str, Any]:
        """
        Enhanced Interactive PM1 evaluation with intelligent guidance system.
        
        Features:
        - Gene-specific hotspot databases
        - Functional domain analysis
        - Structural impact assessment
        - Evidence quality validation
        """
        result = {'applies': False, 'strength': 'Moderate', 'details': '', 'guidance_used': True}
        
        print(f"\n🎯 ACMG PM1 Interactive Guidance System")
        print("=" * 60)
        print(f"🧬 Evaluating: {gene}", end='')
        if position:
            print(f" position {position}")
        else:
            print()
        print(f"📖 Criterion: PM1 - Mutational hotspot/critical functional domain")
        print()
        
        # STEP 1: Gene-specific database recommendations
        print("🗄️  STEP 1: GENE-SPECIFIC DATABASE SEARCH")
        print("─" * 40)
        
        gene_databases = self._get_gene_specific_databases(gene)
        if gene_databases:
            print(f"🎯 Recommended databases for {gene}:")
            for db_name, db_info in gene_databases.items():
                print(f"   🔸 {db_name}: {db_info['url']}")
                print(f"      Purpose: {db_info['purpose']}")
        else:
            print(f"🔍 General databases to search:")
            print(f"   🔸 COSMIC: https://cancer.sanger.ac.uk/cosmic")
            print(f"   🔸 ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/")
            print(f"   🔸 UniProt: https://www.uniprot.org/")
            print(f"   🔸 InterPro: https://www.ebi.ac.uk/interpro/")
        
        print(f"\n📚 Search strategy:")
        print(f"   🔸 Hotspots: '{gene} hotspot mutation'")
        print(f"   🔸 Domains: '{gene} functional domain'")
        if position:
            print(f"   🔸 Position: '{gene} residue {position}'")
        print(f"   🔸 Structure: '{gene} crystal structure'")
        
        # STEP 2: Hotspot analysis
        print(f"\n🔥 STEP 2: MUTATIONAL HOTSPOT ANALYSIS")
        print("─" * 40)
        
        while True:
            print(f"\n❓ Q1: Is this position a documented mutational hotspot?")
            print("   (Check cancer databases, mutation frequency data)")
            q1 = input("   Enter: (y)es / (n)o / (u)nknown / (h)elp: ").strip().lower()
            
            if q1 in ['h', 'help']:
                print(f"\n💡 HELP: Hotspot identification criteria:")
                print(f"   🔥 HOTSPOT INDICATORS:")
                print(f"      • High mutation frequency in disease samples")
                print(f"      • Recurrent mutations in cancer/disease cohorts")
                print(f"      • Documented in COSMIC, TCGA, or disease databases")
                print(f"      • Multiple independent reports of pathogenic variants")
                print(f"   📊 QUANTITATIVE THRESHOLDS:")
                print(f"      • Cancer: ≥10 samples with mutations at position")
                print(f"      • Mendelian: ≥3 independent pathogenic reports")
                print(f"      • Population: Significantly higher than background")
                continue
            elif q1 in ['y', 'yes']:
                # Ask for evidence strength
                while True:
                    print(f"\n❓ Q1b: How strong is the hotspot evidence?")
                    print("   (s)trong - Multiple databases, high frequency")
                    print("   (m)oderate - Some evidence, moderate frequency")
                    strength = input("   Enter: (s)trong / (m)oderate / (h)elp: ").strip().lower()
                    
                    if strength in ['h', 'help']:
                        print(f"\n💡 HELP: Evidence strength for hotspots:")
                        print(f"   🔥 STRONG: COSMIC tier 1, >50 samples, multiple studies")
                        print(f"   🔶 MODERATE: Some database evidence, 10-50 samples")
                        continue
                    elif strength in ['s', 'strong']:
                        result['applies'] = True
                        result['details'] = f"Strong hotspot evidence: Position {position or 'unknown'} in {gene} is a well-documented mutational hotspot (PM1)"
                        result['confidence'] = 'high'
                        break
                    elif strength in ['m', 'moderate']:
                        result['applies'] = True
                        result['details'] = f"Moderate hotspot evidence: Position {position or 'unknown'} in {gene} shows hotspot characteristics (PM1)"
                        result['confidence'] = 'medium'
                        break
                    else:
                        print("❌ Please enter 's', 'm', or 'h'")
                break
            elif q1 in ['n', 'no']:
                break  # Continue to domain analysis
            elif q1 in ['u', 'unknown']:
                print(f"\n⚠️  Hotspot status unclear - continuing to domain analysis")
                break
            else:
                print("❌ Please enter 'y', 'n', 'u', or 'h'")
        
        # STEP 3: Functional domain analysis (if not hotspot)
        if not result['applies']:
            print(f"\n🧬 STEP 3: FUNCTIONAL DOMAIN ANALYSIS")
            print("─" * 40)
            
            while True:
                print(f"\n❓ Q2: Is this variant in a critical functional domain?")
                print("   (Check protein structure, domain databases)")
                q2 = input("   Enter: (y)es / (n)o / (u)nknown / (h)elp: ").strip().lower()
                
                if q2 in ['h', 'help']:
                    print(f"\n💡 HELP: Critical functional domain criteria:")
                    print(f"   🧬 CRITICAL DOMAINS:")
                    print(f"      • Active sites, binding sites, catalytic domains")
                    print(f"      • DNA-binding domains, kinase domains")
                    print(f"      • Structural domains essential for function")
                    print(f"      • Domains with known disease associations")
                    print(f"   ❌ NOT CRITICAL:")
                    print(f"      • Linker regions, flexible loops")
                    print(f"      • Domains with frequent benign variation")
                    print(f"      • Non-conserved regions")
                    continue
                elif q2 in ['y', 'yes']:
                    # Check for benign variation in domain
                    while True:
                        print(f"\n❓ Q2b: Does this domain have frequent benign variation?")
                        print("   (PM1 requires domain WITHOUT benign variation)")
                        benign = input("   Enter: (y)es / (n)o / (u)nknown / (h)elp: ").strip().lower()
                        
                        if benign in ['h', 'help']:
                            print(f"\n💡 HELP: Benign variation assessment:")
                            print(f"   ❌ FREQUENT BENIGN VARIATION (PM1 does not apply):")
                            print(f"      • Multiple benign/likely benign variants in domain")
                            print(f"      • High population frequency variants in domain")
                            print(f"      • Domain tolerates missense changes")
                            print(f"   ✅ NO BENIGN VARIATION (PM1 may apply):")
                            print(f"      • No or very few benign variants in domain")
                            print(f"      • Domain is highly conserved")
                            print(f"      • Missense changes typically pathogenic")
                            continue
                        elif benign in ['n', 'no']:
                            result['applies'] = True
                            result['details'] = f"Variant in critical functional domain without benign variation in {gene} (PM1)"
                            result['confidence'] = 'medium'
                            break
                        elif benign in ['y', 'yes']:
                            result['applies'] = False
                            result['details'] = f"Variant in functional domain with benign variation in {gene} - PM1 does not apply"
                            break
                        elif benign in ['u', 'unknown']:
                            result['applies'] = False
                            result['details'] = f"Functional domain status unclear for {gene} - insufficient evidence for PM1"
                            break
                        else:
                            print("❌ Please enter 'y', 'n', 'u', or 'h'")
                    break
                elif q2 in ['n', 'no']:
                    result['applies'] = False
                    result['details'] = f"Variant not in critical functional domain in {gene} - PM1 does not apply"
                    break
                elif q2 in ['u', 'unknown']:
                    result['applies'] = False
                    result['details'] = f"Functional domain status unclear for {gene} - insufficient evidence for PM1"
                    break
                else:
                    print("❌ Please enter 'y', 'n', 'u', or 'h'")
        
        # STEP 4: Final validation and recommendations
        print(f"\n✅ STEP 4: FINAL VALIDATION")
        print("─" * 40)
        
        if result['applies']:
            print(f"🎉 PM1 APPLIES: {result['details']}")
            print(f"💪 Strength: {result['strength']} (ACMG Standard)")
            print(f"📊 Confidence: {result.get('confidence', 'medium')}")
            
            print(f"\n🔗 SYNERGISTIC CRITERIA TO CONSIDER:")
            print(f"   • PS1: If same amino acid change is pathogenic")
            print(f"   • PM5: If different pathogenic change at same position")
            print(f"   • PP2: If gene has low rate of benign missense")
            print(f"   • PP3: If computational predictors support pathogenicity")
        else:
            print(f"❌ PM1 DOES NOT APPLY: {result['details']}")
            print(f"\n💡 ALTERNATIVE CRITERIA TO CONSIDER:")
            print(f"   • PM5: Different amino acid change at same position")
            print(f"   • PP2: Missense in gene with low benign rate")
            print(f"   • PP3: Computational evidence supports pathogenicity")
        
        print(f"\n📝 DOCUMENTATION REMINDER:")
        print(f"   • Record specific databases and evidence used")
        print(f"   • Note domain boundaries and functional importance")
        print(f"   • Consider structural/functional validation studies")
        
        return result
    
    def _get_gene_specific_databases(self, gene: str) -> Dict[str, Dict[str, str]]:
        """Get gene-specific database recommendations for hotspot analysis."""
        gene = gene.upper()
        
        # Cancer genes
        cancer_genes = {
            'TP53': {
                'IARC TP53 Database': {
                    'url': 'https://p53.iarc.fr/',
                    'purpose': 'Comprehensive TP53 mutation database'
                }
            },
            'BRCA1': {
                'BRCA Exchange': {
                    'url': 'https://brcaexchange.org/',
                    'purpose': 'BRCA1/2 variant classification'
                }
            },
            'BRCA2': {
                'BRCA Exchange': {
                    'url': 'https://brcaexchange.org/',
                    'purpose': 'BRCA1/2 variant classification'
                }
            }
        }
        
        # Mendelian disease genes
        mendelian_genes = {
            'CFTR': {
                'CFTR2 Database': {
                    'url': 'https://cftr2.org/',
                    'purpose': 'Cystic fibrosis variant database'
                }
            },
            'DMD': {
                'Leiden DMD Database': {
                    'url': 'https://www.dmd.nl/',
                    'purpose': 'Duchenne muscular dystrophy variants'
                }
            }
        }
        
        # Return gene-specific databases if available
        if gene in cancer_genes:
            return cancer_genes[gene]
        elif gene in mendelian_genes:
            return mendelian_genes[gene]
        else:
            return {}
    
    def _evaluate_pm2(self, variant_data) -> Dict[str, Any]:
        """
        Evaluate PM2 - Absent from controls.
        
        **ENHANCED: Multi-source population frequency validation**
        - Primary: gnomAD v4 GraphQL API
        - Secondary: Manual input (ExAC, 1000 Genomes)
        - Cross-validation: Multiple population databases
        """
        result = {
            'applies': False, 
            'strength': 'Moderate', 
            'details': '',
            'confidence': 'high',
            'data_source': 'population_database'
        }
        
        # MULTI-SOURCE: Collect frequencies from all available sources
        frequencies = {}
        data_sources = []
        
        # Source 1: gnomAD v4 API
        gnomad_af = None
        if self.api_client and self.api_enabled:
            try:
                # Get genomic coordinates if available
                basic_info = variant_data.basic_info
                chrom = basic_info.get('chromosome')
                pos = basic_info.get('position')
                ref = basic_info.get('ref_allele')
                alt = basic_info.get('alt_allele')
                
                if chrom and pos and ref and alt:
                    freq_result = self.api_client.get_variant_frequency(
                        chrom=chrom, pos=pos, ref=ref, alt=alt
                    )
                    
                    if freq_result and 'allele_frequency' in freq_result:
                        gnomad_af = freq_result.get('allele_frequency')
                        frequencies['gnomad_v4'] = gnomad_af
                        data_sources.append('gnomAD v4 API')
                        
                        # If variant not found in gnomAD
                        if 'note' in freq_result and 'not found' in freq_result['note']:
                            frequencies['gnomad_v4'] = 0.0
            except Exception as e:
                print(f"⚠️  gnomAD API error in PM2: {str(e)}")
        
        # Source 2: Manual population_data (ExAC, 1000 Genomes, etc.)
        pop_data = variant_data.population_data or {}
        
        if gnomad_af is None:
            # Fallback to manual gnomAD input
            gnomad_af = pop_data.get('gnomad_af')
            if gnomad_af is not None:
                frequencies['gnomad_manual'] = gnomad_af
                data_sources.append('gnomAD (manual)')
        
        # Additional population databases from manual input
        if pop_data.get('exac_af') is not None:
            frequencies['exac'] = pop_data.get('exac_af')
            data_sources.append('ExAC')
        
        if pop_data.get('1000g_af') is not None:
            frequencies['1000g'] = pop_data.get('1000g_af')
            data_sources.append('1000 Genomes')
        
        if pop_data.get('topmed_af') is not None:
            frequencies['topmed'] = pop_data.get('topmed_af')
            data_sources.append('TOPMed')
        
        # MULTI-SOURCE VALIDATION: Consensus from multiple databases
        if len(frequencies) > 0:
            max_af = max(frequencies.values())
            min_af = min(frequencies.values())
            avg_af = sum(frequencies.values()) / len(frequencies)
            
            # Check for discrepancies between databases
            discrepancy_threshold = 0.01  # 1% difference
            has_discrepancy = (max_af - min_af) > discrepancy_threshold if len(frequencies) > 1 else False
            
            # CRITICAL CHECK: PM2 should NOT apply if frequency is high enough for BA1/BS1
            gene = variant_data.basic_info.get('gene', '').upper()
            gene_thresholds = GENE_SPECIFIC_THRESHOLDS.get(gene, GENE_SPECIFIC_THRESHOLDS['default'])
            bs1_threshold = gene_thresholds.get('BS1', 0.01)
            
            # Use maximum frequency for safety (conservative approach)
            if max_af >= bs1_threshold:
                result['applies'] = False
                result['details'] = (
                    f"Variant present in population at significant frequency "
                    f"(Max AF: {max_af:.4f} from {len(frequencies)} sources) "
                    f"[{', '.join(data_sources)}]"
                )
                result['confidence'] = 'high'
                result['population_frequencies'] = frequencies
                return result
            
            # PM2 applies if variant is absent or extremely rare across all databases
            if max_af == 0.0:
                result['applies'] = True
                result['details'] = (
                    f"Variant absent from {len(frequencies)} population database(s) "
                    f"(PM2) [{', '.join(data_sources)}]"
                )
                result['confidence'] = 'high' if len(frequencies) >= 2 else 'medium'
            elif max_af < 0.0001:  # Extremely rare (< 0.01%)
                result['applies'] = True
                result['details'] = (
                    f"Variant extremely rare in population "
                    f"(Max AF: {max_af:.6f} from {len(frequencies)} sources) "
                    f"(PM2) [{', '.join(data_sources)}]"
                )
                result['confidence'] = 'high' if len(frequencies) >= 2 else 'medium'
                
                # Warn about discrepancies
                if has_discrepancy:
                    result['details'] += f" ⚠️ Note: Frequency varies across databases (range: {min_af:.6f}-{max_af:.6f})"
            else:
                result['applies'] = False
                result['details'] = (
                    f"Variant present in population databases "
                    f"(Max AF: {max_af:.6f}) [{', '.join(data_sources)}]"
                )
                result['confidence'] = 'medium'
            
            result['population_frequencies'] = frequencies
            result['data_source'] = ', '.join(data_sources)
        else:
            # No frequency data available from any source
            # Check if absent from controls flag is explicitly set
            if pop_data.get('absent_from_controls') is True:
                result['applies'] = True
                result['details'] = "Variant reported as absent from population databases (PM2) [user-provided]"
                result['confidence'] = 'low'
            else:
                result['applies'] = False
                result['details'] = "No population frequency data available from any source"
                result['confidence'] = 'low'
                result['data_source'] = 'none'
        
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
        
        # **QUICK WIN #3: Query ClinVar for different variants at same position**
        if self.api_client and gene and hgvs_p:
            try:
                clinvar_search = self.api_client.search_clinvar_variants_at_position(gene, hgvs_p)
                
                different_aa_pathogenic = clinvar_search.get('different_aa_pathogenic', [])
                if different_aa_pathogenic:
                    result['applies'] = True
                    pathogenic_changes = [v['protein_change'] for v in different_aa_pathogenic[:3]]  # Show first 3
                    result['details'] = f"✓ PM5 applies: Different pathogenic missense at same position: {', '.join(pathogenic_changes)}"
                    result['source'] = 'ClinVar E-utilities'
                    return result
                else:
                    result['details'] = f"No different pathogenic variants found at same position for {gene} {hgvs_p}"
                    return result
            except Exception as e:
                print(f"⚠️ ClinVar search failed: {e}")
        
        # FALLBACK: Check ClinVar data for same residue different AA
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
            
            # **QUICK WIN #1: Use gnomAD missense constraint data**
            # Try to get gene constraint data from API
            constraint_data = None
            if self.api_client and gene:
                try:
                    constraint_data = self.api_client.get_gene_constraint(gene)
                except Exception as e:
                    print(f"⚠️ Could not fetch constraint data: {e}")
            
            # Check missense constraint metrics
            if constraint_data:
                mis_z = constraint_data.get('mis_z')
                oe_mis_upper = constraint_data.get('oe_mis_upper')
                is_mis_constrained = constraint_data.get('is_mis_constrained', False)
                mis_classification = constraint_data.get('mis_classification', 'uncertain')
                
                # Apply PP2 if gene shows significant missense constraint
                if is_mis_constrained and mis_classification == 'missense_constrained':
                    result['applies'] = True
                    if mis_z is not None:
                        result['details'] = f"✓ PP2 applies: {gene} shows missense constraint (mis_z={mis_z:.2f} > 3.09)"
                    elif oe_mis_upper is not None:
                        result['details'] = f"✓ PP2 applies: {gene} shows missense constraint (oe_mis_upper={oe_mis_upper:.3f} < 0.6)"
                    result['source'] = 'gnomAD missense constraint'
                    return result
                elif mis_classification == 'missense_tolerant':
                    result['details'] = f"Gene {gene} is missense tolerant (oe_mis={constraint_data.get('oe_mis', 'N/A'):.3f}), PP2 does not apply"
                    return result
            
            # Fallback: Define genes known to have low benign missense variation
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
                result['details'] = f"Gene {gene} not in low-benign-missense gene list (and no constraint data available)"
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
        """
        Evaluate PP3 - Computational evidence suggests pathogenic impact.
        
        For missense variants, uses the MissenseEvaluator composite score.
        For other variants, falls back to simple predictor counting.
        
        Returns dict with:
            - applies: bool
            - strength: 'Supporting', 'Moderate', or 'Strong'
            - details: str explanation
            - composite_score: float (for missense variants)
        """
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        
        # Use MissenseEvaluator for missense variants
        if variant_type == 'missense':
            return self._evaluate_missense_pp3_bp4(variant_data, direction='pathogenic')
        
        # Fallback for non-missense variants: simple predictor counting
        insilico_data = variant_data.insilico_data or {}
        
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
                    if score < threshold:
                        pathogenic_count += 1
                else:
                    if score > threshold:
                        pathogenic_count += 1
        
        if total_count > 0:
            pathogenic_ratio = pathogenic_count / total_count
            if pathogenic_ratio >= 0.7:
                result['applies'] = True
                result['details'] = f"Computational evidence suggests pathogenic impact ({pathogenic_count}/{total_count} predictors)"
            else:
                result['details'] = f"Computational evidence inconclusive ({pathogenic_count}/{total_count} predictors suggest pathogenic)"
        else:
            result['details'] = "No computational prediction scores available"
        
        return result
    
    def _evaluate_missense_pp3_bp4(self, variant_data, direction: str = 'pathogenic') -> Dict[str, Any]:
        """
        Evaluate PP3/BP4 for missense variants using composite score.
        
        Uses MissenseEvaluator to generate a composite score and maps it
        to ACMG evidence categories.
        
        Args:
            variant_data: VariantData object
            direction: 'pathogenic' for PP3, 'benign' for BP4
            
        Returns:
            Dict with applies, strength, details, and composite_score
        """
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        # Get composite score from MissenseEvaluator
        eval_result = self.missense_evaluator.evaluate_missense_variant(variant_data)
        
        composite_score = eval_result.get('composite_score', 0.5)
        evidence_category = eval_result.get('evidence_category', 'neutral')
        strength = eval_result.get('strength')
        eval_direction = eval_result.get('direction', 'neutral')
        confidence = eval_result.get('confidence', 'low')
        
        result['composite_score'] = composite_score
        result['sub_scores'] = eval_result.get('sub_scores', {})
        result['confidence'] = confidence
        
        # Check if evidence applies in the requested direction
        if direction == 'pathogenic' and eval_direction == 'pathogenic':
            result['applies'] = True
            result['strength'] = strength.capitalize() if strength else 'Supporting'
            result['details'] = (
                f"Missense composite score {composite_score:.3f} → {evidence_category} "
                f"(confidence: {confidence})"
            )
        elif direction == 'benign' and eval_direction == 'benign':
            result['applies'] = True
            result['strength'] = strength.capitalize() if strength else 'Supporting'
            result['details'] = (
                f"Missense composite score {composite_score:.3f} → {evidence_category} "
                f"(confidence: {confidence})"
            )
        else:
            result['details'] = (
                f"Missense composite score {composite_score:.3f} → {evidence_category} "
                f"(does not support {direction} evidence)"
            )
        
        return result
    
    def _evaluate_pp4(self, variant_data) -> Dict[str, Any]:
        """
        Evaluate PP4 - Patient phenotype is highly specific for a disease with a single genetic etiology.
        
        This method uses the PhenotypeMatcher to evaluate phenotype-genotype correlation
        using local gene-phenotype association data and Jaccard similarity scoring.
        
        PP4 Evidence Thresholds:
            - similarity >= 0.8: PP4 applies (strong phenotype match)
            - similarity >= 0.5: PP4_supporting (moderate match)
            - similarity <= 0.2: Consider BP5 instead (phenotype inconsistent)
        
        NOTE: This is an educational approximation using local data files.
        For production use, consider integrating with HPO semantic similarity APIs.
        """
        result = {
            'applies': False, 
            'strength': 'Supporting', 
            'details': '',
            'confidence': 'low',
            'data_source': 'local_phenotype_db'
        }
        
        gene = variant_data.basic_info.get('gene')
        patient_phenotypes = getattr(variant_data, 'patient_phenotypes', None)
        
        # CRITICAL FIX: Check genetic_data for phenotype_specificity (manual override)
        genetic_data = variant_data.genetic_data or {}
        phenotype_specificity = genetic_data.get('phenotype_specificity')
        
        # Manual override: If explicit phenotype specificity is provided by user
        if phenotype_specificity == 'high' and gene:
            result['applies'] = True
            result['confidence'] = 'high'
            result['data_source'] = 'user_provided'
            result['details'] = f"Patient phenotype highly specific for {gene}-related disease (PP4) [user-specified]"
            return result
        
        # Use PhenotypeMatcher for automated evaluation with local data
        if gene:
            # Get phenotypes from patient_phenotypes or genetic_data
            if not patient_phenotypes:
                patient_phenotypes = genetic_data.get('patient_phenotypes') or genetic_data.get('phenotypes')
            
            if patient_phenotypes:
                # Call PhenotypeMatcher with full result structure
                match_result = self.phenotype_matcher.evaluate_phenotype_match(variant_data, patient_phenotypes)
                
                if match_result:
                    similarity = match_result.get('similarity', 0.0)
                    evidence_code = match_result.get('evidence_code')
                    explanation = match_result.get('explanation', '')
                    
                    # Store match details for reporting
                    result['similarity'] = similarity
                    result['patient_terms'] = list(match_result.get('patient_terms', set()))
                    result['gene_terms'] = list(match_result.get('gene_terms', set()))
                    
                    if evidence_code == 'PP4':
                        result['applies'] = True
                        result['confidence'] = 'high' if similarity >= 0.8 else 'medium'
                        result['details'] = explanation + " (PP4)"
                    elif evidence_code == 'PP4_supporting':
                        result['applies'] = True
                        result['confidence'] = 'medium'
                        result['details'] = explanation + " (PP4_supporting)"
                    elif evidence_code == 'BP5':
                        # PP4 doesn't apply, but note that BP5 might
                        result['applies'] = False
                        result['details'] = f"Phenotype does not match {gene}-related disease. See BP5 evaluation."
                    else:
                        result['details'] = explanation
                else:
                    result['details'] = f"Phenotype matching failed for {gene}"
            else:
                result['details'] = "No patient phenotype data provided for PP4 evaluation"
        else:
            result['details'] = "No gene information available for PP4 evaluation"
        
        return result
    
    def _evaluate_pp5(self, variant_data) -> Dict[str, Any]:
        """Evaluate PP5 - Reputable source with ClinVar API integration."""
        from config.constants import REPUTABLE_SOURCE_REQUIREMENTS
        from datetime import datetime
        
        result = {
            'applies': False, 
            'strength': 'Supporting', 
            'details': '',
            'confidence': 'very_low',
            'data_source': 'none'
        }
        
        # Try API-based ClinVar lookup first
        if self.api_client and self.api_enabled:
            try:
                gene = variant_data.gene
                hgvs = variant_data.hgvs_c
                
                clinvar_result = self.api_client.get_clinvar_classification(
                    gene=gene, 
                    hgvs=hgvs
                )
                
                if clinvar_result and 'classification' in clinvar_result:
                    classification = clinvar_result.get('classification', '').lower()
                    review_status = clinvar_result.get('review_status', '').lower()
                    star_rating = clinvar_result.get('star_rating', 0)
                    last_evaluated = clinvar_result.get('date_last_evaluated', '')
                    
                    # Check for non-pathogenicity classifications (drug response, risk factor, etc.)
                    non_pathogenicity_terms = ['drug response', 'risk factor', 'protective', 
                                              'affects', 'association', 'confers sensitivity',
                                              'other', 'not provided']
                    is_non_pathogenicity = any(term in classification for term in non_pathogenicity_terms)
                    
                    if is_non_pathogenicity:
                        result['details'] = (
                            f"ClinVar classification '{classification}' is not a pathogenicity assessment. "
                            f"This is a {classification} annotation and does not apply to PP5/BP6 criteria."
                        )
                        result['data_source'] = 'clinvar_non_pathogenicity'
                        return result
                    
                    # Check if classification is recent (within 5 years)
                    classification_recent = True
                    if last_evaluated and last_evaluated != 'Unknown':
                        try:
                            # Parse date (format: YYYY/MM/DD HH:MM or YYYY/MM/DD)
                            date_str = last_evaluated.split(' ')[0]  # Remove time if present
                            date_parts = date_str.split('/')
                            if len(date_parts) >= 3:
                                year, month, day = int(date_parts[0]), int(date_parts[1]), int(date_parts[2])
                                eval_date = datetime(year, month, day)
                                age_years = (datetime.now() - eval_date).days / 365.25
                                classification_recent = age_years <= REPUTABLE_SOURCE_REQUIREMENTS['max_age_years']
                        except Exception as e:
                            classification_recent = False
                    
                    # Check if meets reputable source requirements
                    meets_requirements = (
                        star_rating >= REPUTABLE_SOURCE_REQUIREMENTS['min_stars'] and
                        classification_recent
                    )
                    
                    if 'pathogenic' in classification and 'benign' not in classification:
                        if meets_requirements:
                            result['applies'] = True
                            result['confidence'] = 'high' if star_rating >= 3 else 'medium'
                            result['data_source'] = 'clinvar_api'
                            result['details'] = (
                                f"ClinVar reports variant as {classification} "
                                f"({star_rating}★, {review_status}) (PP5)"
                            )
                        else:
                            result['confidence'] = 'low'
                            result['data_source'] = 'clinvar_insufficient_review'
                            result['details'] = (
                                f"ClinVar reports {classification} but does not meet PP5 criteria "
                                f"(stars: {star_rating}, recent: {classification_recent})"
                            )
                    elif 'benign' in classification:
                        result['details'] = f"ClinVar reports variant as {classification} - see BP6"
                    else:
                        result['details'] = f"ClinVar classification: {classification}"
                    
                    return result
            except Exception as e:
                print(f"⚠️  ClinVar API error in PP5: {str(e)}")
        
        # Fallback to manual clinvar_data attribute if API unavailable
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
        """Evaluate BA1 - High allele frequency with gnomAD API."""
        result = {'applies': False, 'strength': 'Stand-alone', 'details': ''}
        
        # Try API-based gnomAD frequency lookup first
        gnomad_af = None
        data_source = 'manual'
        
        if self.api_client and self.api_enabled:
            try:
                # Get genomic coordinates if available
                basic_info = variant_data.basic_info
                chrom = basic_info.get('chromosome')
                pos = basic_info.get('position')
                ref = basic_info.get('ref_allele')
                alt = basic_info.get('alt_allele')
                
                if chrom and pos and ref and alt:
                    freq_result = self.api_client.get_variant_frequency(
                        chrom=chrom, pos=pos, ref=ref, alt=alt
                    )
                    
                    if freq_result and 'allele_frequency' in freq_result:
                        gnomad_af = freq_result.get('allele_frequency')
                        data_source = 'gnomad_api'
            except Exception as e:
                print(f"⚠️  gnomAD API error in BA1: {str(e)}")
        
        # Fallback to manual population_data
        if gnomad_af is None:
            pop_data = variant_data.population_data
            gnomad_af = pop_data.get('gnomad_af')
        
        gene = variant_data.basic_info.get('gene', '').upper()
        
        # Get gene-specific threshold or use default
        gene_thresholds = GENE_SPECIFIC_THRESHOLDS.get(gene, GENE_SPECIFIC_THRESHOLDS['default'])
        ba1_threshold = gene_thresholds.get('BA1', 0.05)
        
        # BA1 applies if variant is common in population
        if gnomad_af is not None and gnomad_af > ba1_threshold:
            result['applies'] = True
            result['details'] = f"Variant common in population (gnomAD AF: {gnomad_af:.4f}, threshold: {ba1_threshold}) [{data_source}]"
        else:
            af_display = f"{gnomad_af:.6f}" if gnomad_af is not None else 'N/A'
            result['details'] = f"Variant not common in population (gnomAD AF: {af_display}, threshold: {ba1_threshold}) [{data_source}]"
        
        return result
    
    def _evaluate_bs1(self, variant_data) -> Dict[str, Any]:
        """Evaluate BS1 - Allele frequency higher than expected for disorder with gnomAD API."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        # Try API-based gnomAD frequency lookup first
        gnomad_af = None
        data_source = 'manual'
        
        if self.api_client and self.api_enabled:
            try:
                # Get genomic coordinates if available
                basic_info = variant_data.basic_info
                chrom = basic_info.get('chromosome')
                pos = basic_info.get('position')
                ref = basic_info.get('ref_allele')
                alt = basic_info.get('alt_allele')
                
                if chrom and pos and ref and alt:
                    freq_result = self.api_client.get_variant_frequency(
                        chrom=chrom, pos=pos, ref=ref, alt=alt
                    )
                    
                    if freq_result and 'allele_frequency' in freq_result:
                        gnomad_af = freq_result.get('allele_frequency')
                        data_source = 'gnomad_api'
            except Exception as e:
                print(f"⚠️  gnomAD API error in BS1: {str(e)}")
        
        # Fallback to manual population_data
        if gnomad_af is None:
            pop_data = variant_data.population_data
            gnomad_af = pop_data.get('gnomad_af')
        else:
            pop_data = variant_data.population_data
        
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
        
        # Get gene thresholds for BA1 check
        gene_thresholds = GENE_SPECIFIC_THRESHOLDS.get(gene, GENE_SPECIFIC_THRESHOLDS['default'])
        
        # BS1 applies if variant frequency is higher than expected for disorder
        # Must be significantly higher (not just slightly above)
        if gnomad_af is not None:
            # Check if BA1 threshold is already exceeded (>5%)
            ba1_threshold = gene_thresholds.get('BA1', 0.05) if expected_max_af is None else 0.05
            
            if gnomad_af > ba1_threshold:
                # BA1 takes precedence - don't apply BS1
                result['applies'] = False
                result['details'] = f"Variant frequency exceeds BA1 threshold (gnomAD AF: {gnomad_af:.4f}) [{data_source}]"
            elif gnomad_af > bs1_threshold:
                result['applies'] = True
                result['details'] = f"Variant frequency higher than expected for disorder (gnomAD AF: {gnomad_af:.4f}, expected max: {bs1_threshold}) [{data_source}]"
            else:
                af_display = f"{gnomad_af:.6f}" if gnomad_af > 0 else "0.0"
                result['details'] = f"Variant frequency not higher than expected (gnomAD AF: {af_display}, expected max: {bs1_threshold}) [{data_source}]"
        else:
            result['details'] = f"No population frequency data available [{data_source}]"
        
        return result
    
    def _evaluate_bs2(self, variant_data) -> Dict[str, Any]:
        """Evaluate BS2 - Observed in a healthy adult individual for a recessive disorder."""
        result = {'applies': False, 'strength': 'Strong', 'details': ''}
        
        # **QUICK WIN #2: Check gnomAD homozygote counts for recessive disorders**
        population_data = variant_data.population_data or {}
        ac_hom = population_data.get('ac_hom')
        gnomad_af = population_data.get('gnomad_af')
        
        # If variant observed as homozygous in gnomAD, apply BS2 for recessive disorders
        # Rationale: Healthy homozygotes indicate variant is unlikely causative
        if ac_hom is not None and ac_hom > 0:
            result['applies'] = True
            result['details'] = f"✓ BS2 applies: Variant observed in {ac_hom} healthy homozygotes in gnomAD (recessive disorder)"
            result['source'] = 'gnomAD homozygote count'
            result['manual_review'] = True
            result['guidance'] = "BS2 applies for recessive disorders. For dominant, consider if homozygous state is more severe."
            return result
        
        # CRITICAL FIX: Check genetic_data for observed_in_healthy
        genetic_data = variant_data.genetic_data or {}
        observed_in_healthy = genetic_data.get('observed_in_healthy', False)
        zygosity = genetic_data.get('zygosity', '')
        age_observed = genetic_data.get('age_observed')
        disease_onset_age = genetic_data.get('disease_onset_age')
        penetrance = genetic_data.get('penetrance', 1.0)
        
        # BS2 requires: healthy individual + appropriate zygosity + age considerations
        if not observed_in_healthy:
            result['details'] = "Variant not observed in healthy individuals or homozygotes"
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
            # **QUICK WIN #1: Use gnomAD missense constraint data**
            # Try to get gene constraint data from API
            constraint_data = None
            if self.api_client and gene:
                try:
                    constraint_data = self.api_client.get_gene_constraint(gene)
                except Exception as e:
                    print(f"⚠️ Could not fetch constraint data: {e}")
            
            # Check if gene is missense tolerant (suggests LOF is mechanism)
            if constraint_data:
                oe_mis = constraint_data.get('oe_mis')
                mis_classification = constraint_data.get('mis_classification', 'uncertain')
                is_lof_intolerant = constraint_data.get('is_lof_intolerant', False)
                
                # Apply BP1 if:
                # 1. Gene is missense tolerant (oe_mis > 1.0)
                # 2. Gene is LOF intolerant (suggests LOF is disease mechanism)
                if mis_classification == 'missense_tolerant' and is_lof_intolerant:
                    result['applies'] = True
                    result['details'] = f"✓ BP1 applies: {gene} is missense tolerant (oe_mis={oe_mis:.3f}) and LOF intolerant (truncating is mechanism)"
                    result['source'] = 'gnomAD constraint'
                    return result
                elif mis_classification == 'missense_constrained':
                    result['details'] = f"Gene {gene} is missense constrained, BP1 does not apply"
                    return result
            
            # Fallback: Define genes where primarily truncating variants cause disease
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
                result['details'] = f"Gene {gene} not recognized as primarily truncating mechanism (and no constraint data available)"
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
        """
        Evaluate BP4 - Computational evidence suggests no impact.
        
        For missense variants, uses the MissenseEvaluator composite score.
        For other variants, falls back to simple predictor counting.
        
        Returns dict with:
            - applies: bool
            - strength: 'Supporting', 'Moderate', or 'Strong'
            - details: str explanation
            - composite_score: float (for missense variants)
        """
        result = {'applies': False, 'strength': 'Supporting', 'details': ''}
        
        variant_type = variant_data.basic_info.get('variant_type', '').lower()
        
        # Use MissenseEvaluator for missense variants
        if variant_type == 'missense':
            return self._evaluate_missense_pp3_bp4(variant_data, direction='benign')
        
        # Fallback for non-missense variants: simple predictor counting
        insilico_data = variant_data.insilico_data or {}
        
        benign_count = 0
        total_count = 0
        
        # Check prediction labels first
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
            'revel': 0.5,
            'sift': 0.05,
            'polyphen2': 0.5,
            'mutation_taster': 0.5,
            'dann': 0.5,
            'fathmm': 0.5,
        }
        
        for predictor, threshold in predictor_thresholds.items():
            score = insilico_data.get(predictor)
            if score is not None:
                total_count += 1
                if predictor in ['sift', 'fathmm']:
                    if score > threshold:
                        benign_count += 1
                else:
                    if score < threshold:
                        benign_count += 1
        
        if total_count > 0:
            benign_ratio = benign_count / total_count
            if benign_ratio >= 0.7:
                result['applies'] = True
                result['details'] = f"Computational evidence suggests no impact ({benign_count}/{total_count} predictors)"
            else:
                result['details'] = f"Computational evidence inconclusive ({benign_count}/{total_count} predictors suggest benign)"
        else:
            result['details'] = "No computational prediction scores available"
        
        return result
    
    def _evaluate_bp5(self, variant_data) -> Dict[str, Any]:
        """
        Evaluate BP5 - Variant found in a case with an alternate molecular basis for disease,
        OR phenotype is inconsistent with gene-associated disease.
        
        BP5 can apply in two scenarios:
        1. Traditional: Alternate pathogenic variant explains the phenotype
        2. Phenotype-based: Patient phenotypes show very poor overlap with gene-associated disease
        
        Phenotype-Based BP5 Threshold:
            - similarity <= 0.2: BP5 applies (phenotype clearly inconsistent)
        
        NOTE: This is an educational approximation using local phenotype data.
        """
        result = {
            'applies': False, 
            'strength': 'Supporting', 
            'details': '',
            'confidence': 'low',
            'data_source': 'none'
        }
        
        gene = variant_data.basic_info.get('gene')
        variant_name = variant_data.basic_info.get('variant_name', 'variant')
        
        # SCENARIO 1: Check genetic_data for alternate cause information (traditional BP5)
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
                result['confidence'] = 'high'
                result['data_source'] = 'alternate_variant'
                result['details'] = f"Alternate molecular basis identified ({alternate_gene} {alternate_classification}) explaining phenotype (BP5)"
                return result
        
        # SCENARIO 2: Phenotype-based BP5 using PhenotypeMatcher
        # Check if patient phenotypes are inconsistent with gene-associated disease
        if gene:
            patient_phenotypes = getattr(variant_data, 'patient_phenotypes', None)
            if not patient_phenotypes:
                patient_phenotypes = genetic_data.get('patient_phenotypes') or genetic_data.get('phenotypes')
            
            if patient_phenotypes:
                # Use PhenotypeMatcher to check for phenotype inconsistency
                match_result = self.phenotype_matcher.evaluate_phenotype_match(variant_data, patient_phenotypes)
                
                if match_result:
                    similarity = match_result.get('similarity', 0.0)
                    evidence_code = match_result.get('evidence_code')
                    explanation = match_result.get('explanation', '')
                    
                    # Store match details
                    result['similarity'] = similarity
                    result['patient_terms'] = list(match_result.get('patient_terms', set()))
                    result['gene_terms'] = list(match_result.get('gene_terms', set()))
                    
                    if evidence_code == 'BP5':
                        result['applies'] = True
                        result['confidence'] = 'medium'
                        result['data_source'] = 'phenotype_mismatch'
                        result['details'] = explanation + " (BP5 - phenotype-based)"
                        return result
                    else:
                        # Phenotypes match or inconclusive - BP5 doesn't apply from phenotype
                        result['details'] = f"Phenotype consistent with {gene}-related disease (similarity: {similarity:.0%}). "
        
        # SCENARIO 3: Interactive or test mode for traditional BP5
        if gene:
            if self.test_mode:
                if not result['details']:
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
        """Evaluate BP6 - Reputable source reports variant as benign with ClinVar API."""
        from config.constants import REPUTABLE_SOURCE_REQUIREMENTS
        from datetime import datetime
        
        result = {
            'applies': False, 
            'strength': 'Supporting', 
            'details': '',
            'confidence': 'very_low',
            'data_source': 'none'
        }
        
        # Try API-based ClinVar lookup first
        if self.api_client and self.api_enabled:
            try:
                gene = variant_data.gene
                hgvs = variant_data.hgvs_c
                
                clinvar_result = self.api_client.get_clinvar_classification(
                    gene=gene, 
                    hgvs=hgvs
                )
                
                if clinvar_result and 'classification' in clinvar_result:
                    classification = clinvar_result.get('classification', '').lower()
                    review_status = clinvar_result.get('review_status', '').lower()
                    star_rating = clinvar_result.get('star_rating', 0)
                    last_evaluated = clinvar_result.get('date_last_evaluated', '')
                    
                    # Check for non-pathogenicity classifications (drug response, risk factor, etc.)
                    non_pathogenicity_terms = ['drug response', 'risk factor', 'protective', 
                                              'affects', 'association', 'confers sensitivity',
                                              'other', 'not provided']
                    is_non_pathogenicity = any(term in classification for term in non_pathogenicity_terms)
                    
                    if is_non_pathogenicity:
                        result['details'] = (
                            f"ClinVar classification '{classification}' is not a pathogenicity assessment. "
                            f"This is a {classification} annotation and does not apply to PP5/BP6 criteria."
                        )
                        result['data_source'] = 'clinvar_non_pathogenicity'
                        return result
                    
                    # Check if classification is recent (within 5 years)
                    classification_recent = True
                    if last_evaluated and last_evaluated != 'Unknown':
                        try:
                            # Parse date (format: YYYY/MM/DD HH:MM or YYYY/MM/DD)
                            date_str = last_evaluated.split(' ')[0]  # Remove time if present
                            date_parts = date_str.split('/')
                            if len(date_parts) >= 3:
                                year, month, day = int(date_parts[0]), int(date_parts[1]), int(date_parts[2])
                                eval_date = datetime(year, month, day)
                                age_years = (datetime.now() - eval_date).days / 365.25
                                classification_recent = age_years <= REPUTABLE_SOURCE_REQUIREMENTS['max_age_years']
                        except Exception as e:
                            classification_recent = False
                    
                    # Check if meets reputable source requirements
                    meets_requirements = (
                        star_rating >= REPUTABLE_SOURCE_REQUIREMENTS['min_stars'] and
                        classification_recent
                    )
                    
                    if 'benign' in classification and 'pathogenic' not in classification:
                        if meets_requirements:
                            result['applies'] = True
                            result['confidence'] = 'high' if star_rating >= 3 else 'medium'
                            result['data_source'] = 'clinvar_api'
                            result['details'] = (
                                f"ClinVar reports variant as {classification} "
                                f"({star_rating}★, {review_status}) (BP6)"
                            )
                        else:
                            result['confidence'] = 'low'
                            result['data_source'] = 'clinvar_insufficient_review'
                            result['details'] = (
                                f"ClinVar reports {classification} but does not meet BP6 criteria "
                                f"(stars: {star_rating}, recent: {classification_recent})"
                            )
                    elif 'pathogenic' in classification:
                        result['details'] = f"ClinVar reports variant as {classification} - see PP5"
                    else:
                        result['details'] = f"ClinVar classification: {classification}"
                    
                    return result
            except Exception as e:
                print(f"⚠️  ClinVar API error in BP6: {str(e)}")
        
        # Fallback to manual clinvar_data attribute if API unavailable
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
    
    # =========================================================================
    # Interactive Evidence Collection Methods
    # =========================================================================
    def collect_manual_evidence(self, variant_data=None, criteria: Optional[List[str]] = None) -> 'ManualEvidence':
        """
        Collect literature-based evidence through interactive prompts.
        
        This method initializes the InteractiveEvidenceCollector and guides the user
        through questions about PS3/BS3, PS4, PP1/BS4, PS1/PM5, and PP5/BP6 criteria.
        
        Args:
            variant_data: Optional VariantData object for context
            criteria: Optional list of specific criteria to collect 
                     (e.g., ['PS3_BS3', 'PP1_BS4']). If None, collects all.
        
        Returns:
            ManualEvidence object with collected evidence codes and explanations
        
        Example:
            >>> evaluator = EvidenceEvaluator()
            >>> manual_ev = evaluator.collect_manual_evidence()
            >>> print(manual_ev.codes)
            ['PS3_moderate', 'PP1_supporting']
        """
        from core.interactive_evidence import InteractiveEvidenceCollector, ManualEvidence
        
        # Initialize collector if not already done
        if self._interactive_collector is None:
            self._interactive_collector = InteractiveEvidenceCollector(show_prompts=True)
        
        # Collect evidence
        if criteria:
            self._manual_evidence = self._interactive_collector.collect_selective(criteria, variant_data)
        else:
            self._manual_evidence = self._interactive_collector.collect_all(variant_data)
        
        return self._manual_evidence
    
    def set_manual_evidence(self, manual_evidence: 'ManualEvidence') -> None:
        """
        Set pre-collected manual evidence (useful for testing or batch processing).
        
        Args:
            manual_evidence: ManualEvidence object with codes and explanations
        """
        self._manual_evidence = manual_evidence
    
    def get_manual_evidence(self) -> Optional['ManualEvidence']:
        """
        Get the currently stored manual evidence.
        
        Returns:
            ManualEvidence object or None if not collected
        """
        return self._manual_evidence
    
    def merge_manual_with_automated(self, automated_results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Merge manually collected evidence with automated evaluation results.
        
        This method takes the automated ACMG criteria evaluation results and
        supplements them with user-provided literature-based evidence.
        
        Manual evidence takes precedence for criteria that require human judgment
        (PS3, BS3, PS4, PP1, BS4, PS1, PM5, PP5, BP6).
        
        Args:
            automated_results: Results from evaluate_all_criteria()
        
        Returns:
            Updated results dict with manual evidence integrated
        
        Example:
            >>> evaluator = EvidenceEvaluator()
            >>> auto_results = evaluator.evaluate_all_criteria(variant_data)
            >>> manual_ev = evaluator.collect_manual_evidence()
            >>> merged = evaluator.merge_manual_with_automated(auto_results)
        """
        if self._manual_evidence is None or not self._manual_evidence.has_evidence():
            return automated_results
        
        # Create a copy of results to avoid mutating original
        merged = {
            'pathogenic_criteria': dict(automated_results.get('pathogenic_criteria', {})),
            'benign_criteria': dict(automated_results.get('benign_criteria', {})),
            'applied_criteria': dict(automated_results.get('applied_criteria', {})),
            'evidence_details': dict(automated_results.get('evidence_details', {})),
            'vampp_score': automated_results.get('vampp_score'),
            'statistical_tests': automated_results.get('statistical_tests', {}),
            'manual_evidence': {
                'codes': self._manual_evidence.codes,
                'explanations': self._manual_evidence.explanations
            }
        }
        
        # Map manual evidence codes to ACMG criteria
        code_to_criterion = {
            'PS3': 'PS3', 'PS3_strong': 'PS3', 'PS3_moderate': 'PS3', 'PS3_supporting': 'PS3',
            'BS3': 'BS3', 'BS3_strong': 'BS3', 'BS3_moderate': 'BS3', 'BS3_supporting': 'BS3',
            'PS4': 'PS4', 'PS4_strong': 'PS4', 'PS4_moderate': 'PS4', 'PS4_supporting': 'PS4',
            'PP1': 'PP1', 'PP1_strong': 'PP1', 'PP1_moderate': 'PP1', 'PP1_supporting': 'PP1',
            'BS4': 'BS4', 'BS4_strong': 'BS4', 'BS4_moderate': 'BS4', 'BS4_supporting': 'BS4',
            'PS1': 'PS1', 'PS1_strong': 'PS1', 'PS1_moderate': 'PS1',
            'PM5': 'PM5', 'PM5_moderate': 'PM5', 'PM5_supporting': 'PM5',
            'PP5': 'PP5', 'PP5_supporting': 'PP5',
            'BP6': 'BP6', 'BP6_supporting': 'BP6',
        }
        
        strength_mapping = {
            'strong': 'Strong',
            'moderate': 'Moderate',
            'supporting': 'Supporting',
        }
        
        # Process each manual evidence code
        for code in self._manual_evidence.codes:
            # Determine base criterion and strength
            parts = code.split('_')
            base_criterion = parts[0]
            strength_suffix = parts[1] if len(parts) > 1 else None
            
            # Get the ACMG criterion name
            criterion = code_to_criterion.get(code) or code_to_criterion.get(base_criterion)
            if not criterion:
                continue
            
            # Determine strength
            if strength_suffix and strength_suffix.lower() in strength_mapping:
                strength = strength_mapping[strength_suffix.lower()]
            else:
                # Default strengths per criterion
                default_strengths = {
                    'PS3': 'Strong', 'BS3': 'Strong',
                    'PS4': 'Strong',
                    'PP1': 'Supporting', 'BS4': 'Strong',
                    'PS1': 'Strong', 'PM5': 'Moderate',
                    'PP5': 'Supporting', 'BP6': 'Supporting',
                }
                strength = default_strengths.get(base_criterion, 'Supporting')
            
            # Build the result entry
            explanation = self._manual_evidence.explanations.get(code, 'User-provided evidence')
            
            result_entry = {
                'applies': True,
                'strength': strength,
                'details': explanation,
                'data_source': 'manual_evidence',
                'confidence': 'user_provided'
            }
            
            # Determine if pathogenic or benign criterion
            is_pathogenic = base_criterion.startswith(('PVS', 'PS', 'PM', 'PP'))
            is_benign = base_criterion.startswith(('BA', 'BS', 'BP'))
            
            if is_pathogenic:
                merged['pathogenic_criteria'][criterion] = result_entry
                merged['applied_criteria'][criterion] = result_entry
            elif is_benign:
                merged['benign_criteria'][criterion] = result_entry
                merged['applied_criteria'][criterion] = result_entry
        
        # Add manual evidence details to evidence_details
        merged['evidence_details']['manual_evidence'] = {
            'collected': True,
            'codes': self._manual_evidence.codes,
            'count': len(self._manual_evidence.codes)
        }
        
        return merged
