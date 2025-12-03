"""
Missense Variant Evaluator
==========================

Provides comprehensive evaluation of missense variants by combining multiple
evidence sources into a composite score that can be mapped to ACMG PP3/BP4
evidence categories.

The composite score (0.0-1.0) integrates:
- Conservation scores (phyloP, phastCons, GERP++)
- Functional predictions (REVEL, CADD, AlphaMissense, etc.)
- Structural impact estimates
- Functional domain context
- Population frequency context

IMPORTANT: This is a research/educational implementation. The weights and
thresholds are approximations and should be validated against known variant
datasets before clinical use.

Multi-Source Data Flow (Jan 2025):
- Predictor scores are pre-fetched by EvidenceEvaluator._fetch_external_data()
- This class acts as a "pure interpreter" of pre-fetched data
- Falls back to variant_data.insilico_data for backward compatibility
- Gracefully handles partial data with weight renormalization

Author: Can Sevilmiş
License: MIT License
"""

from typing import Dict, Any, Optional, Tuple
from config.predictors import (
    INSILICO_WEIGHTS,
    INSILICO_THRESHOLDS,
    MISSENSE_COMPOSITE_THRESHOLDS,
    MISSENSE_SCORE_WEIGHTS,
    INVERTED_PREDICTORS,
    MIN_PREDICTORS_FOR_COMPOSITE,
)


class MissenseEvaluator:
    """
    Evaluates missense variants using multiple evidence sources.
    
    Produces a composite score (0.0-1.0) where:
    - Higher scores indicate more likely damaging/pathogenic
    - Lower scores indicate more likely benign/tolerated
    
    The score can be mapped to ACMG evidence categories:
    - PP3 (pathogenic computational evidence) for high scores
    - BP4 (benign computational evidence) for low scores
    """
    
    def __init__(self):
        """Initialize the missense evaluator with data sources."""
        self.domain_regions = self._load_domain_regions()
        self.conservation_scores = self._load_conservation_data()
        self.structural_data = self._load_structural_data()
        self.thresholds = MISSENSE_COMPOSITE_THRESHOLDS
        self.weights = MISSENSE_SCORE_WEIGHTS
    
    def evaluate_missense_variant(self, variant_data) -> Dict[str, Any]:
        """
        Comprehensive missense variant evaluation.
        
        Args:
            variant_data: VariantData object containing variant information
            
        Returns:
            Dict containing:
                - composite_score: float (0.0-1.0)
                - sub_scores: dict of individual component scores
                - evidence_category: str (e.g., 'PP3_moderate', 'BP4_supporting')
                - confidence: str ('high', 'medium', 'low')
                - details: str with human-readable explanation
        """
        # Calculate individual sub-scores
        sub_scores = {
            'conservation': self._calculate_conservation_score(variant_data),
            'functional': self._calculate_functional_score(variant_data),
            'structural': self._calculate_structural_impact(variant_data),
            'domain': self._calculate_domain_impact(variant_data),
            'population': self._calculate_population_context(variant_data)
        }
        
        # Generate composite score and evidence category
        result = self._generate_composite_score(sub_scores)
        result['sub_scores'] = sub_scores
        
        return result
    
    def _calculate_conservation_score(self, variant_data) -> float:
        """
        Calculate conservation score from phyloP, phastCons, and GERP++.
        
        Higher scores indicate more conserved positions (more likely damaging).
        
        Args:
            variant_data: VariantData object
            
        Returns:
            float: Normalized conservation score (0.0-1.0)
        """
        insilico = getattr(variant_data, 'insilico_data', {}) or {}
        
        # Try multiple possible key names for each score
        phylop = (insilico.get('phylop') or 
                  insilico.get('phylop_score') or
                  insilico.get('phylop100way') or
                  getattr(variant_data, 'phylop_score', None))
        
        phastcons = (insilico.get('phastcons') or
                     insilico.get('phastcons_score') or
                     insilico.get('phastcons100way') or
                     getattr(variant_data, 'phastcons_score', None))
        
        gerp = (insilico.get('gerp') or
                insilico.get('gerp_pp') or
                insilico.get('gerp_score') or
                getattr(variant_data, 'gerp_score', None))
        
        # Calculate weighted average of available scores
        total_weight = 0.0
        weighted_sum = 0.0
        
        if phylop is not None:
            normalized = self._normalize_phylop(phylop)
            weighted_sum += 0.4 * normalized
            total_weight += 0.4
        
        if phastcons is not None:
            normalized = self._normalize_phastcons(phastcons)
            weighted_sum += 0.3 * normalized
            total_weight += 0.3
        
        if gerp is not None:
            normalized = self._normalize_gerp(gerp)
            weighted_sum += 0.3 * normalized
            total_weight += 0.3
        
        # Return weighted average, or 0.5 (neutral) if no data
        if total_weight > 0:
            return weighted_sum / total_weight
        return 0.5  # Neutral when no conservation data available
    
    def _calculate_functional_score(self, variant_data) -> float:
        """
        Calculate functional impact score from in silico predictors.
        
        Combines REVEL, CADD, AlphaMissense, SIFT, PolyPhen2, etc.
        Higher scores indicate more likely damaging.
        
        Multi-Source Data Flow:
        - First tries typed predictor_scores (from multi-source API)
        - Falls back to insilico_data dict for backward compatibility
        - Renormalizes weights when predictors are missing
        
        Args:
            variant_data: VariantData object
            
        Returns:
            float: Normalized functional score (0.0-1.0)
        """
        scores = []
        weights = []
        sources = []
        
        # Try to use typed predictor_scores first (multi-source API data)
        predictor_scores = getattr(variant_data, 'predictor_scores', None)
        if predictor_scores:
            for predictor_name, weight in INSILICO_WEIGHTS.items():
                score_obj = predictor_scores.get(predictor_name)
                if score_obj is not None:
                    # Handle PredictorScore dataclass
                    if hasattr(score_obj, 'value') and score_obj.value is not None:
                        normalized = self._normalize_predictor_score(
                            predictor_name, 
                            score_obj.value,
                            is_inverted=getattr(score_obj, 'is_inverted', False)
                        )
                        scores.append(normalized)
                        weights.append(weight)
                        source = getattr(score_obj, 'source', 'unknown')
                        sources.append(f"{predictor_name}({source})")
        
        # Fall back to insilico_data if no typed scores available
        if not scores:
            insilico = getattr(variant_data, 'insilico_data', {}) or {}
            scores, weights, sources = self._extract_from_insilico_dict(insilico)
        
        # Calculate weighted average with weight renormalization
        if scores and weights:
            total_weight = sum(weights)
            weighted_sum = sum(s * w for s, w in zip(scores, weights))
            
            # Log available predictors for debugging
            if len(scores) < MIN_PREDICTORS_FOR_COMPOSITE:
                # Warn about low predictor count
                pass
            
            return weighted_sum / total_weight
        
        return 0.5  # Neutral when no functional data available
    
    def _normalize_predictor_score(
        self, 
        predictor: str, 
        value: float, 
        is_inverted: bool = False
    ) -> float:
        """
        Normalize a predictor score to 0-1 range where higher = more pathogenic.
        
        Args:
            predictor: Predictor name
            value: Raw score value
            is_inverted: Whether the predictor uses inverted scoring
            
        Returns:
            float: Normalized score (0.0-1.0)
        """
        if predictor == 'cadd_phred':
            # CADD phred: 0-99, use 40 as practical max
            return min(float(value) / 40.0, 1.0)
        elif predictor in INVERTED_PREDICTORS or is_inverted:
            if predictor == 'sift':
                # SIFT: 0-1, lower = more damaging → invert
                return 1.0 - float(value)
            elif predictor == 'fathmm':
                # FATHMM: roughly -10 to +10, negative = more damaging
                return max(0.0, min(1.0, (5.0 - float(value)) / 10.0))
            else:
                # Generic inversion for 0-1 scores
                return 1.0 - float(value)
        else:
            # Scores already 0-1 with higher = more damaging
            return min(max(float(value), 0.0), 1.0)
    
    def _extract_from_insilico_dict(self, insilico: dict) -> tuple:
        """
        Extract predictor scores from legacy insilico_data dict.
        
        This provides backward compatibility with the original data format.
        
        Args:
            insilico: Dictionary with predictor scores
            
        Returns:
            Tuple of (scores, weights, sources) lists
        """
        scores = []
        weights = []
        sources = []
        
        # REVEL (0-1, higher = more damaging)
        revel = insilico.get('revel') or insilico.get('revel_score')
        if revel is not None:
            scores.append(float(revel))
            weights.append(INSILICO_WEIGHTS.get('revel', 0.25))
            sources.append('revel(insilico_data)')
        
        # CADD phred (0-99, higher = more damaging)
        cadd = insilico.get('cadd_phred') or insilico.get('cadd')
        if cadd is not None:
            normalized_cadd = min(float(cadd) / 40.0, 1.0)
            scores.append(normalized_cadd)
            weights.append(INSILICO_WEIGHTS.get('cadd_phred', 0.20))
            sources.append('cadd(insilico_data)')
        
        # AlphaMissense (0-1, higher = more damaging)
        alphamissense = insilico.get('alphamissense') or insilico.get('alphamissense_score')
        if alphamissense is not None:
            scores.append(float(alphamissense))
            weights.append(INSILICO_WEIGHTS.get('alphamissense', 0.15))
            sources.append('alphamissense(insilico_data)')
        
        # SIFT (0-1, LOWER = more damaging, so we invert)
        sift = insilico.get('sift') or insilico.get('sift_score')
        if sift is not None:
            inverted_sift = 1.0 - float(sift)
            scores.append(inverted_sift)
            weights.append(INSILICO_WEIGHTS.get('sift', 0.10))
            sources.append('sift(insilico_data)')
        
        # PolyPhen2 (0-1, higher = more damaging)
        polyphen = (insilico.get('polyphen2') or 
                    insilico.get('polyphen2_score') or
                    insilico.get('polyphen'))
        if polyphen is not None:
            scores.append(float(polyphen))
            weights.append(INSILICO_WEIGHTS.get('polyphen2', 0.10))
            sources.append('polyphen2(insilico_data)')
        
        # MetaSVM (0-1, higher = more damaging)
        metasvm = insilico.get('metasvm') or insilico.get('metasvm_score')
        if metasvm is not None:
            scores.append(float(metasvm))
            weights.append(INSILICO_WEIGHTS.get('metasvm', 0.10))
            sources.append('metasvm(insilico_data)')
        
        # VEST4 (0-1, higher = more damaging)
        vest4 = insilico.get('vest4') or insilico.get('vest4_score')
        if vest4 is not None:
            scores.append(float(vest4))
            weights.append(INSILICO_WEIGHTS.get('vest4', 0.05))
            sources.append('vest4(insilico_data)')
        
        # FATHMM (negative = more damaging)
        fathmm = insilico.get('fathmm') or insilico.get('fathmm_score')
        if fathmm is not None:
            # Normalize: roughly -10 to +10 → 0-1 (invert so higher = more damaging)
            normalized_fathmm = max(0.0, min(1.0, (5.0 - float(fathmm)) / 10.0))
            scores.append(normalized_fathmm)
            weights.append(INSILICO_WEIGHTS.get('fathmm', 0.05))
            sources.append('fathmm(insilico_data)')
        
        return scores, weights, sources
    
    def _calculate_structural_impact(self, variant_data) -> float:
        """
        Calculate structural impact score based on protein structure context.
        
        Considers:
        - Amino acid physicochemical property changes (Grantham distance)
        - Secondary structure context
        - Solvent accessibility
        
        Args:
            variant_data: VariantData object
            
        Returns:
            float: Normalized structural impact score (0.0-1.0)
        """
        basic_info = getattr(variant_data, 'basic_info', {}) or {}
        insilico = getattr(variant_data, 'insilico_data', {}) or {}
        
        scores = []
        
        # Check for explicit structural scores in insilico_data
        struct_score = insilico.get('structural_impact') or insilico.get('protein_stability')
        if struct_score is not None:
            return min(max(float(struct_score), 0.0), 1.0)
        
        # Estimate from amino acid change if available
        aa_change = basic_info.get('amino_acid_change') or basic_info.get('hgvs_p')
        if aa_change:
            grantham_score = self._estimate_grantham_impact(aa_change)
            if grantham_score is not None:
                scores.append(grantham_score)
        
        # Check for secondary structure disruption
        secondary_struct = insilico.get('secondary_structure')
        if secondary_struct:
            # Helix/sheet disruptions are more impactful
            if secondary_struct in ['helix', 'sheet']:
                scores.append(0.7)
            else:
                scores.append(0.4)
        
        if scores:
            return sum(scores) / len(scores)
        
        return 0.5  # Neutral when no structural data available
    
    def _calculate_domain_impact(self, variant_data) -> float:
        """
        Calculate domain impact score based on functional domain context.
        
        Variants in critical functional domains are more likely damaging.
        
        Args:
            variant_data: VariantData object
            
        Returns:
            float: Normalized domain impact score (0.0-1.0)
        """
        functional_data = getattr(variant_data, 'functional_data', {}) or {}
        basic_info = getattr(variant_data, 'basic_info', {}) or {}
        
        # Check for explicit domain flags
        in_hotspot = functional_data.get('in_hotspot', False)
        in_functional_domain = functional_data.get('in_functional_domain', False)
        domain_name = functional_data.get('domain_name', '')
        benign_variation = functional_data.get('benign_variation_in_domain', True)
        
        # Hotspot = very high impact
        if in_hotspot:
            return 0.9
        
        # Functional domain without benign variation = high impact
        if in_functional_domain and not benign_variation:
            return 0.8
        
        # Functional domain with benign variation = moderate impact
        if in_functional_domain:
            return 0.6
        
        # Check for critical domain keywords
        critical_domains = ['active_site', 'catalytic', 'dna_binding', 'atp_binding',
                           'kinase', 'ring', 'brct', 'zinc_finger']
        if domain_name:
            domain_lower = domain_name.lower()
            for critical in critical_domains:
                if critical in domain_lower:
                    return 0.75
        
        return 0.5  # Neutral when no domain information
    
    def _calculate_population_context(self, variant_data) -> float:
        """
        Calculate population context score based on allele frequency.
        
        Rare variants are more likely pathogenic; common variants are more
        likely benign. This provides additional context beyond BA1/BS1/PM2.
        
        Multi-Source Data Flow:
        - First tries typed population_stats (from multi-source API)
        - Falls back to population_data dict for backward compatibility
        
        Args:
            variant_data: VariantData object
            
        Returns:
            float: Normalized population score (0.0-1.0)
                   Higher = more likely damaging (rare)
                   Lower = more likely benign (common)
        """
        af = None
        
        # Try to use typed population_stats first (multi-source API data)
        population_stats = getattr(variant_data, 'population_stats', None)
        if population_stats:
            # Use maximum AF across all population sources (conservative)
            for source, stats in population_stats.items():
                if hasattr(stats, 'af') and stats.af is not None:
                    source_af = stats.af
                    if af is None or source_af > af:
                        af = source_af
                elif isinstance(stats, dict) and stats.get('af') is not None:
                    source_af = stats['af']
                    if af is None or source_af > af:
                        af = source_af
        
        # Fall back to population_data dict
        if af is None:
            population_data = getattr(variant_data, 'population_data', {}) or {}
            af = (population_data.get('gnomad_af_popmax') or
                  population_data.get('gnomad_af') or
                  population_data.get('exac_af'))
        
        if af is None:
            return 0.5  # Neutral when no frequency data
        
        af = float(af)
        
        # Convert frequency to pathogenicity score
        # Absent/very rare = high score (likely damaging)
        # Common = low score (likely benign)
        if af == 0:
            return 0.9  # Absent from population
        elif af < 0.00001:  # < 0.001%
            return 0.85
        elif af < 0.0001:   # < 0.01%
            return 0.75
        elif af < 0.001:    # < 0.1%
            return 0.6
        elif af < 0.01:     # < 1%
            return 0.4
        elif af < 0.05:     # < 5%
            return 0.2
        else:               # >= 5%
            return 0.1
    
    def _generate_composite_score(self, sub_scores: Dict[str, float]) -> Dict[str, Any]:
        """
        Generate composite score from sub-scores and map to evidence category.
        
        Uses weighted combination of sub-scores with explicit, documented weights.
        
        Args:
            sub_scores: Dict with keys 'conservation', 'functional', 'structural',
                       'domain', 'population' and float values (0.0-1.0)
        
        Returns:
            Dict containing:
                - composite_score: float (0.0-1.0)
                - evidence_category: str (e.g., 'PP3_moderate', 'BP4_supporting', 'neutral')
                - strength: str ('strong', 'moderate', 'supporting', None)
                - direction: str ('pathogenic', 'benign', 'neutral')
                - confidence: str ('high', 'medium', 'low')
                - details: str with explanation
        """
        # Calculate weighted composite score
        weighted_sum = 0.0
        total_weight = 0.0
        available_components = []
        
        for component, score in sub_scores.items():
            weight = self.weights.get(component, 0.0)
            if score is not None and weight > 0:
                weighted_sum += score * weight
                total_weight += weight
                available_components.append(component)
        
        # Calculate final composite score
        if total_weight > 0:
            composite_score = weighted_sum / total_weight
        else:
            composite_score = 0.5  # Neutral if no data
        
        # Determine confidence based on data availability
        if len(available_components) >= 4:
            confidence = 'high'
        elif len(available_components) >= 2:
            confidence = 'medium'
        else:
            confidence = 'low'
        
        # Map to evidence category using thresholds
        evidence_category, strength, direction = self._map_score_to_evidence(composite_score)
        
        # Generate details string
        component_details = ', '.join([f"{c}={sub_scores[c]:.2f}" for c in available_components])
        details = (f"Composite missense score: {composite_score:.3f} "
                   f"({evidence_category}). Components: {component_details}")
        
        return {
            'composite_score': composite_score,
            'evidence_category': evidence_category,
            'strength': strength,
            'direction': direction,
            'confidence': confidence,
            'details': details
        }
    
    def _map_score_to_evidence(self, score: float) -> Tuple[str, Optional[str], str]:
        """
        Map composite score to ACMG evidence category.
        
        Args:
            score: Composite score (0.0-1.0)
            
        Returns:
            Tuple of (evidence_category, strength, direction)
        """
        t = self.thresholds
        
        if score >= t['PP3_strong']:
            return ('PP3_strong', 'strong', 'pathogenic')
        elif score >= t['PP3_moderate']:
            return ('PP3_moderate', 'moderate', 'pathogenic')
        elif score >= t['PP3_supporting']:
            return ('PP3_supporting', 'supporting', 'pathogenic')
        elif score <= t['BP4_strong']:
            return ('BP4_strong', 'strong', 'benign')
        elif score <= t['BP4_moderate']:
            return ('BP4_moderate', 'moderate', 'benign')
        elif score <= t['BP4_supporting']:
            return ('BP4_supporting', 'supporting', 'benign')
        else:
            return ('neutral', None, 'neutral')
    
    def _estimate_grantham_impact(self, aa_change: str) -> Optional[float]:
        """
        Estimate structural impact from amino acid change using Grantham distance.
        
        Grantham distance measures physicochemical difference between amino acids.
        Higher distance = more likely to be damaging.
        
        Args:
            aa_change: Amino acid change string (e.g., 'p.Arg273His', 'R273H')
            
        Returns:
            float or None: Normalized Grantham-based impact score (0.0-1.0)
        """
        import re
        
        # Simplified Grantham distance categories
        # Based on Grantham R. (1974) Science 185:862-864
        grantham_categories = {
            # Conservative changes (low impact)
            ('D', 'E'): 45, ('N', 'D'): 23, ('Q', 'E'): 29, ('S', 'T'): 58,
            ('I', 'L'): 5, ('I', 'V'): 29, ('L', 'V'): 32, ('F', 'Y'): 22,
            ('K', 'R'): 26, ('A', 'S'): 99, ('A', 'G'): 60,
            # Radical changes (high impact)
            ('R', 'W'): 101, ('C', 'W'): 215, ('G', 'W'): 184, ('P', 'W'): 147,
            ('D', 'W'): 181, ('E', 'W'): 152, ('C', 'F'): 205, ('G', 'R'): 125,
        }
        
        # Parse amino acid change
        # Handle formats: p.Arg273His, R273H, Arg273His
        match = re.search(r'([A-Z])(\d+)([A-Z])', aa_change.upper())
        if not match:
            # Try three-letter code
            match = re.search(r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', aa_change)
            if match:
                aa_map = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                         'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                         'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                         'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}
                ref_aa = aa_map.get(match.group(1).capitalize())
                alt_aa = aa_map.get(match.group(3).capitalize())
                if ref_aa and alt_aa:
                    return self._grantham_to_score(ref_aa, alt_aa, grantham_categories)
            return None
        
        ref_aa = match.group(1)
        alt_aa = match.group(3)
        
        return self._grantham_to_score(ref_aa, alt_aa, grantham_categories)
    
    def _grantham_to_score(self, ref_aa: str, alt_aa: str, 
                           categories: Dict) -> float:
        """Convert Grantham distance to normalized score."""
        if ref_aa == alt_aa:
            return 0.0  # Synonymous
        
        # Look up distance (check both directions)
        distance = categories.get((ref_aa, alt_aa)) or categories.get((alt_aa, ref_aa))
        
        if distance is None:
            # Use average for unknown pairs
            return 0.5
        
        # Normalize: Grantham distances range roughly 5-215
        # Map to 0-1 where higher = more damaging
        normalized = min(distance / 200.0, 1.0)
        return normalized
    
    # Data loading methods (placeholders for external data sources)
    def _load_domain_regions(self) -> Dict:
        """Load functional domain region data."""
        return {}
    
    def _load_conservation_data(self) -> Dict:
        """Load conservation score data."""
        return {}
    
    def _load_structural_data(self) -> Dict:
        """Load protein structural data."""
        return {}
    
    # Normalization helper methods
    def _normalize_phylop(self, phylop: float) -> float:
        """
        Normalize phyloP score to 0-1 range.
        
        phyloP ranges roughly from -20 to +10, where:
        - Positive values indicate conservation (slower evolution)
        - Negative values indicate acceleration (faster evolution)
        
        Args:
            phylop: Raw phyloP score
            
        Returns:
            float: Normalized score (0.0-1.0), higher = more conserved
        """
        if phylop is None:
            return 0.5
        # Map [-5, +10] to [0, 1], clamping outliers
        return min(max((float(phylop) + 5) / 15.0, 0.0), 1.0)
    
    def _normalize_phastcons(self, phastcons: float) -> float:
        """
        Normalize phastCons score to 0-1 range.
        
        phastCons is already 0-1, representing probability of conservation.
        
        Args:
            phastcons: Raw phastCons score (0-1)
            
        Returns:
            float: Normalized score (0.0-1.0)
        """
        if phastcons is None:
            return 0.5
        return min(max(float(phastcons), 0.0), 1.0)
    
    def _normalize_gerp(self, gerp: float) -> float:
        """
        Normalize GERP++ score to 0-1 range.
        
        GERP++ ranges roughly from -12.3 to +6.17, where:
        - Positive values indicate conservation (rejected substitutions)
        - Negative values indicate acceleration
        
        Args:
            gerp: Raw GERP++ score
            
        Returns:
            float: Normalized score (0.0-1.0), higher = more conserved
        """
        if gerp is None:
            return 0.5
        # Map [-12.3, +6.17] to [0, 1]
        return min(max((float(gerp) + 12.3) / 18.47, 0.0), 1.0)
