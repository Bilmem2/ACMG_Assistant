"""
Statistical Utilities Module
============================

Provides statistical analysis functions for evidence evaluation.
"""

import math
from typing import Dict, Optional, Tuple
from scipy import stats


class StatisticalAnalyzer:
    """Performs statistical analyses for ACMG criteria evaluation."""
    
    def __init__(self):
        self.min_lod_supporting = 1.5
        self.min_lod_strong = 3.0
        self.min_families = 3
    
    def calculate_fishers_exact(self, cases_with: int, cases_total: int,
                                controls_with: int, controls_total: int) -> Dict:
        """
        Calculate Fisher's exact test for case-control data (PS4).
        
        Args:
            cases_with: Number of cases with variant
            cases_total: Total number of cases
            controls_with: Number of controls with variant
            controls_total: Total number of controls
            
        Returns:
            Dict with p_value, odds_ratio, significant, and interpretation
        """
        if cases_total <= 0 or controls_total <= 0:
            return {
                'valid': False,
                'error': 'Invalid sample sizes'
            }
        
        # Create 2x2 contingency table
        cases_without = cases_total - cases_with
        controls_without = controls_total - controls_with
        
        # Ensure non-negative values
        if cases_without < 0 or controls_without < 0:
            return {
                'valid': False,
                'error': 'Invalid variant counts'
            }
        
        table = [[cases_with, cases_without],
                 [controls_with, controls_without]]
        
        try:
            # Perform Fisher's exact test
            odds_ratio, p_value = stats.fisher_exact(table, alternative='greater')
            
            # Interpret results
            significant = p_value < 0.05 and odds_ratio >= 2.0
            
            return {
                'valid': True,
                'p_value': p_value,
                'odds_ratio': odds_ratio,
                'significant': significant,
                'confidence': 'high' if significant else 'medium',
                'interpretation': self._interpret_fishers(p_value, odds_ratio),
                'sample_sizes': {
                    'cases_total': cases_total,
                    'controls_total': controls_total,
                    'cases_with_variant': cases_with,
                    'controls_with_variant': controls_with
                }
            }
        except Exception as e:
            return {
                'valid': False,
                'error': f'Statistical calculation failed: {str(e)}'
            }
    
    def _interpret_fishers(self, p_value: float, odds_ratio: float) -> str:
        """Interpret Fisher's exact test results."""
        if p_value < 0.001:
            strength = "very strong"
        elif p_value < 0.01:
            strength = "strong"
        elif p_value < 0.05:
            strength = "moderate"
        else:
            strength = "weak"
        
        if odds_ratio >= 5.0:
            or_desc = "substantially"
        elif odds_ratio >= 2.0:
            or_desc = "significantly"
        else:
            or_desc = "marginally"
        
        if p_value < 0.05 and odds_ratio >= 2.0:
            return f"{strength.capitalize()} evidence: variant is {or_desc} enriched in cases (OR={odds_ratio:.2f}, p={p_value:.2e})"
        else:
            return f"No significant enrichment in cases (OR={odds_ratio:.2f}, p={p_value:.2e})"
    
    def calculate_lod_score(self, families_data: list) -> Dict:
        """
        Calculate LOD score for segregation analysis (PP1/BS4).
        
        Args:
            families_data: List of dicts with 'affected_with', 'affected_total',
                          'unaffected_with', 'unaffected_total'
        
        Returns:
            Dict with lod_score, strength, confidence, and interpretation
        """
        if not families_data or len(families_data) < self.min_families:
            return {
                'valid': False,
                'error': f'Insufficient families (minimum {self.min_families} required)',
                'families_provided': len(families_data) if families_data else 0
            }
        
        total_lod = 0.0
        family_lods = []
        
        for family in families_data:
            aff_with = family.get('affected_with', 0)
            aff_total = family.get('affected_total', 0)
            unaff_with = family.get('unaffected_with', 0)
            unaff_total = family.get('unaffected_total', 0)
            
            # Skip families with insufficient data
            if aff_total == 0:
                continue
            
            # Calculate likelihood ratio for this family
            # Assuming autosomal dominant inheritance (0.5 segregation probability)
            # P(data|variant causal) vs P(data|variant not causal)
            
            # Simplified LOD calculation
            # For proper LOD: log10(likelihood variant segregates / likelihood random)
            if aff_with > 0 and unaff_with == 0:
                # Perfect segregation with affected
                family_lod = aff_with * 0.3  # Simplified score
            elif aff_with == 0 and unaff_with > 0:
                # Perfect non-segregation (benign evidence)
                family_lod = -unaff_with * 0.3
            else:
                # Mixed segregation
                family_lod = (aff_with - unaff_with) * 0.15
            
            family_lods.append(family_lod)
            total_lod += family_lod
        
        # Determine strength based on LOD score
        if total_lod >= self.min_lod_strong:
            strength = 'strong'
            applies = True
            interpretation = f"Strong segregation evidence (LOD={total_lod:.2f} from {len(family_lods)} families)"
        elif total_lod >= self.min_lod_supporting:
            strength = 'supporting'
            applies = True
            interpretation = f"Supporting segregation evidence (LOD={total_lod:.2f} from {len(family_lods)} families)"
        elif total_lod <= -self.min_lod_supporting:
            strength = 'benign_strong'
            applies = True
            interpretation = f"Strong non-segregation evidence (LOD={total_lod:.2f} from {len(family_lods)} families)"
        else:
            strength = 'inconclusive'
            applies = False
            interpretation = f"Inconclusive segregation (LOD={total_lod:.2f} from {len(family_lods)} families)"
        
        return {
            'valid': True,
            'lod_score': total_lod,
            'family_count': len(family_lods),
            'individual_lods': family_lods,
            'strength': strength,
            'applies': applies,
            'confidence': 'high' if len(family_lods) >= 5 else 'medium',
            'interpretation': interpretation
        }
    
    def assess_sample_size_adequacy(self, cases: int, controls: int) -> Dict:
        """Assess if sample sizes are adequate for case-control analysis."""
        min_cases = 20
        min_controls = 50
        min_ratio = 0.2
        
        adequate = (cases >= min_cases and 
                   controls >= min_controls and
                   cases / controls >= min_ratio)
        
        if adequate:
            power = "adequate"
        elif cases >= min_cases // 2:
            power = "moderate"
        else:
            power = "insufficient"
        
        return {
            'adequate': adequate,
            'power': power,
            'cases': cases,
            'controls': controls,
            'recommendation': self._sample_size_recommendation(cases, controls)
        }
    
    def _sample_size_recommendation(self, cases: int, controls: int) -> str:
        """Provide recommendation for sample size improvement."""
        if cases < 20:
            return f"Increase cases to at least 20 (current: {cases})"
        if controls < 50:
            return f"Increase controls to at least 50 (current: {controls})"
        if cases / controls < 0.2:
            return f"Improve case-to-control ratio (current: {cases/controls:.2f})"
        return "Sample sizes adequate for analysis"
