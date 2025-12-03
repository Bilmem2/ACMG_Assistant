# =============================================================================
# Changes in this documentation pass (Dec 2024):
# - Added docstrings to ConflictResolver and WeightedACMGClassifier classes
# - Added type hints to public methods
# - Clarified relationship between weighted and standard classifiers
# =============================================================================
"""
ACMG Classifier Module
====================

This module implements the ACMG classification logic based on
evidence criteria evaluation results.
"""

from typing import Dict, List, Optional, Any, Tuple
from config.constants import CLASSIFICATION_RULES, EVIDENCE_WEIGHTS

# =============================================================================
# Alternative/Experimental Classifiers
# =============================================================================
# NOTE: ConflictResolver and WeightedACMGClassifier below are ALTERNATIVE
# implementations that use weighted scoring. The CANONICAL implementation
# is the ACMGClassifier class further below, which follows standard ACMG rules.
# TODO: Consider consolidating these into a single configurable classifier.
# =============================================================================


class ConflictResolver:
    """
    Resolves conflicts between pathogenic and benign evidence.
    
    Used by WeightedACMGClassifier when both pathogenic and benign
    evidence are present for a variant.
    """
    
    def resolve_conflict(self, pathogenic_evidence: list, benign_evidence: list,
                        pathogenic_score: float, benign_score: float) -> str:
        """
        Resolve classification conflict based on weighted scores.
        
        Args:
            pathogenic_evidence: List of pathogenic criteria codes
            benign_evidence: List of benign criteria codes
            pathogenic_score: Weighted pathogenic score
            benign_score: Weighted benign score
            
        Returns:
            str: Resolved classification with conflict notation
        """
        if pathogenic_score > benign_score:
            return "Likely Pathogenic (conflict resolved)"
        elif benign_score > pathogenic_score:
            return "Likely Benign (conflict resolved)"
        else:
            return "Uncertain Significance (conflict unresolved)"


class WeightedACMGClassifier:
    """
    Alternative ACMG classifier using weighted evidence scoring.
    
    This classifier assigns numerical weights to evidence criteria
    and calculates aggregate scores for classification. This is an
    EXPERIMENTAL approach - use ACMGClassifier for standard ACMG rules.
    
    Attributes:
        evidence_weights: Dict mapping evidence prefixes to weights
        conflict_resolver: ConflictResolver instance for handling conflicts
    """
    
    def __init__(self):
        """Initialize the weighted classifier with default weights."""
        self.evidence_weights = self._load_evidence_weights()
        self.conflict_resolver = ConflictResolver()
    
    def _load_evidence_weights(self) -> dict:
        """
        Load evidence weights for scoring.
        
        Returns:
            Dict mapping evidence type prefixes to numerical weights
        """
        # Example weights, can be loaded from config file
        return {
            'PVS': 8.0,
            'PS': 4.0,
            'PM': 2.0,
            'PP': 1.0,
            'BA': 8.0,
            'BS': 4.0,
            'BP': 1.0
        }
    
    def _get_evidence_confidence(self, evidence: str) -> float:
        """
        Get confidence multiplier for evidence.
        
        Args:
            evidence: Evidence criterion code (e.g., 'PVS1', 'PM2')
            
        Returns:
            float: Confidence multiplier (currently always 1.0)
        """
        # Placeholder: always 1.0, can be extended for real confidence
        return 1.0
    
    def _calculate_weighted_score(self, evidence_list: list) -> float:
        """
        Calculate weighted score for a list of evidence criteria.
        
        Args:
            evidence_list: List of evidence criterion codes
            
        Returns:
            float: Total weighted score
        """
        score = 0
        for evidence in evidence_list:
            base_weight = self.evidence_weights.get(evidence[:3], 1.0)
            confidence = self._get_evidence_confidence(evidence)
            score += base_weight * confidence
        return score
    
    def classify_variant_weighted(self, evidence_list: list) -> str:
        """
        Classify variant using weighted evidence scoring.
        
        Args:
            evidence_list: List of all applied evidence criteria
            
        Returns:
            str: Classification result
        """
        pathogenic_evidence = [e for e in evidence_list if e.startswith(('PVS', 'PS', 'PM', 'PP'))]
        benign_evidence = [e for e in evidence_list if e.startswith(('BA', 'BS', 'BP'))]
        pathogenic_score = self._calculate_weighted_score(pathogenic_evidence)
        benign_score = self._calculate_weighted_score(benign_evidence)
        if pathogenic_score > 0 and benign_score > 0:
            return self.conflict_resolver.resolve_conflict(
                pathogenic_evidence, benign_evidence, pathogenic_score, benign_score
            )
        return self._determine_classification(pathogenic_score, benign_score)
    
    def _determine_classification(self, pathogenic_score: float, benign_score: float) -> str:
        """
        Determine classification based on weighted scores.
        
        Args:
            pathogenic_score: Total pathogenic evidence score
            benign_score: Total benign evidence score
            
        Returns:
            str: Classification result
        """
        if pathogenic_score > benign_score and pathogenic_score >= 4:
            return "Pathogenic"
        elif benign_score > pathogenic_score and benign_score >= 4:
            return "Benign"
        else:
            return "Uncertain Significance"


class ACMGClassifier:
    """
    ACMG variant classifier implementing the official classification rules.
    
    This class takes evidence evaluation results and applies ACMG rules
    to determine the final variant classification.
    """
    
    def __init__(self, use_2023_guidelines: bool = False):
        """
        Initialize the ACMG classifier.
        
        Args:
            use_2023_guidelines (bool): Whether to use ACMG 2023 guidelines
        """
        self.use_2023_guidelines = use_2023_guidelines
        self.classification_rules = CLASSIFICATION_RULES
    
    def classify(self, evidence_results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Classify variant based on evidence results.
        
        Args:
            evidence_results (Dict[str, Any]): Evidence evaluation results
            
        Returns:
            Dict[str, Any]: Classification results
        """
        # Count evidence by strength
        pathogenic_counts = self._count_pathogenic_evidence(evidence_results)
        benign_counts = self._count_benign_evidence(evidence_results)
        
        # Apply classification rules
        classification = self._apply_classification_rules(pathogenic_counts, benign_counts)
        
        # Calculate confidence
        confidence = self._calculate_confidence(evidence_results, classification)
        
        # Generate suggestions
        suggestions = self._generate_suggestions(evidence_results, classification)
        
        # Check for conflicts
        conflicts = self._check_conflicts(pathogenic_counts, benign_counts)
        
        return {
            'classification': classification,
            'confidence': confidence,
            'pathogenic_counts': pathogenic_counts,
            'benign_counts': benign_counts,
            'applied_criteria': evidence_results.get('applied_criteria', {}),
            'suggestions': suggestions,
            'conflicts': conflicts,
            'vampp_score': evidence_results.get('vampp_score'),
            'statistical_tests': evidence_results.get('statistical_tests', {}),
            'guidelines_version': '2023' if self.use_2023_guidelines else '2015'
        }
    
    def _count_pathogenic_evidence(self, evidence_results: Dict[str, Any]) -> Dict[str, int]:
        """Count pathogenic evidence by strength."""
        counts = {
            'very_strong': 0,
            'strong': 0,
            'moderate': 0,
            'supporting': 0
        }
        
        pathogenic_criteria = evidence_results.get('pathogenic_criteria', {})
        applied_criteria = evidence_results.get('applied_criteria', {})
        
        for criterion, result in applied_criteria.items():
            if criterion in pathogenic_criteria:
                strength = result.get('strength', '').lower().replace(' ', '_')
                
                if strength == 'very_strong':
                    counts['very_strong'] += 1
                elif strength == 'strong':
                    counts['strong'] += 1
                elif strength == 'moderate':
                    counts['moderate'] += 1
                elif strength == 'supporting':
                    counts['supporting'] += 1
        
        return counts
    
    def _count_benign_evidence(self, evidence_results: Dict[str, Any]) -> Dict[str, int]:
        """Count benign evidence by strength."""
        counts = {
            'stand_alone': 0,
            'strong': 0,
            'supporting': 0
        }
        
        benign_criteria = evidence_results.get('benign_criteria', {})
        applied_criteria = evidence_results.get('applied_criteria', {})
        
        for criterion, result in applied_criteria.items():
            if criterion in benign_criteria:
                strength = result.get('strength', '').lower().replace(' ', '_').replace('-', '_')
                
                if strength == 'stand_alone':
                    counts['stand_alone'] += 1
                elif strength == 'strong':
                    counts['strong'] += 1
                elif strength == 'supporting':
                    counts['supporting'] += 1
        
        return counts
    
    def _apply_classification_rules(self, pathogenic_counts: Dict[str, int], 
                                  benign_counts: Dict[str, int]) -> str:
        """Apply ACMG classification rules."""
        
        # Check for benign classification first (takes precedence)
        if self._meets_benign_criteria(benign_counts):
            return 'Benign'
        
        if self._meets_likely_benign_criteria(benign_counts):
            return 'Likely Benign'
        
        # Check for pathogenic classification
        if self._meets_pathogenic_criteria(pathogenic_counts):
            return 'Pathogenic'
        
        if self._meets_likely_pathogenic_criteria(pathogenic_counts):
            return 'Likely Pathogenic'
        
        # Default to VUS if no criteria met
        return 'VUS'
    
    def _meets_pathogenic_criteria(self, counts: Dict[str, int]) -> bool:
        """Check if variant meets pathogenic criteria."""
        rules = self.classification_rules['Pathogenic']['rules']
        
        for rule in rules:
            if self._matches_rule(counts, rule):
                return True
        
        return False
    
    def _meets_likely_pathogenic_criteria(self, counts: Dict[str, int]) -> bool:
        """Check if variant meets likely pathogenic criteria."""
        rules = self.classification_rules['Likely Pathogenic']['rules']
        
        for rule in rules:
            if self._matches_rule(counts, rule):
                return True
        
        return False
    
    def _meets_benign_criteria(self, counts: Dict[str, int]) -> bool:
        """Check if variant meets benign criteria."""
        rules = self.classification_rules['Benign']['rules']
        
        for rule in rules:
            if self._matches_rule(counts, rule):
                return True
        
        return False
    
    def _meets_likely_benign_criteria(self, counts: Dict[str, int]) -> bool:
        """Check if variant meets likely benign criteria."""
        rules = self.classification_rules['Likely Benign']['rules']
        
        for rule in rules:
            if self._matches_rule(counts, rule):
                return True
        
        return False
    
    def _matches_rule(self, counts: Dict[str, int], rule: Dict[str, int]) -> bool:
        """Check if evidence counts match a specific rule."""
        for strength, required_count in rule.items():
            if counts.get(strength, 0) < required_count:
                return False
        
        return True
    
    def _calculate_confidence(self, evidence_results: Dict[str, Any], 
                            classification: str) -> str:
        """Calculate confidence level for the classification."""
        
        # Get total evidence counts
        pathogenic_counts = self._count_pathogenic_evidence(evidence_results)
        benign_counts = self._count_benign_evidence(evidence_results)
        
        total_pathogenic = sum(pathogenic_counts.values())
        total_benign = sum(benign_counts.values())
        
        # Calculate confidence based on evidence strength and quantity
        if classification in ['Pathogenic', 'Benign']:
            # Strong classifications
            if pathogenic_counts.get('very_strong', 0) > 0 or benign_counts.get('stand_alone', 0) > 0:
                return 'High'
            elif pathogenic_counts.get('strong', 0) >= 2 or benign_counts.get('strong', 0) >= 2:
                return 'High'
            elif total_pathogenic >= 3 or total_benign >= 2:
                return 'Medium'
            else:
                return 'Low'
        
        elif classification in ['Likely Pathogenic', 'Likely Benign']:
            # Moderate classifications
            if total_pathogenic >= 3 or total_benign >= 2:
                return 'Medium'
            elif total_pathogenic >= 2 or total_benign >= 1:
                return 'Low'
            else:
                return 'Very Low'
        
        else:  # VUS
            return 'N/A'
    
    def _generate_suggestions(self, evidence_results: Dict[str, Any], 
                           classification: str) -> List[str]:
        """Generate suggestions for improving classification."""
        suggestions = []
        
        # Get applied criteria
        applied_pathogenic = evidence_results.get('applied_criteria', {}).get('pathogenic', [])
        applied_benign = evidence_results.get('applied_criteria', {}).get('benign', [])
        
        # Suggest missing key evidence
        if classification == 'VUS':
            if 'PS2' not in applied_pathogenic and 'PS2_Very_Strong' not in applied_pathogenic:
                suggestions.append("Consider de novo status evaluation (PS2)")
            
            if 'PS3' not in applied_pathogenic and 'BS3' not in applied_benign:
                suggestions.append("Consider functional studies (PS3/BS3)")
            
            if 'PP1' not in applied_pathogenic and 'BS4' not in applied_benign:
                suggestions.append("Consider segregation analysis (PP1/BS4)")
            
            if 'PS4' not in applied_pathogenic:
                suggestions.append("Consider case-control studies (PS4)")
        
        # Suggest database checks
        if not evidence_results.get('applied_criteria', {}).get('pathogenic'):
            suggestions.append("Check ClinVar for existing classifications")
            suggestions.append("Verify variant in gnomAD and other population databases")
        
        # Suggest additional in silico tools
        vampp_score = evidence_results.get('vampp_score')
        if vampp_score is None:
            suggestions.append("Obtain additional in silico prediction scores")
        
        # Suggest literature review
        suggestions.append("Review recent literature for this variant or gene")
        
        return suggestions
    
    def _check_conflicts(self, pathogenic_counts: Dict[str, int], 
                        benign_counts: Dict[str, int]) -> List[str]:
        """Check for conflicts between pathogenic and benign evidence."""
        conflicts = []
        
        total_pathogenic = sum(pathogenic_counts.values())
        total_benign = sum(benign_counts.values())
        
        # Check for conflicting evidence
        if total_pathogenic > 0 and total_benign > 0:
            conflicts.append("Conflicting pathogenic and benign evidence detected")
        
        # Check for stand-alone benign with pathogenic evidence
        if benign_counts.get('stand_alone', 0) > 0 and total_pathogenic > 0:
            conflicts.append("Stand-alone benign evidence conflicts with pathogenic evidence")
        
        # Check for very strong pathogenic with benign evidence
        if pathogenic_counts.get('very_strong', 0) > 0 and total_benign > 0:
            conflicts.append("Very strong pathogenic evidence conflicts with benign evidence")
        
        return conflicts
    
    def get_classification_explanation(self, classification_result: Dict[str, Any]) -> str:
        """Get detailed explanation of the classification."""
        
        classification = classification_result['classification']
        pathogenic_counts = classification_result['pathogenic_counts']
        benign_counts = classification_result['benign_counts']
        applied_criteria = classification_result['applied_criteria']
        
        explanation = f"Classification: {classification}\n\n"
        
        # Evidence summary
        explanation += "Evidence Summary:\n"
        explanation += "-" * 20 + "\n"
        
        if applied_criteria.get('pathogenic'):
            explanation += f"Pathogenic criteria: {', '.join(applied_criteria['pathogenic'])}\n"
            explanation += f"  Very Strong: {pathogenic_counts['very_strong']}\n"
            explanation += f"  Strong: {pathogenic_counts['strong']}\n"
            explanation += f"  Moderate: {pathogenic_counts['moderate']}\n"
            explanation += f"  Supporting: {pathogenic_counts['supporting']}\n"
        
        if applied_criteria.get('benign'):
            explanation += f"Benign criteria: {', '.join(applied_criteria['benign'])}\n"
            explanation += f"  Stand-alone: {benign_counts['stand_alone']}\n"
            explanation += f"  Strong: {benign_counts['strong']}\n"
            explanation += f"  Supporting: {benign_counts['supporting']}\n"
        
        # Classification rationale
        explanation += "\nClassification Rationale:\n"
        explanation += "-" * 25 + "\n"
        
        if classification == 'Pathogenic':
            explanation += "Meets criteria for Pathogenic classification based on:\n"
            if pathogenic_counts['very_strong'] >= 1:
                explanation += "- One Very Strong pathogenic criterion\n"
            elif pathogenic_counts['strong'] >= 2:
                explanation += "- Two Strong pathogenic criteria\n"
            else:
                explanation += "- Combination of Strong, Moderate, and Supporting criteria\n"
        
        elif classification == 'Likely Pathogenic':
            explanation += "Meets criteria for Likely Pathogenic classification based on:\n"
            explanation += "- Combination of pathogenic criteria not reaching Pathogenic threshold\n"
        
        elif classification == 'Benign':
            explanation += "Meets criteria for Benign classification based on:\n"
            if benign_counts['stand_alone'] >= 1:
                explanation += "- One Stand-alone benign criterion\n"
            else:
                explanation += "- Combination of Strong and Supporting benign criteria\n"
        
        elif classification == 'Likely Benign':
            explanation += "Meets criteria for Likely Benign classification based on:\n"
            explanation += "- Combination of benign criteria not reaching Benign threshold\n"
        
        else:  # VUS
            explanation += "Variant of Uncertain Significance (VUS):\n"
            explanation += "- Insufficient evidence for pathogenic or benign classification\n"
            explanation += "- Additional data needed for definitive classification\n"
        
        # Additional information
        vampp_score = classification_result.get('vampp_score')
        if vampp_score is not None:
            explanation += f"\nVAMPP-score-like metascore: {vampp_score:.3f}\n"
        
        # Conflicts
        conflicts = classification_result.get('conflicts', [])
        if conflicts:
            explanation += "\nConflicts:\n"
            explanation += "-" * 10 + "\n"
            for conflict in conflicts:
                explanation += f"- {conflict}\n"
        
        # Suggestions
        suggestions = classification_result.get('suggestions', [])
        if suggestions:
            explanation += "\nSuggestions:\n"
            explanation += "-" * 12 + "\n"
            for suggestion in suggestions:
                explanation += f"- {suggestion}\n"
        
        return explanation
    
    def get_acmg_summary(self, classification_result: Dict[str, Any]) -> Dict[str, Any]:
        """Get a concise summary of ACMG classification."""
        
        return {
            'classification': classification_result['classification'],
            'confidence': classification_result['confidence'],
            'pathogenic_criteria': classification_result['applied_criteria'].get('pathogenic', []),
            'benign_criteria': classification_result['applied_criteria'].get('benign', []),
            'total_pathogenic_evidence': sum(classification_result['pathogenic_counts'].values()),
            'total_benign_evidence': sum(classification_result['benign_counts'].values()),
            'vampp_score': classification_result.get('vampp_score'),
            'has_conflicts': bool(classification_result.get('conflicts', [])),
            'guidelines_version': classification_result.get('guidelines_version', '2015')
        }
