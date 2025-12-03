# =============================================================================
# Gene Rules Engine - Core Implementation
# =============================================================================
# Changes in this documentation pass (Dec 2024):
# - Added module and class docstrings
# - Added type hints to public methods
# - Clarified relationship with utils/gene_rules_engine.py
#
# NOTE: This is the CLASS-BASED implementation for gene-specific rule application.
# A FUNCTION-BASED implementation also exists in utils/gene_rules_engine.py which
# provides ClinGen Dosage Sensitivity integration for threshold calculations.
# TODO: Consider merging these two modules in a future refactor to consolidate
#       gene-specific logic in one place.
# =============================================================================
"""
Gene Rules Engine Module
========================

Applies gene-specific rules to modify ACMG evidence criteria based on
known gene characteristics (e.g., BRCA1 RING domain, TP53 DNA binding domain).

See also: utils/gene_rules_engine.py for ClinGen-based threshold calculations.
"""

from typing import Dict, List, Any, Optional


class GeneRulesEngine:
    """
    Engine for applying gene-specific ACMG classification rules.
    
    This class modifies preliminary evidence based on gene-specific knowledge,
    such as critical functional domains, known hotspots, and inheritance patterns.
    
    Attributes:
        gene_rules: Dict mapping gene symbols to lists of modification rules
        inheritance_patterns: Dict of gene inheritance pattern information
    """
    
    def __init__(self):
        """Initialize the gene rules engine with default rules."""
        self.gene_rules = self._load_gene_specific_rules()
        self.inheritance_patterns = self._load_inheritance_patterns()
    
    def _load_gene_specific_rules(self) -> Dict[str, List[Dict[str, Any]]]:
        # Example rules, can be loaded from JSON files
        return {
            'BRCA1': [
                {
                    'condition': 'missense_in_ring_domain',
                    'action': 'upgrade_PM1_to_PM2',
                    'evidence_level': 4
                },
                {
                    'condition': 'nonsense_before_exon_20',
                    'action': 'apply_PVS1',
                    'evidence_level': 5
                }
            ],
            'BRCA2': [
                {
                    'condition': 'missense_in_dna_binding_domain',
                    'action': 'apply_PM1',
                    'evidence_level': 4
                }
            ],
            'TP53': [
                {
                    'condition': 'missense_in_dna_binding_domain',
                    'action': 'upgrade_PM1_to_PS1',
                    'evidence_level': 5
                }
            ]
        }
    
    def _load_inheritance_patterns(self) -> Dict[str, Any]:
        """
        Load gene inheritance pattern data.
        
        Returns:
            Dict mapping gene symbols to inheritance information
        """
        # Placeholder for inheritance patterns
        return {}
    
    def apply_gene_specific_rules(self, variant_data, preliminary_evidence: List[str]) -> List[str]:
        gene = getattr(variant_data, 'gene', None)
        if gene not in self.gene_rules:
            return preliminary_evidence
        rules = self.gene_rules[gene]
        modified_evidence = preliminary_evidence.copy()
        for rule in rules:
            if self._rule_conditions_met(rule, variant_data):
                modified_evidence = self._apply_rule_modifications(rule, modified_evidence)
        return modified_evidence
    
    def _rule_conditions_met(self, rule: Dict[str, Any], variant_data) -> bool:
        """
        Check if rule conditions are met for the given variant.
        
        Args:
            rule: Rule definition with condition and action
            variant_data: VariantData object to evaluate
            
        Returns:
            bool: True if conditions are met
        """
        # Placeholder: always True for demonstration
        return True
    
    def _apply_rule_modifications(self, rule: Dict[str, Any], evidence: List[str]) -> List[str]:
        # Example: upgrade PM1 to PM2 or PS1
        if rule['action'] == 'upgrade_PM1_to_PM2' and 'PM1' in evidence:
            evidence.remove('PM1')
            evidence.append('PM2')
        elif rule['action'] == 'upgrade_PM1_to_PS1' and 'PM1' in evidence:
            evidence.remove('PM1')
            evidence.append('PS1')
        elif rule['action'] == 'apply_PVS1' and 'PVS1' not in evidence:
            evidence.append('PVS1')
        elif rule['action'] == 'apply_PM1' and 'PM1' not in evidence:
            evidence.append('PM1')
        return evidence
