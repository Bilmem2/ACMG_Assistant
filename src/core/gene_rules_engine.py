class GeneRulesEngine:
    def __init__(self):
        self.gene_rules = self._load_gene_specific_rules()
        self.inheritance_patterns = self._load_inheritance_patterns()
    
    def _load_gene_specific_rules(self):
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
    
    def _load_inheritance_patterns(self):
        # Placeholder for inheritance patterns
        return {}
    
    def apply_gene_specific_rules(self, variant_data, preliminary_evidence):
        gene = getattr(variant_data, 'gene', None)
        if gene not in self.gene_rules:
            return preliminary_evidence
        rules = self.gene_rules[gene]
        modified_evidence = preliminary_evidence.copy()
        for rule in rules:
            if self._rule_conditions_met(rule, variant_data):
                modified_evidence = self._apply_rule_modifications(rule, modified_evidence)
        return modified_evidence
    
    def _rule_conditions_met(self, rule, variant_data):
        # Placeholder: always True for demonstration
        return True
    
    def _apply_rule_modifications(self, rule, evidence):
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
