class GeneSpecificRules:
    def __init__(self):
        self.hotspot_regions = self._load_hotspot_regions()
        self.gene_specific_thresholds = self._load_gene_thresholds()
    def evaluate_hotspot_region(self, variant_data):
        gene = getattr(variant_data, 'gene', None)
        position = self._extract_position(getattr(variant_data, 'hgvs_c', None))
        if position is None:
            return None
        if gene in self.hotspot_regions:
            for hotspot in self.hotspot_regions[gene]:
                if hotspot['start'] <= position <= hotspot['end']:
                    if hotspot.get('evidence_level', 0) >= 3:
                        return 'PM1'
        return None
    def _extract_position(self, hgvs_c):
        pass
    def _load_hotspot_regions(self):
        # Örnek veri, ileride dosyadan yüklenecek
        return {
            'KRAS': [
                {'start': 34, 'end': 38, 'evidence_level': 5, 'description': 'Codon 12-13'},
                {'start': 175, 'end': 183, 'evidence_level': 5, 'description': 'Codon 59-61'},
                {'start': 436, 'end': 444, 'evidence_level': 4, 'description': 'Codon 146'}
            ],
            'TP53': [
                {'start': 730, 'end': 750, 'evidence_level': 5, 'description': 'DNA binding domain'},
                {'start': 818, 'end': 825, 'evidence_level': 4, 'description': 'Codon 273'}
            ],
            'BRAF': [
                {'start': 1798, 'end': 1800, 'evidence_level': 5, 'description': 'V600E hotspot'}
            ]
        }
    def _load_gene_thresholds(self):
        pass
