import json

class GnomADClient:
    def get_population_frequencies(self, variant_data):
        key = f"{variant_data.gene}:{variant_data.hgvs_c}"
        try:
            with open('data/population_data/population_frequencies.json', 'r') as f:
                freq_db = json.load(f)
            return freq_db.get(key, {})
        except Exception:
            return {}

class EthnicityAwarePopulationAnalyzer:
    def __init__(self):
        self.gnomad_client = GnomADClient()
        self.ethnicity_thresholds = self._load_ethnicity_thresholds()
    
    def _load_ethnicity_thresholds(self):
        try:
            with open('data/population_data/ethnicity_thresholds.json', 'r') as f:
                return json.load(f)
        except Exception:
            # Default thresholds if file not found
            return {'ALL': 0.01, 'EUR': 0.01, 'AFR': 0.01}
    
    def analyze_population_frequency(self, variant_data, patient_ethnicity=None):
        """Ethnicity-aware population frequency analysis"""
        frequencies = self.gnomad_client.get_population_frequencies(variant_data)
        
        if patient_ethnicity:
            relevant_freq = frequencies.get(patient_ethnicity, frequencies.get('ALL'))
            threshold = self.ethnicity_thresholds.get(patient_ethnicity, 0.01)
        else:
            relevant_freq = frequencies.get('ALL')
            threshold = 0.005
        if relevant_freq is None:
            return None
        if relevant_freq >= threshold:
            return 'BA1'
        if relevant_freq >= threshold * 0.1:
            return 'BS1'
        if relevant_freq <= 0.0001:
            return 'PM2'
        return None
