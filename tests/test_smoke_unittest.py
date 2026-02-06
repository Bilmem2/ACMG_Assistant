import io
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stdout


ROOT_DIR = os.path.dirname(os.path.dirname(__file__))
SRC_DIR = os.path.join(ROOT_DIR, "src")
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from core.variant_data import VariantData
from core.evidence_evaluator import EvidenceEvaluator
from core.population_analyzer import GnomADClient
from core.acmg_classifier import ACMGClassifier
from core.phenotype_matcher import PhenotypeSimilarityCalculator
from utils.report_generator import ReportGenerator
from utils.predictor_api_client import PredictorAPIClient, PopulationAPIClient
from config.predictors import INSILICO_WEIGHTS


class TestSmokeUnittest(unittest.TestCase):
    def test_report_generation_contains_expected_sections(self):
        variant = VariantData(
            basic_info={
                "gene": "BRCA1",
                "chromosome": "17",
                "position": 43071077,
                "ref_allele": "C",
                "alt_allele": "T",
                "variant_type": "missense",
                "consequence": "missense_variant",
                "cdna_change": "c.68_69delAG",
                "protein_change": "p.Glu23Val",
            },
            population_data={"gnomad_af": 0.00001},
            insilico_data={"revel": 0.8, "cadd_phred": 25.0},
        )

        evidence_results = {
            "pathogenic_criteria": {
                "PM2": {"applies": True, "strength": "Moderate", "details": "Low frequency"}
            },
            "benign_criteria": {},
            "applied_criteria": {
                "PM2": {"applies": True, "strength": "Moderate", "details": "Low frequency"}
            },
            "vampp_score": None,
            "statistical_tests": {},
        }

        classifier = ACMGClassifier()
        classification_result = classifier.classify(evidence_results)

        with tempfile.TemporaryDirectory() as temp_dir:
            generator = ReportGenerator(output_dir=temp_dir)
            report_path = generator.generate_report(variant, evidence_results, classification_result)

            self.assertTrue(os.path.exists(report_path))
            with open(report_path, "r", encoding="utf-8") as report_file:
                contents = report_file.read()

            self.assertIn("FINAL CLASSIFICATION", contents)
            self.assertIn("PM2", contents)
            self.assertIn("gnomAD Overall AF", contents)
            self.assertIn("REVEL", contents)

    def test_hgvs_c_fallback_to_cdna_change(self):
        variant = VariantData(basic_info={"cdna_change": "c.123A>G"})
        self.assertEqual(variant.hgvs_c, "c.123A>G")

    def test_manual_insilico_overrides_create_predictor_scores(self):
        evaluator = EvidenceEvaluator(test_mode=True)
        variant = VariantData(insilico_data={"revel": 0.5, "cadd_phred": 20})
        evaluator._apply_manual_predictor_overrides(variant, variant.insilico_data)
        self.assertIsNotNone(variant.predictor_scores)
        self.assertIn("revel", variant.predictor_scores)
        self.assertIn("cadd_phred", variant.predictor_scores)
        self.assertEqual(variant.predictor_scores["revel"].value, 0.5)
        self.assertEqual(variant.predictor_scores["cadd_phred"].value, 20.0)
        self.assertEqual(variant.predictor_scores["revel"].source, "manual_input")

    def test_population_manual_overrides_api_stats(self):
        variant = VariantData(
            population_data={"gnomad_af": 0.02},
            population_stats={"gnomad_v4": {"af": 0.0}}
        )
        client = GnomADClient()
        frequencies = client.get_population_frequencies(variant)
        self.assertEqual(frequencies.get("ALL"), 0.02)

    def test_classifier_applied_criteria_normalization(self):
        classifier = ACMGClassifier()
        evidence_results = {
            "pathogenic_criteria": {"PM2": {}},
            "benign_criteria": {},
            "applied_criteria": {
                "PM2": {"strength": "Moderate", "applies": True}
            },
            "vampp_score": None,
            "statistical_tests": {}
        }
        classification = classifier.classify(evidence_results)
        summary = classifier.get_acmg_summary(classification)
        self.assertIn("PM2", summary["pathogenic_criteria"])
        self.assertEqual(summary["benign_criteria"], [])

    def test_colored_output_contains_classification(self):
        generator = ReportGenerator()
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            generator.print_colored_classification("VUS")
        output = buffer.getvalue()
        self.assertIn("CLASSIFICATION RESULT: VUS", output)

    def test_phenotype_similarity_downweights_generic_terms(self):
        calculator = PhenotypeSimilarityCalculator()
        similarity = calculator.calculate_similarity(
            {"HP:0002664"},
            {"HP:0002664", "HP:0003002"}
        )
        self.assertLess(similarity, 0.4)

    def test_phenotype_similarity_specific_terms_higher(self):
        calculator = PhenotypeSimilarityCalculator()
        generic = calculator.calculate_similarity(
            {"HP:0002664"},
            {"HP:0002664", "HP:0003002"}
        )
        specific = calculator.calculate_similarity(
            {"HP:0003002", "HP:0000137"},
            {"HP:0003002", "HP:0100013", "HP:0000137"}
        )
        self.assertGreater(specific, generic)

    @unittest.skipUnless(
        os.getenv("ACMG_INTEGRATION") == "1",
        "Integration tests disabled. Set ACMG_INTEGRATION=1 to enable."
    )
    def test_predictor_api_client_integration(self):
        client = PredictorAPIClient(api_enabled=True, timeout=20, test_mode=False)
        scores = client.get_predictor_scores(chrom="17", pos=43071077, ref="C", alt="T")
        self.assertIsInstance(scores, dict)
        for predictor in INSILICO_WEIGHTS.keys():
            self.assertIn(predictor, scores)

    @unittest.skipUnless(
        os.getenv("ACMG_INTEGRATION") == "1",
        "Integration tests disabled. Set ACMG_INTEGRATION=1 to enable."
    )
    def test_population_api_client_integration(self):
        client = PopulationAPIClient(api_enabled=True, timeout=20, test_mode=False)
        stats = client.get_population_stats(chrom="17", pos=43071077, ref="C", alt="T")
        self.assertIsInstance(stats, dict)


if __name__ == "__main__":
    unittest.main()
