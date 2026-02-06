import os
import sys
import io
import tempfile
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch


ROOT_DIR = os.path.dirname(os.path.dirname(__file__))
SRC_DIR = os.path.join(ROOT_DIR, "src")
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from core.variant_data import VariantData
from core.acmg_classifier import ACMGClassifier
from core.evidence_evaluator import EvidenceEvaluator
from core.gene_specific_rules import GeneSpecificRules
from core.missense_evaluator import MissenseEvaluator
from core.phenotype_matcher import HPOClient, PhenotypeSimilarityCalculator
from core.interactive_evidence import ManualEvidence, map_functional_studies_to_evidence, MockInputProvider, InteractiveEvidenceCollector
from core.population_analyzer import EthnicityAwarePopulationAnalyzer
from utils.cache import ResultCache, CacheKey
from utils.report_generator import ReportGenerator
from utils.predictor_api_client import PredictorAPIClient
from utils.domain_api_client import HotspotAnnotation
from core.functional_studies_evaluator import FunctionalStudiesEvaluator
from config.constants import TEST_SCENARIOS
from utils.input_handler import InputHandler
from utils.api_client import APIClient


class FakeDomainClient:
    def __init__(self, annotation):
        self.annotation = annotation

    def get_hotspot_annotation(self, gene, position=None, hgvs_p=None):
        return self.annotation


class TestFullStackUnittest(unittest.TestCase):
    def test_cli_flow_test_mode(self):
        cwd = os.getcwd()
        scenario_key = next(iter(TEST_SCENARIOS.keys()))

        with patch("utils.input_handler.InputHandler._select_test_scenario", return_value=scenario_key):
            from acmg_assistant import ACMGAssistant

            assistant = ACMGAssistant(test_mode=True)
            try:
                with tempfile.TemporaryDirectory() as temp_dir:
                    assistant.report_generator = ReportGenerator(output_dir=temp_dir)
                    variant_data = assistant._collect_variant_data()
                    evidence = assistant._evaluate_evidence(variant_data)
                    classification = assistant._classify_variant(evidence)
                    report_path = assistant._generate_report(variant_data, evidence, classification)

                    self.assertTrue(os.path.exists(report_path))
                    with open(report_path, "r", encoding="utf-8") as report_file:
                        contents = report_file.read()
                    self.assertIn("CLASSIFICATION SUMMARY", contents)
            finally:
                os.chdir(cwd)

    def test_cli_interactive_prompt_skip(self):
        cwd = os.getcwd()
        scenario_key = next(iter(TEST_SCENARIOS.keys()))

        with patch("utils.input_handler.InputHandler._select_test_scenario", return_value=scenario_key):
            from acmg_assistant import ACMGAssistant

            assistant = ACMGAssistant(test_mode=True)
            try:
                variant_data = assistant._collect_variant_data()
                with patch("builtins.input", return_value="n"):
                    manual = assistant._offer_interactive_evidence_collection(variant_data)
                self.assertIsNone(manual)
            finally:
                os.chdir(cwd)

    def test_cli_stdout_snapshot(self):
        cwd = os.getcwd()
        scenario_key = next(iter(TEST_SCENARIOS.keys()))

        with patch("utils.input_handler.InputHandler._select_test_scenario", return_value=scenario_key):
            from acmg_assistant import ACMGAssistant

            buffer = io.StringIO()
            with redirect_stdout(buffer):
                assistant = ACMGAssistant(test_mode=True)
            output = buffer.getvalue()

            self.assertIn("ACMG Variant Classification Assistant", output)
            self.assertIn("DISCLAIMER", output)
            self.assertIn("Running in Test Mode", output)

            os.chdir(cwd)

    def test_cli_full_flow_with_mocked_input(self):
        cwd = os.getcwd()
        scenario_key = next(iter(TEST_SCENARIOS.keys()))

        basic_info = {
            "gene": "BRCA1",
            "chromosome": "17",
            "position": 43071077,
            "ref_allele": "C",
            "alt_allele": "T",
            "variant_type": "nonsense",
            "consequence": "stop_gained",
            "cdna_change": "c.68_69delAG",
            "protein_change": "p.Glu23Val",
        }
        population_data = {"gnomad_af": 0.0}
        insilico_data = {"revel": 0.9, "cadd_phred": 30.0}
        genetic_data = {"inheritance": "AD", "zygosity": "heterozygous"}
        functional_data = {}

        with patch("utils.input_handler.InputHandler._select_test_scenario", return_value=scenario_key), \
             patch.object(InputHandler, "collect_basic_info", return_value=basic_info), \
             patch.object(InputHandler, "collect_population_data", return_value=population_data), \
             patch.object(InputHandler, "collect_insilico_data", return_value=insilico_data), \
             patch.object(InputHandler, "collect_genetic_data", return_value=genetic_data), \
             patch.object(InputHandler, "collect_functional_data", return_value=functional_data), \
             patch.object(InputHandler, "collect_patient_phenotypes", return_value=[]), \
             patch.object(APIClient, "get_clinvar_status", return_value={"status": "not_found"}), \
               patch.object(APIClient, "auto_enrich_variant_data", return_value={}), \
               patch.object(EvidenceEvaluator, "_evaluate_pm1", return_value={"applies": False, "strength": "Moderate", "details": "skipped"}), \
               patch.object(EvidenceEvaluator, "_evaluate_pm3", return_value={"applies": False, "strength": "Moderate", "details": "skipped"}), \
                             patch.object(EvidenceEvaluator, "_evaluate_ps4", return_value={"applies": False, "strength": "Strong", "details": "skipped"}), \
                             patch.object(EvidenceEvaluator, "_evaluate_bp5", return_value={"applies": False, "strength": "Supporting", "details": "skipped"}):
            from acmg_assistant import ACMGAssistant

            with patch("builtins.input", side_effect=["", "n"]):
                assistant = ACMGAssistant(test_mode=False)

            try:
                with tempfile.TemporaryDirectory() as temp_dir:
                    assistant.report_generator = ReportGenerator(output_dir=temp_dir)

                    variant_data = assistant._collect_variant_data()
                    with patch("builtins.input", side_effect=["n", "n", "n", "n", "n"]):
                        evidence = assistant._evaluate_evidence(variant_data)
                    classification = assistant._classify_variant(evidence)
                    report_path = assistant._generate_report(variant_data, evidence, classification)

                    self.assertTrue(os.path.exists(report_path))
            finally:
                os.chdir(cwd)
    def test_end_to_end_evaluation_and_report(self):
        variant_data = VariantData(
            basic_info={
                "gene": "TEST",
                "chromosome": "1",
                "position": 1000,
                "ref_allele": "A",
                "alt_allele": "G",
                "variant_type": "missense",
                "consequence": "missense_variant",
                "cdna_change": "c.100A>G",
                "protein_change": "p.Lys34Arg",
            },
            population_data={"gnomad_af": 0.00001},
            insilico_data={"revel": 0.8, "cadd_phred": 25.0},
            genetic_data={"inheritance": "AD", "zygosity": "heterozygous"},
            functional_data={}
        )

        evaluator = EvidenceEvaluator(test_mode=True)
        evaluator._evaluate_pm1 = lambda _: {"applies": False, "strength": "Moderate", "details": "skipped"}
        evidence = evaluator.evaluate_all_criteria(variant_data)

        classifier = ACMGClassifier()
        classification = classifier.classify(evidence)

        with tempfile.TemporaryDirectory() as temp_dir:
            generator = ReportGenerator(output_dir=temp_dir)
            report_path = generator.generate_report(variant_data, evidence, classification)
            self.assertTrue(os.path.exists(report_path))
            with open(report_path, "r", encoding="utf-8") as report_file:
                contents = report_file.read()
            self.assertIn("VARIANT SUMMARY", contents)
            self.assertIn("CLASSIFICATION SUMMARY", contents)

    def test_manual_evidence_merge(self):
        evaluator = EvidenceEvaluator(test_mode=True)
        manual = ManualEvidence()
        manual.add_evidence("PS3_strong", "Functional evidence")
        evaluator.set_manual_evidence(manual)
        merged = evaluator.merge_manual_with_automated({
            "pathogenic_criteria": {},
            "benign_criteria": {},
            "applied_criteria": {},
            "evidence_details": {}
        })
        self.assertIn("PS3", merged["applied_criteria"])

    def test_gene_specific_rules_with_mock_domain(self):
        annotation = HotspotAnnotation(
            is_hotspot=True,
            in_critical_domain=False,
            confidence=0.95,
            source="mock",
            tumor_count=12
        )
        rules = GeneSpecificRules(domain_client=FakeDomainClient(annotation))
        result = rules.evaluate_pm1(gene="TP53", position=248)
        self.assertTrue(result.applies)
        self.assertEqual(result.evidence_code, "PM1")

    def test_missense_evaluator_outputs(self):
        variant = VariantData(
            basic_info={"variant_type": "missense", "amino_acid_change": "R273H"},
            insilico_data={"revel": 0.9, "cadd_phred": 30.0, "sift": 0.01},
            population_data={"gnomad_af": 0.00001}
        )
        evaluator = MissenseEvaluator()
        result = evaluator.evaluate_missense_variant(variant)
        self.assertIn("composite_score", result)
        self.assertGreaterEqual(result["composite_score"], 0.0)
        self.assertLessEqual(result["composite_score"], 1.0)

    def test_population_analyzer_ethnicity_threshold(self):
        variant = VariantData(population_data={"gnomad_af": 0.02})
        analyzer = EthnicityAwarePopulationAnalyzer()
        evidence = analyzer.analyze_population_frequency(variant, patient_ethnicity=None)
        self.assertIn(evidence, ["BA1", "BS1", "PM2", None])

    def test_hpo_client_unknown_term(self):
        client = HPOClient(synonyms_path=os.path.join(ROOT_DIR, "src", "data", "hpo_synonyms.json"))
        terms = client.get_phenotype_terms("unknownphenotype")
        self.assertTrue(any(t.startswith("TEXT:") for t in terms))

    def test_phenotype_similarity_weights(self):
        calc = PhenotypeSimilarityCalculator()
        score = calc.calculate_similarity({"HP:0002664"}, {"HP:0002664", "HP:0003002"})
        self.assertLess(score, 0.5)

    def test_functional_studies_mapping(self):
        code, _ = map_functional_studies_to_evidence(2, 0, "high")
        self.assertEqual(code, "PS3_strong")

    def test_interactive_evidence_collector_no_input(self):
        input_provider = MockInputProvider(["", "", "", "", "", "", "", "", "", "", "", "", ""])
        collector = InteractiveEvidenceCollector(input_provider=input_provider, show_prompts=False)
        evidence = collector.collect_selective(["PS4"])
        self.assertFalse(evidence.has_evidence())

    def test_cache_set_get_invalidate(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cache = ResultCache(cache_dir=temp_dir)
            key = CacheKey(category="predictor", source="test", variant_id="GRCh38:1-1-A-G")
            cache.set(key, {"revel": 0.5})
            self.assertEqual(cache.get(key), {"revel": 0.5})
            cache.invalidate(key)
            self.assertIsNone(cache.get(key))

    def test_predictor_client_mock_scores(self):
        client = PredictorAPIClient(api_enabled=True, test_mode=True)
        scores = client.get_predictor_scores(chrom="1", pos=1, ref="A", alt="G")
        self.assertTrue(any(score.value is not None for score in scores.values()))

    def test_functional_studies_evaluator_placeholder(self):
        evaluator = FunctionalStudiesEvaluator()
        variant = VariantData(basic_info={"gene": "TEST"})
        self.assertIsNone(evaluator.evaluate_functional_evidence(variant))


if __name__ == "__main__":
    unittest.main()
