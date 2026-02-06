import json
import os
import sys
import unittest
from datetime import datetime

import requests

ROOT_DIR = os.path.dirname(os.path.dirname(__file__))
SRC_DIR = os.path.join(ROOT_DIR, "src")
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from utils.api_client import APIClient
from utils.predictor_api_client import PredictorAPIClient, PopulationAPIClient


DATA_PATH = os.path.join(ROOT_DIR, "data", "reference_variants.json")
REPORT_DIR = os.path.join(ROOT_DIR, "test_results")


def _safe_get_status(value):
    if value is None:
        return "no_data"
    if isinstance(value, dict) and value.get("error"):
        return "error"
    return "ok"


def _find_alphamissense_variant() -> dict:
    """Find a variant with AlphaMissense data via MyVariant query."""
    response = requests.get(
        "https://myvariant.info/v1/query",
        params={
            "q": "_exists_:dbnsfp.alphamissense.score AND dbsnp.vartype:snp",
            "size": 10,
            "fields": "_id"
        },
        timeout=20
    )
    response.raise_for_status()
    data = response.json()
    hits = data.get("hits", [])
    if not hits:
        return {}

    import re
    pattern = re.compile(r"^chr(?P<chrom>[^:]+):g\.(?P<pos>\d+)(?P<ref>[ACGT]+)>(?P<alt>[ACGT]+)$")

    for hit in hits:
        hgvs_id = hit.get("_id")
        if not hgvs_id:
            continue
        match = pattern.match(hgvs_id)
        if not match:
            continue
        return {
            "hgvs_id": hgvs_id,
            "chrom": match.group("chrom"),
            "pos": int(match.group("pos")),
            "ref": match.group("ref"),
            "alt": match.group("alt")
        }

    return {}


@unittest.skipUnless(
    os.getenv("ACMG_INTEGRATION") == "1",
    "Integration tests disabled. Set ACMG_INTEGRATION=1 to enable."
)
class TestIntegrationValidationReport(unittest.TestCase):
    def test_build_validation_report(self):
        if not os.path.exists(DATA_PATH):
            self.fail(f"Reference dataset not found: {DATA_PATH}")

        with open(DATA_PATH, "r", encoding="utf-8") as dataset_file:
            dataset = json.load(dataset_file)

        predictor_client = PredictorAPIClient(api_enabled=True, timeout=20, test_mode=False)
        population_client = PopulationAPIClient(api_enabled=True, timeout=20, test_mode=False)
        api_client = APIClient(cache_enabled=True)

        results = {
            "generated_at": datetime.utcnow().isoformat() + "Z",
            "variants": [],
            "hpo_terms": [],
            "ena_accessions": [],
            "summary": {
                "total_variants": 0,
                "variants_with_preferred_data": 0,
                "variants_missing_preferred_data": 0,
                "api_errors": 0
            }
        }

        alphamissense_found = 0
        gnomad_found = 0

        for entry in dataset.get("variants", []):
            variant_result = {
                "id": entry.get("id"),
                "chromosome": entry.get("chromosome"),
                "position": entry.get("position"),
                "ref_allele": entry.get("ref_allele"),
                "alt_allele": entry.get("alt_allele"),
                "gene": entry.get("gene"),
                "preferred_sources": entry.get("preferred_sources", []),
                "sources": {}
            }

            chrom = entry.get("chromosome")
            pos = entry.get("position")
            ref = entry.get("ref_allele")
            alt = entry.get("alt_allele")

            predictor_scores = predictor_client.get_predictor_scores(
                chrom=str(chrom),
                pos=int(pos),
                ref=str(ref),
                alt=str(alt)
            )
            predictor_available = any(
                score.value is not None for score in predictor_scores.values()
            )
            if predictor_scores.get("alphamissense") and predictor_scores["alphamissense"].value is not None:
                alphamissense_found += 1
            variant_result["sources"]["myvariant"] = {
                "status": "ok" if predictor_available else "no_data",
                "available_scores": [
                    name for name, score in predictor_scores.items() if score.value is not None
                ]
            }

            population_stats = population_client.get_population_stats(
                chrom=str(chrom),
                pos=int(pos),
                ref=str(ref),
                alt=str(alt)
            )
            population_available = False
            if population_stats:
                for stats in population_stats.values():
                    if getattr(stats, "af", None) is not None or getattr(stats, "ac", None) is not None:
                        population_available = True
                        break
            if population_available:
                gnomad_found += 1
            variant_result["sources"]["gnomad"] = {
                "status": "ok" if population_available else "no_data",
                "sources": list(population_stats.keys()) if population_stats else []
            }

            clinvar_data = api_client.get_clinvar_status(
                str(chrom), int(pos), str(ref), str(alt)
            )
            clinvar_status = clinvar_data.get("status") if isinstance(clinvar_data, dict) else "error"
            variant_result["sources"]["clinvar"] = {
                "status": clinvar_status,
                "clinvar_id": clinvar_data.get("clinvar_id") if isinstance(clinvar_data, dict) else None,
                "significance": clinvar_data.get("significance") if isinstance(clinvar_data, dict) else None
            }

            conservation_data = api_client.get_conservation_scores(
                str(chrom), int(pos), str(ref), str(alt)
            )
            conservation_scores = conservation_data.get("conservation_scores", {}) if isinstance(conservation_data, dict) else {}
            variant_result["sources"]["ensembl"] = {
                "status": "ok" if conservation_scores else _safe_get_status(conservation_data),
                "scores": conservation_scores
            }

            preferred_sources = set(entry.get("preferred_sources", []))
            preferred_ok = False
            api_errors = 0
            for source_name in preferred_sources:
                source_data = variant_result["sources"].get(source_name, {})
                status = source_data.get("status")
                if status in ("ok", "found"):
                    preferred_ok = True
                if status == "error":
                    api_errors += 1

            results["summary"]["api_errors"] += api_errors
            results["summary"]["total_variants"] += 1
            if preferred_ok:
                results["summary"]["variants_with_preferred_data"] += 1
            else:
                results["summary"]["variants_missing_preferred_data"] += 1

            results["variants"].append(variant_result)

        for term in dataset.get("hpo_terms", []):
            term_id = term.get("id")
            term_label = term.get("label", "")
            response = requests.get(
                "https://www.ebi.ac.uk/ols4/api/ontologies/hp/terms",
                params={"iri": f"http://purl.obolibrary.org/obo/{term_id.replace(':', '_')}"},
                timeout=20
            )
            status = "ok" if response.status_code == 200 else "error"
            label_found = False
            if response.status_code == 200:
                data = response.json()
                embedded = data.get("_embedded", {})
                terms = embedded.get("terms", []) if isinstance(embedded, dict) else []
                if terms:
                    label = terms[0].get("label", "")
                    label_found = term_label.lower() in label.lower()
            results["hpo_terms"].append({
                "id": term_id,
                "expected_label": term_label,
                "status": status,
                "label_match": label_found
            })

        ena_ok = 0
        for accession in dataset.get("ena_accessions", []):
            acc_id = accession.get("id")
            if not acc_id:
                continue
            response = requests.get(
                f"https://www.ebi.ac.uk/ena/browser/api/summary/{acc_id}",
                timeout=20
            )
            ok = response.status_code == 200 and response.text.strip() and response.text.strip() != "[]"
            if ok:
                ena_ok += 1
            results["ena_accessions"].append({
                "id": acc_id,
                "status": "ok" if ok else "no_data",
                "http_status": response.status_code
            })

        # Ensure at least one variant yields AlphaMissense via MyVariant/dbNSFP
        alphamissense_variant = _find_alphamissense_variant()
        if alphamissense_variant:
            scores = predictor_client.get_predictor_scores(
                chrom=str(alphamissense_variant["chrom"]),
                pos=int(alphamissense_variant["pos"]),
                ref=str(alphamissense_variant["ref"]),
                alt=str(alphamissense_variant["alt"])
            )
            if scores.get("alphamissense") and scores["alphamissense"].value is not None:
                alphamissense_found += 1

        os.makedirs(REPORT_DIR, exist_ok=True)
        json_report = os.path.join(REPORT_DIR, "integration_validation_report.json")
        md_report = os.path.join(REPORT_DIR, "integration_validation_report.md")

        with open(json_report, "w", encoding="utf-8") as out_file:
            json.dump(results, out_file, indent=2)

        with open(md_report, "w", encoding="utf-8") as out_file:
            out_file.write("# Integration Validation Report\n\n")
            out_file.write(f"Generated: {results['generated_at']}\n\n")
            out_file.write("## Summary\n")
            out_file.write(f"- Total variants: {results['summary']['total_variants']}\n")
            out_file.write(f"- Variants with preferred data: {results['summary']['variants_with_preferred_data']}\n")
            out_file.write(f"- Variants missing preferred data: {results['summary']['variants_missing_preferred_data']}\n")
            out_file.write(f"- API errors: {results['summary']['api_errors']}\n\n")

            out_file.write("## Variants\n")
            for variant in results["variants"]:
                out_file.write(f"- {variant['id']} ({variant['chromosome']}:{variant['position']} {variant['ref_allele']}>{variant['alt_allele']})\n")
                for source_name, source_data in variant["sources"].items():
                    out_file.write(f"  - {source_name}: {source_data.get('status')}\n")
            out_file.write("\n")

            out_file.write("## HPO Terms\n")
            for term in results["hpo_terms"]:
                out_file.write(f"- {term['id']}: status={term['status']}, label_match={term['label_match']}\n")

        self.assertGreater(
            results["summary"]["variants_with_preferred_data"],
            0,
            "No variants returned data from preferred sources."
        )
        self.assertGreater(
            gnomad_found,
            0,
            "No variants returned population data from gnomAD."
        )
        self.assertGreater(
            alphamissense_found,
            0,
            "No AlphaMissense scores found via MyVariant/dbNSFP."
        )
        if dataset.get("ena_accessions"):
            self.assertGreater(
                ena_ok,
                0,
                "No ENA accessions returned summary data."
            )
