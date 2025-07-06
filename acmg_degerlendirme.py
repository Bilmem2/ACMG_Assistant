# --- Input Collection Functions ---
import os
import sys
import re
import requests
import numpy as np
from typing import Dict, Optional, Tuple

def collect_compound_het_info(variant_number=2):
    print(f"\n{'='*15} COMPOUND HETEROZYGOUS VARIANT INFORMATION {'='*15}")
    print("For recessive diseases, indicate if there is a second heterozygous variant.")
    has_second_variant = validate_choice_input("Is there a second variant? (yes/no)", ["yes", "no"])
    if has_second_variant == "yes":
        return collect_variant_basic_info(variant_number=variant_number)
    return None

def collect_inputs(test_mode=False):
    global first_time
    first_time = True
    inputs: Dict[str, object] = {'variants': []}
    primary_variant = collect_variant_basic_info(test_mode=test_mode, variant_number=1)
    primary_variant.update(collect_population_data())
    primary_variant.update(collect_insilico_inputs(primary_variant['variant_type'], primary_variant, primary_variant['gene']))
    primary_variant.update(collect_genetic_data())
    primary_variant.update(collect_functional_data(primary_variant['variant_type'], primary_variant['gene']))
    if primary_variant['inheritance'] == 'ar' and primary_variant['zygosity'] == 'heterozygous':
        second_variant = collect_compound_het_info()
        if second_variant:
            second_variant['clinvar_status'], second_variant['clinvar_note'] = get_clinvar_status(
                second_variant['chromosome'], second_variant['position'], second_variant['ref_allele'], second_variant['alt_allele']
            )
            second_variant.update(collect_insilico_inputs(second_variant['variant_type'], second_variant, second_variant['gene']))
            if 'variants' not in inputs or not isinstance(inputs['variants'], list):
                inputs['variants'] = []
            inputs['variants'].append(second_variant)
    if 'variants' not in inputs or not isinstance(inputs['variants'], list):
        inputs['variants'] = []
    inputs['variants'].append(primary_variant)
    print(f"\n{'='*15} ADDITIONAL NOTES {'='*15}")
    print("Provide any literature references or additional notes.")
    inputs['notes'] = str(input("Enter notes (press Enter if none): ") or 'None')
    inputs['varsome_url'] = str(generate_varsome_url(primary_variant['chromosome'], primary_variant['position'], primary_variant['ref_allele'], primary_variant['alt_allele']))
    return inputs

# (Tüm fonksiyonlar burada, ana bloktan önce ve sıralı şekilde tanımlı olmalı)

# --- Entry Point ---
def process_variants():
    """Process variants and classify them."""
    print("\n" + "="*60)
    print("      ACMG 2015/2023 VARIANT CLASSIFICATION ASSISTANT")
    print("                 Created by: Can Sevilmiş")
    print("="*60 + "\n")
    while True:
        inputs = collect_inputs()
        evidence = assign_evidence(inputs, use_pp5_bp6=True, use_acmg_2023=True)
        classification = classify_variant(evidence, use_acmg_2023=True)
        suggestions = check_missing_data(inputs)[2]
        save_results(inputs, classification, evidence, suggestions)
        print_results(inputs, classification, evidence, suggestions)
        another_variant = validate_choice_input("Classify another variant? (yes/no)", ["yes", "no"])
        if another_variant == "no":
            break
        clear_cache()

def print_terminal_warnings_and_suggestions(report_data):
    # Print final ACMG classification at the end, colored
    final_class = str(report_data.get('classification', '[MISSING]'))
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    RESET = '\033[0m'
    if 'benign' in final_class.lower():
        color_code = GREEN
    elif 'pathogenic' in final_class.lower():
        color_code = RED
    else:
        color_code = YELLOW
    print(f"\n{color_code}FINAL ACMG CLASSIFICATION: {final_class}{RESET}")
    """Print missing/optional field warnings and suggestions in color, grouped at the end."""
    RED = '\033[91m'
    YELLOW = '\033[93m'
    CYAN = '\033[96m'
    RESET = '\033[0m'

    # Define required and optional fields for each section (with suggestions)
    required_fields = [
        ('gene', 'Gene', 'Obtain from sequencing report or VCF/annotation.'),
        ('chromosome', 'Chromosome', 'Obtain from VCF or annotation tools.'),
        ('position', 'Position', 'Obtain from VCF or annotation tools.'),
        ('ref_allele', 'Reference Allele', 'Obtain from VCF or annotation tools.'),
        ('alt_allele', 'Alternate Allele', 'Obtain from VCF or annotation tools.'),
        ('cdna_change', 'cDNA Change', 'Obtain from sequencing report or annotation.'),
        ('protein_change', 'Protein Change', 'Obtain from sequencing report or annotation.'),
        ('vep_consequence', 'VEP Consequence', 'Use Ensembl VEP or similar annotation tools.'),
        ('variant_type', 'Final Variant Type', 'User-selected or auto-detected.'),
        ('allele_frequency', 'Allele Frequency', 'Check gnomAD: https://gnomad.broadinstitute.org/.'),
        ('inheritance', 'Inheritance', 'Check inheritance pattern in family or literature.'),
        ('zygosity', 'Zygosity', 'Check zygosity from sequencing data.'),
        ('functional_test', 'Functional Test', 'Search for functional test results in PubMed.'),
        ('phenotype_match', 'Phenotype Match', 'Check phenotype match in Orphanet, OMIM, or clinical databases.'),
    ]
    optional_fields = [
        ('date', 'Date/Time', 'Set automatically by the tool.'),
        ('varsome_url', 'Varsome URL', 'Generated by the tool; check https://varsome.com/.'),
        ('clinvar_status', 'ClinVar Status', 'Check https://www.ncbi.nlm.nih.gov/clinvar/.'),
        ('disease_prevalence', 'Disease Prevalence', 'Check Orphanet: https://www.orpha.net/.'),
        ('subpopulation_frequencies', 'Subpopulation Frequencies', 'Check gnomAD subpopulations.'),
        ('allelic_data', 'Allelic State', 'Check for trans/cis status of variants in recessive inheritance.'),
        ('family_history', 'Family History', 'Check family history in clinical records.'),
        ('segregation', 'Segregation', 'Check family segregation data in genetic reports or publications.'),
        ('de_novo', 'De Novo Status', 'Check for de novo status in family studies or trio sequencing.'),
        ('parental_testing', 'Parental Testing', 'Check if parental testing was performed.'),
        ('consanguinity', 'Consanguinity', 'Check for consanguinity in family history.'),
        ('affected_siblings', 'Affected Siblings', 'Check for affected siblings in family.'),
        ('genetic_notes', 'Additional Genetic Notes', 'Add any relevant notes.'),
        ('paternity_confirmed', 'Paternity Confirmed', 'Check for paternity confirmation in de novo cases.'),
    ]

    # Collect missing fields and suggestions
    missing_required = []
    missing_optional = []
    suggestions = []
    for field, label, suggestion in required_fields:
        val = report_data.get(field, 'N/A')
        if val is None or val == '' or val == 'N/A' or val == 'Not specified' or val == 0:
            missing_required.append((label, suggestion))
            suggestions.append(f"{label}: {suggestion}")
    for field, label, suggestion in optional_fields:
        val = report_data.get(field, 'N/A')
        if val is None or val == '' or val == 'N/A' or val == 'Not specified' or val == 0:
            missing_optional.append((label, suggestion))
            suggestions.append(f"{label}: {suggestion}")

    # Print warnings for missing required fields
    if missing_required:
        print(f"{RED}WARNING: Missing required fields:{RESET}")
        for label, suggestion in missing_required:
            print(f"  {RED}[MISSING] {label}: Not provided. Suggestion: {suggestion}{RESET}")

    # Print warnings for missing optional fields
    if missing_optional:
        print(f"{YELLOW}Note: Missing optional fields:{RESET}")
        for label, suggestion in missing_optional:
            print(f"  {YELLOW}[MISSING] {label}: Not provided. Suggestion: {suggestion}{RESET}")

    # Print all suggestions at the end
    if suggestions:
        print(f"\n{CYAN}SUGGESTIONS:{RESET}")
        for s in sorted(set(suggestions)):
            print(f"  {CYAN}{s}{RESET}")
# --- Functional and Clinical Data Collection ---
def collect_functional_data(variant_type, gene):
    """Collect functional and clinical data with user-friendly prompts."""
    inputs = {}
    print(f"\n{'='*15} FUNCTIONAL AND CLINICAL DATA {'='*15}")
    print("Provide functional test results and clinical observations for the variant.")

    # Functional test section
    if variant_type == 'missense':
        inputs['functional_test'] = validate_choice_input(
            "Enter protein activity test result\n - Loss: Reduced or no protein function\n - Normal: Normal protein function",
            ["Loss", "Normal", "unknown"]
        )
    elif variant_type in ['splice_donor', 'splice_acceptor', 'intronic']:
        inputs['functional_test'] = validate_choice_input(
            "Enter splicing alteration result\n - Yes: Splicing is altered\n - No: Splicing is unaffected",
            ["Yes", "No", "unknown"]
        )
    elif variant_type in ['nonsense', 'frameshift']:
        inputs['functional_test'] = validate_choice_input(
            "Enter NMD (Nonsense-Mediated Decay) status\n - Yes: NMD is triggered\n - No: NMD is not triggered",
            ["Yes", "No", "unknown"]
        )
    else:
        inputs['functional_test'] = validate_choice_input(
            "Enter functional test result\n - Loss: Reduced or no function\n - Normal: Normal function",
            ["Loss", "Normal", "unknown"]
        )



    # Phenotype match section
    print(f"\n{'='*15} PHENOTYPE MATCH {'='*15}")
    print("How well does the patient's phenotype match the gene?")
    inputs['phenotype_match'] = validate_choice_input(
        f"Enter phenotype match strength for {gene} gene\n - Strong: Highly specific to gene\n - Moderate: Partially specific\n - Weak: Low specificity",
        ["Strong", "Moderate", "Weak", "unknown"]
    )
    # Remove de novo and parental testing prompts from here (they are already in collect_genetic_data)
    return inputs

"""
MIT License
Copyright (c) 2025 Can Sevilmiş
See LICENSE file for full license text.
"""

import time
import argparse
import requests
import scipy.stats as stats
import numpy as np
import re



from typing import Dict, Optional, Tuple


# --- Imports (ensure all required modules are imported) ---
import os
import sys
import time
import traceback
import inspect
import re
import requests
from typing import Dict, Optional, Tuple

# Remove all logging usage (no logging to file, only print or silent error handling)
def _log_info(msg):
    pass
def _log_warning(msg):
    pass
def _log_error(msg):
    pass

# --- Report Writing ---
def write_structured_report(report_data: dict, filename: str = "variant_classification_report.txt"):
    """Write a structured, human-readable report for the ACMG assessment."""
    import os
    out_path = os.path.join(os.getcwd(), filename)
    try:
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("\n" + "="*70 + "\n")
            f.write(" ACMG VARIANT CLASSIFICATION REPORT\n")
            f.write("="*70 + "\n\n")

            # [VARIANT]
            f.write("[VARIANT]\n")
            f.write(f"  Gene:         {report_data.get('gene', '-'):<12}  cDNA:   {report_data.get('cdna_change', '-'):<18}\n")
            f.write(f"  Chromosome:   {report_data.get('chromosome', '-'):<12}  Protein: {report_data.get('protein_change', '-'):<18}\n")
            f.write(f"  Position:     {report_data.get('position', '-'):<12}  VEP:    {report_data.get('vep_consequence', '-'):<18}\n")
            f.write(f"  Ref/Alt:      {report_data.get('ref_allele', '-')}/{report_data.get('alt_allele', '-'): <10}  Type:   {report_data.get('variant_type', '-'):<18}\n")
            f.write(f"  Date:         {report_data.get('date', '-')}\n\n")

            # [POPULATION]
            f.write("[POPULATION]\n")
            f.write(f"  Allele Freq:        {report_data.get('allele_frequency', '-')}\n")
            f.write(f"  Disease Prevalence: {report_data.get('disease_prevalence', '-')}\n")
            subpops = report_data.get('subpopulation_frequencies', [])
            if isinstance(subpops, list):
                subpop_str = ', '.join(str(x) for x in subpops) if subpops else '-'
            else:
                subpop_str = str(subpops) if subpops else '-'
            f.write(f"  Subpop Freqs:       {subpop_str}\n\n")

            # [IN SILICO PREDICTIONS]
            f.write("[IN SILICO PREDICTIONS]\n")
            f.write("  Tool            Score      Prediction\n")
            f.write("  -----------------------------------------\n")
            insilico_scores = report_data.get('insilico_scores', {})
            insilico_preds = report_data.get('insilico_predictions', {})
            for tool in [
                'revel','cadd','metarnn','clinpred','bayesdel','alphamissense','mutationtaster','polyphen2','sift','fathmm_xf','mutationassessor','provean','mutpred','metolr','esm1b','lrt','phylop_vert','phylop_mamm','phylop_primate','gerp++']:
                score = insilico_scores.get(tool, None)
                pred = insilico_preds.get(f"{tool}_prediction", None)
                if score is not None and (score != 0 or pred not in [None, '-', '', 'unknown']):
                    score_str = str(score) if score not in [None, '', 'N/A'] else '-'
                    pred_str = pred if pred not in [None, '', 'N/A'] else '-'
                    f.write(f"  {tool.upper():<15} {score_str:<10} {pred_str}\n")
            f.write(f"  Combined MetaScore: {report_data.get('metascore', '-')}\n\n")

            # [GENETIC & FAMILY]
            f.write("[GENETIC & FAMILY]\n")
            f.write(f"  Inheritance:      {report_data.get('inheritance', '-'):<12}  Zygosity:      {report_data.get('zygosity', '-'):<12}\n")
            f.write(f"  Allelic State:    {report_data.get('allelic_data', '-'):<12}  Family Hist:   {report_data.get('family_history', '-'):<12}\n")
            f.write(f"  Segregation:      {report_data.get('segregation', '-'):<12}  De Novo:       {report_data.get('de_novo', '-'):<12}\n")
            f.write(f"  Parental Testing: {report_data.get('parental_testing', '-'):<12}  Consanguinity: {report_data.get('consanguinity', '-'):<12}\n")
            f.write(f"  Siblings:         {report_data.get('affected_siblings', '-')}\n")
            notes = report_data.get('genetic_notes', '-')
            f.write(f"  Notes:            {notes if notes else '-'}\n\n")

            # [FUNCTIONAL/CLINICAL]
            f.write("[FUNCTIONAL/CLINICAL]\n")
            f.write(f"  Functional Test:     {report_data.get('functional_test', '-')}\n")
            f.write(f"  Phenotype Match:     {report_data.get('phenotype_match', '-')}\n")
            f.write(f"  Paternity Confirmed: {report_data.get('paternity_confirmed', '-')}\n\n")

            # [ACMG CLASSIFICATION]
            f.write("[ACMG CLASSIFICATION]\n")
            crit = report_data.get('criteria', [])
            if isinstance(crit, list):
                crit_str = ', '.join(crit) if crit else '-'
            else:
                crit_str = str(crit) if crit else '-'
            f.write(f"  Criteria:   {crit_str}\n")
            f.write(f"  Final:      {report_data.get('classification', '-')}\n\n")

            # [REFERENCES]
            f.write("[REFERENCES]\n")
            f.write(f"  Varsome:    {report_data.get('varsome_url', '-')}\n")
            f.write(f"  ClinVar:    {report_data.get('clinvar_status', '-')}\n")
            notes = report_data.get('notes', '-')
            f.write(f"  Notes:      {notes if notes else '-'}\n\n")

            f.write("="*70 + "\n")
            f.flush()
    except Exception as e:
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("[ERROR] Failed to write ACMG report.\n")
            f.write(str(e) + "\n")
        print("[ERROR] Exception in write_structured_report:", str(e))

# --- Caching Mechanisms ---
clinvar_cache = {}
score_cache = {}

# --- Constants ---
EVIDENCE_WEIGHTS = {
    'PVS1': 'Very Strong Pathogenic', 'PS1': 'Strong Pathogenic', 'PS2': 'Strong Pathogenic',
    'PS2_Very_Strong': 'Very Strong Pathogenic', 'PS3': 'Strong Pathogenic', 'PM2': 'Moderate Pathogenic',
    'PM3': 'Moderate Pathogenic', 'PP1': 'Supporting Pathogenic', 'PP3': 'Supporting Pathogenic',
    'PP4': 'Supporting Pathogenic', 'PP5': 'Supporting Pathogenic (Deprecated)', 'BA1': 'Standalone Benign',
    'BS1': 'Strong Benign', 'BS3': 'Strong Benign', 'BS4': 'Strong Benign',
    'BP4': 'Supporting Benign', 'BP6': 'Supporting Benign (Deprecated)'
}
GENE_SPECIFIC_THRESHOLDS = {
    'BRCA1': {'BA1': 0.05, 'BS1': 0.01}, 'TTN': {'BA1': 0.1, 'BS1': 0.05},
    'CAMTA1': {'BA1': 0.05, 'BS1': 0.01}, 'default': {'BA1': 0.05, 'BS1': 0.01}
}
GENE_SPECIFIC_WEIGHTS = {
    'default': {
        'revel': 0.25, 'cadd': 0.20, 'metarnn': 0.15, 'clinpred': 0.15, 'bayesdel': 0.10,
        'alphamissense': 0.10, 'mutationtaster': 0.05, 'polyphen2': 0.05, 'sift': 0.05
    },
    'CAMTA1': {
        'revel': 0.30, 'cadd': 0.25, 'metarnn': 0.20, 'clinpred': 0.15, 'bayesdel': 0.10,
        'alphamissense': 0.10, 'mutationtaster': 0.00, 'polyphen2': 0.00, 'sift': 0.00
    },
    'BRCA1': {
        'revel': 0.35, 'cadd': 0.20, 'metarnn': 0.20, 'clinpred': 0.15, 'bayesdel': 0.10,
        'alphamissense': 0.10, 'mutationtaster': 0.00, 'polyphen2': 0.00, 'sift': 0.00
    }
}
INSILICO_RANGES = {
    'revel': (0.0, 1.0, False), 'cadd': (0.0, 100.0, False), 'metarnn': (0.0, 1.0, False),
    'clinpred': (0.0, 1.0, False), 'bayesdel': (-1.0, 1.0, False), 'alphamissense': (0.0, 1.0, False),
    'mutationtaster': (0.0, 1.0, False), 'polyphen2': (0.0, 1.0, False), 'sift': (0.0, 1.0, True),
    'fathmm_xf': (0.0, 1.0, False), 'mutationassessor': (0.0, 5.0, False), 'provean': (-14.0, 14.0, True),
    'mutpred': (0.0, 1.0, False), 'metolr': (0.0, 1.0, False), 'esm1b': (-100.0, 100.0, False),
    'lrt': (0.0, 1.0, True), 'gerp++': (-12.3, 6.17, False),
    # phyloP conservation scores now use ranked score (0-1)
    'phylop_vert': (0.0, 1.0, False), 'phylop_mamm': (0.0, 1.0, False), 'phylop_primate': (0.0, 1.0, False)
}
LOF_VARIANTS = {'nonsense', 'frameshift', 'splice_donor', 'splice_acceptor'}
first_time = True


def calculate_metascore(inputs: Dict, gene: str) -> Tuple[float, Optional[str], dict]:
    """Calculate a VAMPP-score-like metascore and determine PP3/BP4 criterion. Returns details for reporting."""
    cache_key = tuple((k, inputs.get(k)) for k in sorted(GENE_SPECIFIC_WEIGHTS.get(gene, GENE_SPECIFIC_WEIGHTS['default']).keys()) if inputs.get(f'{k}_score') is not None)
    if cache_key in score_cache:
        metascore, criterion, details = score_cache[cache_key]
        print(f"Info: Combined MetaScore retrieved from cache: {metascore:.4f} ({criterion if criterion else 'No criterion applied'})")
        return metascore, criterion, details

    weights = GENE_SPECIFIC_WEIGHTS.get(gene, GENE_SPECIFIC_WEIGHTS['default'])
    scores = []
    used_weights = []
    used_tools = []
    normalized_scores = {}

    for tool in weights:
        score = inputs.get(f'{tool}_score')
        normalized = normalize_score(tool, score)
        if normalized is not None:
            scores.append(normalized)
            used_weights.append(weights[tool])
            used_tools.append(tool)
            normalized_scores[tool] = round(normalized, 4)

    if not scores:
        return 0.0, None, {"used_tools": [], "normalized_scores": {}, "weights": {}, "weighted_sum": 0.0, "total_weight": 0.0}

    total_weight = sum(used_weights)
    if total_weight == 0:
        return 0.0, None, {"used_tools": used_tools, "normalized_scores": normalized_scores, "weights": used_weights, "weighted_sum": 0.0, "total_weight": 0.0}

    weighted_sum = sum(s * w for s, w in zip(scores, used_weights))
    metascore = weighted_sum / total_weight

    criterion = 'PP3' if metascore > 0.354 else 'BP4' if metascore < 0.226 else None
    details = {
        "used_tools": used_tools,
        "normalized_scores": normalized_scores,
        "weights": {t: w for t, w in zip(used_tools, used_weights)},
        "weighted_sum": round(weighted_sum, 4),
        "total_weight": round(total_weight, 4),
        "final_metascore": round(metascore, 4),
        "criterion": criterion
    }
    score_cache[cache_key] = (metascore, criterion, details)

    return metascore, criterion, details

def get_chromosome_from_ensembl(gene_name):
    """Retrieve chromosome information for a gene using Ensembl API."""
    if not gene_name or gene_name == "Not specified":
        return "Not specified"
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/homo_sapiens/{gene_name}?expand=0"
    try:
        r = requests.get(server + ext, headers={"Content-Type": "application/json"}, timeout=10)
        r.raise_for_status()
        decoded = r.json()
        return decoded.get("seq_region_name", "Not specified")
    except requests.exceptions.RequestException as e:
        _log_warning(f"Ensembl API request failed: {e}. Chromosome information not retrieved.")
        return "Not specified"

def get_clinvar_status(chromosome, position, ref_allele, alt_allele):
    """Fetch ClinVar status for a variant."""
    if chromosome == "Not specified" or position == 0 or ref_allele == "Not specified" or alt_allele == "Not specified":
        return "unknown", "ClinVar data could not be retrieved: Missing information."
    cache_key = f"{chromosome}-{position}-{ref_allele}-{alt_allele}"
    if cache_key in clinvar_cache:
        _log_info("ClinVar status retrieved from cache.")
        return clinvar_cache[cache_key]
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    query = f"{chromosome}:{position}[chrpos38] AND {ref_allele}>{alt_allele}"
    params = {'db': 'clinvar', 'term': query, 'retmode': 'json', 'retmax': 10}
    try:
        r = requests.get(base_url, params=params, timeout=10)
        r.raise_for_status()
        result = r.json()
        ids = result.get('esearchresult', {}).get('idlist', [])
        if not ids:
            result_tuple = ("unknown", "No variant found in ClinVar.")
        else:
            clinvar_status = "unknown"
            clinvar_note = ""
            # Try all returned ids, pick the first with a valid clinical_significance
            for cid in ids:
                detail_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={cid}&retmode=json"
                r_detail = requests.get(detail_url, timeout=10)
                r_detail.raise_for_status()
                detail = r_detail.json()
                docsum = detail.get('result', {}).get(cid, {})
                cs = docsum.get('clinical_significance')
                if cs:
                    # Sometimes cs is a dict, sometimes a string
                    if isinstance(cs, dict):
                        clinvar_status = cs.get('description', 'unknown').lower()
                        submitters = cs.get('submitter', [])
                        submitter_note = " (Single submitter)" if isinstance(submitters, list) and len(submitters) <= 1 else " (Consensus)"
                    else:
                        clinvar_status = str(cs).lower()
                        submitter_note = ""
                    clinvar_note = f"ClinVar status: {clinvar_status}{submitter_note}"
                    break
            if clinvar_status == "unknown":
                clinvar_note = "No clinical significance found in ClinVar."
            result_tuple = (clinvar_status, clinvar_note)
        clinvar_cache[cache_key] = result_tuple
        return result_tuple
    except requests.exceptions.RequestException as e:
        _log_error(f"ClinVar API request failed: {e}")
        return "unknown", f"ClinVar API request failed: {e}"

def validate_float_input(prompt, allow_zero=True, min_value=-100.0, max_value=100.0, reference_url=None, tool=None):
    """Validate float input within a specified range with enhanced user guidance."""
    global first_time
    # Araç bazlı örnek değerler
    tool_examples: Dict[str, str] = {
        'revel': '0.75', 'cadd': '25.0', 'metarnn': '0.6', 'clinpred': '0.6', 'bayesdel': '0.2',
        'alphamissense': '0.6', 'mutationtaster': '0.6', 'polyphen2': '0.6', 'sift': '0.01',
        'fathmm_xf': '0.6', 'mutationassessor': '2.0', 'provean': '-3.0', 'mutpred': '0.6',
        'metolr': '0.6', 'esm1b': '1.0', 'lrt': '0.01', 'gerp++': '3.0',
        'phylop_vert': '0.9', 'phylop_mamm': '0.9', 'phylop_primate': '0.9',
        'spliceai_acceptor_gain': '0.6', 'spliceai_acceptor_loss': '0.6',
        'spliceai_donor_gain': '0.6', 'spliceai_donor_loss': '0.6', 'loftool': '0.6'
    }
    example_value: str = tool_examples[tool] if tool in tool_examples else (f"{min_value + 0.1:.2f}" if min_value != -100.0 else "0.01")
    example_fraction = "1/1000" if max_value <= 1.0 else "5/2"
    reference_text = f" - Check reference: {reference_url}\n" if reference_url else ""
    
    full_prompt = (
        f"{prompt}\n"
        f" - Expected range: {min_value} to {max_value}\n"
        f"{reference_text}"
        f" - Example: {example_value}, {example_fraction}, or 1e-3 (scientific notation)\n"
        f" - Press Enter if unknown{' (sets to 0.0)' if allow_zero else ''}: "
    ) if first_time else (
        f"{prompt}\n"
        f" - Example: {example_value}, {example_fraction}\n"
        f" - Press Enter if unknown{' (sets to 0.0)' if allow_zero else ''}: "
    )
    
    attempt_count = 0
    max_attempts = 3
    
    while attempt_count < max_attempts:
        try:
            value = input(full_prompt).strip()
            if not value:
                if allow_zero:
                    _log_info("No value provided, setting to 0.0 (unknown).")
                    first_time = False
                    return 0.0
                else:
                    _log_error("A value is required for this field.")
                    attempt_count += 1
                    continue
            
            if '/' in value:
                parts = value.split('/')
                if len(parts) != 2:
                    _log_error("Fraction must be in the format 'number/number' (e.g., 1/1000).")
                    attempt_count += 1
                    continue
                num = float(parts[0].replace(',', '').replace('<', ''))
                den = float(parts[1].replace(',', ''))
                if den == 0:
                    _log_error("Denominator cannot be zero.")
                    attempt_count += 1
                    continue
                result = num / den
            else:
                result = float(value.replace(',', ''))
            
            if not (min_value <= result <= max_value):
                _log_error(f"Value must be between {min_value} and {max_value}. Try {example_value}.")
                attempt_count += 1
                continue
            
            first_time = False
            return result
        
        except (ValueError, IndexError):
            _log_error(f"Enter a valid number, fraction, or scientific notation (e.g., {example_value}, {example_fraction}, 1e-3).")
            attempt_count += 1
    
    _log_warning(f"Maximum attempts ({max_attempts}) reached. Setting to {'0.0' if allow_zero else 'minimum value'}.")
    first_time = False
    return 0.0 if allow_zero else min_value

def validate_int_input(prompt, min_value=0):
    """Validate integer input with enhanced user guidance."""
    global first_time
    example_value = min_value + 1
    full_prompt = (
        f"{prompt}\n"
        f" - Expected: Integer >= {min_value}, e.g., {example_value}\n"
        f" - Press Enter if unknown (sets to 0): "
    ) if first_time else (
        f"{prompt}\n"
        f" - Example: {example_value}\n"
        f" - Press Enter if unknown: "
    )
    
    attempt_count = 0
    max_attempts = 3
    
    while attempt_count < max_attempts:
        try:
            value = input(full_prompt).strip()
            if not value:
                _log_info("No value provided, setting to 0 (unknown).")
                first_time = False
                return 0
            
            result = int(value)
            if result < min_value:
                _log_error(f"Value must be {min_value} or greater. Try {example_value}.")
                attempt_count += 1
                continue
            
            first_time = False
            return result
        
        except ValueError:
            _log_error(f"Enter a valid integer, e.g., {example_value}.")
            attempt_count += 1
    
    _log_warning(f"Maximum attempts ({max_attempts}) reached. Setting to 0.")
    first_time = False
    return 0

def validate_choice_input(prompt, valid_choices):
    """Validate choice input from a list of valid options."""
    global first_time
    valid_choices_str = ", ".join(valid_choices)
    full_prompt = f"{prompt}\n - Options: {valid_choices_str}\n - Press Enter if unknown: " if first_time else f"{prompt}\n - Options: {valid_choices_str}: "
    while True:
        value = input(full_prompt).strip().lower()
        if not value:
            return "unknown"
        if value in [choice.lower() for choice in valid_choices]:
            first_time = False
            return value
        _log_error(f"Enter one of: {valid_choices_str}")

def validate_hgvs_input(prompt, prefix):
    """Validate HGVS format input."""
    hgvs_pattern = re.compile(rf"^{prefix}\..*")
    full_prompt = f"{prompt}\n - Example: {prefix}.181T>G or {prefix}.123_124del\n - Press Enter if unknown: "
    while True:
        value = input(full_prompt).strip()
        if not value:
            return "Not specified"
        if hgvs_pattern.match(value):
            return value
        _log_error(f"Enter a valid HGVS format (e.g., {prefix}.181T>G or {prefix}.123_124del).")

def validate_allele(prompt, field_name):
    """Validate allele input (A, C, G, T)."""
    valid_chars = set('ACGT')
    full_prompt = f"{prompt}\n - Example: A, T, AGC\n - Press Enter if unknown: "
    while True:
        value = input(full_prompt).strip().upper()
        if not value:
            return "Not specified"
        if all(char in valid_chars for char in value):
            return value
        _log_error(f"{field_name} must contain only A, C, G, T characters (e.g., A, AGC).")

def validate_chromosome(prompt, auto_chromosome="Not specified"):
    """Validate chromosome input."""
    valid_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
    full_prompt = f"{prompt}\n - Enter chromosome (e.g., 1, X, Y)\n - Press Enter if unknown: "
    while True:
        value = input(full_prompt).strip().upper()
        if not value:
            return auto_chromosome if auto_chromosome != "Not specified" else "Not specified"
        if value in valid_chromosomes:
            return value
        _log_error(f"Enter a valid chromosome (e.g., 1, X, Y).")

def validate_gene_name(prompt):
    """Validate gene name input (basic check for non-empty string)."""
    full_prompt = f"{prompt}\n - Example: BRCA1, CFTR, CAMTA1\n - Press Enter if unknown: "
    while True:
        value = input(full_prompt).strip().upper()
        if not value:
            _log_error("Gene name is required for accurate classification.")
            continue
        if re.match(r'^[A-Z0-9\-]+$', value):
            return value
        _log_error("Enter a valid HGNC gene symbol (e.g., BRCA1, CFTR, CAMTA1).")

def validate_zygosity_input(prompt):
    """Validate zygosity input."""
    valid_choices = ["heterozygous", "homozygous", "hemizygous", "compound het", "mosaic", "unknown", "het", "hom"]
    # User-friendly display, only show main options
    display_choices = [
        "heterozygous: One copy of the variant (e.g., AD, AR)",
        "homozygous: Two copies of the variant (e.g., AR)",
        "hemizygous: One copy on X/Y chromosome (males, X-linked)",
        "compound het: Two different variants in the same gene (AR)",
        "mosaic: Variant present in a subset of cells",
        "unknown: Zygosity not determined"
    ]
    print("Select zygosity:")
    for idx, desc in enumerate(display_choices, 1):
        print(f"  {idx}. {desc}")
    print("  Press Enter if unknown.")
    mapping = {
        "1": "heterozygous", "2": "homozygous", "3": "hemizygous", "4": "compound het", "5": "mosaic", "6": "unknown",
        "heterozygous": "heterozygous", "homozygous": "homozygous", "hemizygous": "hemizygous", "compound het": "compound het", "mosaic": "mosaic", "unknown": "unknown",
        "het": "heterozygous", "hom": "homozygous"
    }
    while True:
        value = input("Enter choice (1-6) or type option: ").strip().lower()
        if not value:
            return "unknown"
        if value in mapping:
            return mapping[value]
        _log_error("Please enter a number 1-6 or a valid option (e.g., heterozygous, homozygous, etc.)")

def detect_lof_variant(cdna_change, vep_consequence):
    """Detect Loss of Function (LoF) variant based on cDNA change and VEP consequence."""
    if vep_consequence.lower() in LOF_VARIANTS:
        _log_info(f"LoF variant detected by VEP consequence: {vep_consequence.lower()}")
        return vep_consequence.lower()
    if cdna_change != "Not specified":
        cdna_lower = cdna_change.lower()
        if any(x in cdna_lower for x in ['del', 'ins', 'dup', 'frameshift']):
            _log_info("LoF variant detected by cDNA change: frameshift")
            return 'frameshift'
        if any(x in cdna_lower for x in ['ter', 'x', '*']):
            _log_info("LoF variant detected by cDNA change: nonsense")
            return 'nonsense'
        if re.search(r'c\.\d+[+-][12]', cdna_lower):
            _log_info("LoF variant detected by cDNA change: splice site")
            return 'splice_donor' if '+' in cdna_change else 'splice_acceptor'
    return None

def calculate_prediction(tool: str, score: float) -> str:
    """Calculate automatic prediction based on in silico tool score."""
    if score is None or score == 0.0:
        return 'unknown'
    thresholds = {
        'revel': lambda x: 'damaging' if x > 0.7 else 'tolerated' if x < 0.3 else 'unknown',
        'cadd': lambda x: 'damaging' if x > 20 else 'tolerated' if x < 10 else 'unknown',
        'metarnn': lambda x: 'damaging' if x > 0.5 else 'unknown',
        'clinpred': lambda x: 'damaging' if x > 0.5 else 'unknown',
        'bayesdel': lambda x: 'deleterious' if x > 0.15 else 'tolerated' if x < -0.15 else 'unknown',
        'alphamissense': lambda x: 'likely_pathogenic' if x > 0.564 else 'likely_benign' if x < 0.34 else 'ambiguous',
        'mutationtaster': lambda x: 'damaging' if x > 0.5 else 'polymorphism' if x < 0.5 else 'unknown',
        'polyphen2': lambda x: 'probably_damaging' if x > 0.85 else 'possibly_damaging' if x > 0.5 else 'benign',
        'sift': lambda x: 'damaging' if x < 0.05 else 'tolerated' if x >= 0.05 else 'unknown',
        'fathmm_xf': lambda x: 'damaging' if x > 0.5 else 'neutral' if x < 0.5 else 'unknown',
        'mutationassessor': lambda x: 'high' if x > 3.5 else 'medium' if x > 1.9 else 'low' if x > 0.8 else 'neutral',
        'provean': lambda x: 'damaging' if x < -2.5 else 'neutral' if x >= -2.5 else 'unknown',
        'mutpred': lambda x: 'damaging' if x > 0.5 else 'neutral' if x < 0.5 else 'unknown',
        'metolr': lambda x: 'damaging' if x > 0.5 else 'neutral' if x < 0.5 else 'unknown',
        'esm1b': lambda x: 'damaging' if x > 0.0 else 'neutral' if x <= 0.0 else 'unknown',
        'lrt': lambda x: 'damaging' if x < 0.05 else 'neutral' if x >= 0.05 else 'unknown',
        'gerp++': lambda x: 'conserved' if x > 2.0 else 'non-conserved' if x <= 2.0 else 'unknown',
        # phyloP rank score (0-1): >=0.9 conserved, <0.5 non-conserved
        'phylop_vert': lambda x: 'conserved' if x >= 0.9 else 'non-conserved' if x < 0.5 else 'unknown',
        'phylop_mamm': lambda x: 'conserved' if x >= 0.9 else 'non-conserved' if x < 0.5 else 'unknown',
        'phylop_primate': lambda x: 'conserved' if x >= 0.9 else 'non-conserved' if x < 0.5 else 'unknown',
        'spliceai': lambda x: 'damaging' if x > 0.5 else 'tolerated' if x < 0.5 else 'unknown',
        'loftool': lambda x: 'damaging' if x > 0.5 else 'tolerated' if x < 0.5 else 'unknown'
    }
    return thresholds.get(tool, lambda x: 'unknown')(score)

def validate_prediction(tool, score, prediction, valid_choices):
    """Validate in silico prediction against expected value."""
    expected = calculate_prediction(tool, score)
    if prediction != 'unknown' and prediction != expected and expected != 'unknown':
        _log_warning(f"Expected prediction for {tool.upper()} score ({score}) is '{expected}', but '{prediction}' was entered.")
        _log_info(f"Verify with Varsome: https://varsome.com/")
        override = validate_choice_input("Use entered prediction? (yes/no)", ["yes", "no"])
        return prediction if override == "yes" else expected
    return prediction if prediction != 'unknown' else expected

def normalize_score(tool, score):
    """Normalize in silico scores to a 0-1 range."""
    if tool not in INSILICO_RANGES or score is None:
        return None
    min_val, max_val, reverse = INSILICO_RANGES[tool]
    if max_val == min_val:
        return 0.5
    normalized = (score - min_val) / (max_val - min_val)
    return 1.0 - normalized if reverse else normalized

def generate_varsome_url(chromosome, position, ref_allele, alt_allele):
    """Generate a Varsome URL for the variant."""
    if chromosome != "Not specified" and position != 0 and ref_allele != "Not specified" and alt_allele != "Not specified":
        return f"https://varsome.com/variant/hg38/{chromosome}-{position}-{ref_allele}-{alt_allele}"
    return "Varsome URL could not be generated (Missing data)"

def collect_variant_basic_info(test_mode=False, variant_number=1):
    """Collect basic variant information with user-friendly prompts."""
    global first_time
    inputs = {}
    if test_mode:
        # BRCA1 c.5096G>A (p.Arg1699Gln), chr17:43070940 G>A (GRCh38)
        inputs['gene'] = "BRCA1"
        auto_chromosome = get_chromosome_from_ensembl(inputs['gene'])
        inputs['chromosome'] = auto_chromosome if auto_chromosome != "Not specified" else "17"
        inputs['position'] = 43070940
        inputs['ref_allele'] = "G"
        inputs['alt_allele'] = "A"
        inputs['cdna_change'] = "c.5096G>A"
        inputs['protein_change'] = "p.Arg1699Gln"
        inputs['vep_consequence'] = "missense"
        return inputs



    # Gene name
    inputs['gene'] = validate_gene_name("Enter gene name")
    print()
    # Chromosome
    auto_chromosome = get_chromosome_from_ensembl(inputs['gene'])
    if auto_chromosome != "Not specified":
        print(f"Info: Gene {inputs['gene']} found on chromosome {auto_chromosome}.")
        print()
        use_auto = validate_choice_input("Use this chromosome? (yes/no)", ["yes", "no"])
        print()
        inputs['chromosome'] = auto_chromosome if use_auto == "yes" else validate_chromosome("Enter chromosome")
        print()
    else:
        print("Info: Could not retrieve chromosome automatically from Ensembl.")
        print()
        inputs['chromosome'] = validate_chromosome("Enter chromosome")
        print()
    # Position
    inputs['position'] = validate_int_input("Enter genomic position\n - Example: 7249507 (GRCh38)")
    print()
    # Alleles
    inputs['ref_allele'] = validate_allele("Enter reference allele", "Reference allele")
    print()
    inputs['alt_allele'] = validate_allele("Enter alternate allele", "Alternate allele")
    print()
    if inputs['ref_allele'] == inputs['alt_allele'] and inputs['ref_allele'] != "Not specified":
        print("Error: Reference and alternate alleles cannot be the same.")
        print()
        inputs['alt_allele'] = validate_allele("Enter alternate allele again", "Alternate allele")
        print()
    # cDNA/protein
    inputs['cdna_change'] = validate_hgvs_input("Enter cDNA change", "c")
    print()
    inputs['protein_change'] = validate_hgvs_input("Enter protein change", "p")
    print()

    # VEP consequence
    print("\n" + "="*40)
    print(f"{'VEP CONSEQUENCE':^40}")
    print("="*40)
    print("VEP Consequence describes the effect of the variant on the gene.")
    inputs['vep_consequence'] = validate_choice_input(
        "Enter VEP consequence",
        ["missense", "nonsense", "frameshift", "splice_donor", "splice_acceptor", "intronic", "synonymous", "unknown"]
    )
    print()

    # --- Always ask for final variant consequence/type ---
    lof_type = detect_lof_variant(inputs.get('cdna_change', ''), inputs['vep_consequence'])
    print("="*40)
    print(f"{'VARIANT CONSEQUENCE/TYPE':^40}")
    print("="*40)
    if lof_type:
        print(f"Automatic LoF detection: Based on cDNA change and VEP consequence, suggested type is '{lof_type}'.")
        use_lof = validate_choice_input("Use this as the final variant type? (yes/no)", ["yes", "no"])
        if use_lof == "yes":
            inputs['variant_type'] = lof_type
        else:
            inputs['variant_type'] = validate_choice_input(
                "Select the final variant consequence/type",
                ["missense", "nonsense", "frameshift", "splice_donor", "splice_acceptor", "intronic", "synonymous", "unknown"]
            )
    else:
        inputs['variant_type'] = validate_choice_input(
            "Select the final variant consequence/type",
            ["missense", "nonsense", "frameshift", "splice_donor", "splice_acceptor", "intronic", "synonymous", "unknown"]
        )
    print()

    # ClinVar status (moved here for visibility)
    print("\n" + "="*40)
    print(f"{'CLINVAR STATUS':^40}")
    print("="*40)
    status, note = get_clinvar_status(inputs['chromosome'], inputs['position'], inputs['ref_allele'], inputs['alt_allele'])
    inputs['clinvar_status'] = status
    inputs['clinvar_note'] = note
    print(f"Info: {note}")
    print("ClinVar status indicates the known pathogenic/benign classification of the variant.")
    print()
    override_clinvar = validate_choice_input("Manually override ClinVar status? (yes/no)", ["yes", "no"])
    print()
    if override_clinvar == "yes":
        inputs['clinvar_status'] = validate_choice_input(
            "Enter ClinVar status",
            ["pathogenic", "likely pathogenic", "benign", "likely benign", "vus", "unknown"]
        )
        print()
    return inputs

def collect_population_data():
    """Collect population frequency data with user-friendly prompts."""
    inputs = {}
    print(f"\n{'='*15} POPULATION DATA {'='*15}")
    print("Provide frequency data for the variant in the general population (e.g., from gnomAD).")
    inputs['allele_frequency'] = validate_float_input(
        "Enter allele frequency",
        max_value=1.0,
        reference_url="https://gnomad.broadinstitute.org/"
    )
    print(f"\n{'='*15} SUBPOPULATION FREQUENCIES {'='*15}")
    print("Enter gnomAD subpopulation frequencies (e.g., African, European).")
    print("Enter one frequency at a time; press Enter twice to finish.")
    subpopulation_frequencies: list[float] = []
    while True:
        freq = validate_float_input(
            "Enter subpopulation frequency",
            max_value=1.0,
            reference_url="https://gnomad.broadinstitute.org/"
        )
        if freq == 0.0:
            break
        subpopulation_frequencies.append(freq)
    inputs['subpopulation_frequencies'] = subpopulation_frequencies
    print(f"\n{'='*15} DISEASE PREVALENCE {'='*15}")
    print("Provide the prevalence of the associated disease in the population.")
    inputs['disease_prevalence'] = validate_float_input(
        "Enter disease prevalence",
        max_value=1.0,
        reference_url="https://www.orpha.net/"
    )
    print(f"\n{'='*15} CASE-CONTROL DATA FOR FISHER'S EXACT TEST {'='*15}")
    print("Case-control data compares allele frequencies between affected and unaffected populations.")
    run_fisher_test = validate_choice_input("Include case-control data? (yes/no)", ["yes", "no"])
    if run_fisher_test == "yes":
        print(f"\n{'='*15} CASE-CONTROL DATA INPUT {'='*15}")
        print("Provide allele counts and totals for case and control cohorts.")
        inputs['cohort_allele_count'] = validate_int_input("Enter allele count in case cohort")
        inputs['cohort_total'] = validate_int_input("Enter total alleles in case cohort")
        inputs['control_allele_count'] = validate_int_input("Enter allele count in control population")
        inputs['control_total'] = validate_int_input("Enter total alleles in control population")
    else:
        inputs.update({k: 0 for k in ['cohort_allele_count', 'cohort_total', 'control_allele_count', 'control_total']})
    return inputs

def collect_insilico_inputs(variant_type, variant, gene):
    """Collect in silico predictions. All scores must be entered manually by the user. No local dataset search is performed."""
    global first_time
    print(f"\n{'='*15} IN SILICO PREDICTIONS {'='*15}")
    print("You will now enter in silico prediction scores. Each tool is grouped by category. After each entry, an automatic prediction will be shown; you can confirm or override it.")
    inputs = {}
    # Missense
    if variant_type == 'missense':
        tools = [
            ('revel', 'Enter REVEL score\n - Range: 0.0 to 1.0\n - Thresholds: >0.7 (damaging), <0.3 (tolerated)', 0.0, 1.0, ['damaging', 'tolerated', 'unknown']),
            ('cadd', 'Enter CADD Phred score\n - Range: 0.0 to 100.0\n - Threshold: >20 (damaging)', 0.0, 100.0, ['damaging', 'tolerated', 'unknown']),
            ('alphamissense', 'Enter AlphaMissense score\n - Range: 0.0 to 1.0\n - Thresholds: >0.564 (likely_pathogenic), <0.34 (likely_benign)', 0.0, 1.0, ['likely_pathogenic', 'likely_benign', 'ambiguous', 'unknown']),
            ('sift', 'Enter SIFT score\n - Range: 0.0 to 1.0\n - Threshold: <0.05 (damaging)', 0.0, 1.0, ['damaging', 'tolerated', 'unknown']),
            ('bayesdel', 'Enter BayesDel score\n - Range: -1.0 to 1.0\n - Threshold: >0.15 (deleterious)', -1.0, 1.0, ['deleterious', 'tolerated', 'unknown']),
            ('fathmm_xf', 'Enter FATHMM-XF score\n - Range: 0.0 to 1.0\n - Threshold: >0.5 (damaging)', 0.0, 1.0, ['damaging', 'neutral', 'unknown']),
            ('mutationtaster', 'Enter MutationTaster score\n - Range: 0.0 to 1.0\n - Threshold: >0.5 (damaging)', 0.0, 1.0, ['damaging', 'polymorphism', 'unknown']),
            ('mutationassessor', 'Enter MutationAssessor score\n - Range: 0.0 to 5.0\n - Thresholds: >3.5 (high), >1.9 (medium), >0.8 (low)', 0.0, 5.0, ['high', 'medium', 'low', 'neutral', 'unknown']),
            ('provean', 'Enter PROVEAN score\n - Range: -14.0 to 14.0\n - Threshold: <-2.5 (damaging)', -14.0, 14.0, ['damaging', 'neutral', 'unknown']),
            ('mutpred', 'Enter MutPred score\n - Range: 0.0 to 1.0\n - Threshold: >0.5 (damaging)', 0.0, 1.0, ['damaging', 'neutral', 'unknown']),
            ('metolr', 'Enter MetaLR score\n - Range: 0.0 to 1.0\n - Threshold: >0.5 (damaging)', 0.0, 1.0, ['damaging', 'neutral', 'unknown']),
            ('esm1b', 'Enter ESM1b LLR score\n - Range: -100.0 to 100.0', -100.0, 100.0, ['damaging', 'neutral', 'unknown']),
            ('lrt', 'Enter LRT p-value\n - Range: 0.0 to 1.0\n - Threshold: <0.05 (damaging)', 0.0, 1.0, ['damaging', 'neutral', 'unknown']),
        ]
        for tool, prompt, min_val, max_val, choices in tools:
            inputs[f'{tool}_score'] = validate_float_input(prompt, min_value=min_val, max_value=max_val, tool=tool)
            auto_pred = calculate_prediction(tool, inputs[f'{tool}_score'])
            print(f"Info: Automatic prediction for {tool.upper()}: {auto_pred}")
            override = validate_choice_input(f"Manually enter prediction for {tool.upper()}? (yes/no)", ["yes", "no"])
            if override == "yes":
                inputs[f'{tool}_prediction'] = validate_choice_input(f"Enter {tool.upper()} prediction", choices)
            else:
                inputs[f'{tool}_prediction'] = auto_pred
            inputs[f'{tool}_prediction'] = validate_prediction(tool, inputs[f'{tool}_score'], inputs[f'{tool}_prediction'], choices)
    elif variant_type in ['splice_donor', 'splice_acceptor', 'intronic']:
        spliceai_scores = ['acceptor_gain', 'acceptor_loss', 'donor_gain', 'donor_loss']
        for score_type in spliceai_scores:
            inputs[f'spliceai_{score_type}_score'] = validate_float_input(
                f"Enter SpliceAI {score_type.replace('_', ' ').title()} score\n - Range: 0.0 to 1.0\n - Threshold: >0.5 (damaging)",
                min_value=0.0, max_value=1.0, tool='spliceai'
            )
            auto_pred = calculate_prediction('spliceai', inputs[f'spliceai_{score_type}_score'])
            print(f"Info: Automatic prediction for SpliceAI {score_type.replace('_', ' ').title()}: {auto_pred}")
            override = validate_choice_input(f"Manually enter prediction for SpliceAI {score_type.replace('_', ' ').title()}? (yes/no)", ["yes", "no"])
            if override == "yes":
                inputs[f'spliceai_{score_type}_prediction'] = validate_choice_input(
                    f"Enter SpliceAI {score_type.replace('_', ' ').title()} prediction",
                    ["damaging", "tolerated", "unknown"]
                )
            else:
                inputs[f'spliceai_{score_type}_prediction'] = auto_pred
            inputs[f'spliceai_{score_type}_prediction'] = validate_prediction('spliceai', inputs[f'spliceai_{score_type}_score'], inputs[f'spliceai_{score_type}_prediction'], ["damaging", "tolerated", "unknown"])
    elif variant_type in ['nonsense', 'frameshift']:
        inputs['loftool_score'] = validate_float_input(
            "Enter LoFtool score\n - Range: 0.0 to 1.0\n - Threshold: >0.9 (damaging)",
            min_value=0.0, max_value=1.0, tool='loftool'
        )
        auto_pred = calculate_prediction('loftool', inputs['loftool_score'])
        print(f"Info: Automatic prediction for LoFtool: {auto_pred}")
        override = validate_choice_input("Manually enter LoFtool prediction? (yes/no)", ["yes", "no"])
        if override == "yes":
            inputs['loftool_prediction'] = validate_choice_input("Enter LoFtool prediction", ["damaging", "tolerated", "unknown"])
        else:
            inputs['loftool_prediction'] = auto_pred
        inputs['loftool_prediction'] = validate_prediction('loftool', inputs['loftool_score'], inputs['loftool_prediction'], ["damaging", "tolerated", "unknown"])
    # Conservation tools (always ask)
    for tool, prompt, min_val, max_val, choices in [
        ('phylop_vert', 'phyloP Vertebrate score\n - Range: 0.0 to 1.0\n - Threshold: >=0.9 (conserved), <0.5 (non-conserved)', 0.0, 1.0, ['conserved', 'non-conserved', 'unknown']),
        ('phylop_mamm', 'phyloP Mammalian score\n - Range: 0.0 to 1.0\n - Threshold: >=0.9 (conserved), <0.5 (non-conserved)', 0.0, 1.0, ['conserved', 'non-conserved', 'unknown']),
        ('phylop_primate', 'phyloP Primate score\n - Range: 0.0 to 1.0\n - Threshold: >=0.9 (conserved), <0.5 (non-conserved)', 0.0, 1.0, ['conserved', 'non-conserved', 'unknown']),
        ('gerp++', 'GERP++ RS score\n - Range: -12.3 to 6.17\n - Threshold: >2.0 (conserved)', -12.3, 6.17, ['conserved', 'non-conserved', 'unknown'])
    ]:
        inputs[f'{tool}_score'] = validate_float_input(prompt, min_value=min_val, max_value=max_val, tool=tool)
        auto_pred = calculate_prediction(tool, inputs[f'{tool}_score'])
        print(f"Info: Automatic prediction for {tool.upper()}: {auto_pred}")
        override = validate_choice_input(f"Manually enter prediction for {tool.upper()}? (yes/no)", ["yes", "no"])
        if override == "yes":
            inputs[f'{tool}_prediction'] = validate_choice_input(f"Enter {tool.upper()} prediction", choices)
        else:
            inputs[f'{tool}_prediction'] = auto_pred
        inputs[f'{tool}_prediction'] = validate_prediction(tool, inputs[f'{tool}_score'], inputs[f'{tool}_prediction'], choices)
    # Return in silico scores and predictions in the expected structure
    insilico_scores = {k.replace('_score',''): v for k,v in inputs.items() if k.endswith('_score')}
    insilico_predictions = {k: v for k,v in inputs.items() if k.endswith('_prediction')}
    for k,v in insilico_scores.items():
        variant[f'{k}_score'] = v
    for k,v in insilico_predictions.items():
        variant[k] = v
    return {'insilico_scores': insilico_scores, 'insilico_predictions': insilico_predictions}

def collect_genetic_data():
    """Collect genetic data such as inheritance and zygosity with user-friendly prompts."""
    inputs = {}
    print(f"\n{'='*15} GENETIC DATA {'='*15}")
    print("Provide information about inheritance, zygosity, allelic state, and family data.")

    # Inheritance pattern (ask only once)
    inheritance = validate_choice_input(
        "Enter inheritance pattern\n - AD: Autosomal Dominant\n - AR: Autosomal Recessive\n - X-linked: X-linked inheritance\n - Mitochondrial: Mitochondrial\n - Unknown: Unknown",
        ["AD", "AR", "X-linked", "Mitochondrial", "unknown", "xlr", "xld", "mito"]
    ).lower()
    if inheritance in ["xlr", "xld"]:
        inheritance = "x-linked"
    if inheritance == "mito":
        inheritance = "mitochondrial"
    inputs['inheritance'] = inheritance

    # Zygosity (ask only once)
    zygosity = validate_zygosity_input(
        "Enter zygosity\n - heterozygous (het): One copy of variant\n - homozygous (hom): Two copies of variant\n - hemizygous: One copy (X/Y chromosome in males)\n - compound het: Compound heterozygous\n - mosaic: Mosaic\n - unknown: Unknown"
    ).lower()
    if zygosity == "het":
        zygosity = "heterozygous"
    if zygosity == "hom":
        zygosity = "homozygous"
    inputs['zygosity'] = zygosity

    # Allelic state (for AR, compound het, etc.)
    allelic_data = "unknown"
    if inheritance == "ar":
        if zygosity == "heterozygous" or zygosity == "compound het":
            print(f"\n{'='*15} ALLELIC DATA {'='*15}")
            print("For recessive diseases, indicate if the variant is in trans with another pathogenic variant.")
            allelic_data = validate_choice_input(
                "Is the variant in trans with a pathogenic variant? (yes/no/unknown)",
                ["yes", "no", "unknown"]
            )
        elif zygosity == "homozygous":
            allelic_data = "homozygous"
    elif inheritance == "x-linked":
        if zygosity == "hemizygous":
            allelic_data = "hemizygous"
    elif inheritance == "mitochondrial":
        # For most clinical ACMG workflows, detailed mitochondrial heteroplasmy is not required for variant classification.
        # If needed, uncomment and clarify:
        # print("Mitochondrial variants can be present in all (homoplasmic) or some (heteroplasmic) mitochondria.")
        # allelic_data = validate_choice_input(
        #     "Is the variant homoplasmic (present in all mitochondria) or heteroplasmic (present in some mitochondria)? (homoplasmic/heteroplasmic/unknown)",
        #     ["homoplasmic", "heteroplasmic", "unknown"]
        # )
        # For most users, this is not required, so default to unknown:
        allelic_data = "unknown"
    inputs['allelic_data'] = allelic_data


    # FAMILY DATA
    print(f"\n{'='*15} FAMILY DATA {'='*15}")
    print("Provide family history and segregation information.")
    inputs['family_history'] = validate_choice_input(
        "Is there a family history of the disease? (yes/no/unknown)",
        ["yes", "no", "unknown"]
    )
    inputs['segregation'] = validate_choice_input(
        "Does the variant segregate with disease in the family? (yes/no/unknown)",
        ["yes", "no", "unknown"]
    )

    # Only ask de novo, parental testing, consanguinity, siblings if not already present in inputs (avoid double ask)
    if 'de_novo' not in inputs:
        print(f"\n{'='*15} DE NOVO STATUS {'='*15}")
        print("Is the variant de novo (not inherited from parents)?")
        inputs['de_novo'] = validate_choice_input(
            "Is the variant de novo? (yes/no/unknown)",
            ["yes", "no", "unknown"]
        )
    if 'parental_testing' not in inputs:
        print(f"\n{'='*15} PARENTAL TESTING {'='*15}")
        print("Has parental testing been performed?")
        inputs['parental_testing'] = validate_choice_input(
            "Has parental testing been performed? (yes/no/unknown)",
            ["yes", "no", "unknown"]
        )
    if 'consanguinity' not in inputs:
        print(f"\n{'='*15} CONSANGUINITY {'='*15}")
        print("Is there consanguinity in the family?")
        inputs['consanguinity'] = validate_choice_input(
            "Is there consanguinity in the family? (yes/no/unknown)",
            ["yes", "no", "unknown"]
        )
    if 'affected_siblings' not in inputs:
        print(f"\n{'='*15} SIBLINGS {'='*15}")
        print("Are there affected siblings?")
        inputs['affected_siblings'] = validate_choice_input(
            "Are there affected siblings? (yes/no/unknown)",
            ["yes", "no", "unknown"]
        )

    # Additional notes
    print(f"\n{'='*15} ADDITIONAL NOTES {'='*15}")
    print("Any additional notes or comments? (Press Enter to skip)")
    notes = input("Notes: ").strip()
    inputs['genetic_notes'] = notes if notes else ""

    return inputs


# --- ACMG Evidence Assignment and Classification Functions ---


def run_acmg_assessment():
    """Main entry point for ACMG assessment and structured report saving."""
    # Collect all user inputs (basic info, population, in silico, then clinical/genetic/functional)
    basic_info = collect_variant_basic_info()
    pop_data = collect_population_data()
    variant = {**basic_info, **pop_data}

    # In silico scores (organized and user-friendly)
    print("\n" + "="*15 + " IN SILICO PREDICTIONS " + "="*15)
    print("You will now enter in silico prediction scores. Each tool is grouped by category. After each entry, an automatic prediction will be shown; you can confirm or override it.")

    # Define tool categories and explanations
    variant_level_tools = [
        ('revel', 'REVEL: Ensemble missense pathogenicity score (0-1)'),
        ('cadd', 'CADD: Combined Annotation Dependent Depletion (0-100)'),
        ('metarnn', 'MetaRNN: Missense pathogenicity (0-1)'),
        ('clinpred', 'ClinPred: Missense pathogenicity (0-1)'),
        ('bayesdel', 'BayesDel: Missense pathogenicity (-1 to 1)'),
        ('alphamissense', 'AlphaMissense: Missense pathogenicity (0-1)'),
        ('mutationtaster', 'MutationTaster: Disease-causing prediction (0-1)'),
        ('polyphen2', 'PolyPhen-2: Missense effect (0-1)'),
        ('sift', 'SIFT: Missense effect (0-1, lower is more damaging)'),
        ('fathmm_xf', 'FATHMM-XF: Functional effect (0-1)'),
        ('mutationassessor', 'MutationAssessor: Functional impact (0-5)'),
        ('provean', 'PROVEAN: Protein variation effect (-14 to 14)'),
        ('mutpred', 'MutPred: Missense pathogenicity (0-1)'),
        ('metolr', 'MetaLR: Missense pathogenicity (0-1)'),
        ('esm1b', 'ESM-1b: Protein language model (-100 to 100)'),
        ('lrt', 'LRT: Likelihood ratio test (0-1)'),
    ]
    conservation_tools = [
        ('phylop_vert', 'phyloP Vertebrate: Conservation ranked score (0-1). High = conserved.'),
        ('phylop_mamm', 'phyloP Mammal: Conservation ranked score (0-1). High = conserved.'),
        ('phylop_primate', 'phyloP Primate: Conservation ranked score (0-1). High = conserved.'),
        ('gerp++', 'GERP++: Evolutionary constraint score (-12.3 to 6.17). High = conserved. See https://mendel.stanford.edu/SidowLab/downloads/gerp/index.html'),
    ]
# --- Utility: Population Outlier Check ---

    insilico_scores = {}
    predictions = {}

    print("\n--- VARIANT-LEVEL PREDICTORS ---")
    for tool, explanation in variant_level_tools:
        if tool not in INSILICO_RANGES:
            continue
        print(f"{explanation}")
        score = validate_float_input(f"Enter {tool.upper()} score", min_value=INSILICO_RANGES[tool][0], max_value=INSILICO_RANGES[tool][1], tool=tool)
        insilico_scores[tool] = score
        variant[f"{tool}_score"] = score
        auto_pred = calculate_prediction(tool, score)
        print(f"  {tool.upper()} automatic prediction: {auto_pred}")
        # Prediction choices
        if tool in ['revel', 'cadd', 'metarnn', 'clinpred', 'alphamissense', 'mutationtaster', 'polyphen2', 'sift', 'fathmm_xf', 'mutationassessor', 'provean', 'mutpred', 'metolr', 'esm1b', 'lrt', 'bayesdel']:
            pred_choices = ['damaging', 'tolerated', 'unknown']
            if tool == 'bayesdel':
                pred_choices = ['deleterious', 'tolerated', 'unknown']
            elif tool == 'alphamissense':
                pred_choices = ['likely_pathogenic', 'likely_benign', 'ambiguous', 'unknown']
            elif tool == 'mutationtaster':
                pred_choices = ['damaging', 'polymorphism', 'unknown']
            elif tool == 'polyphen2':
                pred_choices = ['probably_damaging', 'possibly_damaging', 'benign', 'unknown']
            elif tool == 'mutationassessor':
                pred_choices = ['high', 'medium', 'low', 'neutral', 'unknown']
            elif tool == 'provean':
                pred_choices = ['damaging', 'neutral', 'unknown']
        else:
            pred_choices = ['unknown']
        user_pred = validate_choice_input(f"Prediction for {tool.upper()} (auto: {auto_pred})?\n(Press Enter to continue)", pred_choices)
        final_pred = validate_prediction(tool, score, user_pred, pred_choices)
        predictions[f"{tool}_prediction"] = final_pred
        variant[f"{tool}_prediction"] = final_pred

    print("\n--- CONSERVATION/PROTEIN-LEVEL PREDICTORS ---")
    for tool, explanation in conservation_tools:
        if tool not in INSILICO_RANGES:
            continue
        print(f"{explanation}")
        score = validate_float_input(f"Enter {tool.upper()} score", min_value=INSILICO_RANGES[tool][0], max_value=INSILICO_RANGES[tool][1], tool=tool)
        insilico_scores[tool] = score
        variant[f"{tool}_score"] = score
        auto_pred = calculate_prediction(tool, score)
        print(f"  {tool.upper()} automatic prediction: {auto_pred}")
        pred_choices = ['conserved', 'non-conserved', 'unknown']
        user_pred = validate_choice_input(f"Prediction for {tool.upper()} (auto: {auto_pred})?\n(Press Enter to continue)", pred_choices)
        final_pred = validate_prediction(tool, score, user_pred, pred_choices)
        predictions[f"{tool}_prediction"] = final_pred
        variant[f"{tool}_prediction"] = final_pred

    variant['insilico_scores'] = insilico_scores
    variant['insilico_predictions'] = predictions

    # --- Restore clinical/functional/genetic questions after in silico input ---
    # Inheritance, zygosity, allelic data
    # (functions are now defined above)
    genetic_data = collect_genetic_data()
    variant.update(genetic_data)

    # Functional and clinical data (segregation, de novo, paternity, phenotype match, functional test)
    functional_data = collect_functional_data(variant.get('vep_consequence', 'unknown'), variant.get('gene', ''))
    variant.update(functional_data)

    # Metascore calculation (for PP3/BP4 only, not final classification)
    metascore, criterion, meta_details = calculate_metascore(variant, variant['gene'])
    variant['metascore'] = round(metascore, 4)
    variant['metascore_criterion'] = criterion if criterion else ''




# --- Entry Point ---
if __name__ == "__main__":
    # (Açılış mesajı sadece process_variants içinde gösterilecek)
    try:
        # Move process_variants definition above this block if needed
        process_variants()
        print("Varsome: https://varsome.com/")
        print("Info: All in silico scores must be entered manually. No local dataset search is performed.")
    except Exception as e:
        print("\nAn error occurred during execution:")
        print(str(e))

    print("\nPress Enter to exit the program.")
    print("Or type 'R' and press Enter to start a new ACMG assessment.")
    choice = input("").strip().lower()
    if choice == "r":
        print("\nRestarting ACMG assessment...\n")
        import os, sys
        os.execv(sys.executable, [sys.executable] + sys.argv)
    else:
        print("Exiting. Goodbye!")
        import sys
        sys.exit(0)








def calculate_ensemble_score(inputs, variant_type):
    """Calculate ensemble score based on in silico predictions."""
    damaging_count = 0
    benign_count = 0
    
    if variant_type == 'missense':
        tools = ['revel', 'cadd', 'metarnn', 'clinpred', 'bayesdel', 'alphamissense', 'mutationtaster', 'polyphen2', 'sift', 'fathmm_xf', 'mutationassessor', 'provean', 'mutpred', 'metolr', 'esm1b', 'lrt', 'gerp']
        for tool in tools:
            prediction = inputs.get(f'{tool}_prediction', 'unknown')
            if prediction in ['damaging', 'likely_pathogenic', 'deleterious', 'conserved', 'high', 'medium']:
                damaging_count += 1
            elif prediction in ['tolerated', 'likely_benign', 'neutral', 'non-conserved', 'low', 'polymorphism']:
                benign_count += 1
    elif variant_type in ['splice_donor', 'splice_acceptor', 'intronic']:
        for score_type in ['acceptor_gain', 'acceptor_loss', 'donor_gain', 'donor_loss']:
            prediction = inputs.get(f'spliceai_{score_type}_prediction', 'unknown')
            if prediction == 'damaging':
                damaging_count += 1
            elif prediction == 'tolerated':
                benign_count += 1
    elif variant_type in ['nonsense', 'frameshift']:
        prediction = inputs.get('loftool_prediction', 'unknown')
        if prediction == 'damaging':
            damaging_count += 1
        elif prediction == 'tolerated':
            benign_count += 1
    
    for tool in ['phylop_vert', 'phylop_mamm', 'phylop_primate']:
        prediction = inputs.get(f'{tool}_prediction', 'unknown')
        if prediction == 'conserved':
            damaging_count += 1
        elif prediction == 'non-conserved':
            benign_count += 1
    
    return damaging_count, benign_count

def fisher_exact_test(cohort_allele_count, cohort_total, control_allele_count, control_total):
    """Perform Fisher's Exact Test for case-control data."""
    if cohort_total == 0 or control_total == 0:
        return 1.0, "Invalid: Total allele count is zero."
    
    table = [
        [cohort_allele_count, cohort_total - cohort_allele_count],
        [control_allele_count, control_total - control_allele_count]
    ]
    
    try:
        oddsratio, p_value = stats.fisher_exact(table)
        return p_value, f"Fisher's Exact Test p-value: {p_value:.4f}, Odds Ratio: {oddsratio:.4f}"
    except Exception as e:
        return 1.0, f"Fisher's Exact Test failed: {str(e)}"

def check_population_outlier(allele_frequency, subpopulation_frequencies):
    """Check if allele frequency is an outlier compared to subpopulations."""
    if not subpopulation_frequencies:
        return False, "No subpopulation data provided."
    
    mean_freq = np.mean(subpopulation_frequencies)
    std_freq = np.std(subpopulation_frequencies) if len(subpopulation_frequencies) > 1 else 0.01
    z_score = (allele_frequency - mean_freq) / std_freq if std_freq > 0 else 0
    is_outlier = abs(z_score) > 2
    return is_outlier, f"Z-score: {z_score:.4f}, Mean: {mean_freq:.6f}, Std: {std_freq:.6f}"

def assign_evidence(inputs, use_pp5_bp6=False, use_acmg_2023=False):
    """Assign ACMG evidence criteria based on input data."""
    evidence = []
    primary_variant = inputs['variants'][0]
    gene_thresholds = GENE_SPECIFIC_THRESHOLDS.get(primary_variant['gene'], GENE_SPECIFIC_THRESHOLDS['default'])
    
    ba1_threshold = gene_thresholds['BA1']
    bs1_threshold = gene_thresholds['BS1']
    if primary_variant['allele_frequency'] > ba1_threshold:
        evidence.append('BA1')
    elif primary_variant['allele_frequency'] > bs1_threshold and primary_variant['disease_prevalence'] > 0:
        evidence.append('BS1')
    elif primary_variant['allele_frequency'] > 0 and primary_variant['allele_frequency'] < 0.001:
        evidence.append('PM2')
    
    is_outlier, outlier_result = check_population_outlier(
        primary_variant['allele_frequency'], primary_variant.get('subpopulation_frequencies', [])
    )
    primary_variant['outlier_result'] = outlier_result
    if is_outlier and primary_variant['allele_frequency'] < 0.001:
        evidence.append('PM2')
    
    if primary_variant['variant_type'] in LOF_VARIANTS:
        evidence.append('PVS1')
    
    p_value, fisher_result = fisher_exact_test(
        primary_variant['cohort_allele_count'], primary_variant['cohort_total'],
        primary_variant['control_allele_count'], primary_variant['control_total']
    )
    primary_variant['fisher_p_value'] = p_value
    primary_variant['fisher_result'] = fisher_result
    # p_value bazen tuple olabiliyor, sadece float ise karşılaştır
    if isinstance(p_value, float) and p_value < 0.05:
        _log_info(f"Statistically significant difference: {fisher_result}")
        if primary_variant['cohort_allele_count'] / primary_variant['cohort_total'] > primary_variant['control_allele_count'] / primary_variant['control_total']:
            evidence.append('PS1')
        else:
            evidence.append('BS1')
    
    if 'single submitter' not in primary_variant['clinvar_note'].lower():
        if primary_variant['clinvar_status'] == 'pathogenic':
            evidence.append('PS1')
        elif primary_variant['clinvar_status'] == 'benign' and use_pp5_bp6:
            evidence.append('BP6')
        elif primary_variant['clinvar_status'] == 'likely pathogenic' and use_pp5_bp6:
            evidence.append('PP5')
    
    if primary_variant['metascore_criterion'] == 'PP3':
        evidence.append('PP3')
        primary_variant['insilico_prediction'] = 'damaging'
    elif primary_variant['metascore_criterion'] == 'BP4':
        evidence.append('BP4')
        primary_variant['insilico_prediction'] = 'benign'
    else:
        primary_variant['insilico_prediction'] = 'unknown'
    
    damaging_count, benign_count = calculate_ensemble_score(primary_variant, primary_variant['variant_type'])
    if damaging_count >= 3 and 'PP3' not in evidence:
        evidence.append('PP3')
    elif benign_count >= 3 and 'BP4' not in evidence:
        evidence.append('BP4')
    
    if primary_variant['variant_type'] == 'missense' and primary_variant['functional_test'] == 'loss':
        evidence.append('PS3')
    elif primary_variant['variant_type'] == 'missense' and primary_variant['functional_test'] == 'normal':
        evidence.append('BS3')
    elif primary_variant['variant_type'] in ['splice_donor', 'splice_acceptor', 'intronic'] and primary_variant['functional_test'] == 'yes':
        evidence.append('PS3')
    elif primary_variant['variant_type'] in ['splice_donor', 'splice_acceptor', 'intronic'] and primary_variant['functional_test'] == 'no':
        evidence.append('BS3')
    elif primary_variant['variant_type'] in ['nonsense', 'frameshift'] and primary_variant['functional_test'] == 'yes':
        evidence.append('PS3')
    elif primary_variant['variant_type'] in ['nonsense', 'frameshift'] and primary_variant['functional_test'] == 'no':
        evidence.append('BS3')
    
    if primary_variant['segregation'] == 'yes':
        evidence.append('PP1')
    elif primary_variant['segregation'] == 'no':
        evidence.append('BS4')
    
    if use_acmg_2023 and primary_variant['de_novo'] == 'yes' and primary_variant['paternity_confirmed'] == 'yes':
        evidence.append('PS2_Very_Strong')
    elif primary_variant['de_novo'] == 'yes' and primary_variant['inheritance'] in ['ad', 'x-linked']:
        evidence.append('PS2')
    
    if primary_variant['phenotype_match'] in ['strong', 'moderate']:
        evidence.append('PP4')
    
    if primary_variant['inheritance'] == 'ar' and primary_variant['zygosity'] == 'heterozygous' and len(inputs['variants']) > 1:
        second_variant = inputs['variants'][1]
        if second_variant['clinvar_status'] in ['pathogenic', 'likely pathogenic'] or detect_lof_variant(second_variant['cdna_change'], second_variant['vep_consequence']):
            evidence.append('PM3')
    
    return evidence

def format_evidence(evidence):
    """Format evidence list for display."""
    if not evidence:
        return "No evidence provided"
    return ", ".join(
        f"{e} ({EVIDENCE_WEIGHTS.get(e, 'Undefined')})" if e in EVIDENCE_WEIGHTS else e
        for e in evidence
    )

def classify_variant(evidence, use_acmg_2023=False):
    """Classify the variant based on ACMG criteria."""
    pvs_count = sum(1 for e in evidence if e == 'PVS1')
    ps_count = sum(1 for e in evidence if e.startswith('PS') and e != 'PS2_Very_Strong')
    ps_very_strong_count = sum(1 for e in evidence if e == 'PS2_Very_Strong')
    pm_count = sum(1 for e in evidence if e.startswith('PM'))
    pp_count = sum(1 for e in evidence if e.startswith('PP'))
    ba_count = sum(1 for e in evidence if e == 'BA1')
    bs_count = sum(1 for e in evidence if e.startswith('BS'))
    bp_count = sum(1 for e in evidence if e.startswith('BP'))

    if use_acmg_2023:
        if ba_count >= 1:
            return "Benign (BA1)"
        elif bs_count >= 2:
            return "Likely Benign (>=2 BS)"
        elif bs_count == 1 and bp_count >= 1:
            return "Likely Benign (BS+BP)"
        elif (pvs_count >= 1 and (ps_count >= 1 or ps_very_strong_count >= 1)) or \
             (ps_very_strong_count >= 1 and pm_count >= 1) or \
             (ps_count >= 2) or \
             (ps_count >= 1 and pm_count >= 3) or \
             (pm_count >= 4):
            return "Pathogenic"
        elif (ps_count == 1 and pm_count >= 1) or \
             (ps_count == 1 and pp_count >= 2) or \
             (pm_count >= 3) or \
             (pm_count >= 2 and pp_count >= 2):
            return "Likely Pathogenic"
        elif bs_count == 1 or bp_count >= 2:
            return "Likely Benign (BS or >=2 BP)"
        else:
            return "Uncertain Significance"
    else:
        if ba_count >= 1:
            return "Benign (BA1)"
        elif bs_count >= 2:
            return "Likely Benign (>=2 BS)"
        elif bs_count == 1 and bp_count >= 1:
            return "Likely Benign (BS+BP)"
        elif (pvs_count >= 1 and ps_count >= 1) or \
             (pvs_count >= 1 and pm_count >= 2) or \
             (ps_count >= 2) or \
             (ps_count >= 1 and pm_count >= 3) or \
             (pm_count >= 4):
            return "Pathogenic"
        elif (ps_count == 1 and pm_count >= 1) or \
             (ps_count == 1 and pp_count >= 2) or \
             (pm_count >= 3) or \
             (pm_count >= 2 and pp_count >= 2):
            return "Likely Pathogenic"
        elif bs_count == 1 or bp_count >= 2:
            return "Likely Benign (BS or >=2 BP)"
        else:
            return "Uncertain Significance"

def check_missing_data(inputs):
    """Check for missing data and provide suggestions."""
    missing_count = 0
    missing_details = []
    suggestions = {'varsome': [], 'orphanet_omim': [], 'clinvar': [], 'pubmed': [], 'genetic_reports': []}
    primary_variant = inputs['variants'][0]

    # Check for missing key data fields and provide suggestions
    if primary_variant.get('clinvar_status', "unknown") == "unknown":
        missing_count += 1
        missing_details.append('clinvar_status')
        suggestions['clinvar'].append("Check the variant's pathogenic/benign status in ClinVar.")

    if primary_variant.get('functional_test', "unknown") == "unknown":
        missing_count += 1
        missing_details.append('functional_test')
        suggestions['pubmed'].append("Search for functional test results in PubMed.")

    if primary_variant.get('segregation', "unknown") == "unknown":
        missing_count += 1
        missing_details.append('segregation')
        suggestions['genetic_reports'].append("Check family segregation data in genetic reports or publications.")

    if primary_variant.get('de_novo', "unknown") == "unknown":
        missing_count += 1
        missing_details.append('de_novo')
        suggestions['genetic_reports'].append("Check for de novo status in family studies or trio sequencing.")

    if primary_variant.get('paternity_confirmed', "unknown") == "unknown" and primary_variant.get('de_novo', '') == 'yes':
        missing_count += 1
        missing_details.append('paternity_confirmed')
        suggestions['genetic_reports'].append("Check for paternity confirmation in de novo cases.")

    if primary_variant.get('phenotype_match', "unknown") == "unknown":
        missing_count += 1
        missing_details.append('phenotype_match')
        suggestions['orphanet_omim'].append("Check phenotype match in Orphanet, OMIM, or clinical databases.")

    if primary_variant.get('allelic_data', "unknown") == "unknown" and primary_variant.get('inheritance', '') == 'ar':
        missing_count += 1
        missing_details.append('allelic_data')
        suggestions['genetic_reports'].append("Check for trans/cis status of variants in recessive inheritance.")

    if primary_variant.get('inheritance', "unknown") == "unknown":
        missing_count += 1
        missing_details.append('inheritance')
        suggestions['genetic_reports'].append("Check inheritance pattern in family or literature.")

    if primary_variant.get('allele_frequency', 0.0) == 0.0:
        missing_count += 1
        missing_details.append('allele_frequency')
        suggestions['varsome'].append("Check population frequency in gnomAD or Varsome.")

    if primary_variant.get('disease_prevalence', 0.0) == 0.0:
        missing_count += 1
        missing_details.append('disease_prevalence')
        suggestions['orphanet_omim'].append("Check disease prevalence in Orphanet or OMIM.")

    if primary_variant.get('variant_type', '') == 'missense':
        # Check for missing in silico scores
        for tool in ['revel', 'cadd', 'metarnn', 'clinpred', 'bayesdel', 'alphamissense', 'mutationtaster', 'polyphen2', 'sift']:
            if primary_variant.get(f'{tool}_score', 0.0) == 0.0:
                missing_count += 1
                missing_details.append(f'{tool}_score')
                suggestions['varsome'].append(f"Check {tool.upper()} score in Varsome or local datasets.")

    # Print summary of missing data and suggestions
    print(f"\n{'='*15} MISSING DATA SUMMARY {'='*15}")
    if missing_count == 0:
        print("All key data fields are present.")
    else:
        print(f"Missing {missing_count} key data fields: {', '.join(missing_details)}")
        print("Suggestions for data sources:")
        for source, items in suggestions.items():
            if items:
                print(f"- {source.title()}:\n  - " + "\n  - ".join(set(items)))
    return missing_count, missing_details, suggestions


def print_results(inputs, classification, evidence, suggestions):
    """Print classification results to the console with detailed explanations."""
    print(f"\n{'='*60}")
    print(f"VARIANT CLASSIFICATION REPORT")
    print(f"{'='*60}\n")
    print("SUMMARY")
    print(f"{'-'*40}")
    for i, variant in enumerate(inputs['variants'], 1):
        print(f"Variant {i}:")
        for k, v in variant.items():
            print(f"  {k}: {v}")
    print(f"\nClassification: {classification}")
    print(f"Evidence: {format_evidence(evidence)}")
    print(f"Varsome URL: {inputs.get('varsome_url', 'N/A')}")
    print(f"Combined MetaScore: {inputs['variants'][0].get('combined_pathogenicity_score', 0.0):.4f}")
    print(f"Fisher's Exact Test: {inputs['variants'][0].get('fisher_result', 'N/A')}")
    print(f"Subpopulation Outlier: {inputs['variants'][0].get('outlier_result', 'N/A')}")
    
    print(f"\n{'-'*40}")
    print("ACMG EVIDENCE DESCRIPTIONS")
    print(f"{'-'*40}")
    evidence_descriptions = {
        'BA1': f"Standalone benign: Population frequency >{GENE_SPECIFIC_THRESHOLDS.get(inputs['variants'][0]['gene'], GENE_SPECIFIC_THRESHOLDS['default'])['BA1']} ({inputs['variants'][0]['allele_frequency']}).",
        'BS1': f"Strong benign: Frequency ({inputs['variants'][0]['allele_frequency']}) above gene-specific threshold ({GENE_SPECIFIC_THRESHOLDS.get(inputs['variants'][0]['gene'], GENE_SPECIFIC_THRESHOLDS['default'])['BS1']}).",
        'PM2': f"Moderate pathogenic: Very rare variant (frequency {inputs['variants'][0]['allele_frequency']} <0.001, outlier: {inputs['variants'][0]['outlier_result']}).",
        'PVS1': f"Very strong pathogenic: LoF variant ({inputs['variants'][0]['variant_type']}).",
        'PS1': f"Strong pathogenic: ClinVar pathogenic or Fisher's Exact Test significant ({inputs['variants'][0]['clinvar_status']}, {inputs['variants'][0]['fisher_result']}).",
        'PS2': f"Strong pathogenic: De novo variant (Confirmed: {inputs['variants'][0]['paternity_confirmed']}).",
        'PS2_Very_Strong': f"Very strong pathogenic: De novo variant, parentage confirmed ({inputs['variants'][0]['paternity_confirmed']}).",
        'PS3': f"Strong pathogenic: Functional test shows loss of function ({inputs['variants'][0]['functional_test']}).",
        'BS3': f"Strong benign: Functional test shows no effect ({inputs['variants'][0]['functional_test']}).",
        'PP1': f"Supporting pathogenic: Variant segregates with disease in family ({inputs['variants'][0]['segregation']}).",
        'BS4': f"Strong benign: Variant present in healthy family members ({inputs['variants'][0]['segregation']}).",
        'PP4': f"Supporting pathogenic: Phenotype match ({inputs['variants'][0]['phenotype_match']}).",
        'PM3': f"Moderate pathogenic: Compound heterozygous with pathogenic variant ({inputs['variants'][0]['allelic_data']}).",
        'PP3': f"Supporting pathogenic: Combined MetaScore {inputs['variants'][0]['combined_pathogenicity_score']:.4f} (>0.354) or ≥3 damaging predictions.",
        'BP4': f"Supporting benign: Combined MetaScore {inputs['variants'][0]['combined_pathogenicity_score']:.4f} (<0.226) or ≥3 tolerated predictions.",
        'BP6': f"Supporting benign: ClinVar benign ({inputs['variants'][0]['clinvar_status']}).",
        'PP5': f"Supporting pathogenic: ClinVar likely pathogenic ({inputs['variants'][0]['clinvar_status']})."
    }
    for e in evidence:
        print(f"{e} ({EVIDENCE_WEIGHTS.get(e, 'Undefined')}): {evidence_descriptions.get(e, 'Undefined evidence')}")
    
    is_missing_data, missing_details, suggestions = check_missing_data(inputs)
    if is_missing_data:
        print(f"\n{'='*40}")
        print("WARNING: Critical data missing, classification may be unreliable")
        print(f"{'='*40}")
        print(f"Number of missing fields: {len(missing_details)}")
        print(f"Missing fields: {', '.join(missing_details)}")
        print(f"\nSUGGESTIONS")
        print(f"{'-'*40}")
        if suggestions['varsome']:
            print("Varsome Suggestions:")
            for suggestion in suggestions['varsome']:
                print(f"  - {suggestion} ({inputs['varsome_url']})")
        if suggestions['orphanet_omim']:
            print("Orphanet/OMIM Suggestions:")
            for suggestion in suggestions['orphanet_omim']:
                print(f"  - {suggestion} (https://www.orpha.net/, https://www.omim.org/)")
        if suggestions['clinvar']:
            print("ClinVar Suggestions:")
            for suggestion in suggestions['clinvar']:
                print(f"  - {suggestion} (https://www.ncbi.nlm.nih.gov/clinvar/)")
        if suggestions['pubmed']:
            print("PubMed Suggestions:")
            for suggestion in suggestions['pubmed']:
                print(f"  - {suggestion} (https://pubmed.ncbi.nlm.nih.gov/)")
        if suggestions['genetic_reports']:
            print("Genetic Report Suggestions:")
            for suggestion in suggestions['genetic_reports']:
                print(f"  - {suggestion}")
    
    print(f"\n{'-'*40}")
    print("ADDITIONAL NOTES")
    print(f"{'-'*40}")
    print(f"{inputs['notes']}\n")

def save_results(inputs, classification, evidence, suggestions):
    """Save classification results to a file with detailed explanations."""
    evidence_descriptions = {
        'BA1': f"Standalone benign: Population frequency >{GENE_SPECIFIC_THRESHOLDS.get(inputs['variants'][0]['gene'], GENE_SPECIFIC_THRESHOLDS['default'])['BA1']} ({inputs['variants'][0]['allele_frequency']}).",
        'BS1': f"Strong benign: Frequency ({inputs['variants'][0]['allele_frequency']}) above gene-specific threshold ({GENE_SPECIFIC_THRESHOLDS.get(inputs['variants'][0]['gene'], GENE_SPECIFIC_THRESHOLDS['default'])['BS1']}).",
        'PM2': f"Moderate pathogenic: Very rare variant (frequency {inputs['variants'][0]['allele_frequency']} <0.001, outlier: {inputs['variants'][0]['outlier_result']}).",
        'PVS1': f"Very strong pathogenic: LoF variant ({inputs['variants'][0]['variant_type']}).",
        'PS1': f"Strong pathogenic: ClinVar pathogenic or Fisher's Exact Test significant ({inputs['variants'][0]['clinvar_status']}, {inputs['variants'][0]['fisher_result']}).",
        'PS2': f"Strong pathogenic: De novo variant (Confirmed: {inputs['variants'][0]['paternity_confirmed']}).",
        'PS2_Very_Strong': f"Very strong pathogenic: De novo variant, parentage confirmed ({inputs['variants'][0]['paternity_confirmed']}).",
        'PS3': f"Strong pathogenic: Functional test shows loss of function ({inputs['variants'][0]['functional_test']}).",
        'BS3': f"Strong benign: Functional test shows no effect ({inputs['variants'][0]['functional_test']}).",
        'PP1': f"Supporting pathogenic: Variant segregates with disease in family ({inputs['variants'][0]['segregation']}).",
        'BS4': f"Strong benign: Variant present in healthy family members ({inputs['variants'][0]['segregation']}).",
        'PP4': f"Supporting pathogenic: Phenotype match ({inputs['variants'][0]['phenotype_match']}).",
        'PM3': f"Moderate pathogenic: Compound heterozygous with pathogenic variant ({inputs['variants'][0]['allelic_data']}).",
        'PP3': f"Supporting pathogenic: Combined MetaScore {inputs['variants'][0]['combined_pathogenicity_score']:.4f} (>0.354) or ≥3 damaging predictions.",
        'BP4': f"Supporting benign: Combined MetaScore {inputs['variants'][0]['combined_pathogenicity_score']:.4f} (<0.226) or ≥3 tolerated predictions.",
        'BP6': f"Supporting benign: ClinVar benign ({inputs['variants'][0]['clinvar_status']}).",
        'PP5': f"Supporting pathogenic: ClinVar likely pathogenic ({inputs['variants'][0]['clinvar_status']})."
    }
    
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    log_file = 'variant_classification_log.txt'
    # Ensure the log file is created if it does not exist
    import os
    if not os.path.exists(log_file):
        with open(log_file, 'w', encoding='utf-8') as f_init:
            f_init.write('')
    with open(log_file, 'a', encoding='utf-8') as f:
        f.write(f"{'='*60}\n")
        f.write(f"VARIANT CLASSIFICATION REPORT ({timestamp})\n")
        f.write(f"{'='*60}\n\n")
        f.write("SUMMARY\n")
        f.write(f"{'-'*40}\n")
        for i, variant in enumerate(inputs['variants'], 1):
            f.write(f"Variant {i}:\n")
            f.write(f"  Gene: {variant['gene']}\n")
            f.write(f"  Chromosomal Location: {variant['chromosome']}-{variant['position']}-{variant['ref_allele']}-{variant['alt_allele']}\n")
            f.write(f"  cDNA Change: {variant['cdna_change']}\n")
            f.write(f"  Protein Change: {variant['protein_change']}\n")
        f.write(f"\nClassification: {classification}\n")
        f.write(f"Evidence: {format_evidence(evidence)}\n")
        f.write(f"Varsome URL: {inputs['varsome_url']}\n")
        f.write(f"Combined MetaScore: {inputs['variants'][0]['combined_pathogenicity_score']:.4f}\n")
        f.write(f"Fisher's Exact Test: {inputs['variants'][0]['fisher_result']}\n")
        f.write(f"Subpopulation Outlier: {inputs['variants'][0]['outlier_result']}\n")
        
        f.write(f"\n{'-'*40}\n")
        f.write("ACMG EVIDENCE DESCRIPTIONS\n")
        f.write(f"{'-'*40}\n")
        for e in evidence:
            f.write(f"{e} ({EVIDENCE_WEIGHTS.get(e, 'Undefined')}): {evidence_descriptions.get(e, 'Undefined evidence')}\n")
        
        f.write(f"\n{'-'*40}\n")
        f.write("INPUTS\n")
        f.write(f"{'-'*40}\n")
        for i, variant in enumerate(inputs['variants'], 1):
            f.write(f"\nVariant {i}:\n")
            f.write(f"  Gene: {variant['gene']}\n")
            f.write(f"  Chromosome: {variant['chromosome']}\n")
            f.write(f"  Position: {variant['position']}\n")
            f.write(f"  Reference Allele: {variant['ref_allele']}\n")
            f.write(f"  Alternate Allele: {variant['alt_allele']}\n")
            f.write(f"  cDNA Change: {variant['cdna_change']}\n")
            f.write(f"  Protein Change: {variant['protein_change']}\n")
            f.write(f"  VEP Consequence: {variant['vep_consequence']}\n")
            f.write(f"  Population Frequency: {variant['allele_frequency']}\n")
            f.write(f"  Subpopulation Frequencies: {variant.get('subpopulation_frequencies', [])})\n")
            f.write(f"  Disease Prevalence: {variant['disease_prevalence']}\n")
            f.write(f"  Case Allele Count: {variant['cohort_allele_count']}\n")
            f.write(f"  Case Total Alleles: {variant['cohort_total']}\n")
            f.write(f"  Control Allele Count: {variant['control_allele_count']}\n")
            f.write(f"  Control Total Alleles: {variant['control_total']}\n")
            f.write(f"  Inheritance Pattern: {variant['inheritance']}\n")
            f.write(f"  Variant Type: {variant['variant_type']}\n")
            f.write(f"  ClinVar Status: {variant['clinvar_status']} ({variant['clinvar_note']})\n")
            f.write(f"  In silico Prediction: {variant['insilico_prediction']}\n")
            f.write(f"  In silico Tools Used: {variant['insilico_tools']}\n")
            for tool in GENE_SPECIFIC_WEIGHTS.get(variant['gene'], GENE_SPECIFIC_WEIGHTS['default']):
                if variant.get(f'{tool}_score', 0.0) != 0.0:
                    f.write(f"  {tool.upper()}: {variant.get(f'{tool}_prediction', 'unknown')}/{variant.get(f'{tool}_score', 0.0)}\n")
            if variant['variant_type'] in ['splice_donor', 'splice_acceptor', 'intronic']:
                for score_type in ['acceptor_gain', 'acceptor_loss', 'donor_gain', 'donor_loss']:
                    if variant.get(f'spliceai_{score_type}_score', 0.0) != 0.0:
                        f.write(f"  SpliceAI {score_type.replace('_', ' ').title()}: {variant.get(f'spliceai_{score_type}_prediction', 'unknown')}/{variant.get(f'spliceai_{score_type}_score', 0.0)}\n")
            if variant['variant_type'] in ['nonsense', 'frameshift']:
                if variant.get('loftool_score', 0.0) != 0.0:
                    f.write(f"  LoFtool: {variant.get('loftool_prediction', 'unknown')}/{variant.get('loftool_score', 0.0)}\n")
            for tool in ['phylop_vert', 'phylop_mamm', 'phylop_primate']:
                if variant.get(f'{tool}_score', 0.0) != 0.0:
                    f.write(f"  {tool.upper()}: {variant.get(f'{tool}_prediction', 'unknown')}/{variant.get(f'{tool}_score', 0.0)}\n")
            f.write(f"  Zygosity: {variant['zygosity']}\n")
            f.write(f"  Allelic Data: {variant['allelic_data']}\n")
            f.write(f"  Functional Test: {variant['functional_test']}\n")
            f.write(f"  Segregation: {variant['segregation']}\n")
            f.write(f"  De Novo: {variant['de_novo']}\n")
            f.write(f"  Paternity Confirmed: {variant['paternity_confirmed']}\n")
            f.write(f"  Phenotype Match: {variant['phenotype_match']}\n")
        
        f.write(f"\n{'-'*40}\n")
        f.write("ADDITIONAL NOTES\n")
        f.write(f"{'-'*40}\n")
        f.write(f"{inputs['notes']}\n")
        
        is_missing_data, missing_details, suggestions = check_missing_data(inputs)
        if is_missing_data:
            f.write(f"\n{'='*40}\n")
            f.write("WARNING: Critical data missing, classification may be unreliable\n")
            f.write(f"{'='*40}\n")
            f.write(f"Number of missing fields: {len(missing_details)}\n")
            f.write(f"Missing fields: {', '.join(missing_details)}\n")
            f.write(f"\nSUGGESTIONS\n")
            f.write(f"{'-'*40}\n")
            if suggestions['varsome']:
                f.write("Varsome Suggestions:\n")
                for suggestion in suggestions['varsome']:
                    f.write(f"  - {suggestion} ({inputs['varsome_url']})\n")
            if suggestions['orphanet_omim']:
                f.write("Orphanet/OMIM Suggestions:\n")
                for suggestion in suggestions['orphanet_omim']:
                    f.write(f"  - {suggestion} (https://www.orpha.net/, https://www.omim.org/)\n")
            if suggestions['clinvar']:
                f.write("ClinVar Suggestions:\n")
                for suggestion in suggestions['clinvar']:
                    f.write(f"  - {suggestion} (https://www.ncbi.nlm.nih.gov/clinvar/)\n")
            if suggestions['pubmed']:
                f.write("PubMed Suggestions:\n")
                for suggestion in suggestions['pubmed']:
                    f.write(f"  - {suggestion} (https://pubmed.ncbi.nlm.nih.gov/)\n")
            if suggestions['genetic_reports']:
                f.write("Genetic Report Suggestions:\n")
                for suggestion in suggestions['genetic_reports']:
                    f.write(f"  - {suggestion}\n")
    
    print(f"\nInfo: Results saved to '{log_file}'.\n")

def clear_cache():
    """Clear ClinVar and score caches."""
    clinvar_cache.clear()
    score_cache.clear()
    print("Info: Caches cleared.")


# --- Entry Point ---
def process_variants():
    """Process variants and classify them."""
    print("\n" + "="*60)
    print("      ACMG 2015/2023 VARIANT CLASSIFICATION ASSISTANT")
    print("                 Created by: Can Sevilmiş")
    print("="*60 + "\n")
    while True:
        inputs = collect_inputs()
        evidence = assign_evidence(inputs, use_pp5_bp6=True, use_acmg_2023=True)
        classification = classify_variant(evidence, use_acmg_2023=True)
        suggestions = check_missing_data(inputs)[2]
        save_results(inputs, classification, evidence, suggestions)
        print_results(inputs, classification, evidence, suggestions)
        another_variant = validate_choice_input("Classify another variant? (yes/no)", ["yes", "no"])
        if another_variant == "no":
            break
        clear_cache()

def print_terminal_warnings_and_suggestions(report_data):
    # Print final ACMG classification at the end, colored
    final_class = str(report_data.get('classification', '[MISSING]'))
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    RESET = '\033[0m'
    if 'benign' in final_class.lower():
        color_code = GREEN
    elif 'pathogenic' in final_class.lower():
        color_code = RED
    else:
        color_code = YELLOW
    print(f"\n{color_code}FINAL ACMG CLASSIFICATION: {final_class}{RESET}")
    """Print missing/optional field warnings and suggestions in color, grouped at the end."""
    RED = '\033[91m'
    YELLOW = '\033[93m'
    CYAN = '\033[96m'
    RESET = '\033[0m'

    # Define required and optional fields for each section (with suggestions)
    required_fields = [
        ('gene', 'Gene', 'Obtain from sequencing report or VCF/annotation.'),
        ('chromosome', 'Chromosome', 'Obtain from VCF or annotation tools.'),
        ('position', 'Position', 'Obtain from VCF or annotation tools.'),
        ('ref_allele', 'Reference Allele', 'Obtain from VCF or annotation tools.'),
        ('alt_allele', 'Alternate Allele', 'Obtain from VCF or annotation tools.'),
        ('cdna_change', 'cDNA Change', 'Obtain from sequencing report or annotation.'),
        ('protein_change', 'Protein Change', 'Obtain from sequencing report or annotation.'),
        ('vep_consequence', 'VEP Consequence', 'Use Ensembl VEP or similar annotation tools.'),
        ('variant_type', 'Final Variant Type', 'User-selected or auto-detected.'),
        ('allele_frequency', 'Allele Frequency', 'Check gnomAD: https://gnomad.broadinstitute.org/.'),
        ('inheritance', 'Inheritance', 'Check inheritance pattern in family or literature.'),
        ('zygosity', 'Zygosity', 'Check zygosity from sequencing data.'),
        ('functional_test', 'Functional Test', 'Search for functional test results in PubMed.'),
        ('phenotype_match', 'Phenotype Match', 'Check phenotype match in Orphanet, OMIM, or clinical databases.'),
    ]
    optional_fields = [
        ('date', 'Date/Time', 'Set automatically by the tool.'),
        ('varsome_url', 'Varsome URL', 'Generated by the tool; check https://varsome.com/.'),
        ('clinvar_status', 'ClinVar Status', 'Check https://www.ncbi.nlm.nih.gov/clinvar/.'),
        ('disease_prevalence', 'Disease Prevalence', 'Check Orphanet: https://www.orpha.net/.'),
        ('subpopulation_frequencies', 'Subpopulation Frequencies', 'Check gnomAD subpopulations.'),
        ('allelic_data', 'Allelic State', 'Check for trans/cis status of variants in recessive inheritance.'),
        ('family_history', 'Family History', 'Check family history in clinical records.'),
        ('segregation', 'Segregation', 'Check family segregation data in genetic reports or publications.'),
        ('de_novo', 'De Novo Status', 'Check for de novo status in family studies or trio sequencing.'),
        ('parental_testing', 'Parental Testing', 'Check if parental testing was performed.'),
        ('consanguinity', 'Consanguinity', 'Check for consanguinity in family history.'),
        ('affected_siblings', 'Affected Siblings', 'Check for affected siblings in family.'),
        ('genetic_notes', 'Additional Genetic Notes', 'Add any relevant notes.'),
        ('paternity_confirmed', 'Paternity Confirmed', 'Check for paternity confirmation in de novo cases.'),
    ]

    # Collect missing fields and suggestions
    missing_required = []
    missing_optional = []
    suggestions = []
    for field, label, suggestion in required_fields:
        val = report_data.get(field, 'N/A')
        if val is None or val == '' or val == 'N/A' or val == 'Not specified' or val == 0:
            missing_required.append((label, suggestion))
            suggestions.append(f"{label}: {suggestion}")
    for field, label, suggestion in optional_fields:
        val = report_data.get(field, 'N/A')
        if val is None or val == '' or val == 'N/A' or val == 'Not specified' or val == 0:
            missing_optional.append((label, suggestion))
            suggestions.append(f"{label}: {suggestion}")

    # Print warnings for missing required fields
    if missing_required:
        print(f"{RED}WARNING: Missing required fields:{RESET}")
        for label, suggestion in missing_required:
            print(f"  {RED}[MISSING] {label}: Not provided. Suggestion: {suggestion}{RESET}")

    # Print warnings for missing optional fields
    if missing_optional:
        print(f"{YELLOW}Note: Missing optional fields:{RESET}")
        for label, suggestion in missing_optional:
            print(f"  {YELLOW}[MISSING] {label}: Not provided. Suggestion: {suggestion}{RESET}")

    # Print all suggestions at the end
    if suggestions:
        print(f"\n{CYAN}SUGGESTIONS:{RESET}")
        for s in sorted(set(suggestions)):
            print(f"  {CYAN}{s}{RESET}")
# --- Functional and Clinical Data Collection ---
def collect_functional_data(variant_type, gene):
    """Collect functional and clinical data with user-friendly prompts."""
    inputs = {}
    print(f"\n{'='*15} FUNCTIONAL AND CLINICAL DATA {'='*15}")
    print("Provide functional test results and clinical observations for the variant.")

    # Functional test section
    if variant_type == 'missense':
        inputs['functional_test'] = validate_choice_input(
            "Enter protein activity test result\n - Loss: Reduced or no protein function\n - Normal: Normal protein function",
            ["Loss", "Normal", "unknown"]
        )
    elif variant_type in ['splice_donor', 'splice_acceptor', 'intronic']:
        inputs['functional_test'] = validate_choice_input(
            "Enter splicing alteration result\n - Yes: Splicing is altered\n - No: Splicing is unaffected",
            ["Yes", "No", "unknown"]
        )
    elif variant_type in ['nonsense', 'frameshift']:
        inputs['functional_test'] = validate_choice_input(
            "Enter NMD (Nonsense-Mediated Decay) status\n - Yes: NMD is triggered\n - No: NMD is not triggered",
            ["Yes", "No", "unknown"]
        )
    else:
        inputs['functional_test'] = validate_choice_input(
            "Enter functional test result\n - Loss: Reduced or no function\n - Normal: Normal function",
            ["Loss", "Normal", "unknown"]
        )



    # Phenotype match section
    print(f"\n{'='*15} PHENOTYPE MATCH {'='*15}")
    print("How well does the patient's phenotype match the gene?")
    inputs['phenotype_match'] = validate_choice_input(
        f"Enter phenotype match strength for {gene} gene\n - Strong: Highly specific to gene\n - Moderate: Partially specific\n - Weak: Low specificity",
        ["Strong", "Moderate", "Weak", "unknown"]
    )
    # Remove de novo and parental testing prompts from here (they are already in collect_genetic_data)
    return inputs

"""
MIT License
Copyright (c) 2025 Can Sevilmiş
See LICENSE file for full license text.
"""

import time
import argparse
import requests
import scipy.stats as stats
import numpy as np
import re



from typing import Dict, Optional, Tuple


# --- Imports (ensure all required modules are imported) ---
import os
import sys
import time
import traceback
import inspect
import re
import requests
from typing import Dict, Optional, Tuple

# Remove all logging usage (no logging to file, only print or silent error handling)
def _log_info(msg):
    pass
def _log_warning(msg):
    pass
def _log_error(msg):
    pass

# --- Report Writing ---
def write_structured_report(report_data: dict, filename: str = "variant_classification_report.txt"):
    """Write a structured, human-readable report for the ACMG assessment."""
    import os
    out_path = os.path.join(os.getcwd(), filename)
    try:
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("\n" + "="*70 + "\n")
            f.write(" ACMG VARIANT CLASSIFICATION REPORT\n")
            f.write("="*70 + "\n\n")

            # [VARIANT]
            f.write("[VARIANT]\n")
            f.write(f"  Gene:         {report_data.get('gene', '-'):<12}  cDNA:   {report_data.get('cdna_change', '-'):<18}\n")
            f.write(f"  Chromosome:   {report_data.get('chromosome', '-'):<12}  Protein: {report_data.get('protein_change', '-'):<18}\n")
            f.write(f"  Position:     {report_data.get('position', '-'):<12}  VEP:    {report_data.get('vep_consequence', '-'):<18}\n")
            f.write(f"  Ref/Alt:      {report_data.get('ref_allele', '-')}/{report_data.get('alt_allele', '-'): <10}  Type:   {report_data.get('variant_type', '-'):<18}\n")
            f.write(f"  Date:         {report_data.get('date', '-')}\n\n")

            # [POPULATION]
            f.write("[POPULATION]\n")
            f.write(f"  Allele Freq:        {report_data.get('allele_frequency', '-')}\n")
            f.write(f"  Disease Prevalence: {report_data.get('disease_prevalence', '-')}\n")
            subpops = report_data.get('subpopulation_frequencies', [])
            if isinstance(subpops, list):
                subpop_str = ', '.join(str(x) for x in subpops) if subpops else '-'
            else:
                subpop_str = str(subpops) if subpops else '-'
            f.write(f"  Subpop Freqs:       {subpop_str}\n\n")

            # [IN SILICO PREDICTIONS]
            f.write("[IN SILICO PREDICTIONS]\n")
            f.write("  Tool            Score      Prediction\n")
            f.write("  -----------------------------------------\n")
            insilico_scores = report_data.get('insilico_scores', {})
            insilico_preds = report_data.get('insilico_predictions', {})
            for tool in [
                'revel','cadd','metarnn','clinpred','bayesdel','alphamissense','mutationtaster','polyphen2','sift','fathmm_xf','mutationassessor','provean','mutpred','metolr','esm1b','lrt','phylop_vert','phylop_mamm','phylop_primate','gerp++']:
                score = insilico_scores.get(tool, None)
                pred = insilico_preds.get(f"{tool}_prediction", None)
                if score is not None and (score != 0 or pred not in [None, '-', '', 'unknown']):
                    score_str = str(score) if score not in [None, '', 'N/A'] else '-'
                    pred_str = pred if pred not in [None, '', 'N/A'] else '-'
                    f.write(f"  {tool.upper():<15} {score_str:<10} {pred_str}\n")
            f.write(f"  Combined MetaScore: {report_data.get('metascore', '-')}\n\n")

            # [GENETIC & FAMILY]
            f.write("[GENETIC & FAMILY]\n")
            f.write(f"  Inheritance:      {report_data.get('inheritance', '-'):<12}  Zygosity:      {report_data.get('zygosity', '-'):<12}\n")
            f.write(f"  Allelic State:    {report_data.get('allelic_data', '-'):<12}  Family Hist:   {report_data.get('family_history', '-'):<12}\n")
            f.write(f"  Segregation:      {report_data.get('segregation', '-'):<12}  De Novo:       {report_data.get('de_novo', '-'):<12}\n")
            f.write(f"  Parental Testing: {report_data.get('parental_testing', '-'):<12}  Consanguinity: {report_data.get('consanguinity', '-'):<12}\n")
            f.write(f"  Siblings:         {report_data.get('affected_siblings', '-')}\n")
            notes = report_data.get('genetic_notes', '-')
            f.write(f"  Notes:            {notes if notes else '-'}\n\n")

            # [FUNCTIONAL/CLINICAL]
            f.write("[FUNCTIONAL/CLINICAL]\n")
            f.write(f"  Functional Test:     {report_data.get('functional_test', '-')}\n")
            f.write(f"  Phenotype Match:     {report_data.get('phenotype_match', '-')}\n")
            f.write(f"  Paternity Confirmed: {report_data.get('paternity_confirmed', '-')}\n\n")

            # [ACMG CLASSIFICATION]
            f.write("[ACMG CLASSIFICATION]\n")
            crit = report_data.get('criteria', [])
            if isinstance(crit, list):
                crit_str = ', '.join(crit) if crit else '-'
            else:
                crit_str = str(crit) if crit else '-'
            f.write(f"  Criteria:   {crit_str}\n")
            f.write(f"  Final:      {report_data.get('classification', '-')}\n\n")

            # [REFERENCES]
            f.write("[REFERENCES]\n")
            f.write(f"  Varsome:    {report_data.get('varsome_url', '-')}\n")
            f.write(f"  ClinVar:    {report_data.get('clinvar_status', '-')}\n")
            notes = report_data.get('notes', '-')
            f.write(f"  Notes:      {notes if notes else '-'}\n\n")

            f.write("="*70 + "\n")
            f.flush()
    except Exception as e:
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("[ERROR] Failed to write ACMG report.\n")
            f.write(str(e) + "\n")
        print("[ERROR] Exception in write_structured_report:", str(e))

# --- Caching Mechanisms ---
clinvar_cache = {}
score_cache = {}

# --- Constants ---
EVIDENCE_WEIGHTS = {
    'PVS1': 'Very Strong Pathogenic', 'PS1': 'Strong Pathogenic', 'PS2': 'Strong Pathogenic',
    'PS2_Very_Strong': 'Very Strong Pathogenic', 'PS3': 'Strong Pathogenic', 'PM2': 'Moderate Pathogenic',
    'PM3': 'Moderate Pathogenic', 'PP1': 'Supporting Pathogenic', 'PP3': 'Supporting Pathogenic',
    'PP4': 'Supporting Pathogenic', 'PP5': 'Supporting Pathogenic (Deprecated)', 'BA1': 'Standalone Benign',
    'BS1': 'Strong Benign', 'BS3': 'Strong Benign', 'BS4': 'Strong Benign',
    'BP4': 'Supporting Benign', 'BP6': 'Supporting Benign (Deprecated)'
}
GENE_SPECIFIC_THRESHOLDS = {
    'BRCA1': {'BA1': 0.05, 'BS1': 0.01}, 'TTN': {'BA1': 0.1, 'BS1': 0.05},
    'CAMTA1': {'BA1': 0.05, 'BS1': 0.01}, 'default': {'BA1': 0.05, 'BS1': 0.01}
}
GENE_SPECIFIC_WEIGHTS = {
    'default': {
        'revel': 0.25, 'cadd': 0.20, 'metarnn': 0.15, 'clinpred': 0.15, 'bayesdel': 0.10,
        'alphamissense': 0.10, 'mutationtaster': 0.05, 'polyphen2': 0.05, 'sift': 0.05
    },
    'CAMTA1': {
        'revel': 0.30, 'cadd': 0.25, 'metarnn': 0.20, 'clinpred': 0.15, 'bayesdel': 0.10,
        'alphamissense': 0.10, 'mutationtaster': 0.00, 'polyphen2': 0.00, 'sift': 0.00
    },
    'BRCA1': {
        'revel': 0.35, 'cadd': 0.20, 'metarnn': 0.20, 'clinpred': 0.15, 'bayesdel': 0.10,
        'alphamissense': 0.10, 'mutationtaster': 0.00, 'polyphen2': 0.00, 'sift': 0.00
    }
}
INSILICO_RANGES = {
    'revel': (0.0, 1.0, False), 'cadd': (0.0, 100.0, False), 'metarnn': (0.0, 1.0, False),
    'clinpred': (0.0, 1.0, False), 'bayesdel': (-1.0, 1.0, False), 'alphamissense': (0.0, 1.0, False),
    'mutationtaster': (0.0, 1.0, False), 'polyphen2': (0.0, 1.0, False), 'sift': (0.0, 1.0, True),
    'fathmm_xf': (0.0, 1.0, False), 'mutationassessor': (0.0, 5.0, False), 'provean': (-14.0, 14.0, True),
    'mutpred': (0.0, 1.0, False), 'metolr': (0.0, 1.0, False), 'esm1b': (-100.0, 100.0, False),
    'lrt': (0.0, 1.0, True), 'gerp++': (-12.3, 6.17, False),
    # phyloP conservation scores now use ranked score (0-1)
    'phylop_vert': (0.0, 1.0, False), 'phylop_mamm': (0.0, 1.0, False), 'phylop_primate': (0.0, 1.0, False)
}
LOF_VARIANTS = {'nonsense', 'frameshift', 'splice_donor', 'splice_acceptor'}
first_time = True


def calculate_metascore(inputs: Dict, gene: str) -> Tuple[float, Optional[str], dict]:
    """Calculate a VAMPP-score-like metascore and determine PP3/BP4 criterion. Returns details for reporting."""
    cache_key = tuple((k, inputs.get(k)) for k in sorted(GENE_SPECIFIC_WEIGHTS.get(gene, GENE_SPECIFIC_WEIGHTS['default']).keys()) if inputs.get(f'{k}_score') is not None)
    if cache_key in score_cache:
        metascore, criterion, details = score_cache[cache_key]
        print(f"Info: Combined MetaScore retrieved from cache: {metascore:.4f} ({criterion if criterion else 'No criterion applied'})")
        return metascore, criterion, details

    weights = GENE_SPECIFIC_WEIGHTS.get(gene, GENE_SPECIFIC_WEIGHTS['default'])
    scores = []
    used_weights = []
    used_tools = []
    normalized_scores = {}

    for tool in weights:
        score = inputs.get(f'{tool}_score')
        normalized = normalize_score(tool, score)
        if normalized is not None:
            scores.append(normalized)
            used_weights.append(weights[tool])
            used_tools.append(tool)
            normalized_scores[tool] = round(normalized, 4)

    if not scores:
        return 0.0, None, {"used_tools": [], "normalized_scores": {}, "weights": {}, "weighted_sum": 0.0, "total_weight": 0.0}

    total_weight = sum(used_weights)
    if total_weight == 0:
        return 0.0, None, {"used_tools": used_tools, "normalized_scores": normalized_scores, "weights": used_weights, "weighted_sum": 0.0, "total_weight": 0.0}

    weighted_sum = sum(s * w for s, w in zip(scores, used_weights))
    metascore = weighted_sum / total_weight

    criterion = 'PP3' if metascore > 0.354 else 'BP4' if metascore < 0.226 else None
    details = {
        "used_tools": used_tools,
        "normalized_scores": normalized_scores,
        "weights": {t: w for t, w in zip(used_tools, used_weights)},
        "weighted_sum": round(weighted_sum, 4),
        "total_weight": round(total_weight, 4),
        "final_metascore": round(metascore, 4),
        "criterion": criterion
    }
    score_cache[cache_key] = (metascore, criterion, details)

    return metascore, criterion, details

def get_chromosome_from_ensembl(gene_name):
    """Retrieve chromosome information for a gene using Ensembl API."""
    if not gene_name or gene_name == "Not specified":
        return "Not specified"
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/homo_sapiens/{gene_name}?expand=0"
    try:
        r = requests.get(server + ext, headers={"Content-Type": "application/json"}, timeout=10)
        r.raise_for_status()
        decoded = r.json()
        return decoded.get("seq_region_name", "Not specified")
    except requests.exceptions.RequestException as e:
        _log_warning(f"Ensembl API request failed: {e}. Chromosome information not retrieved.")
        return "Not specified"

def get_clinvar_status(chromosome, position, ref_allele, alt_allele):
    """Fetch ClinVar status for a variant."""
    if chromosome == "Not specified" or position == 0 or ref_allele == "Not specified" or alt_allele == "Not specified":
        return "unknown", "ClinVar data could not be retrieved: Missing information."
    cache_key = f"{chromosome}-{position}-{ref_allele}-{alt_allele}"
    if cache_key in clinvar_cache:
        _log_info("ClinVar status retrieved from cache.")
        return clinvar_cache[cache_key]
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    query = f"{chromosome}:{position}[chrpos38] AND {ref_allele}>{alt_allele}"
    params = {'db': 'clinvar', 'term': query, 'retmode': 'json', 'retmax': 10}
    try:
        r = requests.get(base_url, params=params, timeout=10)
        r.raise_for_status()
        result = r.json()
        ids = result.get('esearchresult', {}).get('idlist', [])
        if not ids:
            result_tuple = ("unknown", "No variant found in ClinVar.")
        else:
            clinvar_status = "unknown"
            clinvar_note = ""
            # Try all returned ids, pick the first with a valid clinical_significance
            for cid in ids:
                detail_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={cid}&retmode=json"
                r_detail = requests.get(detail_url, timeout=10)
                r_detail.raise_for_status()
                detail = r_detail.json()
                docsum = detail.get('result', {}).get(cid, {})
                cs = docsum.get('clinical_significance')
                if cs:
                    # Sometimes cs is a dict, sometimes a string
                    if isinstance(cs, dict):
                        clinvar_status = cs.get('description', 'unknown').lower()
                        submitters = cs.get('submitter', [])
                        submitter_note = " (Single submitter)" if isinstance(submitters, list) and len(submitters) <= 1 else " (Consensus)"
                    else:
                        clinvar_status = str(cs).lower()
                        submitter_note = ""
                    clinvar_note = f"ClinVar status: {clinvar_status}{submitter_note}"
                    break
            if clinvar_status == "unknown":
                clinvar_note = "No clinical significance found in ClinVar."
            result_tuple = (clinvar_status, clinvar_note)
        clinvar_cache[cache_key] = result_tuple
        return result_tuple
    except requests.exceptions.RequestException as e:
        _log_error(f"ClinVar API request failed: {e}")
        return "unknown", f"ClinVar API request failed: {e}"

def validate_float_input(prompt, allow_zero=True, min_value=-100.0, max_value=100.0, reference_url=None, tool=None):
    """Validate float input within a specified range with enhanced user guidance."""
    global first_time
    # Araç bazlı örnek değerler
    tool_examples: Dict[str, str] = {
        'revel': '0.75', 'cadd': '25.0', 'metarnn': '0.6', 'clinpred': '0.6', 'bayesdel': '0.2',
        'alphamissense': '0.6', 'mutationtaster': '0.6', 'polyphen2': '0.6', 'sift': '0.01',
        'fathmm_xf': '0.6', 'mutationassessor': '2.0', 'provean': '-3.0', 'mutpred': '0.6',
        'metolr': '0.6', 'esm1b': '1.0', 'lrt': '0.01', 'gerp++': '3.0',
        'phylop_vert': '0.9', 'phylop_mamm': '0.9', 'phylop_primate': '0.9',
        'spliceai_acceptor_gain': '0.6', 'spliceai_acceptor_loss': '0.6',
        'spliceai_donor_gain': '0.6', 'spliceai_donor_loss': '0.6', 'loftool': '0.6'
    }
    example_value: str = tool_examples[tool] if tool in tool_examples else (f"{min_value + 0.1:.2f}" if min_value != -100.0 else "0.01")
    example_fraction = "1/1000" if max_value <= 1.0 else "5/2"
    reference_text = f" - Check reference: {reference_url}\n" if reference_url else ""
    
    full_prompt = (
        f"{prompt}\n"
        f" - Expected range: {min_value} to {max_value}\n"
        f"{reference_text}"
        f" - Example: {example_value}, {example_fraction}, or 1e-3 (scientific notation)\n"
        f" - Press Enter if unknown{' (sets to 0.0)' if allow_zero else ''}: "
    ) if first_time else (
        f"{prompt}\n"
        f" - Example: {example_value}, {example_fraction}\n"
        f" - Press Enter if unknown{' (sets to 0.0)' if allow_zero else ''}: "
    )
    
    attempt_count = 0
    max_attempts = 3
    
    while attempt_count < max_attempts:
        try:
            value = input(full_prompt).strip()
            if not value:
                if allow_zero:
                    _log_info("No value provided, setting to 0.0 (unknown).")
                    first_time = False
                    return 0.0
                else:
                    _log_error("A value is required for this field.")
                    attempt_count += 1
                    continue
            
            if '/' in value:
                parts = value.split('/')
                if len(parts) != 2:
                    _log_error("Fraction must be in the format 'number/number' (e.g., 1/1000).")
                    attempt_count += 1
                    continue
                num = float(parts[0].replace(',', '').replace('<', ''))
                den = float(parts[1].replace(',', ''))
                if den == 0:
                    _log_error("Denominator cannot be zero.")
                    attempt_count += 1
                    continue
                result = num / den
            else:
                result = float(value.replace(',', ''))
            
            if not (min_value <= result <= max_value):
                _log_error(f"Value must be between {min_value} and {max_value}. Try {example_value}.")
                attempt_count += 1
                continue
            
            first_time = False
            return result
        
        except (ValueError, IndexError):
            _log_error(f"Enter a valid number, fraction, or scientific notation (e.g., {example_value}, {example_fraction}, 1e-3).")
            attempt_count += 1
    
    _log_warning(f"Maximum attempts ({max_attempts}) reached. Setting to {'0.0' if allow_zero else 'minimum value'}.")
    first_time = False
    return 0.0 if allow_zero else min_value

def validate_int_input(prompt, min_value=0):
    """Validate integer input with enhanced user guidance."""
    global first_time
    example_value = min_value + 1
    full_prompt = (
        f"{prompt}\n"
        f" - Expected: Integer >= {min_value}, e.g., {example_value}\n"
        f" - Press Enter if unknown (sets to 0): "
    ) if first_time else (
        f"{prompt}\n"
        f" - Example: {example_value}\n"
        f" - Press Enter if unknown: "
    )
    
    attempt_count = 0
    max_attempts = 3
    
    while attempt_count < max_attempts:
        try:
            value = input(full_prompt).strip()
            if not value:
                _log_info("No value provided, setting to 0 (unknown).")
                first_time = False
                return 0
            
            result = int(value)
            if result < min_value:
                _log_error(f"Value must be {min_value} or greater. Try {example_value}.")
                attempt_count += 1
                continue
            
            first_time = False
            return result
        
        except ValueError:
            _log_error(f"Enter a valid integer, e.g., {example_value}.")
            attempt_count += 1
    
    _log_warning(f"Maximum attempts ({max_attempts}) reached. Setting to 0.")
    first_time = False
    return 0

def validate_choice_input(prompt, valid_choices):
    """Validate choice input from a list of valid options."""
    global first_time
    valid_choices_str = ", ".join(valid_choices)
    full_prompt = f"{prompt}\n - Options: {valid_choices_str}\n - Press Enter if unknown: " if first_time else f"{prompt}\n - Options: {valid_choices_str}: "
    while True:
        value = input(full_prompt).strip().lower()
        if not value:
            return "unknown"
        if value in [choice.lower() for choice in valid_choices]:
            first_time = False
            return value
        _log_error(f"Enter one of: {valid_choices_str}")

def validate_hgvs_input(prompt, prefix):
    """Validate HGVS format input."""
    hgvs_pattern = re.compile(rf"^{prefix}\..*")
    full_prompt = f"{prompt}\n - Example: {prefix}.181T>G or {prefix}.123_124del\n - Press Enter if unknown: "
    while True:
        value = input(full_prompt).strip()
        if not value:
            return "Not specified"
        if hgvs_pattern.match(value):
            return value
        _log_error(f"Enter a valid HGVS format (e.g., {prefix}.181T>G or {prefix}.123_124del).")

def validate_allele(prompt, field_name):
    """Validate allele input (A, C, G, T)."""
    valid_chars = set('ACGT')
    full_prompt = f"{prompt}\n - Example: A, T, AGC\n - Press Enter if unknown: "
    while True:
        value = input(full_prompt).strip().upper()
        if not value:
            return "Not specified"
        if all(char in valid_chars for char in value):
            return value
        _log_error(f"{field_name} must contain only A, C, G, T characters (e.g., A, AGC).")

def validate_chromosome(prompt, auto_chromosome="Not specified"):
    """Validate chromosome input."""
    valid_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
    full_prompt = f"{prompt}\n - Enter chromosome (e.g., 1, X, Y)\n - Press Enter if unknown: "
    while True:
        value = input(full_prompt).strip().upper()
        if not value:
            return auto_chromosome if auto_chromosome != "Not specified" else "Not specified"
        if value in valid_chromosomes:
            return value
        _log_error(f"Enter a valid chromosome (e.g., 1, X, Y).")

def validate_gene_name(prompt):
    """Validate gene name input (basic check for non-empty string)."""
    full_prompt = f"{prompt}\n - Example: BRCA1, CFTR, CAMTA1\n - Press Enter if unknown: "
    while True:
        value = input(full_prompt).strip().upper()
        if not value:
            _log_error("Gene name is required for accurate classification.")
            continue
        if re.match(r'^[A-Z0-9\-]+$', value):
            return value
        _log_error("Enter a valid HGNC gene symbol (e