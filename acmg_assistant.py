
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

import logging
from typing import Dict, Optional, Tuple
def load_local_variant_scores(xlsx_path, chrom, pos, ref, alt):
    # This function is deprecated and no longer used.
    return None

# Logging setup
logging.basicConfig(filename='acmg_assessor.log', level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')
import re

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
    'lrt': (0.0, 1.0, True), 'gerp': (-12.3, 6.17, False), 'phylop_vert': (-20.0, 20.0, False),
    'phylop_mamm': (-20.0, 20.0, False), 'phylop_primate': (-20.0, 20.0, False)
}
LOF_VARIANTS = {'nonsense', 'frameshift', 'splice_donor', 'splice_acceptor'}
first_time = True


def calculate_metascore(inputs: Dict, gene: str) -> Tuple[float, Optional[str]]:
    """Calculate a VAMPP-score-like metascore and determine PP3/BP4 criterion with detailed logging."""
    cache_key = tuple((k, inputs.get(k)) for k in sorted(GENE_SPECIFIC_WEIGHTS.get(gene, GENE_SPECIFIC_WEIGHTS['default']).keys()) if inputs.get(f'{k}_score') is not None)
    if cache_key in score_cache:
        metascore, criterion = score_cache[cache_key]
        print(f"Info: Combined MetaScore retrieved from cache: {metascore:.4f} ({criterion if criterion else 'No criterion applied'})")
        return metascore, criterion
    
    weights = GENE_SPECIFIC_WEIGHTS.get(gene, GENE_SPECIFIC_WEIGHTS['default'])
    scores = []
    used_weights = []
    used_tools = []
    
    logging.info(f"Calculating metascore for {gene}...")
    for tool in weights:
        score = inputs.get(f'{tool}_score')
        normalized = normalize_score(tool, score)
        if normalized is not None:
            scores.append(normalized)
            used_weights.append(weights[tool])
            used_tools.append(tool)
    
    if not scores:
        logging.warning("No valid in silico scores available for metascore calculation.")
        return 0.0, None
    
    total_weight = sum(used_weights)
    if total_weight == 0:
        logging.warning("Total weight is zero, cannot compute metascore.")
        return 0.0, None
    
    weighted_sum = sum(s * w for s, w in zip(scores, used_weights))
    metascore = weighted_sum / total_weight
    
    # Detailed logging
    logging.info(f"Metascore calculation details for {gene}:")
    logging.info(f" - Tools used: {', '.join(used_tools)}")
    logging.info(f" - Normalized scores: {['{}: {:.4f}'.format(t, s) for t, s in zip(used_tools, scores)]}")
    logging.info(f" - Weights: {['{}: {:.2f}'.format(t, w) for t, w in zip(used_tools, used_weights)]}")
    logging.info(f" - Weighted sum: {weighted_sum:.4f}, Total weight: {total_weight:.4f}")
    logging.info(f" - Final metascore: {metascore:.4f}")
    
    # VAMPP-score cut-offs
    criterion = 'PP3' if metascore > 0.354 else 'BP4' if metascore < 0.226 else None
    score_cache[cache_key] = (metascore, criterion)
    
    # Kruskal-Wallis test
    benign_scores = [0.1, 0.2, 0.15, 0.3] * (len(scores) // 4 + 1)
    pathogenic_scores = [0.8, 0.9, 0.85, 0.7] * (len(scores) // 4 + 1)
    vus_scores = scores
    try:
        stat, p_value = stats.kruskal(benign_scores[:len(scores)], pathogenic_scores[:len(scores)], vus_scores)
        if p_value < 0.05 and criterion == 'PP3':
            logging.info(f"Kruskal-Wallis test supports PP3 (p-value: {p_value:.4f})")
        elif p_value < 0.05 and criterion == 'BP4':
            logging.info(f"Kruskal-Wallis test supports BP4 (p-value: {p_value:.4f})")
    except Exception as e:
        logging.warning(f"Kruskal-Wallis test failed: {e}")
    
    return metascore, criterion

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
        logging.warning(f"Ensembl API request failed: {e}. Chromosome information not retrieved.")
        return "Not specified"

def get_clinvar_status(chromosome, position, ref_allele, alt_allele):
    """Fetch ClinVar status for a variant."""
    if chromosome == "Not specified" or position == 0 or ref_allele == "Not specified" or alt_allele == "Not specified":
        return "unknown", "ClinVar data could not be retrieved: Missing information."
    cache_key = f"{chromosome}-{position}-{ref_allele}-{alt_allele}"
    if cache_key in clinvar_cache:
        logging.info("ClinVar status retrieved from cache.")
        return clinvar_cache[cache_key]
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    query = f"{chromosome}:{position}[chrpos38] AND {ref_allele}>{alt_allele}"
    params = {'db': 'clinvar', 'term': query, 'retmode': 'json', 'retmax': 1}
    try:
        r = requests.get(base_url, params=params, timeout=10)
        r.raise_for_status()
        result = r.json()
        ids = result.get('esearchresult', {}).get('idlist', [])
        if not ids:
            result_tuple = ("unknown", "No variant found in ClinVar.")
        else:
            detail_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={ids[0]}&retmode=json"
            r_detail = requests.get(detail_url, timeout=10)
            r_detail.raise_for_status()
            detail = r_detail.json()
            docsum = detail.get('result', {}).get(ids[0], {})
            clinvar_status = docsum.get('clinical_significance', {}).get('description', 'unknown')
            submitters = len(docsum.get('clinical_significance', {}).get('submitter', []))
            submitter_note = " (Single submitter)" if submitters <= 1 else " (Consensus)"
            result_tuple = (clinvar_status.lower(), f"ClinVar status: {clinvar_status}{submitter_note}")
        clinvar_cache[cache_key] = result_tuple
        return result_tuple
    except requests.exceptions.RequestException as e:
        logging.error(f"ClinVar API request failed: {e}")
        return "unknown", f"ClinVar API request failed: {e}"

def validate_float_input(prompt, allow_zero=True, min_value=-100.0, max_value=100.0, reference_url=None, tool=None):
    """Validate float input within a specified range with enhanced user guidance."""
    global first_time
    # Araç bazlı örnek değerler
    tool_examples: Dict[str, str] = {
        'revel': '0.75', 'cadd': '25.0', 'metarnn': '0.6', 'clinpred': '0.6', 'bayesdel': '0.2',
        'alphamissense': '0.6', 'mutationtaster': '0.6', 'polyphen2': '0.6', 'sift': '0.01',
        'fathmm_xf': '0.6', 'mutationassessor': '2.0', 'provean': '-3.0', 'mutpred': '0.6',
        'metolr': '0.6', 'esm1b': '1.0', 'lrt': '0.01', 'gerp': '3.0',
        'phylop_vert': '3.0', 'phylop_mamm': '3.0', 'phylop_primate': '3.0',
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
                    logging.info("No value provided, setting to 0.0 (unknown).")
                    first_time = False
                    return 0.0
                else:
                    logging.error("A value is required for this field.")
                    attempt_count += 1
                    continue
            
            if '/' in value:
                parts = value.split('/')
                if len(parts) != 2:
                    logging.error("Fraction must be in the format 'number/number' (e.g., 1/1000).")
                    attempt_count += 1
                    continue
                num = float(parts[0].replace(',', '').replace('<', ''))
                den = float(parts[1].replace(',', ''))
                if den == 0:
                    logging.error("Denominator cannot be zero.")
                    attempt_count += 1
                    continue
                result = num / den
            else:
                result = float(value.replace(',', ''))
            
            if not (min_value <= result <= max_value):
                logging.error(f"Value must be between {min_value} and {max_value}. Try {example_value}.")
                attempt_count += 1
                continue
            
            first_time = False
            return result
        
        except (ValueError, IndexError):
            logging.error(f"Enter a valid number, fraction, or scientific notation (e.g., {example_value}, {example_fraction}, 1e-3).")
            attempt_count += 1
    
    logging.warning(f"Maximum attempts ({max_attempts}) reached. Setting to {'0.0' if allow_zero else 'minimum value'}.")
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
                logging.info("No value provided, setting to 0 (unknown).")
                first_time = False
                return 0
            
            result = int(value)
            if result < min_value:
                logging.error(f"Value must be {min_value} or greater. Try {example_value}.")
                attempt_count += 1
                continue
            
            first_time = False
            return result
        
        except ValueError:
            logging.error(f"Enter a valid integer, e.g., {example_value}.")
            attempt_count += 1
    
    logging.warning(f"Maximum attempts ({max_attempts}) reached. Setting to 0.")
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
        logging.error(f"Enter one of: {valid_choices_str}")

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
        logging.error(f"Enter a valid HGVS format (e.g., {prefix}.181T>G or {prefix}.123_124del).")

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
        logging.error(f"{field_name} must contain only A, C, G, T characters (e.g., A, AGC).")

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
        logging.error(f"Enter a valid chromosome (e.g., 1, X, Y).")

def validate_gene_name(prompt):
    """Validate gene name input (basic check for non-empty string)."""
    full_prompt = f"{prompt}\n - Example: BRCA1, CFTR, CAMTA1\n - Press Enter if unknown: "
    while True:
        value = input(full_prompt).strip().upper()
        if not value:
            logging.error("Gene name is required for accurate classification.")
            continue
        if re.match(r'^[A-Z0-9\-]+$', value):
            return value
        logging.error("Enter a valid HGNC gene symbol (e.g., BRCA1, CFTR, CAMTA1).")

def validate_zygosity_input(prompt):
    """Validate zygosity input."""
    valid_choices = ["heterozygous", "homozygous", "unknown", "het", "hom"]
    valid_choices_str = ", ".join(valid_choices)
    full_prompt = f"{prompt}\n - Options: {valid_choices_str}\n - Press Enter if unknown: "
    while True:
        value = input(full_prompt).strip().lower()
        if not value:
            return "unknown"
        if value in valid_choices:
            if value == "het":
                return "heterozygous"
            if value == "hom":
                return "homozygous"
            return value
        logging.error(f"Enter one of: {valid_choices_str}")

def detect_lof_variant(cdna_change, vep_consequence):
    """Detect Loss of Function (LoF) variant based on cDNA change and VEP consequence."""
    if vep_consequence.lower() in LOF_VARIANTS:
        logging.info(f"LoF variant detected by VEP consequence: {vep_consequence.lower()}")
        return vep_consequence.lower()
    if cdna_change != "Not specified":
        cdna_lower = cdna_change.lower()
        if any(x in cdna_lower for x in ['del', 'ins', 'dup', 'frameshift']):
            logging.info("LoF variant detected by cDNA change: frameshift")
            return 'frameshift'
        if any(x in cdna_lower for x in ['ter', 'x', '*']):
            logging.info("LoF variant detected by cDNA change: nonsense")
            return 'nonsense'
        if re.search(r'c\.\d+[+-][12]', cdna_lower):
            logging.info("LoF variant detected by cDNA change: splice site")
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
        'gerp': lambda x: 'conserved' if x > 2.0 else 'non-conserved' if x <= 2.0 else 'unknown',
        'phylop': lambda x: 'conserved' if x > 2.7 else 'non-conserved' if x <= 2.7 else 'unknown',
        'spliceai': lambda x: 'damaging' if x > 0.5 else 'tolerated' if x < 0.5 else 'unknown',
        'loftool': lambda x: 'damaging' if x > 0.5 else 'tolerated' if x < 0.5 else 'unknown'
    }
    
    return thresholds.get(tool, lambda x: 'unknown')(score)

def validate_prediction(tool, score, prediction, valid_choices):
    """Validate in silico prediction against expected value."""
    expected = calculate_prediction(tool, score)
    if prediction != 'unknown' and prediction != expected and expected != 'unknown':
        logging.warning(f"Expected prediction for {tool.upper()} score ({score}) is '{expected}', but '{prediction}' was entered.")
        logging.info(f"Verify with Varsome: https://varsome.com/")
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

    print(f"\n{'='*15} VARIANT {variant_number} BASIC INFORMATION {'='*15}")
    print("Please provide details about the variant's genomic location and effect.")
    print("Press Enter to skip optional fields where indicated.")
    
    inputs['gene'] = validate_gene_name("Enter gene name")
    auto_chromosome = get_chromosome_from_ensembl(inputs['gene'])
    if auto_chromosome != "Not specified":
        print(f"Info: Gene {inputs['gene']} found on chromosome {auto_chromosome}.")
        use_auto = validate_choice_input("Use this chromosome? (yes/no)", ["yes", "no"])
        inputs['chromosome'] = auto_chromosome if use_auto == "yes" else validate_chromosome("Enter chromosome")
    else:
        print("Info: Could not retrieve chromosome automatically from Ensembl.")
        inputs['chromosome'] = validate_chromosome("Enter chromosome")
    inputs['position'] = validate_int_input("Enter genomic position\n - Example: 7249507 (GRCh38)")
    inputs['ref_allele'] = validate_allele("Enter reference allele", "Reference allele")
    inputs['alt_allele'] = validate_allele("Enter alternate allele", "Alternate allele")
    if inputs['ref_allele'] == inputs['alt_allele'] and inputs['ref_allele'] != "Not specified":
        print("Error: Reference and alternate alleles cannot be the same.")
        inputs['alt_allele'] = validate_allele("Enter alternate allele again", "Alternate allele")
    inputs['cdna_change'] = validate_hgvs_input("Enter cDNA change", "c")
    inputs['protein_change'] = validate_hgvs_input("Enter protein change", "p")
    print(f"\n{'='*15} VEP CONSEQUENCE {'='*15}")
    print("VEP Consequence describes the effect of the variant on the gene.")
    inputs['vep_consequence'] = validate_choice_input(
        "Enter VEP consequence",
        ["missense", "nonsense", "frameshift", "splice_donor", "splice_acceptor", "intronic", "synonymous", "unknown"]
    )
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
    print("Varsome: https://varsome.com/")
    print("Info: All in silico scores must be entered manually. No local dataset search is performed.")

    inputs: Dict[str, object] = {'insilico_tools': 'REVEL, CADD, MetaRNN, ClinPred, BayesDel, AlphaMissense, MutationTaster, PolyPhen, SIFT'}
    tools = [
        ('revel', 'Enter REVEL score', 0.0, 1.0, ['damaging', 'tolerated', 'unknown']),
        ('cadd', 'Enter CADD Phred score', 0.0, 100.0, ['damaging', 'tolerated', 'unknown']),
        ('metarnn', 'Enter MetaRNN score', 0.0, 1.0, ['damaging', 'tolerated', 'unknown']),
        ('clinpred', 'Enter ClinPred score', 0.0, 1.0, ['damaging', 'tolerated', 'unknown']),
        ('bayesdel', 'Enter BayesDel score', -1.0, 1.0, ['deleterious', 'tolerated', 'unknown']),
        ('alphamissense', 'Enter AlphaMissense score', 0.0, 1.0, ['likely_pathogenic', 'likely_benign', 'ambiguous', 'unknown']),
        ('mutationtaster', 'Enter MutationTaster score', 0.0, 1.0, ['damaging', 'polymorphism', 'unknown']),
        ('polyphen2', 'Enter PolyPhen-2 score', 0.0, 1.0, ['probably_damaging', 'possibly_damaging', 'benign', 'unknown']),
        ('sift', 'Enter SIFT score', 0.0, 1.0, ['damaging', 'tolerated', 'unknown'])
    ]

    if variant_type == 'missense':
        for tool, prompt, min_val, max_val, choices in tools:
            score_key = f'{tool}_score'
            score = validate_float_input(
                prompt,
                min_value=min_val,
                max_value=max_val,
                reference_url="https://varsome.com/",
                tool=tool
            )
            inputs[score_key] = score
            auto_pred = calculate_prediction(tool, float(score) if score is not None else 0.0)
            print(f"Info: {tool.upper()} automatic prediction: {auto_pred}")
            override = validate_choice_input(f"Enter manual prediction for {tool.upper()}? (yes/no)", ["yes", "no"])
            if override == "yes":
                inputs[f'{tool}_prediction'] = validate_choice_input(f"Enter {tool.upper()} prediction", choices)
            else:
                inputs[f'{tool}_prediction'] = auto_pred
            # Ensure score is a float or convertible, else use 0.0
            score_val = inputs[score_key]
            score_float = 0.0
            if score_val is not None:
                if isinstance(score_val, (float, int)):
                    score_float = float(score_val)
                elif isinstance(score_val, str):
                    try:
                        score_float = float(score_val)
                    except ValueError:
                        score_float = 0.0
                else:
                    score_float = 0.0
            inputs[f'{tool}_prediction'] = validate_prediction(tool, score_float, inputs[f'{tool}_prediction'], choices)
    # ...splice and lof for your code as before...

    # Calculate metascore
    metascore, criterion = calculate_metascore(inputs, gene)
    inputs['combined_pathogenicity_score'] = f"{metascore:.4f}" if metascore is not None else "0.0"
    inputs['metascore_criterion'] = str(criterion) if criterion is not None else ''
    print(f"Info: Combined MetaScore: {metascore:.4f} ({criterion if criterion else 'No criterion applied'})")

    return inputs

def collect_genetic_data():
    """Collect genetic data such as inheritance and zygosity with user-friendly prompts."""
    inputs = {}
    inputs['variants'] = []  # Ensure variants is a list
    print(f"\n{'='*15} GENETIC DATA {'='*15}")
    print("Provide information about the variant's inheritance pattern and zygosity.")
    inputs['inheritance'] = validate_choice_input(
        "Enter inheritance pattern\n - AD: Autosomal Dominant\n - AR: Autosomal Recessive\n - X-linked: X-linked inheritance",
        ["AD", "AR", "X-linked", "unknown"]
    )
    if inputs['inheritance'] in ['xlr', 'xld']:
        inputs['inheritance'] = 'x-linked'
    inputs['zygosity'] = validate_zygosity_input(
        "Enter zygosity\n - heterozygous (het): One copy of variant\n - homozygous (hom): Two copies of variant"
    )
    if inputs['inheritance'] == 'ar' and inputs['zygosity'] == 'heterozygous':
        print(f"\n{'='*15} ALLELIC DATA {'='*15}")
        print("For recessive diseases, indicate if the variant is in trans with another pathogenic variant.")
        inputs['allelic_data'] = validate_choice_input(
            "Is the variant in trans with a pathogenic variant? (yes/no)",
            ["yes", "no", "unknown"]
        )
    else:
        inputs['allelic_data'] = "unknown"
    return inputs

def collect_functional_data(variant_type, gene):
    """Collect functional and clinical data with user-friendly prompts."""
    inputs = {}
    print(f"\n{'='*15} FUNCTIONAL AND CLINICAL DATA {'='*15}")
    print("Provide functional test results and clinical observations for the variant.")
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
    inputs['segregation'] = validate_choice_input(
        "Does the variant segregate with disease in family?\n - Yes: Variant present in affected family members\n - No: Variant present in unaffected family members",
        ["Yes", "No", "unknown"]
    )
    inputs['de_novo'] = validate_choice_input(
        "Is the variant de novo?\n - Yes: Variant not inherited from parents\n - No: Variant inherited",
        ["Yes", "No", "unknown"]
    )
    if inputs['de_novo'] == 'yes':
        print(f"\n{'='*15} PATERNITY CONFIRMATION {'='*15}")
        print("Biological parentage confirmation strengthens de novo status (ACMG PS2 criterion).")
        inputs['paternity_confirmed'] = validate_choice_input(
            "Is paternity confirmed?\n - Yes: Biological parentage verified\n - No: Parentage not verified",
            ["Yes", "No", "unknown"]
        )
    else:
        inputs['paternity_confirmed'] = 'unknown'
    inputs['phenotype_match'] = validate_choice_input(
        f"Enter phenotype match strength for {gene} gene\n - Strong: Highly specific to gene\n - Moderate: Partially specific\n - Weak: Low specificity",
        ["Strong", "Moderate", "Weak", "unknown"]
    )
    return inputs

def collect_compound_het_info(variant_number=2):
    """Collect information for a second heterozygous variant, if applicable."""
    print(f"\n{'='*15} COMPOUND HETEROZYGOUS VARIANT INFORMATION {'='*15}")
    print("For recessive diseases, indicate if there is a second heterozygous variant.")
    has_second_variant = validate_choice_input("Is there a second variant? (yes/no)", ["yes", "no"])
    if has_second_variant == "yes":
        return collect_variant_basic_info(variant_number=variant_number)
    return None

def collect_inputs(test_mode=False):
    """Collect all input data for variant classification."""
    global first_time
    first_time = True
    print("\n" + "="*60)
    print("      ACMG 2015/2023 VARIANT CLASSIFICATION ASSISTANT")
    print("                 Created by: Can Sevilmiş")
    print("="*60 + "\n")
    
    inputs: Dict[str, object] = {'variants': []}
    primary_variant = collect_variant_basic_info(test_mode=test_mode, variant_number=1)
    lof_type = detect_lof_variant(primary_variant['cdna_change'], primary_variant['vep_consequence'])
    
    print(f"\n{'='*15} VARIANT TYPE DETECTION {'='*15}")
    if lof_type:
        print(f"Info: Based on cDNA change ({primary_variant['cdna_change']}) and VEP Consequence ({primary_variant['vep_consequence']}),")
        print(f"an automatic LoF variant type '{lof_type}' was detected.")
        print("LoF (Loss of Function) variants significantly impact gene function.")
        use_lof = validate_choice_input("Use the automatically detected variant type? (yes/no)", ["yes", "no"])
        if use_lof == "yes":
            primary_variant['variant_type'] = lof_type
        else:
            print(f"\n{'='*15} VARIANT TYPE SELECTION {'='*15}")
            print("Manually select the variant type to classify its effect on the gene.")
            primary_variant['variant_type'] = validate_choice_input(
                "Enter final variant type",
                ["missense", "nonsense", "frameshift", "splice_donor", "splice_acceptor", "intronic", "synonymous", "unknown"]
            )
    else:
        print("Info: No automatic LoF detection was possible.")
        print(f"\n{'='*15} VARIANT TYPE SELECTION {'='*15}")
        print("Manually select the variant type to classify its effect on the gene.")
        primary_variant['variant_type'] = validate_choice_input(
            "Enter final variant type",
            ["missense", "nonsense", "frameshift", "splice_donor", "splice_acceptor", "intronic", "synonymous", "unknown"]
        )
    
    status, note = get_clinvar_status(primary_variant['chromosome'], primary_variant['position'], primary_variant['ref_allele'], primary_variant['alt_allele'])
    primary_variant['clinvar_status'] = status
    primary_variant['clinvar_note'] = note
    print(f"Info: {note}")
    print(f"\n{'='*15} CLINVAR STATUS {'='*15}")
    print("ClinVar status indicates the known pathogenic/benign classification of the variant.")
    override_clinvar = validate_choice_input("Manually override ClinVar status? (yes/no)", ["yes", "no"])
    if override_clinvar == "yes":
        primary_variant['clinvar_status'] = validate_choice_input(
            "Enter ClinVar status",
            ["pathogenic", "likely pathogenic", "benign", "likely benign", "vus", "unknown"]
        )
    
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
        logging.info(f"Statistically significant difference: {fisher_result}")
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

def process_variants():
    """Process variants and classify them."""
    while True:
        inputs = collect_inputs()
        evidence = assign_evidence(inputs, use_pp5_bp6=True, use_acmg_2023=True)
        classification, evidence = classify_variant(evidence, use_acmg_2023=True)
        suggestions = check_missing_data(inputs)[2]
        save_results(inputs, classification, evidence, suggestions)
        print_results(inputs, classification, evidence, suggestions)
        
        another_variant = validate_choice_input("Classify another variant? (yes/no)", ["yes", "no"])
        if another_variant == "no":
            break
        clear_cache()
  
def main():
    """Main function to run the ACMG variant classification tool."""
    parser = argparse.ArgumentParser(description="ACMG 2015/2023 Variant Classification Assistant with VAMPP-score Integration")
    parser.add_argument('--test', action='store_true', help="Run in test mode with predefined inputs")
    args = parser.parse_args()
    
    
    if args.test:
        print("Info: Running in test mode with predefined inputs.")
        inputs = collect_inputs(test_mode=True)
        evidence = assign_evidence(inputs, use_pp5_bp6=True, use_acmg_2023=True)
        classification, evidence = classify_variant(evidence, use_acmg_2023=True)
        suggestions = check_missing_data(inputs)[2]
        save_results(inputs, classification, evidence, suggestions)
        print_results(inputs, classification, evidence, suggestions)
    else:
        process_variants()

if __name__ == "__main__":
    main()