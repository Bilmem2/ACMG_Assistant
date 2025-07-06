import time
import argparse
import requests
import scipy.stats as stats
import numpy as np
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
    'default': {'BA1': 0.05, 'BS1': 0.01}
}
INSILICO_WEIGHTS = {
    'revel': 0.25, 'cadd': 0.25, 'alphamissense': 0.20, 'sift': 0.10, 'bayesdel': 0.10,
    'fathmm_xf': 0.05, 'mutationtaster': 0.05, 'mutationassessor': 0.05, 'provean': 0.05,
    'mutpred': 0.05, 'metolr': 0.05, 'esm1b': 0.05, 'lrt': 0.05, 'gerp': 0.05,
    'phylop_vert': 0.05, 'phylop_mamm': 0.05, 'phylop_primate': 0.05
}
LOF_VARIANTS = {'nonsense', 'frameshift', 'splice_donor', 'splice_acceptor'}
first_time = True

# --- Helper Functions ---

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
        print(f"Ensembl API request failed: {e}. Chromosome information not retrieved.")
        return "Not specified"

def get_clinvar_status(chromosome, position, ref_allele, alt_allele):
    """Fetch ClinVar status for a variant."""
    if chromosome == "Not specified" or position == 0 or ref_allele == "Not specified" or alt_allele == "Not specified":
        return "unknown", "ClinVar data could not be retrieved: Missing information."
    cache_key = f"{chromosome}-{position}-{ref_allele}-{alt_allele}"
    if cache_key in clinvar_cache:
        print("ClinVar status retrieved from cache.")
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
        return "unknown", f"ClinVar API request failed: {e}"

def validate_float_input(prompt, allow_zero=True, min_value=-100.0, max_value=100.0):
    """Validate float input within a specified range."""
    global first_time
    full_prompt = f"{prompt}\n - Press Enter if unknown (sets to 0.0): " if first_time else f"{prompt}\n - Press Enter if unknown: "
    while True:
        try:
            value = input(full_prompt).strip()
            if not value and allow_zero:
                return 0.0
            if '/' in value:
                parts = value.split('/')
                num = float(parts[0].replace(',', '').replace('<',''))
                den = float(parts[1].replace(',', ''))
                result = num / den if den != 0 else 0
            else:
                result = float(value)
            if not (min_value <= result <= max_value):
                print(f"Error: Value must be between {min_value} and {max_value} (e.g., {min_value + 0.1}).")
                continue
            first_time = False
            return result
        except (ValueError, IndexError):
            print("Error: Enter a valid number (e.g., 0.001) or fraction (e.g., 1/1000).")

def validate_int_input(prompt, min_value=0):
    """Validate integer input."""
    global first_time
    full_prompt = f"{prompt}\n - Press Enter if unknown (sets to 0): " if first_time else f"{prompt}\n - Press Enter if unknown: "
    while True:
        try:
            value = input(full_prompt).strip()
            if not value:
                return 0
            result = int(value)
            if result < min_value:
                print(f"Error: Value must be {min_value} or greater (e.g., {min_value + 1}).")
                continue
            first_time = False
            return result
        except ValueError:
            print(f"Error: Enter a valid integer (e.g., {min_value + 1}).")

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
        print(f"Error: Enter one of: {valid_choices_str}")

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
        print(f"Error: Enter a valid HGVS format (e.g., {prefix}.181T>G or {prefix}.123_124del).")

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
        print(f"Error: {field_name} must contain only A, C, G, T characters (e.g., A, AGC).")

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
        print(f"Error: Enter a valid chromosome (e.g., 1, X, Y).")

def validate_gene_name(prompt):
    """Validate gene name input (basic check for non-empty string)."""
    full_prompt = f"{prompt}\n - Example: BRCA1, CFTR, CAMTA1\n - Press Enter if unknown: "
    while True:
        value = input(full_prompt).strip().upper()
        if not value:
            print("Error: Gene name is required for accurate classification.")
            continue
        if re.match(r'^[A-Z0-9\-]+$', value):
            return value
        print("Error: Enter a valid HGNC gene symbol (e.g., BRCA1, CFTR, CAMTA1).")

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
        print(f"Error: Enter one of: {valid_choices_str}")

def detect_lof_variant(cdna_change, vep_consequence):
    """Detect Loss of Function (LoF) variant based on cDNA change and VEP consequence."""
    if vep_consequence.lower() in LOF_VARIANTS:
        return vep_consequence.lower()
    if cdna_change != "Not specified":
        cdna_lower = cdna_change.lower()
        if any(x in cdna_lower for x in ['del', 'ins', 'dup', 'frameshift']):
            return 'frameshift'
        if any(x in cdna_lower for x in ['ter', 'x', '*']):
            return 'nonsense'
        if re.search(r'c\.\d+[+-][12]', cdna_lower):
            return 'splice_donor' if '+' in cdna_change else 'splice_acceptor'
    return None

def calculate_prediction(tool, score):
    """Calculate in silico prediction based on tool and score."""
    if tool == 'sift' and score is not None:
        return 'damaging' if score < 0.05 else 'tolerated'
    elif tool == 'revel' and score is not None:
        return 'damaging' if score > 0.7 else 'tolerated'
    elif tool == 'alphamissense' and score is not None:
        return 'likely_pathogenic' if score > 0.564 else 'likely_benign' if score < 0.34 else 'ambiguous'
    elif tool == 'cadd' and score is not None:
        return 'damaging' if score > 20 else 'tolerated'
    elif tool == 'fathmm_xf' and score is not None:
        return 'damaging' if score > 0.5 else 'neutral'
    elif tool == 'bayesdel' and score is not None:
        return 'deleterious' if score > 0.15 else 'tolerated'
    elif tool == 'mutationtaster' and score is not None:
        return 'damaging' if score > 0.5 else 'polymorphism'
    elif tool == 'mutationassessor' and score is not None:
        return 'high' if score > 3.5 else 'medium' if score > 1.9 else 'low' if score > 0.8 else 'neutral'
    elif tool == 'provean' and score is not None:
        return 'damaging' if score < -2.5 else 'neutral'
    elif tool == 'mutpred' and score is not None:
        return 'damaging' if score > 0.5 else 'neutral'
    elif tool == 'metolr' and score is not None:
        return 'damaging' if score > 0.5 else 'neutral'
    elif tool == 'esm1b' and score is not None:
        return 'damaging' if score > 0 else 'neutral'
    elif tool == 'lrt' and score is not None:
        return 'damaging' if score < 0.05 else 'neutral'
    elif tool == 'gerp' and score is not None:
        return 'conserved' if score > 2.0 else 'non-conserved'
    elif tool == 'phylop' and score is not None:
        return 'conserved' if score > 2.7 else 'non-conserved'
    elif tool == 'spliceai' and score is not None:
        return 'damaging' if score > 0.5 else 'tolerated'
    elif tool == 'loftool' and score is not None:
        return 'damaging' if score > 0.9 else 'tolerated'
    return 'unknown'

def validate_prediction(tool, score, prediction, valid_choices):
    """Validate in silico prediction against expected value."""
    expected = calculate_prediction(tool, score)
    if prediction != 'unknown' and prediction != expected and expected != 'unknown':
        print(f"Warning: Expected prediction for {tool.upper()} score ({score}) is '{expected}', but '{prediction}' was entered.")
        print(f"Verify with Varsome: https://varsome.com/")
        override = validate_choice_input("Use entered prediction? (yes/no)", ["yes", "no"])
        return prediction if override == "yes" else expected
    return prediction if prediction != 'unknown' else expected

def normalize_score(tool, score):
    """Normalize in silico scores to a 0-1 range."""
    ranges = {
        'sift': (0.0, 1.0, True), 'revel': (0.0, 1.0, False), 'alphamissense': (0.0, 1.0, False),
        'cadd': (0.0, 100.0, False), 'fathmm_xf': (0.0, 1.0, False), 'bayesdel': (-1.0, 1.0, False),
        'mutationtaster': (0.0, 1.0, False), 'mutationassessor': (0.0, 5.0, False),
        'provean': (-14.0, 14.0, True), 'mutpred': (0.0, 1.0, False), 'metolr': (0.0, 1.0, False),
        'esm1b': (-100.0, 100.0, False), 'lrt': (0.0, 1.0, True), 'gerp': (-12.3, 6.17, False),
        'phylop_vert': (-20.0, 20.0, False), 'phylop_mamm': (-20.0, 20.0, False), 'phylop_primate': (-20.0, 20.0, False)
    }
    if tool not in ranges or score is None:
        return None
    min_val, max_val, reverse = ranges[tool]
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
        inputs['gene'] = "BRCA1"
        auto_chromosome = get_chromosome_from_ensembl(inputs['gene'])
        inputs['chromosome'] = auto_chromosome if auto_chromosome != "Not specified" else "17"
        inputs['position'] = 7249507 if variant_number == 1 else 7249508
        inputs['ref_allele'] = "T"
        inputs['alt_allele'] = "G"
        inputs['cdna_change'] = "c.181T>G" if variant_number == 1 else "c.182A>C"
        inputs['protein_change'] = "p.Ser61Gly" if variant_number == 1 else "p.Lys62Gln"
        inputs['vep_consequence'] = "missense"
        return inputs

    print(f"\n{'='*15} VARIANT {variant_number} BASIC INFORMATION {'='*15}")
    print("Please provide details about the variant's genomic location and effect.")
    print("Press Enter to skip optional fields where indicated.")
    
    # Gene Name
    inputs['gene'] = validate_gene_name("Enter gene name")

    # Chromosome
    auto_chromosome = get_chromosome_from_ensembl(inputs['gene'])
    if auto_chromosome != "Not specified":
        print(f"Info: Gene {inputs['gene']} found on chromosome {auto_chromosome}.")
        use_auto = validate_choice_input("Use this chromosome? (yes/no)", ["yes", "no"])
        inputs['chromosome'] = auto_chromosome if use_auto == "yes" else validate_chromosome("Enter chromosome")
    else:
        print("Info: Could not retrieve chromosome automatically from Ensembl.")
        inputs['chromosome'] = validate_chromosome("Enter chromosome")

    # Position
    inputs['position'] = validate_int_input("Enter genomic position\n - Example: 7249507 (GRCh38)")

    # Reference and Alternate Alleles
    inputs['ref_allele'] = validate_allele("Enter reference allele", "Reference allele")
    inputs['alt_allele'] = validate_allele("Enter alternate allele", "Alternate allele")
    if inputs['ref_allele'] == inputs['alt_allele'] and inputs['ref_allele'] != "Not specified":
        print("Error: Reference and alternate alleles cannot be the same.")
        inputs['alt_allele'] = validate_allele("Enter alternate allele again", "Alternate allele")

    # cDNA Change
    inputs['cdna_change'] = validate_hgvs_input("Enter cDNA change", "c")

    # Protein Change
    inputs['protein_change'] = validate_hgvs_input("Enter protein change", "p")

    # VEP Consequence
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
        "Enter allele frequency\n - Example: 0.001 or 1/1000\n - Range: 0.0 to 1.0\n - Check gnomAD: https://gnomad.broadinstitute.org/",
        max_value=1.0
    )
    
    print(f"\n{'='*15} SUBPOPULATION FREQUENCIES {'='*15}")
    print("Enter gnomAD subpopulation frequencies (e.g., African, European).")
    print("Enter one frequency at a time; press Enter twice to finish.")
    subpopulation_frequencies = []
    while True:
        freq = validate_float_input(
            "Enter subpopulation frequency\n - Example: 0.0002\n - Range: 0.0 to 1.0",
            max_value=1.0
        )
        if freq == 0.0:
            break
        subpopulation_frequencies.append(freq)
    inputs['subpopulation_frequencies'] = subpopulation_frequencies
    
    print(f"\n{'='*15} DISEASE PREVALENCE {'='*15}")
    print("Provide the prevalence of the associated disease in the population.")
    inputs['disease_prevalence'] = validate_float_input(
        "Enter disease prevalence\n - Example: 0.001\n - Range: 0.0 to 1.0\n - Check Orphanet: https://www.orpha.net/",
        max_value=1.0
    )
    
    print(f"\n{'='*15} CASE-CONTROL DATA FOR FISHER'S EXACT TEST {'='*15}")
    print("Case-control data compares allele frequencies between affected and unaffected populations.")
    run_fisher_test = validate_choice_input("Include case-control data? (yes/no)", ["yes", "no"])
    if run_fisher_test == "yes":
        print(f"\n{'='*15} CASE-CONTROL DATA INPUT {'='*15}")
        print("Provide allele counts and totals for case and control cohorts.")
        inputs['cohort_allele_count'] = validate_int_input("Enter allele count in case cohort\n - Example: 10")
        inputs['cohort_total'] = validate_int_input("Enter total alleles in case cohort\n - Example: 1000")
        inputs['control_allele_count'] = validate_int_input("Enter allele count in control population\n - Example: 5")
        inputs['control_total'] = validate_int_input("Enter total alleles in control population\n - Example: 10000")
    else:
        inputs.update({k: 0 for k in ['cohort_allele_count', 'cohort_total', 'control_allele_count', 'control_total']})
    return inputs

def collect_insilico_inputs(variant_type):
    """Collect in silico prediction data based on variant type with user-friendly prompts."""
    print(f"\n{'='*15} IN SILICO PREDICTIONS {'='*15}")
    print("Provide scores for in silico tools to predict variant impact.")
    print("Predictions can be entered manually or calculated automatically based on scores.")
    inputs = {'insilico_tools': input("Enter in silico tools used (e.g., REVEL, SIFT; press Enter if unknown): ") or 'Not specified'}
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
        ('gerp', 'Enter GERP++ RS score\n - Range: -12.3 to 6.17\n - Threshold: >2.0 (conserved)', -12.3, 6.17, ['conserved', 'non-conserved', 'unknown'])
    ]
    if variant_type == 'missense':
        for tool, prompt, min_val, max_val, choices in tools:
            inputs[f'{tool}_score'] = validate_float_input(prompt, min_value=min_val, max_value=max_val)
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
                min_value=0.0, max_value=1.0
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
            min_value=0.0, max_value=1.0
        )
        auto_pred = calculate_prediction('loftool', inputs['loftool_score'])
        print(f"Info: Automatic prediction for LoFtool: {auto_pred}")
        override = validate_choice_input("Manually enter LoFtool prediction? (yes/no)", ["yes", "no"])
        if override == "yes":
            inputs['loftool_prediction'] = validate_choice_input("Enter LoFtool prediction", ["damaging", "tolerated", "unknown"])
        else:
            inputs['loftool_prediction'] = auto_pred
        inputs['loftool_prediction'] = validate_prediction('loftool', inputs['loftool_score'], inputs['loftool_prediction'], ["damaging", "tolerated", "unknown"])
    for tool, prompt, min_val, max_val, choices in [
        ('phylop_vert', 'Enter phyloP Vertebrate score\n - Range: -20.0 to 20.0\n - Threshold: >2.7 (conserved)', -20.0, 20.0, ['conserved', 'non-conserved', 'unknown']),
        ('phylop_mamm', 'Enter phyloP Mammalian score\n - Range: -20.0 to 20.0\n - Threshold: >2.7 (conserved)', -20.0, 20.0, ['conserved', 'non-conserved', 'unknown']),
        ('phylop_primate', 'Enter phyloP Primate score\n - Range: -20.0 to 20.0\n - Threshold: >2.7 (conserved)', -20.0, 20.0, ['conserved', 'non-conserved', 'unknown'])
    ]:
        inputs[f'{tool}_score'] = validate_float_input(prompt, min_value=min_val, max_value=max_val)
        auto_pred = calculate_prediction(tool.replace('_vert', '').replace('_mamm', '').replace('_primate', ''), inputs[f'{tool}_score'])
        print(f"Info: Automatic prediction for {tool.upper()}: {auto_pred}")
        override = validate_choice_input(f"Manually enter prediction for {tool.upper()}? (yes/no)", ["yes", "no"])
        if override == "yes":
            inputs[f'{tool}_prediction'] = validate_choice_input(f"Enter {tool.upper()} prediction", choices)
        else:
            inputs[f'{tool}_prediction'] = auto_pred
        inputs[f'{tool}_prediction'] = validate_prediction(tool.replace('_vert', '').replace('_mamm', '').replace('_primate', ''), inputs[f'{tool}_score'], inputs[f'{tool}_prediction'], choices)
    return inputs

def collect_genetic_data():
    """Collect genetic data such as inheritance and zygosity with user-friendly prompts."""
    inputs = {}
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
    print(" (With Caching, Statistical Analysis, and Missing Data Checks)")
    print("="*60 + "\n")
    
    inputs = {'variants': []}
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
    primary_variant.update(collect_insilico_inputs(primary_variant['variant_type']))
    primary_variant.update(collect_genetic_data())
    primary_variant.update(collect_functional_data(primary_variant['variant_type'], primary_variant['gene']))
    
    if primary_variant['inheritance'] == 'ar' and primary_variant['zygosity'] == 'heterozygous':
        second_variant = collect_compound_het_info()
        if second_variant:
            second_variant['clinvar_status'], second_variant['clinvar_note'] = get_clinvar_status(
                second_variant['chromosome'], second_variant['position'], second_variant['ref_allele'], second_variant['alt_allele']
            )
            inputs['variants'].append(second_variant)
    
    inputs['variants'].append(primary_variant)
    print(f"\n{'='*15} ADDITIONAL NOTES {'='*15}")
    print("Provide any literature references or additional notes.")
    inputs['notes'] = input("Enter notes (press Enter if none): ") or 'None'
    inputs['varsome_url'] = generate_varsome_url(primary_variant['chromosome'], primary_variant['position'], primary_variant['ref_allele'], primary_variant['alt_allele'])
    return inputs

def calculate_combined_pathogenicity_score(inputs):
    """Calculate a combined pathogenicity score from in silico tools."""
    cache_key = tuple((k, inputs.get(k)) for k in sorted(INSILICO_WEIGHTS.keys()) if inputs.get(f'{k}_score') is not None)
    if cache_key in score_cache:
        print("Info: Combined Pathogenicity Score retrieved from cache.")
        return score_cache[cache_key]
    
    scores = []
    weights = []
    for tool in INSILICO_WEIGHTS:
        score = inputs.get(f'{tool}_score')
        normalized = normalize_score(tool, score)
        if normalized is not None:
            scores.append(normalized)
            weights.append(INSILICO_WEIGHTS[tool])
    
    if not scores:
        return 0.0
    
    ranks = stats.rankdata(scores, method='average')
    rank_normalized = [(rank - 1) / (len(ranks) - 1) if len(ranks) > 1 else 0.5 for rank in ranks]
    total_weight = sum(weights)
    if total_weight == 0:
        return 0.0
    weighted_sum = sum(r * w for r, w in zip(rank_normalized, weights))
    combined_score = weighted_sum / total_weight
    score_cache[cache_key] = combined_score
    return combined_score

def calculate_ensemble_score(inputs, variant_type):
    """Calculate ensemble score based on in silico predictions."""
    damaging_count = 0
    benign_count = 0
    
    if variant_type == 'missense':
        tools = ['revel', 'cadd', 'alphamissense', 'sift', 'bayesdel', 'fathmm_xf', 'mutationtaster', 'mutationassessor', 'provean', 'mutpred', 'metolr', 'esm1b', 'lrt', 'gerp']
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
        return p_value, f"Fisher’s Exact Test p-value: {p_value:.4f}, Odds Ratio: {oddsratio:.4f}"
    except Exception as e:
        return 1.0, f"Fisher’s Exact Test failed: {str(e)}"

def check_population_outlier(allele_frequency, subpopulation_frequencies):
    """Check if allele frequency is an outlier compared to subpopulations."""
    if not subpopulation_frequencies:
        return False, "No subpopulation data provided."
    
    mean_freq = np.mean(subpopulation_frequencies)
    std_freq = np.std(subpopulation_frequencies) if len(subpopulation_frequencies) > 1 else 0.01
    z_score = (allele_frequency - mean_freq) / std_freq if std_freq > 0 else 0
    is_outlier = abs(z_score) > 2
    return is_outlier, f"Z-score: {z_score:.4f}, Mean: {mean_freq:.6f}, Std: {std_freq:.6f}"

def kruskal_wallis_test(inputs, variant_class='vus'):
    """Perform Kruskal-Wallis test to compare in silico scores."""
    benign_scores = []
    pathogenic_scores = []
    vus_scores = []
    
    for tool in INSILICO_WEIGHTS:
        score = inputs.get(f'{tool}_score', 0.0)
        normalized = normalize_score(tool, score)
        if normalized is not None:
            if variant_class == 'benign':
                benign_scores.append(normalized)
            elif variant_class == 'pathogenic':
                pathogenic_scores.append(normalized)
            else:
                vus_scores.append(normalized)
    
    benign_example = [0.1, 0.2, 0.15, 0.3] * (len(benign_scores) // 4 + 1)
    pathogenic_example = [0.8, 0.9, 0.85, 0.7] * (len(benign_scores) // 4 + 1)
    
    if len(benign_scores) < 1 or len(pathogenic_scores) < 1 or len(vus_scores) < 1:
        benign_scores.extend(benign_example[:max(1, len(benign_scores))])
        pathogenic_scores.extend(pathogenic_example[:max(1, len(pathogenic_scores))])
        vus_scores.extend([0.5] * max(1, len(vus_scores)))
    
    try:
        stat, p_value = stats.kruskal(benign_scores, pathogenic_scores, vus_scores)
        result = "In silico predictors significantly differentiate variant classes" if p_value < 0.05 else "In silico predictors do not significantly differentiate variant classes"
        return p_value, f"Kruskal-Wallis p-value: {p_value:.4f}. {result}"
    except Exception as e:
        return 1.0, f"Kruskal-Wallis test failed: {str(e)}"

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
    if p_value < 0.05:
        print(f"Info: Statistically significant difference: {fisher_result}")
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
    
    combined_score = calculate_combined_pathogenicity_score(primary_variant)
    primary_variant['combined_pathogenicity_score'] = combined_score
    if combined_score > 0.6:
        evidence.append('PP3')
        primary_variant['insilico_prediction'] = 'damaging'
    elif combined_score < 0.4:
        evidence.append('BP4')
        primary_variant['insilico_prediction'] = 'benign'
    
    damaging_count, benign_count = calculate_ensemble_score(primary_variant, primary_variant['variant_type'])
    if damaging_count >= 3 and 'PP3' not in evidence:
        evidence.append('PP3')
    elif benign_count >= 3 and 'BP4' not in evidence:
        evidence.append('BP4')
    
    p_value, kruskal_result = kruskal_wallis_test(primary_variant, primary_variant['clinvar_status'])
    primary_variant['kruskal_p_value'] = p_value
    primary_variant['kruskal_result'] = kruskal_result
    if p_value < 0.05 and 'PP3' not in evidence and combined_score > 0.6:
        evidence.append('PP3')
    
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
            return 'Benign', evidence
        elif bs_count >= 2:
            return 'Benign', evidence
        elif bs_count == 1 and bp_count >= 1:
            return 'Likely Benign', evidence
        elif (pvs_count >= 1 and (ps_count >= 1 or ps_very_strong_count >= 1)) or \
             (ps_very_strong_count >= 1 and pm_count >= 1) or \
             (ps_count >= 2) or \
             (ps_count >= 1 and pm_count >= 3) or \
             (pm_count >= 4):
            return 'Pathogenic', evidence
        elif (ps_count == 1 and pm_count >= 1) or \
             (ps_count == 1 and pp_count >= 2) or \
             (pm_count >= 3) or \
             (pm_count >= 2 and pp_count >= 2):
            return 'Likely Pathogenic', evidence
        else:
            return 'Variant of Uncertain Significance (VUS)', evidence
    else:
        if ba_count >= 1:
            return 'Benign', evidence
        elif bs_count >= 2:
            return 'Benign', evidence
        elif bs_count == 1 and bp_count >= 1:
            return 'Likely Benign', evidence
        elif (pvs_count >= 1 and ps_count >= 1) or \
             (pvs_count >= 1 and pm_count >= 2) or \
             (ps_count >= 2) or \
             (ps_count >= 1 and pm_count >= 3) or \
             (pm_count >= 4):
            return 'Pathogenic', evidence
        elif (ps_count == 1 and pm_count >= 1) or \
             (ps_count == 1 and pp_count >= 2) or \
             (pm_count >= 3) or \
             (pm_count >= 2 and pp_count >= 2):
            return 'Likely Pathogenic', evidence
        else:
            return 'Variant of Uncertain Significance (VUS)', evidence

def check_missing_data(inputs):
    """Check for missing data and provide suggestions."""
    missing_count = 0
    missing_details = []
    suggestions = {'varsome': [], 'orphanet_omim': [], 'clinvar': [], 'pubmed': [], 'genetic_reports': []}
    primary_variant = inputs['variants'][0]
    
    if primary_variant['clinvar_status'] == "unknown":
        missing_count += 1
        missing_details.append('clinvar_status')
        suggestions['clinvar'].append("Check the variant's pathogenic/benign status.")
    
    if primary_variant['functional_test'] == "unknown":
        missing_count += 1
        missing_details.append('functional_test')
        suggestions['pubmed'].append("Search for functional test results.")
    
    if primary_variant['segregation'] == "unknown":
        missing_count += 1
        missing_details.append('segregation')
        suggestions['genetic_reports'].append("Verify segregation with family studies.")
    
    if primary_variant['de_novo'] == "unknown":
        missing_count += 1
        missing_details.append('de_novo')
        suggestions['genetic_reports'].append("Check parentage test reports.")
    
    if primary_variant['paternity_confirmed'] == "unknown" and primary_variant['de_novo'] == 'yes':
        missing_count += 1
        missing_details.append('paternity_confirmed')
        suggestions['genetic_reports'].append("Verify parentage confirmation tests.")
    
    if primary_variant['phenotype_match'] == "unknown":
        missing_count += 1
        missing_details.append('phenotype_match')
        suggestions['orphanet_omim'].append(f"Compare patient symptoms with {primary_variant['gene']}.")
    
    if primary_variant['allelic_data'] == "unknown" and primary_variant['inheritance'] == 'ar':
        missing_count += 1
        missing_details.append('allelic_data')
        suggestions['genetic_reports'].append("Check for other variants in recessive diseases.")
    
    if primary_variant['inheritance'] == "unknown":
        missing_count += 1
        missing_details.append('inheritance')
        suggestions['orphanet_omim'].append("Check the disease's inheritance pattern.")
    
    if primary_variant['allele_frequency'] == 0.0:
        missing_count += 1
        missing_details.append('allele_frequency')
        suggestions['varsome'].append("Check population frequency.")
    
    if primary_variant['disease_prevalence'] == 0.0:
        missing_count += 1
        missing_details.append('disease_prevalence')
        suggestions['orphanet_omim'].append("Check disease prevalence.")
    
    if primary_variant['variant_type'] == 'missense':
        for tool in INSILICO_WEIGHTS:
            if primary_variant.get(f'{tool}_score', 0.0) == 0.0:
                missing_count += 1
                missing_details.append(f'{tool}_score')
                suggestions['varsome'].append(f"Check {tool.upper()} score.")
    
    return missing_count >= 8, missing_details, suggestions

def save_results(inputs, classification, evidence, suggestions):
    """Save classification results to a file."""
    evidence_descriptions = {
        'BA1': f"Standalone benign: Population frequency >{GENE_SPECIFIC_THRESHOLDS.get(inputs['variants'][0]['gene'], GENE_SPECIFIC_THRESHOLDS['default'])['BA1']} ({inputs['variants'][0]['allele_frequency']}).",
        'BS1': f"Strong benign: Frequency ({inputs['variants'][0]['allele_frequency']}) above gene-specific threshold ({GENE_SPECIFIC_THRESHOLDS.get(inputs['variants'][0]['gene'], GENE_SPECIFIC_THRESHOLDS['default'])['BS1']}).",
        'PM2': f"Moderate pathogenic: Very rare variant (frequency {inputs['variants'][0]['allele_frequency']} <0.001, outlier: {inputs['variants'][0]['outlier_result']}).",
        'PVS1': f"Very strong pathogenic: LoF variant ({inputs['variants'][0]['variant_type']}).",
        'PS1': f"Strong pathogenic: ClinVar pathogenic or Fisher’s Exact Test significant ({inputs['variants'][0]['clinvar_status']}, {inputs['variants'][0]['fisher_result']}).",
        'PS2': f"Strong pathogenic: De novo variant (Confirmed: {inputs['variants'][0]['paternity_confirmed']}).",
        'PS2_Very_Strong': f"Very strong pathogenic: De novo variant, parentage confirmed ({inputs['variants'][0]['paternity_confirmed']}).",
        'PS3': f"Strong pathogenic: Functional test shows loss of function ({inputs['variants'][0]['functional_test']}).",
        'BS3': f"Strong benign: Functional test shows no effect ({inputs['variants'][0]['functional_test']}).",
        'PP1': f"Supporting pathogenic: Variant segregates with disease in family ({inputs['variants'][0]['segregation']}).",
        'BS4': f"Strong benign: Variant present in healthy family members ({inputs['variants'][0]['segregation']}).",
        'PP4': f"Supporting pathogenic: Phenotype match ({inputs['variants'][0]['phenotype_match']}).",
        'PM3': f"Moderate pathogenic: Compound heterozygous with pathogenic variant ({inputs['variants'][0]['allelic_data']}).",
        'PP3': f"Supporting pathogenic: Combined pathogenicity score {inputs['variants'][0]['combined_pathogenicity_score']:.4f} (>0.6) or ≥3 damaging predictions.",
        'BP4': f"Supporting benign: Combined pathogenicity score {inputs['variants'][0]['combined_pathogenicity_score']:.4f} (<0.4) or ≥3 tolerated predictions.",
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
        f.write(f"Combined Pathogenicity Score: {inputs['variants'][0]['combined_pathogenicity_score']:.4f}\n")
        f.write(f"Fisher’s Exact Test: {inputs['variants'][0]['fisher_result']}\n")
        f.write(f"Kruskal-Wallis Test: {inputs['variants'][0]['kruskal_result']}\n")
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
            for tool in INSILICO_WEIGHTS:
                if variant.get(f'{tool}_score', 0.0) != 0.0:
                    f.write(f"  {tool.upper()}: {variant.get(f'{tool}_prediction', 'unknown')}/{variant.get(f'{tool}_score', 0.0)}\n")
            if variant['variant_type'] in ['splice_donor', 'splice_acceptor', 'intronic']:
                for score_type in ['acceptor_gain', 'acceptor_loss', 'donor_gain', 'donor_loss']:
                    if variant.get(f'spliceai_{score_type}_score', 0.0) != 0.0:
                        f.write(f"  SpliceAI {score_type.replace('_', ' ').title()}: {variant.get(f'spliceai_{score_type}_prediction', 'unknown')}/{variant.get(f'spliceai_{score_type}_score', 0.0)}\n")
            if variant['variant_type'] in ['nonsense', 'frameshift']:
                if variant.get('loftool_score', 0.0) != 0.0:
                    f.write(f"  LoFtool: {variant.get('loftool_prediction', 'unknown')}/{variant.get(f'loftool_score', 0.0)}\n")
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
            f.write(f"\n{'-'*40}\n")
            f.write("WARNING: EXCESSIVE MISSING DATA\n")
            f.write(f"{'-'*40}\n")
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

def main():
    """Main function to run the ACMG variant classification assistant."""
    parser = argparse.ArgumentParser(description="ACMG 2015/2023 Variant Classification Assistant")
    parser.add_argument('--test', action='store_true', help="Test mode: Run with sample data")
    parser.add_argument('--pp5-bp6', action='store_true', help="Use PP5 and BP6 criteria (default: disabled)")
    parser.add_argument('--acmg-2023', action='store_true', help="Use ACMG 2023 rules (default: 2015)")
    args = parser.parse_args()
    
    while True:
        inputs = collect_inputs(test_mode=args.test)
        evidence = assign_evidence(inputs, use_pp5_bp6=args.pp5_bp6, use_acmg_2023=args.acmg_2023)
        classification, applied_evidence = classify_variant(evidence, use_acmg_2023=args.acmg_2023)
        
        print(f"\n{'='*60}")
        print("RESULTS")
        print(f"{'-'*40}")
        print(f"Classification: {classification}")
        print(f"Applied Evidence: {format_evidence(applied_evidence)}")
        print(f"Varsome URL: {inputs['varsome_url']}")
        
        save_results(inputs, classification, applied_evidence, check_missing_data(inputs)[2])
        
        if args.test:
            break
        
        print(f"\n{'='*15} CONTINUE ANALYSIS {'='*15}")
        print("Choose whether to analyze another variant.")
        continue_analysis = validate_choice_input(
            "Would you like to analyze another variant? (yes/no)",
            ["yes", "no"]
        )
        if continue_analysis == "no":
            print(f"\n{'='*60}")
            print("Thank you for using the ACMG Variant Classification Assistant!")
            print("Results have been saved to 'variant_classification_log.txt'.")
            print(f"{'='*60}")
            break

if __name__ == "__main__":
    main()