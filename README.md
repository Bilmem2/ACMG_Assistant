# 🧬 ACMG Variant Classification Assistant v4.0.0

> **A research-oriented, transparent, and extensible variant interpretation pipeline implementing ACMG/AMP 2015 & 2023 guidelines.**

## 📥 Quick Start

> **[Download Latest Windows Executable from Google Drive](https://drive.google.com/file/d/1Xmtoix2V55DXWxRykt7X2v8UNJ7m7xMZ/view?usp=sharing)**
>
> Ready-to-use standalone `.exe` — No Python installation required!

---

## 📋 Table of Contents

1. [Overview](#overview)
2. [Philosophy](#philosophy-critical)
3. [Data Flow Architecture](#data-flow-architecture)
4. [Caching Model](#caching-model)
5. [Evidence Types Supported](#evidence-types-supported)
6. [Limitations & Safety Notes](#limitations--safety-notes)
7. [Versioning & Changelog](#versioning--changelog)
8. [Installation & Usage](#installation--usage)
9. [Binary Releases](#binary-releases--executable-distribution)
10. [Project Structure](#project-structure)
11. [Citation & References](#citation--references)
12. [Contact](#contact)

---

## Overview

**ACMG Assistant** is a comprehensive variant classification tool that implements the ACMG/AMP 2015 and 2023 guidelines for interpreting sequence variants. It provides:

- **Automatic evaluation** of population frequency criteria (BA1, BS1, PM2)
- **Automatic evaluation** of computational/in-silico criteria (PP3, BP4)
- **Automatic evaluation** of functional domain criteria (PM1)
- **Automatic evaluation** of phenotype matching (PP4, BP5)
- **Interactive evaluation** of literature-based criteria (PS3/BS3, PS4, PP1/BS4, PS1/PM5, PP5/BP6)

### Purpose

This tool is designed for:
- **Educational use**: Understanding ACMG classification logic
- **Research pipelines**: Reproducible, transparent variant interpretation
- **Workflow augmentation**: Pre-screening variants before expert review

> ⚠️ **NOT for direct clinical decision-making** — All results require validation by qualified professionals.

---

## Philosophy (CRITICAL)

### 🚫 No Hardcoded Biological Truth

The core design principle of ACMG Assistant is:

> **"Local code contains ONLY interpretive logic — all factual biological data must come from external sources."**

This means:
- ✅ **Thresholds, weights, scoring formulas** → Defined locally
- ✅ **ACMG evidence combination rules** → Defined locally
- ❌ **Predictor scores (REVEL, CADD, etc.)** → Must be fetched from APIs
- ❌ **Population allele frequencies** → Must be fetched from APIs
- ❌ **Functional domains, hotspots** → Must be fetched from APIs
- ❌ **Gene-specific rules** → Must be fetched from APIs or entered by user

### Data Sources (in priority order)

All factual variant-level data must come from one of:

1. **External APIs** (gnomAD, ClinVar, UniProt, CancerHotspots, myvariant.info, etc.)
2. **User input** (interactive evidence collection for literature-based criteria)
3. **Validated cache** (previously fetched and validated API responses)

The local codebase **NEVER fabricates** biological values — it only interprets them.

### Interactive Evidence for Literature-Based Criteria

Criteria that require literature review are handled through **structured interactive prompts**:

| Criterion | What User Provides |
|-----------|-------------------|
| PS3 / BS3 | Functional study details, assay type, quality level |
| PS4 | Case-control counts, odds ratio data |
| PP1 / BS4 | Segregation data, LOD scores, family structure |
| PS1 / PM5 | Prior variant at same codon, ClinVar status |
| PP5 / BP6 | External lab assertions, submission quality |

### Phenotype Matching

The `PhenotypeMatcher` uses:
- Local HPO ontology for term similarity
- Reproducible Jaccard/IC-based similarity scores
- **NOT a replacement for clinical phenotyping**

---

## Data Flow Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              CLI / Entry Point                               │
│                            (acmg_assistant.py)                               │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            EvidenceEvaluator                                 │
│                    (Central orchestration engine)                            │
│                                                                              │
│  ┌──────────────────────────────────────────────────────────────────────┐   │
│  │                    _fetch_external_data()                             │   │
│  │         "Fetch once, interpret many" — pre-loads all data            │   │
│  └──────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘
          │                    │                    │                    │
          ▼                    ▼                    ▼                    ▼
┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐
│ PredictorAPI    │  │ PopulationAPI   │  │ GeneSpecific    │  │ Interactive     │
│ Client          │  │ Client          │  │ Rules           │  │ Evidence        │
│                 │  │                 │  │                 │  │ Collector       │
│ • myvariant.info│  │ • gnomAD GraphQL│  │ • CancerHotspots│  │                 │
│ • AlphaMissense │  │ • ExAC REST     │  │ • UniProt       │  │ • PS3/BS3       │
│ • CADD API      │  │ • TOPMed        │  │ • ClinGen       │  │ • PS4           │
│ • VEP           │  │                 │  │                 │  │ • PP1/BS4       │
│                 │  │                 │  │                 │  │ • PS1/PM5       │
│ Multi-source    │  │ Multi-source    │  │ PM1 external    │  │ • PP5/BP6       │
│ priority-based  │  │ gnomAD-first    │  │ hotspots/domains│  │                 │
└────────┬────────┘  └────────┬────────┘  └────────┬────────┘  └────────┬────────┘
         │                    │                    │                    │
         └───────────────┬────┴────────────────────┴────────────────────┘
                         │
                         ▼
          ┌─────────────────────────────────────┐
          │          ResultCache                 │
          │   (Strict validation layer)          │
          │                                      │
          │  • Validates all entries before use  │
          │  • Rejects invalid/corrupted data    │
          │  • TTL-based expiration              │
          │  • Thread-safe file storage          │
          │                                      │
          │  "Cache is optimization, NOT truth"  │
          └─────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                           Pure Interpretation Layer                          │
│                                                                              │
│  ┌────────────────┐  ┌────────────────┐  ┌────────────────┐                 │
│  │ MissenseEval   │  │ PopulationAnal │  │ PhenotypeMatcher│                 │
│  │                │  │                │  │                │                 │
│  │ Composite score│  │ AF thresholds  │  │ HPO similarity │                 │
│  │ → PP3 / BP4    │  │ → BA1/BS1/PM2  │  │ → PP4 / BP5    │                 │
│  └────────────────┘  └────────────────┘  └────────────────┘                 │
│                                                                              │
│  "Data missing → safely degrade evidence strength, never fabricate"         │
└─────────────────────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            ACMGClassifier                                    │
│                                                                              │
│  Merges all evidence (automatic + interactive) → Final classification       │
│  Pathogenic / Likely Pathogenic / VUS / Likely Benign / Benign              │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Key Principles

1. **Fetch Once, Interpret Many**: External data is fetched at the start; evaluators are pure interpreters
2. **Cache is Optimization, Not Truth**: Invalid cache entries are rejected, not trusted
3. **Graceful Degradation**: Missing data → reduced evidence strength, never fabricated values

---

## Caching Model

ACMG Assistant uses a **strict, validated caching layer** to minimize redundant API calls while ensuring data integrity.

### CacheKey Structure

```python
CacheKey(
    category='predictor' | 'population',
    source='dbNSFP' | 'gnomAD_GraphQL' | ...,
    variant_id='GRCh38:17-7674234-G-A',
    version='v4.0'  # optional
)
```

### TTL (Time-To-Live) Rules

| Data Type | Default TTL | Rationale |
|-----------|-------------|-----------|
| Predictor scores | 7 days | Scores rarely change |
| Population data | 30 days | gnomAD updates infrequently |

### Validation Before Trust

All cached entries are validated before use:

- **Predictor scores**: Must be within valid ranges (REVEL ∈ [0,1], CADD ∈ [0,60], etc.)
- **Population stats**: AF ∈ [0,1], AC ≤ AN, no negative values
- **Invalid entries**: Automatically invalidated and re-fetched

### Corruption Handling

- JSON decode errors → Cache miss (file removed)
- Hash mismatch → Cache miss (file removed)
- Expired entries → Cache miss (file removed)

### Thread Safety

`ResultCache` uses `threading.RLock` for safe concurrent access.

### Cache Location

Default: `src/api_cache/` (organized by category/source)

---

## Evidence Types Supported

| Evidence | Source | Type | Details |
|----------|--------|------|---------|
| **BA1 / BS1 / PM2** | Population data | ✅ Automatic | gnomAD, ExAC allele frequencies |
| **PP3 / BP4** | In-silico predictors | ✅ Automatic | Multi-source: REVEL, CADD, AlphaMissense, SIFT, PolyPhen2, etc. |
| **PM1** | Functional domains | ✅ Automatic | CancerHotspots API + UniProt domains |
| **PP4 / BP5** | Phenotype matching | ✅ Automatic | HPO ontology similarity scoring |
| **PS3 / BS3** | Functional studies | 🔄 Interactive | User enters assay details, quality |
| **PS4** | Case-control data | 🔄 Interactive | User enters case/control counts |
| **PP1 / BS4** | Segregation | 🔄 Interactive | User enters pedigree, LOD scores |
| **PS1 / PM5** | Prior variants | 🔄 Interactive | User enters codon-based prior data |
| **PP5 / BP6** | External assertions | 🔄 Interactive | User enters lab submissions |
| **PVS1** | Null variants | ✅ Automatic | LOF in haploinsufficient genes |
| **PS2 / PM6** | De novo | 🔄 Interactive | User confirms parental testing |
| **PM3 / BP2** | In-trans/cis | 🔄 Interactive | User enters phase data |

---

## Limitations & Safety Notes

### ⚠️ Critical Disclaimers

1. **NOT a Clinical Decision-Making System**
   - All evidence must be reviewed by a qualified clinical geneticist
   - Classifications are suggestions, not diagnoses

2. **API Dependency**
   - Network failures or API downtime may limit automatic criteria
   - Missing data results in reduced evidence, never forced calls

3. **Interactive Evidence Accuracy**
   - PS3, PS4, PP1, etc. rely on truthful user input
   - Garbage in → garbage out

4. **Missense Composite Score**
   - The PP3/BP4 composite score is a **research approximation**
   - It is NOT the validated VAMPP score — only inspired by similar methodology
   - Clinical labs should use their own validated thresholds

5. **Phenotype Matching**
   - HPO similarity is approximate and algorithmic
   - NOT a replacement for clinical phenotyping by experts
   - Low sensitivity for rare/novel phenotypes

6. **Cache Validity**
   - Cache is an optimization layer, not a source of truth
   - Stale cache may return outdated scores if TTL not managed

### What This Tool Does NOT Do

- ❌ Replace expert clinical judgment
- ❌ Guarantee 100% accuracy
- ❌ Provide legally binding classifications
- ❌ Automatically retrieve all possible data sources
- ❌ Handle structural variants (current focus: SNVs, indels)

---

## Versioning & Changelog

### Current Version: **4.0.0**

Released: December 2025

### Major Changes in v4.0.0

| Feature | Description |
|---------|-------------|
| **Multi-source predictor system** | Fetches from myvariant.info, AlphaMissense API, CADD API with source priority |
| **Multi-source population AF** | gnomAD GraphQL (primary), ExAC, TOPMed fallbacks |
| **Strict validated caching** | `ResultCache` with validation, TTL, thread safety |
| **PM1 via external hotspots** | CancerHotspots API + UniProt functional domains |
| **Phenotype matcher overhaul** | HPO-based similarity with IC weighting |
| **Interactive evidence subsystem** | Structured prompts for PS3/BS3, PS4, PP1/BS4, PS1/PM5, PP5/BP6 |
| **198 tests passing** | Comprehensive test coverage |

### Previous Versions

- **v3.5.0**: Gene-specific rules, enhanced API integration
- **v3.3.0**: Statistical framework (Fisher's exact, LOD scoring)
- **v3.0.0**: Initial ACMG 2023 support
- **v2.x**: Core ACMG 2015 implementation

---

## Installation & Usage

### Python Installation (Recommended for Development)

```bash
# Clone the repository
git clone https://github.com/Bilmem2/ACMG_Assistant
cd ACMG_Assistant

# Create virtual environment (optional but recommended)
python -m venv venv
venv\Scripts\activate  # Windows
# source venv/bin/activate  # Linux/Mac

# Install dependencies
pip install -r requirements.txt

# Run from src directory
cd src
python acmg_assistant.py
```

### Command Options

```bash
# Standard mode (ACMG 2015)
python acmg_assistant.py

# ACMG 2023 mode
python acmg_assistant.py --acmg-2023

# Test mode (mock data, no API calls)
python acmg_assistant.py --test

# Show version
python acmg_assistant.py --version
```

### API Dependencies

The tool requires internet access for:
- **gnomAD** (population frequencies)
- **myvariant.info** (predictor scores via dbNSFP)
- **CancerHotspots** (PM1 hotspot detection)
- **UniProt** (functional domains)
- **ClinVar** (external assertions)

Offline mode uses cached data only.

---

## Binary Releases / Executable Distribution

### Windows Standalone Executable

ACMG Assistant is distributed as a standalone Windows `.exe` for users without Python:

> **[Download from Google Drive](https://drive.google.com/file/d/1Xmtoix2V55DXWxRykt7X2v8UNJ7m7xMZ/view?usp=sharing)**

### Usage

1. Download and extract the ZIP file
2. Run `ACMG_Assistant.exe`
3. Follow the interactive prompts

### Notes on Executable

- **First startup may be slow** (PyInstaller unpacking)
- **Internet required** for API calls
- **Cache persists** between runs in the same directory
- **No Python installation needed**

### Rebuilding the Executable

```bash
# Install build dependencies
pip install -r requirements_build.txt

# Build executable
python build_executable_new.py

# Output: dist/ACMG_Assistant.exe
```

---

## Project Structure

```
ACMG_Assistant/
├── src/
│   ├── acmg_assistant.py              # Main CLI entry point
│   ├── config/
│   │   ├── __init__.py
│   │   ├── constants.py               # Thresholds, API settings
│   │   ├── predictors.py              # PredictorScore, PopulationStats dataclasses
│   │   └── version.py                 # Version metadata
│   ├── core/
│   │   ├── __init__.py
│   │   ├── acmg_classifier.py         # Final classification engine
│   │   ├── evidence_evaluator.py      # Central orchestration
│   │   ├── variant_data.py            # VariantData dataclass
│   │   ├── missense_evaluator.py      # PP3/BP4 composite scoring
│   │   ├── population_analyzer.py     # BA1/BS1/PM2 evaluation
│   │   ├── gene_specific_rules.py     # Gene-specific PM1, thresholds
│   │   ├── phenotype_matcher.py       # PP4/BP5 HPO matching
│   │   └── functional_studies_evaluator.py  # PS3/BS3 evaluation
│   └── utils/
│       ├── __init__.py
│       ├── api_client.py              # ClinVar, Ensembl clients
│       ├── predictor_api_client.py    # Multi-source predictor/population
│       ├── cache.py                   # ResultCache with validation
│       ├── input_handler.py           # Interactive evidence collection
│       ├── report_generator.py        # Report output
│       └── validators.py              # Input validation
├── tests/
│   ├── test_acmg_classifier.py        # 20 tests
│   ├── test_gene_specific_pm1.py      # 31 tests
│   ├── test_interactive_evidence.py  # 56 tests
│   ├── test_predictor_population_api.py  # 40 tests
│   └── test_cache_and_validation.py  # 51 tests
├── data/
│   ├── gene_rules/                    # Gene-specific configuration
│   └── domain_annotations/            # Functional domain data
├── requirements.txt                   # Runtime dependencies
├── requirements_build.txt             # Build dependencies
├── pyproject.toml                     # Project configuration
└── README.md                          # This file
```

---

## Citation & References

### Citing ACMG Assistant

If you use this tool in your research, please cite:

> **ACMG Variant Classification Assistant**  
> https://doi.org/10.5281/zenodo.15831866

### Related Work

This tool uses a VAMPP-score-inspired metascore approach. If you use this methodology, please also cite:

> Eylul Aydin, Berk Ergun, et al. "A New Era in Missense Variant Analysis: Statistical Insights and the Introduction of VAMPP-Score for Pathogenicity Assessment." *bioRxiv* (2024). [DOI: 10.1101/2024.07.11.602867](https://doi.org/10.1101/2024.07.11.602867)

### ACMG/AMP Guidelines

> Richards S, et al. "Standards and guidelines for the interpretation of sequence variants." *Genet Med.* 2015;17(5):405-424.

> Plon SE, et al. "Sequence variant classification and reporting: recommendations for improving the interpretation of cancer susceptibility genetic test results." *Hum Mutat.* 2008.

---

## Contact

- **Author**: Can Sevilmiş
- **Email**: cansevilmiss@gmail.com
- **LinkedIn**: [cansevilmiss](https://linkedin.com/in/cansevilmiss)
- **GitHub**: [Bilmem2/ACMG_Assistant](https://github.com/Bilmem2/ACMG_Assistant)

---

<p align="center">
  <strong>ACMG Assistant v4.0.0</strong><br>
</p>
