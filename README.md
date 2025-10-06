# 🧬 ACMG Variant Classification Assistant

**Version 3.4.0**

## **Quick Start**

> **📥 [Download Latest Version (v3.4.0) from Google Drive](https://drive.google.com/file/d/10pU9QmJkHky23B39JefU8RUsdNSmOAJX/view?usp=sharing)**  
> 
> **Ready-to-use standalone executable - No Python installation required!**  
> 1. Download and extract the zip file  
> 2. Run `ACMG_Assistant.exe`  
> 3. Start classifying variants immediately  

---

![v3 3 0 screenshot](https://github.com/user-attachments/assets/dad9b471-65cd-4ac9-b5da-e5da3ddb9c21)
Screenshot of the CLI welcome page

## ⚙️ Key Features

- **Complete ACMG/AMP Implementation**: All 28 criteria (2015 & 2023 guidelines)
- **Statistical Framework**: Fisher's Exact Test, LOD scoring, confidence tracking
- **30+ In Silico Predictors**: REVEL, CADD, AlphaMissense, VEST4, ESM1b, SpliceAI, MetaSVM...
- **Gene-Specific Thresholds**: Custom population frequencies for high-variability genes
- **Interactive & Automated**: Mix of data-driven automation and literature-guided evaluation

## 💻 Installation & Usage

### Standalone Executable (Recommended)
```bash
# Download and extract ACMG_Assistant.zip
# Run: ACMG_Assistant.exe
```

### Python Installation
```bash
git clone https://github.com/Bilmem2/ACMG_Assistant
cd src
pip install -r ../requirements.txt
python acmg_assistant.py
```

### Command Options
```bash
# Normal mode
ACMG_Assistant.exe             # Executable
python acmg_assistant.py               # Python (from src/ directory)

# ACMG 2023 guidelines
ACMG_Assistant.exe --acmg-2023
python acmg_assistant.py --acmg-2023

# Test mode (Python only)
python acmg_assistant.py --test
```

## 🔧 Building Executable

```bash
pip install -r requirements_build.txt
python build_executable_new.py
```

## 📊 Latest Release & Version History

For version history, release notes, and previous versions, visit:
**[GitHub Releases](https://github.com/Bilmem2/ACMG_Assistant/releases)**

Each release includes:
- Source code archives

## 📊 In Silico Predictors

This algorithm integrates **30+ computational prediction tools** across multiple categories for comprehensive variant pathogenicity assessment:

### 🎯 **Primary Metascores & Ensemble Methods**
- **REVEL** - Rare Exome Variant Ensemble Learner
- **CADD** - Combined Annotation Dependent Depletion  
- **AlphaMissense** - DeepMind's protein structure-based predictor
- **MetaRNN** - Recurrent neural network metapredictor
- **ClinPred** - Clinical significance predictor
- **BayesDel** - Bayesian deleteriousness score
- **MetaSVM/MetaLR** - SVM/Logistic regression ensemble methods

### 🧬 **Missense Variant Predictors**
- **SIFT** - Sorting Intolerant From Tolerant
- **PolyPhen-2** - Polymorphism Phenotyping v2 (HDiv/HVar)
- **PROVEAN** - Protein Variation Effect Analyzer
- **VEST4** - Variant Effect Scoring Tool v4
- **ESM1b** - Evolutionary Scale Modeling (protein language model)
- **MutationTaster** - Disease-causing potential predictor
- **FATHMM** - Functional Analysis through Hidden Markov Models
- **MutationAssessor** - Functional impact assessment
- **MutPred** - Pathogenicity prediction with structural features
- **LRT** - Likelihood Ratio Test

### 🧮 **Conservation & Evolutionary Analysis**
- **PhyloP** - Multiple alignments (100/30/17-way vertebrates/mammals/primates)
- **phastCons** - Phylogenetic conservation (100/30/17-way)
- **GERP++** - Genomic Evolutionary Rate Profiling
- **SiPhy** - Site-specific phylogenetic analysis

### ✂️ **Splice Site Prediction**
- **SpliceAI** - Deep learning splice predictor (AG/AL/DG/DL scores)
- **Ada** - Adaptive boosting splice predictor
- **RF** - Random Forest splice predictor

### 🔬 **Functional & Regulatory Predictors**
- **FITCONS** - Functional information content (multiple cell types)
- **Combined Metascore** - Custom VAMPP-like weighted combination

### 📚 **Data Sources & Integration**
- **Varsome** - Comprehensive variant annotation platform
- **ClinVar** - NCBI clinical significance database
- **dbNSFP** - Database of human non-synonymous SNPs
- **Manual Entry** - Custom predictor score input interface

> **Note**: All predictor scores require manual entry - no automatic database retrieval is performed to ensure data accuracy and user control.


## 🏗️ Project Structure

```
src/
├── acmg_assistant.py          # Main application entry point  
├── config/
│   ├── __init__.py
│   └── constants.py           # ACMG criteria thresholds, predictor configs (v3.3.0)
├── core/
│   ├── __init__.py
│   ├── acmg_classifier.py            # Main classification engine
│   ├── evidence_evaluator.py         # Evidence scoring with confidence tracking (v3.3.0)
│   ├── variant_data.py               # Variant data structures
│   ├── functional_studies_evaluator.py  # PS3/BS3 functional studies (v3.3.0)
│   ├── gene_rules_engine.py          # Gene-specific rule engine (v3.3.0)
│   ├── gene_specific_rules.py        # BRCA1/2, TTN, MUC16 rules (v3.3.0)
│   ├── missense_evaluator.py         # PS1/PM5/PP2/BP1 missense analysis (v3.3.0)
│   ├── phenotype_matcher.py          # PP4/BP5 phenotype matching (v3.3.0)
│   └── population_analyzer.py        # BA1/BS1/PM2/BS2 population analysis (v3.3.0)
└── utils/
    ├── __init__.py
    ├── api_client.py          # ClinVar/Ensembl API integrations
    ├── input_handler.py       # User input processing
    ├── report_generator.py    # Classification report output
    ├── statistical_utils.py   # Fisher's exact, LOD scoring (v3.3.0)
    └── validators.py          # Input validation functions
```

## 📈 ACMG Criteria Implementation & Algorithm

### Core Classification Logic
The algorithm implements a **multi-layered evidence evaluation system** that processes variants through:

1. **Variant Type Detection**: Automatically categorizes variants (missense, nonsense, frameshift, splice, etc.)
2. **Evidence Collection**: Gathers population, functional, computational, and segregation data
3. **Criterion-Specific Evaluation**: Applies ACMG/AMP criteria with variant-type-specific logic
4. **Statistical Integration**: Combines evidence using Fisher's Exact Test and prevalence thresholds
5. **Final Classification**: Determines pathogenicity based on evidence strength combinations

### ACMG/AMP Criteria Support

**All 28 Criteria Fully Implemented**:

**Pathogenic Evidence**:
- **PVS1**: Null variants in haploinsufficient genes with splice prediction
- **PS1-4**: Functional studies, segregation, prevalence analysis with statistical validation  
- **PM1-6**: Domain analysis, gene-specific frequencies, computational predictions
- **PP1-5**: Segregation support, functional confirmation, computational evidence

**Benign Evidence**:
- **BA1**: High-frequency variants with population-specific thresholds
- **BS1-4**: Functional studies, segregation analysis, computational predictions
- **BP1-7**: Comprehensive computational and frequency-based evidence

**Interactive Criteria**: PS1, PS4, PM3, PM5, PP4, BP2, BP5, BP6 feature user-guided literature review
**Automated Criteria**: All others use data-driven logic with manual review options

### Special Algorithm Features
- **Confidence Tracking (v3.3.0)**: Every criterion reports confidence level (high/medium/low/very_low) and data source
- **Statistical Automation (v3.3.0)**: Automated Fisher's Exact Test (PS4) and LOD scoring (PP1/BS4)
- **Enhanced Validation (v3.3.0)**: Stricter PS2 2023 upgrade (both parents required), PP5/BP6 expert panel requirements
- **Gene-Specific PM2**: Custom population frequency thresholds per gene
- **Statistical PS4**: Fisher's Exact Test for case-control prevalence analysis
- **Conservation Analysis**: Multi-species phylogenetic conservation scoring
- **2023 Guidelines**: PP5/BP6 support and PS2_Very_Strong implementation

## ⚠️ Disclaimer & Important Notes

**⚠️ RESEARCH USE ONLY**: This tool is intended for research and educational purposes. Variant classifications should NOT be used for clinical decision-making without validation by qualified medical professionals and certified genetic counselors.

**Technical Notes**:
- Manual predictor score entry required (no automatic database retrieval)
- Internet required for API calls (ClinVar, Ensembl)
- Not validated for direct clinical use

##   Citation & References

**If you use this algorithm in your research or clinical pipeline, please cite:**

[https://doi.org/10.5281/zenodo.15831866](https://doi.org/10.5281/zenodo.15831866)

---

This tool uses a VAMPP-score-like metascore approach for in silico pathogenicity prediction, inspired by the original VAMPP-score framework. If you use this tool, the VAMPP-score, or any component of their statistical framework in your work, please cite the original VAMPP-score publication:

> Eylul Aydin, Berk Ergun, Ozlem Akgun-Dogan, Yasemin Alanay, Ozden Hatirnaz Ng, Ozkan Ozdemir. "A New Era in Missense Variant Analysis: Statistical Insights and the Introduction of VAMPP-Score for Pathogenicity Assessment." *bioRxiv* (2024.07.11.602867). [DOI: 10.1101/2024.07.11.602867](https://doi.org/10.1101/2024.07.11.602867)

For in silico and molecular analysis methodology, see:

> Ozdemir, O., Bychkovsky, B. L., Unal, B., Onder, G., Amanvermez, U., Aydin, E., Ergun, B., Sahin, I., Gokbayrak, M., Ugurtas, C., Koroglu, M. N., Cakir, B., Kalay, I., Cine, N., Ozbek, U., Rana, H. Q., Hatirnaz Ng, O., & Agaoglu, N. B. (2024). "Molecular and In Silico Analysis of the CHEK2 Gene in Individuals with High Risk of Cancer Predisposition from Türkiye." *Cancers*, 16(22), 3876. [DOI: 10.3390/cancers16223876](https://doi.org/10.3390/cancers16223876)

**VAMPP-score** is a registered framework. For more information and access to the original VAMPP-score data and scripts, see: [vamppscore.com](https://vamppscore.com/) and the [VAMPP-score-data Google Drive](https://drive.google.com/drive/folders/1emkHcTlxgjH6G-2Yl4wQQnKi5Wsip4IY?usp=drive_link).

## 👤 Contact

- **Author**: Can Sevilmiş
- **Email**: cansevilmiss@gmail.com
- **LinkedIn**: [cansevilmiss](https://linkedin.com/in/cansevilmiss)
