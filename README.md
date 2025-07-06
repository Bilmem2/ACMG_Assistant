# 🧬 ACMG Variant Classification Assistant

**ACMG/AMP Variant Classification with Statistical Framework**  
*Last Updated: July 2025*

## 🚀 **Quick Start - Download Executable**

> **📥 [Download ACMG_Assistant.zip from Google Drive](https://drive.google.com/drive/folders/1emkHcTlxgjH6G-2Yl4wQQnKi5Wsip4IY?usp=drive_link)**  
> 
> **Ready-to-use standalone executable - No Python installation required!**  
> 1. Download and extract the zip file  
> 2. Run `ACMG_Assistant.exe`  
> 3. Start classifying variants immediately  

---

A comprehensive tool for classifying genetic variants according to ACMG/AMP 2015 and 2023 guidelines. Features VAMPP-score implementation, comprehensive in silico predictor integration, and evidence evaluation algorithms.

## ⚙️ Key Features

- **Complete ACMG/AMP Guidelines**: 2015 & 2023 standards with PP5/BP6 and PS2_Very_Strong
- **30+ In Silico Predictors**: REVEL, CADD, AlphaMissense, VEST4, PrimateAI, ESM1b, SpliceAI, MMSplice
- **VAMPP-Score Integration**: Metascore with weighted predictor combination
- **API Integration**: ClinVar and Ensembl with intelligent caching
- **Statistical Framework**: Fisher's Exact Test, prevalence-based thresholds, conservation analysis

## 💻 Installation & Usage

### Standalone Executable (Recommended)
```bash
# Download and extract ACMG_Assistant.zip
# Run: ACMG_Assistant.exe
```

### Python Installation
```bash
git clone https://github.com/Bilmem2/acmg-assessor.git
cd acmg-assessor
pip install -r requirements.txt
python acmg_assistant.py
```

### Command Options
```bash
# Normal mode
acmg_assistant.exe                    # Executable
python acmg_assistant.py              # Python

# ACMG 2023 guidelines
acmg_assistant.exe --acmg-2023
python acmg_assistant.py --acmg-2023

# Test mode (Python only)
python acmg_assistant.py --test
```

## 🔧 Building Executable

```bash
pip install -r requirements_build.txt
python build_executable.py
```

## 📊 In Silico Predictors

**Primary Metascores**: REVEL, CADD, AlphaMissense, MetaRNN, ClinPred, BayesDel  
**High Priority**: VEST4, PrimateAI, ESM1b, PROVEAN  
**Conservation**: PhyloP (100/30/17-way), GERP++  
**Splice**: SpliceAI, MMSplice, Ada, RF, dbscSNV  
**Individual**: SIFT, PolyPhen-2, MutationTaster, FATHMM

**Score Sources**: Varsome, ClinVar, dbNSFP (manual entry required)


## 🏗️ Project Structure

```
acmg_assessor/
├── acmg_assistant.py              # Main application entry point
├── build_executable.py            # PyInstaller build script
├── requirements.txt               # Python dependencies
├── requirements_build.txt         # Build-specific dependencies
├── config/
│   ├── __init__.py
│   └── constants.py              # ACMG criteria thresholds, predictor configs
├── core/
│   ├── __init__.py
│   ├── acmg_classifier.py        # Main classification engine
│   ├── evidence_evaluator.py     # Evidence scoring logic
│   └── variant_data.py           # Variant data structures
├── utils/
│   ├── __init__.py
│   ├── api_client.py             # ClinVar/Ensembl API integrations
│   ├── input_handler.py          # User input processing
│   ├── report_generator.py       # Classification report output
│   └── validators.py             # Input validation functions
└── tests/
    ├── test_*.py                 # Comprehensive test suites
    └── __pycache__/              # Python bytecode cache
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

**Pathogenic Evidence**:
- **PVS1**: Null variants in haploinsufficient genes with splice prediction
- **PS1-4**: Functional studies, segregation, prevalence analysis with statistical validation
- **PM1-6**: Domain analysis, gene-specific frequencies, computational predictions
- **PP1-5**: Segregation support, functional confirmation, computational evidence

**Benign Evidence**:
- **BA1**: High-frequency variants with population-specific thresholds
- **BS1-4**: Functional studies, segregation analysis, computational predictions
- **BP1-7**: Comprehensive computational and frequency-based evidence

### Special Algorithm Features
- **VAMPP-Score Integration**: Metascore combining 30+ predictors for PP3/BP4
- **Gene-Specific PM2**: Custom population frequency thresholds per gene
- **Statistical PS4**: Fisher's Exact Test for case-control prevalence analysis
- **Conservation Analysis**: Multi-species phylogenetic conservation scoring
- **2023 Guidelines**: PP5/BP6 support and PS2_Very_Strong implementation

## ⚠️ Disclaimer & Important Notes

**⚠️ MEDICAL DISCLAIMER**: This tool is intended for research and educational purposes only. It is **NOT** a substitute for professional medical advice, diagnosis, or treatment. Variant classifications provided by this tool should not be used for clinical decision-making without proper validation by qualified medical professionals and certified genetic counselors. Always consult with healthcare providers for clinical interpretation of genetic variants.

**Important Technical Notes**:
- **Manual Score Entry**: No automatic retrieval from databases - all predictor scores must be manually entered
- **Research Use Only**: Not validated for direct clinical use without additional confirmation
- **Internet Required**: For API calls (ClinVar, Ensembl) and database queries
- **Test Mode**: Available only in Python installation, not in standalone executable

##   Citation & References
This tool uses a VAMPP-score-like metascore approach for in silico pathogenicity prediction, inspired by the original VAMPP-score framework. If you use this tool, the VAMPP-score, or any component of their statistical framework in your work, please cite the original VAMPP-score publication:

> Eylul Aydin, Berk Ergun, Ozlem Akgun-Dogan, Yasemin Alanay, Ozden Hatirnaz Ng, Ozkan Ozdemir. "A New Era in Missense Variant Analysis: Statistical Insights and the Introduction of VAMPP-Score for Pathogenicity Assessment." *bioRxiv* (2024.07.11.602867). [DOI: 10.1101/2024.07.11.602867](https://doi.org/10.1101/2024.07.11.602867)

For in silico and molecular analysis methodology, see:

> Ozdemir, O., Bychkovsky, B. L., Unal, B., Onder, G., Amanvermez, U., Aydin, E., Ergun, B., Sahin, I., Gokbayrak, M., Ugurtas, C., Koroglu, M. N., Cakir, B., Kalay, I., Cine, N., Ozbek, U., Rana, H. Q., Hatirnaz Ng, O., & Agaoglu, N. B. (2024). "Molecular and In Silico Analysis of the CHEK2 Gene in Individuals with High Risk of Cancer Predisposition from Türkiye." *Cancers*, 16(22), 3876. [DOI: 10.3390/cancers16223876](https://doi.org/10.3390/cancers16223876)

**VAMPP-score** is a registered framework. For more information and access to the original VAMPP-score data and scripts, see: [vamppscore.com](https://vamppscore.com/) and the [VAMPP-score-data Google Drive](https://drive.google.com/drive/folders/1emkHcTlxgjH6G-2Yl4wQQnKi5Wsip4IY?usp=drive_link).

## 👤 Contact

- **Author**: Can Sevilmiş
- **LinkedIn**: [cansevilmiss](https://linkedin.com/in/cansevilmiss)
