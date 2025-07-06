# 🧬 ACMG Variant Classification Assistant

**Advanced ACMG/AMP Variant Classification with Enhanced Statistical Framework**  
*Last Updated: July 2025*

The ACMG Variant Classification Assistant is a comprehensive, state-of-the-art tool for classifying genetic variants according to ACMG/AMP 2015 and 2023 guidelines. Designed for clinical geneticists, bioinformaticians, researchers, and students, this tool integrates advanced statistical frameworks, population genetics, conservation analysis, and machine learning-based scoring systems to provide evidence-based, standardized variant classification. The tool features an enhanced VAMPP-score implementation, prevalence-based threshold adjustment, comprehensive in silico predictor integration, and sophisticated evidence evaluation algorithms.
3. **Double-click `ACMG_Assistant.exe`** or use the provided batch file
4. **Follow the interactive prompts** for variant analysis

**No Python installation required!** The executable includes everything needed.

**Note:** Test mode is only available through Python installation (Option 2).

### Option 2: Python Installation (For Developers)

### Prerequisites
- **Python 3.8 or higher**
- **Internet connection** (for API calls; optional for test mode)
- **Operating System**: Windows, macOS, or Linux

### Installation

```bash
# Clone the repository
git clone https://github.com/Bilmem2/acmg-assessor.git
cd acmg-assessor

# Install dependencies
pip install -r requirements.txt
```rk**  
*Last Updated: January 2025*

The ACMG Variant Classification Assistant is a comprehensive, state-of-the-art tool for classifying genetic variants according to ACMG/AMP 2015 and 2023 guidelines. Designed for clinical geneti## 📖 Citation & References

### ACMG Guidelines Referencesformaticians, researchers, and students, this tool integrates advanced statistical frameworks, population genetics, conservation analysis, and machine learning-based scoring systems to provide evidence-based, standardized variant classification. The tool features an enhanced VAMPP-score implementation, prevalence-based threshold adjustment, comprehensive in silico predictor integration, and sophisticated evidence evaluation algorithms.

**⚠️ Important:** All in silico scores must be entered manually by the user. This tool does not automatically retrieve scores from any database or local files. Use resources like Varsome, ClinVar, or dbNSFP for score collection.

**License:** This project is licensed under the MIT License. See LICENSE file for details.

## ⚙️ Key Features

### 🔬 **Advanced ACMG Implementation**
- **Complete ACMG/AMP 2015 & 2023 Guidelines**: Full implementation including optional PP5/BP6 and enhanced de novo (PS2_Very_Strong) criteria
- **All Variant Types Supported**: Missense, nonsense, frameshift, splice variants, synonymous, intronic, regulatory, and structural variants
- **Evidence Integration**: Sophisticated algorithms for combining population, computational, functional, and genetic evidence

### 🧮 **Enhanced Statistical Framework**
- **VAMPP-Score Integration**: Advanced metascore algorithm based on VAMPP-score methodology with weighted predictor integration
- **Conservation Analysis**: Multi-species phylogenetic conservation scoring (PhyloP, GERP++, SiPhy)
- **Prevalence-Based Thresholds**: Dynamic threshold adjustment based on disease prevalence and population frequency
- **Statistical Tests**: Fisher's Exact Test for case-control data analysis
- **Gene-Specific Calibration**: Customized thresholds for different gene classes

### 🔍 **Comprehensive In Silico Integration**
- **Primary Predictors**: REVEL, CADD, AlphaMissense, MetaRNN, ClinPred, BayesDel
- **High Priority Predictors**: VEST4, PrimateAI, ESM1b, PROVEAN (enhanced accuracy)
- **Conservation Scores**: PhyloP (100-way, 30-way, 17-way), GERP++
- **Functional Predictors**: SIFT, PolyPhen-2, MutationTaster, FATHMM
- **Splice Prediction**: SpliceAI, MMSplice, Ada, RF, dbscSNV (for intronic and splice variants)
- **Meta-Predictors**: MetaSVM, MetaLR, MutPred, LRT

### 🌐 **API Integration & Automation**
- **ClinVar Integration**: Automatic variant status retrieval with user override options
- **Ensembl API**: Automatic chromosome and gene information fetching
- **Intelligent Caching**: API response caching for improved performance
- **Fallback Mechanisms**: Graceful handling of API failures with manual input options

### 📊 **Advanced Evidence Evaluation**
- **PP3/BP4 Enhancement**: Sophisticated in silico evidence evaluation with metascore thresholds
- **PM2 Refinement**: Gene-specific population frequency thresholds
- **PS2/PM6 Enhancement**: Comprehensive de novo variant evaluation
- **PM1 Implementation**: Functional domain and hotspot analysis
- **PS4 Analysis**: Case-control statistical evaluation

### 🎯 **User Experience**
- **Interactive Interface**: Step-by-step prompts with clear instructions and examples
- **Input Validation**: Comprehensive validation with helpful error messages
- **Alias Support**: User-friendly shortcuts and abbreviations
- **Colored Output**: Enhanced readability with color-coded information
- **Test Mode**: Built-in test scenarios for demonstration and validation

## 🧠 Intended Users

- **Clinical Geneticists**: Standardized, reproducible variant classification for clinical practice
- **Molecular Diagnostics Professionals**: Efficient variant curation workflow with comprehensive evidence integration
- **Bioinformatics Researchers**: Advanced computational framework for variant analysis research
- **Graduate Students**: Educational tool for learning ACMG guidelines and variant interpretation principles
- **Variant Curation Teams**: Collaborative variant interpretation with detailed reporting
- **Genetic Counselors**: Evidence-based variant assessment for patient counseling

## � Quick Start Guide

### 1. Installation
```bash
# Clone the repository
git clone https://github.com/Bilmem2/acmg-assessor.git
cd acmg-assessor

# Install dependencies
pip install -r requirements.txt
```

### 2. First Run (Test Mode)
```bash
# Try the test mode first to see how the tool works
python acmg_assistant.py --test
```

### 3. Real Variant Analysis
```bash
# Run with real variant data
python acmg_assistant.py

# Use ACMG 2023 guidelines
python acmg_assistant.py --acmg-2023
```

### 4. Collect In Silico Scores
Before starting, gather scores from:
- **Varsome** (https://varsome.com) - Primary resource
- **ClinVar** (https://www.ncbi.nlm.nih.gov/clinvar/)
- **dbNSFP** (for academic users)

### 5. Follow Interactive Prompts
The tool will guide you through:
- Basic variant information
- Population frequencies
- In silico predictions
- Functional and genetic data
- Generate comprehensive report

## �💻 Installation & Usage

### Prerequisites
- **Python 3.8 or higher**
- **Internet connection** (for API calls; optional for test mode)
- **Operating System**: Windows, macOS, or Linux

### Installation

```bash
# Clone the repository
git clone https://github.com/Bilmem2/acmg-assessor.git
cd acmg-assessor

# Install dependencies
pip install -r requirements.txt
```

### Required Dependencies
```bash
pip install requests scipy numpy colorama
```

### Basic Usage

**For Standalone Executable:**
```bash
# Normal mode (standard usage)
acmg_assistant.exe

# ACMG 2023 guidelines
acmg_assistant.exe --acmg-2023
```

**For Python Installation:**
```bash
# Standard mode with ACMG 2015 guidelines
python acmg_assistant.py

# Use ACMG 2023 guidelines
python acmg_assistant.py --acmg-2023

# Test mode with sample data (Python only)
python acmg_assistant.py --test

# Show version information
python acmg_assistant.py --version
```

### Interactive Workflow

The tool guides users through comprehensive data collection:

#### 1. **📋 Basic Variant Information**
- Gene symbol (with automatic validation)
- Chromosome and genomic position (optional for API integration)
- Reference and alternate alleles
- HGVS nomenclature (cDNA and protein)
- Variant consequence (with comprehensive SO term support)

#### 2. **👥 Population Data**
- gnomAD frequencies (global and population-specific)
- ExAC, 1000 Genomes, ESP frequencies
- Homozygous/heterozygous counts
- Disease prevalence information

#### 3. **🔬 In Silico Scores**
- **All scores entered manually** (no automatic retrieval)
- Primary predictors: REVEL, CADD, AlphaMissense, etc.
- Conservation scores: PhyloP, GERP++, PhastCons
- Functional predictors: SIFT, PolyPhen-2, MutationTaster
- Splice predictors: SpliceAI (for relevant variants)

#### 4. **🧬 Genetic Data**
- Inheritance pattern and penetrance
- Zygosity and allelic state
- Compound heterozygous information
- De novo status with parental confirmation

#### 5. **⚡ Functional Data**
- Segregation analysis results
- Functional study outcomes
- Phenotype matching assessment
- Case-control data (if available)

## 🔧 Building Standalone Executable

**For developers who want to create the executable:**

### Prerequisites
- Python 3.8+ installed
- All project dependencies installed

### Build Process
```bash
# Install build requirements
pip install -r requirements_build.txt

# Build executable (Windows)
python build_executable.py

# Or use the batch script
build.bat
```

### Distribution
1. **Test the executable** in `ACMG_Assistant_Portable` folder
2. **Zip the entire folder** for distribution
3. **Share with end users** who don't have Python installed

**Executable features:**
- **Self-contained**: No Python installation needed
- **~50-100 MB size**: Includes Python runtime and all dependencies
- **Console interface**: Easy-to-use interactive prompts
- **Portable**: Can be run from any folder or USB drive
- **Normal mode only**: Test mode requires Python installation

## 🏗️ Project Structure

```
acmg_assessor/
├── acmg_assistant.py              # Main application entry point
├── requirements.txt               # Python dependencies
├── requirements_build.txt         # Build dependencies (PyInstaller)
├── build_executable.py            # Executable build script
├── build.bat                      # Windows build batch script
├── README.md                     # Documentation
├── LICENSE                       # MIT License
├── config/
│   ├── __init__.py
│   └── constants.py              # Configuration constants and thresholds
├── core/
│   ├── __init__.py
│   ├── variant_data.py           # Variant data structure and validation
│   ├── acmg_classifier.py        # ACMG classification logic
│   └── evidence_evaluator.py     # Evidence evaluation algorithms
├── utils/
│   ├── __init__.py
│   ├── input_handler.py          # Interactive user input handling
│   ├── api_client.py             # External API integration
│   ├── validators.py             # Data validation utilities
│   └── report_generator.py       # Comprehensive report generation
└── tests/
    ├── test_camta1.py            # CAMTA1 variant test case
    ├── test_vus_scenarios.py     # VUS classification scenarios
    ├── test_vus_upgrade_scenarios.py # VUS upgrade test cases
    └── test_comprehensive_scenarios.py # Comprehensive test suite
```

## 📈 Algorithm Enhancements

### VAMPP-Score Integration
- **Weighted Predictor Combination**: Sophisticated weighting scheme based on predictor performance
- **Normalization Framework**: Min-max normalization with pathogenic direction alignment
- **Gene-Specific Calibration**: Adaptive scoring based on gene characteristics
- **Sigmoid Transformation**: Enhanced score distribution for better discrimination

### Conservation Analysis
- **Multi-Species Integration**: PhyloP 100-way, 30-way, and 17-way vertebrate conservation
- **Evolutionary Constraint**: GERP++ and SiPhy constraint scores
- **Functional Conservation**: PhastCons conserved element detection

### Prevalence-Based Thresholds
- **Dynamic Adjustment**: Threshold modification based on disease prevalence
- **Population Stratification**: Ancestry-specific frequency analysis
- **Penetrance Consideration**: Incomplete penetrance adjustment factors

## 📊 Output Files

- **`variant_classification_report.txt`**: Comprehensive classification report with evidence summary
- **`variant_classification_log.txt`**: Detailed processing log for reproducibility
- **`api_cache.json`**: Cached API responses for performance optimization

## 🧪 Test Scenarios

The tool includes built-in test cases for validation:
- **CAMTA1 c.319A>G**: De novo missense variant test
- **VUS Scenarios**: Diverse VUS cases with metascore algorithm
- **Comprehensive Scenarios**: Real-world variant classification tests

```bash
# Run individual test cases
python test_camta1.py
python test_vus_scenarios.py

# Run comprehensive test suite
python test_comprehensive_scenarios.py
```

## 🔬 ACMG Criteria Implementation

### Enhanced Pathogenic Criteria
- **PVS1**: Loss-of-function with mechanism validation
- **PS1**: Same amino acid change with enhanced matching
- **PS2/PS2_Very_Strong**: De novo with parental confirmation levels
- **PS3**: Functional studies with evidence strength assessment
- **PS4**: Case-control with Fisher's Exact Test
- **PM1**: Functional domain analysis with constraint scoring
- **PM2**: Gene-specific population frequency thresholds
- **PP3**: Enhanced in silico evidence with VAMPP-score integration

### Enhanced Benign Criteria
- **BA1/BS1**: Dynamic frequency thresholds based on disease prevalence
- **BP4**: Sophisticated computational evidence integration
- **BP7**: Comprehensive splice impact prediction

## 🌐 API Integration

### ClinVar Integration
```python
# Automatic variant lookup with user confirmation
clinvar_status = api_client.get_clinvar_status(chr, pos, ref, alt)
```

### Ensembl Integration
```python
# Gene and chromosome information retrieval
gene_info = api_client.get_gene_info(gene_symbol)
```

## 📚 Statistical Framework

### Fisher's Exact Test
```python
# Case-control association analysis
p_value = fisher_exact_test(cases_with, total_cases, controls_with, total_controls)
```

### VAMPP-Score Calculation
```python
# Weighted metascore with prevalence adjustment
vampp_score = calculate_vampp_score(predictors, weights, prevalence)
```

## ⚠️ Important Notes

### Manual Score Entry

⚠️ **CRITICAL: All in silico scores must be entered manually**

**Score Collection Process:**
1. **Use Varsome** (https://varsome.com) - Comprehensive in silico score aggregation
2. **Check ClinVar** (https://www.ncbi.nlm.nih.gov/clinvar/) - Variant annotations
3. **dbNSFP Database** - Academic users can access comprehensive scores
4. **Individual Predictor Websites** - For specific score verification

**Supported Predictors:**
- **Primary Metascores**: REVEL, CADD, AlphaMissense, MetaRNN, ClinPred, BayesDel
- **High Priority Predictors**: VEST4, PrimateAI, ESM1b, PROVEAN (enhanced classification)
- **Conservation Scores**: PhyloP (100-way, 30-way, 17-way), GERP++
- **Individual Predictors**: SIFT, PolyPhen-2, MutationTaster, FATHMM
- **Splice Predictors**: SpliceAI, MMSplice, Ada, RF, dbscSNV (for splice-relevant variants)

**Important Notes:**
- Tool performs NO automatic score retrieval
- NO local database or Excel file searches
- Users must collect scores from external resources
- Enter "unknown" or "not available" if scores cannot be obtained

### Genomic Position
- **Optional but recommended** for ClinVar integration
- **Can be omitted if not available**
- **Manual ClinVar status entry as fallback**

### Clinical Use Disclaimer

⚠️ **IMPORTANT CLINICAL DISCLAIMER** ⚠️

This tool is for **research and educational purposes only**. 

**Clinical Limitations:**
- This tool is NOT intended for direct clinical use without professional oversight
- All results must be validated by certified professionals following institutional guidelines
- Clinical variant interpretations require comprehensive manual review
- The tool cannot replace expert clinical judgment and laboratory protocols

**Validation Requirements:**
- Always cross-check results with primary sources (ClinVar, Varsome, LOVD, etc.)
- Verify all in silico scores independently
- Confirm population frequencies with multiple databases
- Review functional studies and segregation data independently

**Professional Responsibility:**
- The author is not responsible for clinical decisions made using this tool
- Users must follow their institutional and national guidelines
- Clinical reporting must include appropriate disclaimers and limitations
- Consider tool results as ONE piece of evidence, not definitive classification

**Recommended Use:**
- Educational training for ACMG guidelines
- Research variant analysis and method development
- Preliminary analysis prior to clinical review
- Quality control and consistency checking

##  Citation & References

### Primary Citation
If you use this tool in your research, please cite:
```
ACMG Variant Classification Assistant v1.0.0
Author: Can Sevilmiş
GitHub: https://github.com/cansevilmis/acmg-assessor
Year: 2025
```

### ACMG Guidelines References
This tool implements the following guidelines:
```
Richards, S., Aziz, N., Bale, S., et al. (2015). 
"Standards and guidelines for the interpretation of sequence variants: a joint consensus 
recommendation of the American College of Medical Genetics and Genomics and the Association 
for Molecular Pathology." 
Genetics in Medicine, 17(5), 405-423. DOI: 10.1038/gim.2015.30

Brnich, S.E., Abou Tayoun, A.N., Couch, F.J., et al. (2019). 
"Recommendations for application of the functional evidence PS3/BS3 criterion using the 
ACMG/AMP sequence variant interpretation framework." 
Genome Medicine, 12, 3. DOI: 10.1186/s13073-019-0690-2
```

### Statistical Framework References
This tool implements advanced statistical methods. Key references:
```
Aydin, E., Ergun, B., Akgun-Dogan, O., Alanay, Y., Hatirnaz Ng, O., & Ozdemir, O. (2024).
"A New Era in Missense Variant Analysis: Statistical Insights and the Introduction of 
VAMPP-Score for Pathogenicity Assessment." 
bioRxiv, 2024.07.11.602867. DOI: 10.1101/2024.07.11.602867

Ozdemir, O., et al. (2024). "Molecular and In Silico Analysis of the CHEK2 Gene in 
Individuals with High Risk of Cancer Predisposition from Türkiye." 
Cancers, 16(22), 3876. DOI: 10.3390/cancers16223876
```

## 👤 Contact & Support

- **Author**: Can Sevilmiş
- **LinkedIn**: [cansevilmiss](https://linkedin.com/in/cansevilmiss)
- **GitHub**: [Bilmem2](https://github.com/Bilmem2)

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
