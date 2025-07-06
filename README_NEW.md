# 🧬 ACMG Variant Classification Assistant

**Advanced ACMG/AMP Variant Classification with Enhanced Statistical Framework**  
*Last Updated: July 2025*

A comprehensive tool for classifying genetic variants according to ACMG/AMP 2015 and 2023 guidelines. Features enhanced VAMPP-score implementation, comprehensive in silico predictor integration, and sophisticated evidence evaluation algorithms.

## ⚙️ Key Features

- **Complete ACMG/AMP Guidelines**: 2015 & 2023 standards with PP5/BP6 and enhanced PS2_Very_Strong
- **30+ In Silico Predictors**: REVEL, CADD, AlphaMissense, VEST4, PrimateAI, ESM1b, SpliceAI, MMSplice
- **VAMPP-Score Integration**: Advanced metascore with weighted predictor combination
- **API Integration**: ClinVar and Ensembl with intelligent caching
- **Statistical Framework**: Fisher's Exact Test, prevalence-based thresholds, conservation analysis

## 💻 Installation & Usage

### Standalone Executable (Recommended)
```bash
# Download and extract ACMG_Assistant_Portable.zip
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

## 🧪 Test Cases

```bash
python test_camta1.py                    # CAMTA1 variant
python test_comprehensive_scenarios.py  # Multiple scenarios
python test_vus_scenarios.py            # VUS classification
```

## 🏗️ Project Structure

```
acmg_assessor/
├── acmg_assistant.py              # Main application
├── build_executable.py            # Executable builder
├── requirements.txt               # Dependencies
├── requirements_build.txt         # Build dependencies
├── config/constants.py            # Configuration
├── core/                          # Classification logic
├── utils/                         # Input/output handling
└── tests/                         # Test scenarios
```

## 📈 ACMG Criteria Implementation

**Enhanced Pathogenic**: PVS1, PS1-4, PM1-6, PP1-5  
**Enhanced Benign**: BA1, BS1-4, BP1-7  
**Special Features**: VAMPP-score PP3/BP4, gene-specific PM2, statistical PS4

## ⚠️ Important Notes

- **Manual Score Entry**: No automatic retrieval from databases
- **Research Use**: Not for direct clinical decisions without validation
- **Internet Required**: For API calls (ClinVar, Ensembl)
- **Test Mode**: Python installation only

## 📖 References

```
Richards, S., et al. (2015). Standards and guidelines for sequence variant interpretation. 
Genetics in Medicine, 17(5), 405-423.

Aydin, E., et al. (2024). VAMPP-Score for Pathogenicity Assessment. 
bioRxiv, 2024.07.11.602867.
```

## 👤 Contact

- **Author**: Can Sevilmiş
- **LinkedIn**: [cansevilmiss](https://linkedin.com/in/cansevilmiss)
- **GitHub**: [Bilmem2](https://github.com/Bilmem2)

## 📄 License

MIT License - see [LICENSE](LICENSE) file.
