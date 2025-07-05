# 🧬 ACMG Variant Classification Assistant

**Last Updated:** July 5, 2025

The ACMG Variant Classification Assistant is a robust, user-friendly tool for classifying genetic variants according to the ACMG/AMP 2015 and 2023 guidelines. It is designed for clinical geneticists, bioinformaticians, researchers, and students, and is distributed as a portable Windows executable (`.exe`) that requires no Python installation or as a Python script for cross-platform use. The tool guides users through all required data entry, integrates population and in silico data, and applies advanced statistical frameworks (including a VAMPP-score-like metascore) for evidence-based, standardized variant classification. **All in silico scores are entered manually by the user; there is no local Excel or dataset dependency.**

**License:** This project is licensed under the MIT License. See LICENSE file for details.

Download the ZIP file  `acmg_assistant.zip` [HERE](https://drive.google.com/file/d/1vaJ3ue0sQYPycaRSKM1XRxZB4tllVyle/view?usp=drive_link).

## ⚙️ Key Features
- Implements ACMG/AMP 2015 and 2023 guidelines, including optional PP5/BP6 and enhanced de novo (PS2_Very_Strong) criteria.
- Supports all major variant types: missense, nonsense, frameshift, splice, synonymous, and more.
- Integrates in silico predictors: REVEL, CADD, SIFT, AlphaMissense, MetaRNN, ClinPred, BayesDel, MutationTaster, PolyPhen-2, and others.
- **All in silico scores must be entered manually by the user. No local Excel or dataset search, and no automatic in silico score retrieval.** Use resources like Varsome for reference.
- Automatically retrieves chromosome information (Ensembl API) and ClinVar status (with caching for speed).
- Advanced statistical analyses: Fisher’s Exact Test for case-control data and a VAMPP-score-like metascore for missense variants.
- Step-by-step command-line prompts with input validation, error handling, and clear English instructions.
- Identifies missing or incomplete data and suggests resources (Varsome, ClinVar, Orphanet, PubMed) for further investigation.
- Saves detailed results to `variant_classification_report.txt`, including classification, applied criteria, and suggestions for missing data.
- Test mode with sample data for rapid demonstration or educational use.
- Fully portable: runs as a standalone `.exe` for Windows or as a Python script for all platforms, no installation or local data files required.

## 🧠 Intended Users
- Clinical and medical geneticists for standardized, reproducible variant classification.
- Molecular diagnostics professionals for efficient, evidence-based variant curation.
- Bioinformatics researchers for integrating population, in silico, and functional data.
- MSc/PhD students in genetics, genomics, or bioinformatics for learning ACMG guidelines and variant interpretation.
- Variant curation trainees and educators for hands-on practice and teaching.

## 💻 Installation & Usage

### Prerequisites
- **Windows Users**:
  - Windows OS (the `.exe` version is for Windows only).
  - Internet connection (required for Ensembl and ClinVar API calls; optional for test mode).
  - **Download**: Get `acmg_assistant.zip` [HERE](https://drive.google.com/file/d/1vaJ3ue0sQYPycaRSKM1XRxZB4tllVyle/view?usp=drive_link).
- **macOS/Linux or Python Users**:
  - Python 3.x must be installed ([python.org](https://www.python.org/downloads/)).
  - Required dependencies:
    ```bash
    pip install requests scipy numpy colorama
    ```
- Internet connection (required for Ensembl and ClinVar API calls; optional for test mode).

### Installation (Windows Executable)
1. Download the file `acmg_assistant.zip` and extract it.
2. Double-click `acmg_assistant.exe` to launch the command-line interface.

### Installation (Python Script)
1. Download the `acmg_assistant.py` file from the repository: [GitHub Repo](https://github.com/Bilmem2/acmg_assistant).
2. Install dependencies:
   ```bash
   pip install requests scipy numpy colorama
   ```
3. Run the script:
   ```bash
   python src/acmg_assistant.py
   ```

### Usage
1. Launch the tool (`acmg_assistant.exe` or `python src/acmg_assistant.py`) and follow the interactive prompts.
2. Enter the required information:
   - **Basic Information**: Gene, chromosome, position, reference/alternate alleles, cDNA/protein changes, VEP consequence.
   - **Population Data**: Allele frequency (e.g., gnomAD), subpopulation frequencies, disease prevalence.
   - **In Silico Scores**: All scores and predictions (e.g., REVEL, CADD, MetaRNN, AlphaMissense, SIFT, PolyPhen-2, etc.) must be entered manually. **No local Excel or dataset lookup is performed. No automatic in silico score retrieval.** Use resources like Varsome for reference.
   - **Genetic Data**: Inheritance pattern, zygosity, allelic data (for recessive diseases).
   - **Functional Data**: Functional test results, segregation, de novo status, phenotype match.
3. The tool validates all inputs (e.g., HGVS format, valid chromosomes) and provides automatic predictions and suggestions where possible.
4. Review the results:
   - **Classification**: Pathogenic, Likely Pathogenic, VUS, Likely Benign, or Benign.
   - **Applied Criteria**: ACMG evidence codes (e.g., PVS1, PS3, PP4) with descriptions.
   - **Varsome URL**: Direct link to verify variant details.
5. Results are saved to `variant_classification_log.txt` for record-keeping and reproducibility.

**Command-line Options**:
- `--test`: Run in test mode with sample data (no internet required).
- `--acmg-2023`: Use ACMG/AMP 2023 guidelines.

## ⚠️ Disclaimer
This tool is for research and educational purposes only. Clinical variant interpretations must be validated by certified professionals following institutional and national guidelines. Always cross-check results with primary sources (e.g., ClinVar, Varsome). The author is not responsible for clinical decisions made using this tool.

## 🔧 Planned & Future Improvements
- Import variant data from VCF files for batch processing.
- Visualize in silico score distributions with graphical outputs.
- Generate PDF or HTML reports for clinical or research use.
- Improve performance with asynchronous API calls.

## 📚 Citation & References
This tool uses a VAMPP-score-like metascore approach for in silico pathogenicity prediction, inspired by the original VAMPP-score framework. If you use this tool, the VAMPP-score, or any component of their statistical framework in your work, please cite the original VAMPP-score publication:

> Eylul Aydin, Berk Ergun, Ozlem Akgun-Dogan, Yasemin Alanay, Ozden Hatirnaz Ng, Ozkan Ozdemir. "A New Era in Missense Variant Analysis: Statistical Insights and the Introduction of VAMPP-Score for Pathogenicity Assessment." *bioRxiv* (2024.07.11.602867). [DOI: 10.1101/2024.07.11.602867](https://doi.org/10.1101/2024.07.11.602867)

For in silico and molecular analysis methodology, see:

> Ozdemir, O., Bychkovsky, B. L., Unal, B., Onder, G., Amanvermez, U., Aydin, E., Ergun, B., Sahin, I., Gokbayrak, M., Ugurtas, C., Koroglu, M. N., Cakir, B., Kalay, I., Cine, N., Ozbek, U., Rana, H. Q., Hatirnaz Ng, O., & Agaoglu, N. B. (2024). "Molecular and In Silico Analysis of the CHEK2 Gene in Individuals with High Risk of Cancer Predisposition from Türkiye." *Cancers*, 16(22), 3876. [DOI: 10.3390/cancers16223876](https://doi.org/10.3390/cancers16223876)

**VAMPP-score** is a registered framework. For more information and access to the original VAMPP-score data and scripts, see: [vamppscore.com](https://vamppscore.com/) and the [VAMPP-score-data Google Drive](https://drive.google.com/drive/folders/1emkHcTlxgjH6G-2Yl4wQQnKi5Wsip4IY?usp=drive_link).

## 👤 Contact
- **Author**: Can Sevilmiş  
- **LinkedIn**: [cansevilmiss](https://www.linkedin.com/in/cansevilmiss/)
