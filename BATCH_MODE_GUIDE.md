# Batch Mode Usage Guide

## What is Batch Mode?
Batch Mode allows you to analyze multiple genetic variants automatically by providing a CSV file containing all the necessary variant information. The program will read each row, process the data, and generate classification results and reports for each variant without manual input.

---

## How to Use Batch Mode

1. **Prepare your CSV file**
   - The CSV file should contain one variant per row.
   - Each column should represent a specific data field required for analysis.

2. **Required Columns**
   - `chromosome`: Chromosome number or name (e.g., 17)
   - `position`: Genomic position (e.g., 43093464)
   - `ref_allele`: Reference allele (e.g., C)
   - `alt_allele`: Alternate allele (e.g., T)
   - `variant_type`: Type of variant (e.g., missense, frameshift)
   - `consequence`: Variant consequence (e.g., missense_variant)
   - `gene`: Gene name (e.g., TP53)


3. **Optional & Supported Columns**
   - You can add the following optional columns to your CSV file. The program will automatically use these fields in the analysis if present:
     - **Population data:** `gnomad_af`, `exac_af`, `topmed_af`
     - **In silico scores:** `cadd_phred`, `revel`, `sift`, `polyphen2`, `mutation_taster`, `dann`, `fathmm`, `spliceai_ag_score`, `spliceai_al_score`, `spliceai_dg_score`, `spliceai_dl_score`
     - **Genetic data:** `segregation`, `de_novo`, `denovo`, `maternity_confirmed`, `paternity_confirmed`, `inheritance`
     - **Functional data:** `case_control`, `functional_study`
     - **Clinical data:** `clinical_significance`
     - **Patient phenotypes:** `patient_phenotypes`

   - If these fields are present, fill in the relevant columns for each variant. Missing fields will be treated as "no data" in the analysis.



4. **Example CSV Format**

```csv
chromosome,position,ref_allele,alt_allele,variant_type,consequence,gene,gnomad_af,cadd_phred,revel,case_control,clinical_significance,patient_phenotypes
17,43093464,C,T,missense,missense_variant,TP53,0.00001,25.6,0.82,not_observed_in_healthy,pathogenic,"intellectual disability, microcephaly"
13,32906641,G,A,frameshift,frameshift_variant,BRCA2,0.00002,30.1,0.91,observed_in_healthy,benign,"breast cancer"
```


5. **Running Batch Mode**
   - Place your CSV file in the same directory as the program, or provide the full path.
   - Run the following command in your terminal:

```bash
python acmg_assistant.py --batch your_file.csv
```

- If you are using the executable version, select Batch Mode and provide the CSV file path when prompted.

---

## Notes
- All variant data (required and optional fields) should be present in the CSV file. No manual input will be requested for each variant.
- If optional fields (in silico, population, functional, genetic, clinical, etc.) are present, they will be used automatically in the analysis.
- Reports for each variant will be generated automatically and saved in the output directory.
- If the file is not found, check the path and filename.

---

## Troubleshooting
- **File Not Found:** Ensure the CSV file is in the correct directory or provide the full path.
- **Missing Columns:** Make sure all required columns are present and correctly named.
- **Encoding Issues:** Save your CSV file with UTF-8 encoding.

---

For further assistance, refer to the README or contact the developer.
