# VAMPP-score
**Variant Analysis with Multiple Pathogenicity Predictors**

## Description
VAMPP-score is a metascore for missense variant pathogenicity prediction. It utilizes 56 available in silico pathogenicity predictors (ISPPs) from dbNSFP(v4.9A) based on their performance in distinguishing and identifying three main variant groups (Pathogenic, Benign, and Unknown). The tool is trained using ClinVar data.

ClinVar: [Landrum et al., 2020](https://doi.org/10.1093/nar/gkz972)

dbNSFP v4: [Liu et al., 2020](https://doi.org/10.1186/s13073-020-00803-9)

dbNSFP: [Liu et al., 2011](https://doi.org/10.1002/humu.21517)



## Features
* Pathogenicity score: VAMPP-score provides a weighted pathogenicity score for missense variants, highlighting the best performing ISPPs for each gene.
* Gene-based ISPP performance: To calculate the VAMPP-score, the best-performing ISPPs are determined in a gene-based manner for weighting. This gene-based performance is defined as the "Gene coefficient (_M_)." M reflects the ISPP performance based on its ability to distinguish variant groups and accurately assign scores to variants in the order of Benign, Unknown, and Pathogenic in a 0-1 ranking. It focuses specifically on the ISPP's performance in detecting Pathogenic vs. Benign variants, while ensuring the Unknown group is ranked in the middle.


## Usage and Access
**VAMPP-score** is not available for installation. However, you can navigate to the web interface to access scores. VAMPP-score is planned to manually be updated on the first Friday of each month, with recalculations available on the web platform. Each month's annotation data is stored in a folder named YYYY-MM, accessable through Google Drive. The annotation data - calculated VAMPP-scores- is uploaded each month in the same format and can be accessed through [VAMPP-score-data](https://drive.google.com/drive/folders/1emkHcTlxgjH6G-2Yl4wQQnKi5Wsip4IY?usp=drive_link).


VAMPP-score calculation scripts are available in the `scripts` directory, organized into parts:

- `10X`: Data curation and database construction
- `20X`: Calculate VAMPP-score components and the raw VAMPP-score
- `30X`: Calculate VAMPP-ranked and obtain annotation data for dbNSFP variants

The scripts can be run from source after installing the dependencies. The provided scripts illustrate the workflow from creating the background DuckDB database with all input files to calculating the VAMPP-score.


### Dependencies
- Python 3.x
- DuckDB 0.9.3
- pandas
- request
- PyArrow 
- scipy

## Web Interface
VAMPP-score will soon be available at [vamppscore.com ](https://vamppscore.com/) with a dynamic interface. You can access the VAMPP-score of a variant, a gene-score for each ISPP based on their prediction performance, and the best-performing ISPPs for your gene of interest.

## Contact

If you have any questions, please contact us at:
* [Ozkan Ozdemir](mailto:Ozkan.Ozdemir@acibadem.edu.tr)
* [Eylul Aydin](mailto:Eylul.Aydin@live.acibadem.edu.tr)

## Citation
If you use VAMPP-score or any component of our statistical framework in your work, please cite us:
**A New Era in Missense Variant Analysis: Statistical Insights and the Introduction of VAMPP-Score for Pathogenicity Assessment**
Eylul Aydin, Berk Ergun, Ozlem Akgun-Dogan, Yasemin Alanay, Ozden Hatirnaz Ng, Ozkan Ozdemir
bioRxiv (2024.07.11.602867; doi: https://doi.org/10.1101/2024.07.11.602867)
