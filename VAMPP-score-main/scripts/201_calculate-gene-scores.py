# import necessary modules

import pandas as pd
import numpy as np
from scipy import stats

# Load variants for calculation - if you used another file format, adjust accordingly
df = pd.read_parquet('/path/to/your/parquet/with/pre-filtered/missense/clinvar/data/[file_name].parquet')

# Replace unsupported characters with NaN 
df = df.replace(['.', '-', ' ', 'None'], np.nan)

# Optionally rename your dataframe before statistical analysis
df_kruskal = df

# Define the columns with ISPP-rankscores
column_names = list(df_kruskal.keys())
rankscore_cols = list(filter(lambda x: 'rankscore' in x, column_names))

# Define variant groups according to the ClinVar submissions
pathogenic_rows = df_kruskal[df_kruskal['ClinicalSignificance'].isin(['Pathogenic','Likely pathogenic','Pathogenic/Likely pathogenic'])]
benign_rows = df_kruskal[df_kruskal['ClinicalSignificance'].isin(['Benign','Likely benign','Benign/Likely benign'])]
unknown_rows = df_kruskal[df_kruskal['ClinicalSignificance'].isin(['Uncertain significance'])]


# Create a set to store unique consecutive values
unique_genes = set()

# Iterate over each row and extract unique consecutive values
for row in df_kruskal['genename']:
    values = row.split(';')
    cleaned_values = []
    for value in values:
        if value not in cleaned_values:
            cleaned_values.append(value)
    unique_genes.add(','.join(cleaned_values))

# Execute Kruskal-Wallis for the multiple comparison between variant groups
def calc_kruskal(path_df, beni_df, unkn_df) -> dict:
    kruskal_dict = {}
    for rankscore in rankscore_cols:
        try:
            kruskal_dict[rankscore] = stats.kruskal(path_df[rankscore].dropna(), beni_df[rankscore].dropna(), unkn_df[rankscore].dropna(), nan_policy='omit')
        except ValueError:
            continue
    return dict(map(lambda x: (x, kruskal_dict[x]), sorted(kruskal_dict, key=lambda x: kruskal_dict[x].pvalue)))


# Execute Wilcoxon-Rank-Sum for the pairwise comparison between Pathogenic vs. Benign variant groups

def calc_pairwise_comp(path_df, beni_df, unkn_df) -> dict:
    wilcoxon_path_ben_dict = {}
    pairwise_comp_results = {}
    for rankscore in path_df.columns:
        tool_p_scores = path_df[rankscore].dropna()
        tool_b_scores = beni_df[rankscore].dropna()
        wilcoxon_path_ben_dict[rankscore] = stats.ranksums(tool_p_scores, tool_b_scores)
        pairwise_comp_results['p-b'] = dict(map(lambda x: (x, wilcoxon_path_ben_dict[x]), sorted(wilcoxon_path_ben_dict, key=lambda x: wilcoxon_path_ben_dict[x].pvalue)))
    return pairwise_comp_results

# Get the gene scores - 'M' as the gene coefficient metric

def get_gene_score(genename):

    variant_score_multiplier = {}

    p_rows = pathogenic_rows[pathogenic_rows['genename'].str.contains(genename)]
    b_rows = benign_rows[benign_rows['genename'].str.contains(genename)]
    u_rows = unknown_rows[unknown_rows['genename'].str.contains(genename)]

    kruskal_dict = calc_kruskal(p_rows, b_rows, u_rows)
  
  # Include the ISPPs with statistically significant results (p<0.05) from Kruskal-Wallis for further analysis
    chosen_tools_list = list(filter(lambda y: kruskal_dict[y].pvalue < 0.05, kruskal_dict))

    pairwise_comp_dict = calc_pairwise_comp(p_rows[chosen_tools_list], b_rows[chosen_tools_list], u_rows[chosen_tools_list])

    for tool in chosen_tools_list:
        tool_p_scores, tool_b_scores, tool_u_scores = (
            p_rows[tool].dropna().astype(float), b_rows[tool].dropna().astype(float), u_rows[tool].dropna().astype(float)
        )

        kruskal_mult = 1 - kruskal_dict[tool].pvalue
        pairwise_comp_mult = 1 - pairwise_comp_dict['p-b'][tool].pvalue
        theorical_u_mean = abs((tool_p_scores.mean() + tool_b_scores.mean()) / 2)
        u_max_dist = abs(tool_p_scores.mean() - theorical_u_mean)
        proximity_to_unk_mean_mult = 1 - (abs(tool_u_scores.mean() - theorical_u_mean) / u_max_dist)
        score_diff_mult = tool_p_scores.mean() - tool_b_scores.mean()

        s = kruskal_mult * pairwise_comp_mult * proximity_to_unk_mean_mult * score_diff_mult
        if s > 0:
            variant_score_multiplier[tool] = s

    return variant_score_multiplier

all_genenames = set()

all_genenames = all_genenames.union(df_kruskal['genename'])

from concurrent.futures import ThreadPoolExecutor

with ThreadPoolExecutor(52) as executor:
    all_genescores = dict(executor.map(lambda x: (x, get_gene_score(x)), all_genenames))

# write the output gene scores to a JSON file
import json

with open('/path/to/your/output/directory/[gene_scores].json', 'w') as f:
    json.dump(all_genescores, f)
