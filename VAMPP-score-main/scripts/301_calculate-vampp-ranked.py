#import necessary modules

import pandas as pd
import numpy as np
import glob
import os

# Directory where the VAMPP-score by chromosome TSV files are stored
directory = "/file/path/to/your/chromosome/files/with/calculated/raw-vampp"

# Pattern to match the files, adjust as needed - here it is stored for each chromosome in separate tsv files
file_pattern = os.path.join(directory, "vamppscore_dbNSFP4.7a_variant.chr*.tsv.gz")

# List to store scores from all files
all_scores = []

# List to store DataFrames for later concatenation
dataframes = []

for filename in glob.glob(file_pattern):
    # Read the specific columns to save memory and enforce data type for chromosome column, adjust as needed
    df = pd.read_csv(filename, compression='gzip', header=0, delimiter='\t', 
                     usecols=[0, 1, 2, 3, 4, 5], 
                     names=['#chr', 'pos', 'ref', 'alt', 'Gene', 'VAMPP_score'],
                     dtype={'#chr': str})  # Ensure chromosome column is read as string
    
    # Display unique chromosome values to verify consistency (optional: can be removed later)
    print(f"Unique chromosome values in {filename}: {df['#chr'].unique()}")
    
    # Store scores and DataFrame
    all_scores.extend(df['VAMPP_score'].tolist())
    dataframes.append(df)

# Compute global Min and Max for normalization
min_score = min(all_scores)
max_score = max(all_scores)

# Apply Min-Max scaling
for df in dataframes:
    df['VAMPP_score_ranked'] = ((df['VAMPP_score'] - min_score) / (max_score - min_score)) if max_score > min_score else 0

# Concatenate all DataFrames
final_df = pd.concat(dataframes)

# Save the final DataFrame to a new gzipped TSV file, not separated for each chromosome, readily annotated for each variant
final_df.to_csv('/file/path/to/your/output/file/[VAMPP-score_file_name].tsv.gz', sep='\t', index=False, compression='gzip')
