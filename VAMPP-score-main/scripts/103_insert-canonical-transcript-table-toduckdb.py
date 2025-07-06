# Import necessary modules

import pandas as pd
import duckdb

# Open a connection to DuckDB
conn = duckdb.connect(database='/file/path/to/your/duckdb/[your_database].duckdb')

# Ensure the target table exists in DuckDB, adjust table names and column names as needed
conn.execute("""
    CREATE TABLE IF NOT EXISTS [canonical_transcript_table] (
        chromosome_name VARCHAR,
        ensembl_gene_id VARCHAR,
        transcript_id VARCHAR,
        is_canonical VARCHAR
    )
""")

# Base path for chromosome files with Ensembl transcript information
base_path = '/file/path/to/your/ensembl/files/'
# List of chromosome names in the files
chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
# Define chunk size for batch processing
chunksize = 5000  # Adjust based on your system's memory capacity and file size

# Iterate over all transcript ids for chromosomes

for chr_name in chromosomes:
    # Update the file path for each chromosome
    file_path = f'{base_path}chromosome_{chr_name}.tsv'
    
    print(f"Processing chromosome {chr_name}...")  # Indicate the start of processing a new chromosome

    # Process the TSV file in chunks
    for chunk_id, df_chunk in enumerate(pd.read_csv(file_path, sep='\t', chunksize=chunksize)):
        # Filter rows where 'is_canonical' equals 1 to choose the canonical transcript IDs
        canonical_df = df_chunk[df_chunk['is_canonical'] == 1].copy()

        # Select only the necessary columns (assuming these are the column names in your TSV)
        canonical_df = canonical_df[['chromosome_name', 'ensembl_gene_id', 'transcript_id', 'is_canonical']]

        # Register the DataFrame as a temporary view in DuckDB
        temp_view_name = f"temp_canonical_{chr_name}_{chunk_id}"
        conn.register(temp_view_name, canonical_df)
        
        # Insert the data from the temporary view into the [canonical_transcripts_table], adjust table name is needed
        conn.execute(f"""
            INSERT INTO [canonical_transcript_table]
            SELECT * FROM {temp_view_name}
        """)
        
        # Unregister the temporary view to free up resources
        conn.unregister(temp_view_name)
    
    # After successfully processing and inserting each chromosome's data
    print(f"Successfully inserted canonical transcripts for chromosome {chr_name}.")

# Close the DuckDB connection
conn.close()

print("Completed inserting canonical transcripts for all chromosomes.")
