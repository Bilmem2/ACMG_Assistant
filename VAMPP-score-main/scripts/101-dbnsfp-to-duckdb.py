# import necessary modules

import duckdb
import gzip
import pandas as pd
import numpy as np

# Create and connect to your DuckDB, you can manually change settings for DuckDB execution commands accordingly to your system

conn = duckdb.connect('file/path/to/your/database/yourdatabase.duckdb')
conn.execute('SET enable_progress_bar=FALSE')
conn.execute("SET GLOBAL pandas_analyze_sample=100000")
duckdb.default_connection.execute("SET GLOBAL pandas_analyze_sample=100000")

# Create a variable for further use in to convert potentially unsupported characters

dot_to_none = lambda z: None if z == '.' else z

# Create a function to split Ensembl transcript IDs (which will be further used in filtering variants for ISPP performance detection)
# by ';' in the specific column (idx=14 in the dbNSFP4.7A version),
# creating new lines for each with the specific information

def split_by_transcripts(line_data, transcript_col_idx=14):

    transcript_count = len(line_data[transcript_col_idx].split(';'))

    return zip(*map(lambda y: y*transcript_count if len(y) < transcript_count else y[:transcript_count], map(lambda x: x.split(';'), line_data)))

# Read the header to from any of the chrosome files from dbNSFP for table information

with gzip.open('file/path/to/extracted/dbNSFP/chromosome/files/dbNSFP4.7a_variant.chr[chromosome_number].gz', 'rt') as f:
    header = f.readline().strip().split('\t')

# Generate the CREATE TABLE SQL statement with all columns as VARCHAR
column_definitions = ", ".join([f'"{col}" VARCHAR' for col in header])
create_table_sql = f"CREATE TABLE IF NOT EXISTS [table_name] ({column_definitions});"
conn.execute(create_table_sql)

# Adjusted batch sizes, you can set according to your system
first_commit = True
batch_size = 5000  # New batch size for subsequent commits
initial_batch_size = 10000  # Initial batch size for the first commit

# Insert chromosome information from dbNSFP to the DuckDB table

for _chr in list(range(1,23))+["X", "Y"]:

    with gzip.open(f'file/path/to/extracted/dbNSFP/chromosome/files/dbNSFP4.7a_variant.chr{_chr}.gz', 'rt') as dbnsfp_gz_f:
        header = dbnsfp_gz_f.readline().strip().split('\t')
        transcript_col_idx = header.index("Ensembl_transcriptid")
        processed_lines = []

        while line := dbnsfp_gz_f.readline():
           # Adjust condition to check against new batch sizes
            if not len(processed_lines) > (initial_batch_size if first_commit else batch_size):
                line = line.strip().split('\t')
                if ';' in line[14]:
                    processed_lines += list(map(lambda y: tuple(map(dot_to_none, y)), split_by_transcripts(line)))
                else:
                    processed_lines.append(tuple(map(lambda z: None if z == '.' else z, line)))
            else:
                processed_df = pd.DataFrame(processed_lines, columns=header)
                processed_df = processed_df.replace('.', np.nan)
                if first_commit:
                    conn.execute('INSERT INTO [table_name] SELECT * FROM processed_df')
                    first_commit = False
                else:
                    conn.execute('INSERT INTO [table_name] SELECT * FROM processed_df')
                processed_lines = []
                processed_df = None
        processed_df = pd.DataFrame(processed_lines, columns=header)
        processed_df = processed_df.replace('.', np.nan)
        conn.execute('INSERT INTO [table_name] SELECT * FROM processed_df')
