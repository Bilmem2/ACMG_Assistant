# Import necessary modules

import pandas as pd
import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

# Define the list of chromosomes to iterate over
chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']

# Read desired column names from a file, adjusting for special characters
wanted_columns = []
with open('file/path/to/your/wanted/column/info/[wanted_columns].txt', 'rt') as file: #adjust as needed
    for line in file:
        column_name = line.strip()
        if '-' in column_name or '#' in column_name or '+' in column_name:
            column_name = f'"{column_name}"'
        wanted_columns.append(column_name)

# Establish a database connection
conn = duckdb.connect(database='file/path/to/your/duckdb/[database_name].duckdb')

# Create a temporary table for ClinVar variants - can be accessed through out 'clinvar-data' link (only involves missense variants on GRCh38)
conn.execute("""
    CREATE TEMPORARY TABLE temp_clinvar AS
    SELECT * FROM read_csv_auto('file/path/to/your/filtered/missense-variants/from/clinvar/[variant_file].tsv');
""")

# Initialize the output file path for the Parquet writer, you can choose a different file format, then adjust accordingly
output_parquet_file = 'file/path/to/your/output/file/[file_name].parquet'

# Placeholder for schema definition (will be set upon fetching the first batch)
schema = None

for chromosome in chromosomes:
    print(f"Processing chromosome {chromosome}")
    table_name = f"dbNSFP4_7a_chr{chromosome}"

    # Generate the SQL query for the current chromosome, adjusting column names as necessary
    cast_columns = ", ".join([f"CAST({table_name}.{col} AS VARCHAR) AS {col}" for col in wanted_columns])
    query = f"""
        SELECT {cast_columns}, temp_clinvar.*
        FROM {table_name}
        INNER JOIN [inserted_canonical_transcriptinfo_table] ON {table_name}.Ensembl_transcriptid = [inserted_canonical_transcriptinfo_table].transcript_id
        INNER JOIN temp_clinvar ON {table_name}."#chr" = temp_clinvar.Chromosome
        AND {table_name}."pos(1-based)" = temp_clinvar.PositionVCF
        AND {table_name}.ref = temp_clinvar.ReferenceAlleleVCF
        AND {table_name}.alt = temp_clinvar.AlternateAlleleVCF
    """
    
    # Execute the query and fetch the data
    df = conn.execute(query).fetchdf()
    
    if df.empty:
        print(f"No data found for chromosome {chromosome}. Continuing to the next chromosome.")
        continue

    if schema is None:
        # Define the schema based on the first batch of data fetched
        schema = pa.Table.from_pandas(df).schema
        writer = pq.ParquetWriter(output_parquet_file, schema=schema, compression='snappy')
    
    # Convert the DataFrame to a PyArrow Table and write (append) it to the Parquet file
    table = pa.Table.from_pandas(df, schema=schema)
    writer.write_table(table)
    print(f"Data for chromosome {chromosome} has been written to the Parquet file.")

# Close the ParquetWriter to finalize the Parquet file
writer.close()

# Close the database connection
conn.close()

print("Successfully processed all chromosomes and written to the Parquet file.")
