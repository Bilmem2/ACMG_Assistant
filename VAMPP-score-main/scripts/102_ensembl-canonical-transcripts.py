
# Transcript information is fetched from Ensembl for each human gene with chromosome, gene, position and transcript info.

#Import the necessary module for data fetching

import requests

def fetch_chromosome_lengths():
    url = "https://rest.ensembl.org/info/assembly/homo_sapiens?content-type=application/json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        chromosome_lengths = {item['name']: item['length'] for item in data['top_level_region'] if item['name'] in [str(i) for i in range(1, 23)] + ['X', 'Y']}
        return chromosome_lengths
    else:
        print("Failed to retrieve chromosome lengths.")
        return {}

def fetch_genes_transcripts(chromosome, start, end):
    genes_url = f"https://rest.ensembl.org/overlap/region/human/{chromosome}:{start}-{end}?feature=gene;content-type=application/json"
    genes_response = requests.get(genes_url, headers={"Content-Type": "application/json"})
    if genes_response.status_code == 200:
        genes_data = genes_response.json()
        for gene in genes_data:
            transcripts_url = f"https://rest.ensembl.org/lookup/id/{gene['id']}?expand=1;content-type=application/json"
            transcripts_response = requests.get(transcripts_url, headers={"Content-Type": "application/json"})
            if transcripts_response.status_code == 200:
                transcripts_data = transcripts_response.json().get('Transcript', [])
                gene['transcripts'] = [{
                    'transcript_id': transcript['id'],
                    'is_canonical': transcript.get('is_canonical', 0)
                } for transcript in transcripts_data]
        return genes_data
    else:
        print(f"Failed to fetch data for chromosome {chromosome} segment {start}-{end}.")
        return []

def write_genes_to_file(chromosome, genes_data):
    filename = f"path/to/your/output/file/chromosome_{chromosome}.tsv" # adjust as needed
    with open(filename, 'w') as file:
        file.write("ensembl_gene_id\thgnc_symbol\tchromosome_name\tstart_position\tend_position\ttranscript_id\tis_canonical\n")
        for gene in genes_data:
            for transcript in gene.get('transcripts', []):
                file.write(f"{gene.get('id', 'N/A')}\t{gene.get('external_name', 'N/A')}\t{chromosome}\t{gene.get('start', 'N/A')}\t{gene.get('end', 'N/A')}\t{transcript.get('transcript_id', 'N/A')}\t{transcript.get('is_canonical', 'N/A')}\n")

chromosome_lengths = fetch_chromosome_lengths()
segment_size = 5_000_000

for chromosome, length in chromosome_lengths.items():
    print(f"Processing chromosome: {chromosome}")
    all_genes = []
    for start in range(1, length, segment_size):
        end = min(start + segment_size - 1, length)
        segment_genes = fetch_genes_transcripts(chromosome, start, end)
        all_genes.extend(segment_genes)
    
    write_genes_to_file(chromosome, all_genes)
    print(f"Data for chromosome {chromosome} written to file.")
