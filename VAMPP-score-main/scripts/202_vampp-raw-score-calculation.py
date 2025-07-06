# Import necessary modules

import json
import gzip
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

# Create a function to load the calculated gene scores for ISPP performance

def load_genescores_json(genescores_f):
    print("Loading gene scores from JSON file...")
    with open(genescores_f, "r") as f:
        gene_scores = json.load(f)
    print("Gene scores loaded successfully.")
    return gene_scores

# Create a function to calculate the VAMPP-score, utilizing the loaded gene scores iterating over dbNSFP chromosome data

def metascore_dbnsfp(dbnsfp_f, metascore_tsv, all_genescores):

    with gzip.open(dbnsfp_f, 'rt', encoding='utf-8') as f:
        header = f.readline().strip('\n').split('\t')
        indexed_header = list(enumerate(header))
        rankscore_indices_cols = list(filter(lambda x: '_rankscore' in x[1], indexed_header))
        vamppscore_indices_cols = list(map(lambda y: (y[0], y[1].rsplit('_', 1)[0] + '_vampp') , rankscore_indices_cols))
        genename_col = list(filter(lambda x: 'genename' == x[1], indexed_header))[0]
        rank_and_vamp_scores_header_str = '\t'.join([a + '\t' + b  for a, b in zip([i[1] for i in rankscore_indices_cols], [i[1] for i in vamppscore_indices_cols])])

        def get_score(l_score, genescore):
            if not (l_score and genescore) or l_score == '.':
                return None
            else:
                # try:
                metascore = float(l_score) * genescore
                # except:
                #     print((l_score, genescore))
                return str(metascore)

        with gzip.open(metascore_tsv, 'wt', encoding='utf-8') as o:
            o.write(f'#chr\tpos\tref\talt\tGene\tVAMPP_score\t{rank_and_vamp_scores_header_str}\n')

            while l := f.readline():
                l = l.strip('\n').split('\t')
                try:
                    curr_gene = l[genename_col[0]].split(';')[0]
                    #print(curr_gene)
                    curr_gene_scores = all_genescores[curr_gene]
                    # print(curr_gene_scores)
                    positive_curr_gene_scores = dict(filter(lambda g: (g[0],g[1]) if g[1] and g[1] > 0 else None, curr_gene_scores.items()))
                except (IndexError, KeyError):
                    continue
                curr_l_scores = tuple(map(lambda s: (s[1], l[s[0]]), rankscore_indices_cols)) # (scorename, score)
                curr_l_metascores = tuple(map(lambda m: get_score(m[1], positive_curr_gene_scores.get(m[0], None)), curr_l_scores))
                # curr_l_metascores_str = map(lambda s: s if s else '', curr_l_metascores)
                curr_l_rank_and_vampp_scores_str = '\t'. join(map(lambda y: '\t'.join((str(y[0]), str(y[1]) if y[1] else '.')), zip([i[1] for i in curr_l_scores], curr_l_metascores)))
                exist_curr_l_metascores = tuple(map(lambda z: float(z), filter(lambda s: s if s else None, curr_l_metascores)))
                try:
                    curr_l_metascore = sum(exist_curr_l_metascores)/len(exist_curr_l_metascores)
                except ZeroDivisionError:
                    curr_l_metascore = ''
                if curr_l_metascore == '':
                    continue
                out_l = '\t'.join(l[0:4]+[curr_gene]+[str(curr_l_metascore)]+[curr_l_rank_and_vampp_scores_str])+'\n'
                o.write(out_l)
    print(f"Finished processing for {dbnsfp_f}. Data written to {metascore_tsv}")  # Log finish

def execute_dbnsfp(input_dir, output_dir, all_genescores):
    chr_s = list(range(1, 23)) + ["X", "Y"]
    with ThreadPoolExecutor(24) as executor:
        futures = []
        for c in chr_s:
            future = executor.submit(
                metascore_dbnsfp,
                str(Path(input_dir) / f"dbNSFP4.7a_variant.chr{c}.gz"),
                str(Path(output_dir) / f"vamppscore_dbNSFP4.7a_variant.chr{c}.tsv.gz"), #adjust the name if needed with another prefix
                all_genescores
            )
            futures.append(future)

        # Wait for all futures and log as they complete
        for future in futures:
            future.result()  # This will raise exceptions from the threads if any and ensure completion

# Define the needed input and output file paths

def main():
    gene_scores_file = "/path/to/genescore/files/[all_genescores].json"  # Adjust the path as needed
    input_directory = "/path/to/extracted/dbnsfp/chromosome/files"  # Adjust as necessary
    output_directory = "/path/to/output/files/vamppscore/by/chromosomes"  # Adjust as necessary

    # Load gene scores
    gene_scores_dict = load_genescores_json(gene_scores_file)

    # Process each chromosome's DBNSFP file
    execute_dbnsfp(input_directory, output_directory, gene_scores_dict)

    print("Processing complete. All metascores have been calculated and saved.")

if __name__ == "__main__":
    main()
