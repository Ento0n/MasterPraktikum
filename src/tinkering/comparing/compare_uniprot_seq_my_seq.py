import pandas as pd
import numpy as np

df = pd.read_csv("output/uniprot_genebank_homo_sapiens.tsv", sep="\t", index_col=0)

wrong_sequences = set()
counter = 0
gene_names_correct = set()
first_three_bases_correct = set()

for i, row in df.iterrows():
    correct_flag = False
    first_3_flag = False
    uniprot_seq = row["uniprot_sequence"]

    if row["protein_sequences"] is np.nan:
        continue

    strand = row["strand"]

    counter += 1

    for sequence in row["protein_sequences"].split(","):
        if sequence == uniprot_seq:
            correct_flag = True
        elif sequence[:3] == uniprot_seq[:3] and strand == "+":
            first_3_flag = True

    if correct_flag:
        gene_names_correct.add(i)
    elif first_3_flag:
        first_three_bases_correct.add(i)
    else:
        wrong_sequences.add(i)

print(f"sequences in general: {counter}")
print(f"correct sequences: {len(gene_names_correct)}")
# print(f"gene names correct sequences: {gene_names_correct}")
print(f"wrong sequences: {len(wrong_sequences)}")
print(f"wrong sequences: {list(wrong_sequences)[:10]}")
print(f"only first 3 match, rest uncertain: {first_three_bases_correct}")
print(f"length of first 3 match: {len(first_three_bases_correct)}")
