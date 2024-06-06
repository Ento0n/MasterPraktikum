import pandas as pd
import numpy as np

df = pd.read_csv("output/uniprot_genebank_homo_sapiens.tsv", sep="\t", index_col=0)

wrong_sequences = {}

for i, row in df.iterrows():
    flag = False
    uniprot_seq = row["uniprot_sequence"]

    if row["protein_sequences"] is np.nan:
        continue

    for sequence in row["protein_sequences"].split(","):
        if sequence[:-1] == uniprot_seq:
            flag = True

    if not flag:
        wrong_sequences[uniprot_seq] = row["protein_sequences"]

print(len(wrong_sequences))
