import pandas as pd
from tqdm import tqdm
import numpy as np

# get gene names from uniprot

data = {}
with open("data/uniprot_gene3d_human.tsv") as f:
    for line in tqdm(f):
        if line.startswith("Entry"):
            continue

        uniprot_gene_name = line.split("\t")[4].split(" ")[0]
        cath_superfamily = line.split("\t")[7]
        pdb_id = line.split("\t")[8].strip()

        if uniprot_gene_name.endswith(";"):
            uniprot_gene_name = uniprot_gene_name[:-1]

        data[uniprot_gene_name] = {"cath_superfamily": cath_superfamily, "pdb_id": pdb_id}

# extract chromosomes from gene_info
gene_info_homo_sapiens_gene_ids = set()
with open("data/gene_info") as f:
    for line in tqdm(f):
        if not line.startswith("9606"):
            continue

        gene_name = line.split("\t")[2]
        chromosome = line.split("\t")[6]

        if gene_name in data.keys():
            data[gene_name]["chromosome"] = chromosome

# extract genomic locations from gene2accession

with open("data/gene2accession") as f:
    genome_accessions = {}
    old_gene_name = ""
    for line in tqdm(f):
        if not line.startswith("9606"):
            continue

        status = line.split("\t")[2]
        if status not in ["VALIDATED", "REVIEWED"]:
            continue

        gene_name = line.split("\t")[15].strip()

        # skip entries that cannot be mapped
        if gene_name not in data.keys():
            continue

        # clear the record of genomic accesions in case new gene is handled
        if gene_name != old_gene_name:
            data[old_gene_name]["genome_accessions"] = genome_accessions.copy()
            genome_accessions.clear()

        genome_accession = line.split("\t")[12]
        start = line.split("\t")[9]
        end = line.split("\t")[10]

        genome_accessions[genome_accession] = f"{start}:{end}"

        # remember old gene name for saving all genome accesions in 1 dictionary
        old_gene_name = gene_name

# convert to pandas dataframe
df = pd.DataFrame.from_dict(data, orient="index")
df.replace(to_replace="", value="NA")
df.index.name = "gene_name"
df.to_csv("output/uniprot_genebank_homo_sapiens.tsv", sep="\t")
