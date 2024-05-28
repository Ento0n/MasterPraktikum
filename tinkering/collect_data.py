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
gene_id_gene_name_mapping = {}
with open("data/gene_info") as f:
    for line in tqdm(f):
        if not line.startswith("9606"):
            continue

        gene_name = line.split("\t")[2]
        chromosome = line.split("\t")[6]
        gene_id = line.split("\t")[1]

        if gene_name in data.keys():
            data[gene_name]["chromosome"] = chromosome
            data[gene_name]["gene_id"] = gene_id
            gene_id_gene_name_mapping[gene_id] = gene_name

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

# extract coding sequences (CDS) from .gff-file


def extract_gene_id_gff(line: str) -> str:
    gene_id = None

    attributes = line.split("\t")[8]

    for item in attributes.split(";"):
        if item.split("=")[0] == "Dbxref":
            references = item.split("=")[1]

            for item2 in references.split(","):
                if item2.split(":")[0] == "GeneID":
                    gene_id = item2.split(":")[1]
                    break

            # break outer for loop
            break
    return gene_id

with open("data/GCF_000001405.40_GRCh38.p14_genomic.gff") as f:
    CDSs = {}
    as_events = {}
    old_gene_id = ""
    old_gene_name = ""
    old_start = ""
    i = 0
    j = 0
    for line in tqdm(f):
        if line.startswith("#"):
            continue

        if line.split("\t")[2] == "pseudogene":
            gene_id = extract_gene_id_gff(line)

            if gene_id not in gene_id_gene_name_mapping.keys():
                continue

            gene_name = gene_id_gene_name_mapping[gene_id]

            data[gene_name]["pseudogene"] = True

        if line.split("\t")[2] == "CDS":
            gene_id = extract_gene_id_gff(line)

            if gene_id not in gene_id_gene_name_mapping.keys():
                continue

            gene_name = gene_id_gene_name_mapping[gene_id]

            data[gene_name]["pseudogene"] = False

            # clear the record of CDSs in case new gene is handled
            if gene_id != old_gene_id:
                data[old_gene_name]["CDSs"] = as_events.copy()
                as_events.clear()
                i = 0

            # extract start
            start = line.split("\t")[3]

            # new alternative splicing event, new protein
            if gene_id == old_gene_id and start < old_start:
                as_events[i] = CDSs.copy()
                CDSs.clear()
                j = 0
                i += 1

            # extract information
            stop = line.split("\t")[4]
            strand = line.split("\t")[6]
            phase = line.split("\t")[7]

            CDSs[j] = {"start": start, "stop": stop, "strand": strand, "phase": phase}

            # remember old gene id and increme nt CDS counter
            old_gene_id = gene_id
            old_gene_name = gene_name
            old_start = start
            j += 1

# convert to pandas dataframe
df = pd.DataFrame.from_dict(data, orient="index")
df.replace(to_replace="", value="NA")
df.index.name = "gene_name"
df.to_csv("output/uniprot_genebank_homo_sapiens.tsv", sep="\t")
