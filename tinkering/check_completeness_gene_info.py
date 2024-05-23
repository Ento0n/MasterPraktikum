
gene_info_homo_sapiens_gene_ids = set()
with open("data/gene_info") as f:
    for line in f:
        if not line.startswith("9606"):
            continue

        gene_info_homo_sapiens_gene_ids.add(line.split("\t")[2])

uniprot_gene_names = set()
with open("data/uniprot_gene3d_human.tsv") as f:
    for line in f:
        uniprot_gene_name = line.split("\t")[4].split(" ")[0]

        if uniprot_gene_name.endswith(";"):
            uniprot_gene_name = uniprot_gene_name[:-1]

        uniprot_gene_names.add(uniprot_gene_name)

print(uniprot_gene_names.difference(gene_info_homo_sapiens_gene_ids))
print(len(uniprot_gene_names.difference(gene_info_homo_sapiens_gene_ids)))

