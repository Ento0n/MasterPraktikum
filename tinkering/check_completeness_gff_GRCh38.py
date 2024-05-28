
uniprot2genbank_gene_ids = set()
with open("output/uniprot_genebank_homo_sapiens.tsv") as f:
    for line in f:
        gene_id = line.split("\t")[4]
        uniprot2genbank_gene_ids.add(gene_id)

types = set()
refseq_gene_ids = set()
with open("data/GCF_000001405.40_GRCh38.p14_genomic.gff") as f:
    for line in f:
        if line.startswith("#"):
            continue

        types.add(line.split("\t")[2])

        if line.split("\t")[2] == "CDS" or line.split("\t")[2] == "pseudogene":
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

        refseq_gene_ids.add(gene_id)

print(uniprot2genbank_gene_ids.difference(refseq_gene_ids))
print(len(uniprot2genbank_gene_ids.difference(refseq_gene_ids)))
print(types)

# check whether gene names is a better connection
"""
uniprot2genbank_gene_names = set()
with open("output/uniprot_genebank_homo_sapiens.tsv") as f:
    for line in f:
        gene_name = line.split("\t")[0]
        uniprot2genbank_gene_names.add(gene_name)

refseq_gene_names = set()
with open("data/GCF_000001405.40_GRCh38.p14_genomic.gff") as f:
    for line in f:
        if line.startswith("#"):
            continue

        if line.split("\t")[2] == "CDS" or line.split("\t")[2] == "pseudogene":
            attributes = line.split("\t")[8]

            for item in attributes.split(";"):
                if item.split("=")[0] == "gene":
                    gene_name = item.split("=")[1].strip()

        refseq_gene_names.add(gene_name)

print(uniprot2genbank_gene_names.difference(refseq_gene_names))
print(len(uniprot2genbank_gene_names.difference(refseq_gene_names)))
"""
# gene id is still a better connection between the files, 46 compared to 64
