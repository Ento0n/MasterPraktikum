
uniprot2genbank_gene_ids = set()
with open("output/uniprot_genebank_homo_sapiens.tsv") as f:
    for line in f:
        gene_id = line.split("\t")[4]
        uniprot2genbank_gene_ids.add(gene_id)

refseq_gene_ids = set()
with open("data/GCF_000001405.40_GRCh38.p14_genomic.gff") as f:
    for line in f:
        if line.startswith("#"):
            continue

        if line.split("\t")[2] == "CDS":
            attributes = line.split("\t")[8]

            for item in line.split(";"):
                if item.split("=")[0] == "Dbxref":
                    references = item.split("=")[1]

                    for item2 in references.split(","):
                        if item2.split(":")[0] == "GeneID":
                            gene_id = item2.split(":")[1]
                            break

        refseq_gene_ids.add(gene_id)

print(uniprot2genbank_gene_ids.difference(refseq_gene_ids))
print(len(uniprot2genbank_gene_ids.difference(refseq_gene_ids)))
