
gene2accession_homo_sapiens_gene_ids = []
tmp_id = ""
with open("data/gene2accession") as f:
    for line in f:
        if line.startswith("#"):
            print(line)

        if not line.startswith("9606"):
            continue

        tmp_id = line.split("\t")[1]

        print(line)

        line = f.readline()

        while line.split("\t")[1] == tmp_id:
            print(line)
            line = f.readline()

        # gene2accession_homo_sapiens_gene_ids.append(line.split("\t")[1])
        break


gene_info_homo_sapiens_gene_ids = []
with open("data/gene_info") as f:
    for line in f:
        if line.startswith("#"):
            print(line)

        if not line.startswith("9606"):
            continue

        gene_info_homo_sapiens_gene_ids.append(line.split("\t")[1])
        print(line)
        break
