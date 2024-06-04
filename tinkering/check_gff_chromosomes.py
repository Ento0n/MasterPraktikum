

with open("data/GCF_000001405.40_GRCh38.p14_genomic.gff") as f:
    for line in f:
        if line.startswith("##sequence-region"):
            print(line)