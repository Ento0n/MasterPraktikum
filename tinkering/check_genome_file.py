

with open("data/GCA_000001405.29_GRCh38.p14_genomic.fna") as f:
    counter = 0
    for line in f:
        if line.startswith(">"):
            print(counter)
            print(line.strip())
            counter = 0
        else:
            for character in line.strip():
                counter += 1
