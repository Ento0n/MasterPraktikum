from Bio import SeqIO

sequences = list(SeqIO.parse("data/GCA_000001405.29_GRCh38.p14_genomic.fna", "fasta"))

for record in sequences:
    if record.id == "CM000663.2":
        region_seq = record.seq[46760:46760]
        print(f"print: -{region_seq}-")
        break


"""
with open("data/GCA_000001405.29_GRCh38.p14_genomic.fna") as f:
    seq = ""
    flag = False
    start = 22922594
    for line in f:
        if line.startswith(">CM000684.2"):
            while True:
                line = f.readline().strip()

                for character in line:
                    start = start - 1

                    if start == 0:
                        flag = True

                    if flag:
                        seq = seq + character

                if flag:
                    break

        if flag:
            break

    print(seq)

"""


