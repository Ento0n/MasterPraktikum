from Bio import SeqIO

sequences = list(SeqIO.parse("data/GCA_000001405.29_GRCh38.p14_genomic.fna", "fasta"))

for record in sequences:
    if record.id == "CM000663.2":
        region_seq = record.seq[46758:46760]
        print(region_seq)
        break
