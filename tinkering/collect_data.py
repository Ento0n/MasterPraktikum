import pandas as pd
from tqdm import tqdm
from Bio import SeqIO

# get gene names from uniprot
data = {}
with open("data/uniprot_gene3d_human.tsv") as f:
    for line in tqdm(f):
        if line.startswith("Entry"):
            continue

        uniprot_gene_name = line.split("\t")[4].split(" ")[0]
        cath_superfamily = line.split("\t")[7]
        pdb_id = line.split("\t")[8].strip()
        uniprot_sequence = line.split("\t")[9].strip()

        if uniprot_gene_name.endswith(";"):
            uniprot_gene_name = uniprot_gene_name[:-1]

        data[uniprot_gene_name] = {"cath_superfamily": cath_superfamily, "pdb_id": pdb_id, "uniprot_sequence": uniprot_sequence}

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

        # remember old gene name for saving all genome accessions in 1 dictionary
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


def convert_nuc2aa(dna: str, frame: int):
    # size of slices
    n = 3

    # make sure all nucleotides are upper case letters
    dna = dna.upper()

    slices = [dna[i:i + n] for i in range(frame, len(dna), n)]

    prot_seq = ""
    for piece in slices:
        if piece in codon_table.keys():
            prot_seq = prot_seq + codon_table[piece]
        else:
            prot_seq = prot_seq + "-"

    return prot_seq


def reverse_complement(dna: str):
    reverse_mapping = {"A": "T", "T": "A", "G": "C", "C": "G",
                       "a": "t", "t": "a", "g": "c", "c": "g", "N": "N", "n": "n"}
    reverse_dna = ""

    # turn around DNA
    dna = dna[::-1]

    for character in dna:
        if character in reverse_mapping.keys():
            reverse_dna = reverse_dna + reverse_mapping[character]
        else:
            raise Exception(f"Bad character in DNA!: {character}")

    return reverse_dna


# translation of chromosome number to RefSeq chromosome ID
chromosome_to_head_mapping = {"1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12", "4": "NC_000004.12",
                              "5": "NC_000005.10", "6": "NC_000006.12", "7": "NC_000007.14", "8": "NC_000008.11",
                              "9": "NC_000009.12", "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
                              "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10", "16": "NC_000016.10",
                              "17": "NC_000017.11", "18": "NC_000018.10", "19": "NC_000019.10", "20": "NC_000020.11",
                              "21": "NC_000021.9", "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10",
                              "MT": "NC_012920.1"}

sequences = list(SeqIO.parse("data/genomes/homo_sapiens/GCF_000001405.40_GRCh38.p14_genomic.fna", "fasta"))
codon_table = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': 'U',    # Stop (Opal), But also a U in some cases
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}


with open("data/genomes/homo_sapiens/GCF_000001405.40_GRCh38.p14_genomic.gff") as f:
    CDSs = {}
    cds_chain = False
    extracted_sequences = {}
    as_events = {}
    old_gene_id = ""
    old_gene_name = ""
    old_start = ""
    old_stop = ""
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

            cds_chain = False

        elif line.split("\t")[2] == "CDS":
            gene_id = extract_gene_id_gff(line)

            if gene_id not in gene_id_gene_name_mapping.keys():
                continue

            gene_name = gene_id_gene_name_mapping[gene_id]

            # in case an entry is already present, skip this one...
            if "CDSs" in data[gene_name].keys():
                continue

            data[gene_name]["pseudogene"] = False
            strand = line.split("\t")[6]
            data[gene_name]["strand"] = strand
            sequence_region = line.split("\t")[0]
            data[gene_name]["gff_sequence_region"] = sequence_region

            # clear the record of CDSs in case new gene is handled
            if gene_id != old_gene_id:
                # add last collected CDSs
                as_events[i] = CDSs.copy()
                as_events[str(i) + "a"] = extracted_sequences.copy()
                extracted_sequences.clear()
                CDSs.clear()
                j = 0

                # add counter for coding sequences
                cds_counter = []
                for i in range(int(len(as_events) / 2)):  # divide by 2 because of 0a and so on
                    cds_counter.append(len(as_events[i]))
                data[old_gene_name]["cds_counter"] = cds_counter.copy()
                cds_counter.clear()

                # Add protein sequences as new column of df
                pro_seqs = []
                for i in range(int(len(as_events) / 2)):
                    full_seq = ""
                    for seq in as_events[str(i) + "a"].values():
                        seq = str(seq)  # not necessary needed, but pycharm gives warning otherwise
                        full_seq = full_seq + seq

                    # check whether there is even 1 CDS entry
                    if full_seq != "":
                        pro_seq = convert_nuc2aa(full_seq, as_events[i][0]["frame"])

                        if pro_seq.endswith("*") or pro_seq.endswith("U") or pro_seq.endswith("-"):
                            pro_seq = pro_seq[:-1]

                        pro_seqs.append(pro_seq)

                    else:
                        pro_seqs.append("-")

                if len(pro_seqs) > 1:
                    data[old_gene_name]["protein_sequences"] = ",".join(pro_seqs)
                elif len(pro_seqs) == 1:
                    data[old_gene_name]["protein_sequences"] = pro_seqs[0]
                else:
                    data[old_gene_name]["protein_sequences"] = "-"

                # add all collected alternative splicing events to the collected data
                data[old_gene_name]["CDSs"] = as_events.copy()
                as_events.clear()
                i = 0

            # extract start
            start = int(line.split("\t")[3])
            stop = int(line.split("\t")[4])

            # in case of multiple alternative sequencing events, add last CDSs collection to as_events
            if not cds_chain and gene_id == old_gene_id:
                as_events[i] = CDSs.copy()
                as_events[str(i) + "a"] = extracted_sequences.copy()
                extracted_sequences.clear()
                CDSs.clear()
                j = 0
                i += 1

            # extract information
            frame = int(line.split("\t")[7])

            CDSs[j] = {"start": start, "stop": stop, "frame": frame}

            # filter out chromosome names I cannot map
            if data[gene_name]["chromosome"] != "-":
                if data[gene_name]["chromosome"] != "X|Y":
                    # extract sequence out of genome fasta file
                    for record in sequences:
                        if record.id == chromosome_to_head_mapping[data[gene_name]["chromosome"]]:
                            region_seq = record.seq[start-1:stop]

                            # Turn around sequence for reverse strand
                            if strand == "+":
                                extracted_sequences[j] = str(region_seq)
                            else:
                                extracted_sequences[j] = reverse_complement(str(region_seq))

                            # record id found end loop
                            break

            # remember old gene id and increment CDS counter, also indicated that cds has been processed
            cds_chain = True
            old_gene_id = gene_id
            old_gene_name = gene_name
            old_start = start
            old_stop = stop
            j += 1

        else:
            cds_chain = False


# convert to pandas dataframe
df = pd.DataFrame.from_dict(data, orient="index")
df.replace(to_replace="", value="NA")
df.index.name = "gene_name"
df.to_csv("output/uniprot_genebank_homo_sapiens.tsv", sep="\t")
