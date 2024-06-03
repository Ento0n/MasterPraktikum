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


chromosome_to_head_mapping = {"1": "CM000663.2", "2": "CM000664.2", "3": "CM000665.2", "4": "CM000666.2", "5": "CM000667.2",
                              "6": "CM000668.2", "7": "CM000668.2", "8": "CM000670.2", "9": "CM000671.2", "10": "CM000672.2",
                              "11": "CM000673.2", "12": "CM000674.2", "13": "CM000675.2", "14": "CM000676.2",
                              "15": "CM000677.2", "16": "CM000678.2", "17": "CM000679.2", "18": "CM000680.2",
                              "19": "CM000680.2", "20": "CM000682.2", "21": "CM000683.2", "22": "CM000684.2",
                              "X": "CM000685.2", "Y": "CM000686.2", "MT": "GL000209.2"}
sequences = list(SeqIO.parse("data/GCA_000001405.29_GRCh38.p14_genomic.fna", "fasta"))
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
    'TGA': '*',    # Stop
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


with open("data/GCF_000001405.40_GRCh38.p14_genomic.gff") as f:
    CDSs = {}
    extracted_sequences = {}
    as_events = {}
    aa_sequences = {}
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
            strand = line.split("\t")[6]
            data[gene_name]["strand"] = strand

            # clear the record of CDSs in case new gene is handled
            if gene_id != old_gene_id:
                # add last collected CDSs
                as_events[i] = CDSs.copy()
                as_events[str(i) + "a"] = extracted_sequences.copy()
                as_events[str(i) + "b"] = aa_sequences.copy()
                extracted_sequences.clear()
                aa_sequences.clear()
                CDSs.clear()
                j = 0

                # add counter for coding sequences
                cds_counter = []
                for i in range(int(len(as_events) / 3)): # divide by 3 because of 0a & 0b and so on
                    cds_counter.append(len(as_events[i]))
                data[old_gene_name]["cds_counter"] = cds_counter.copy()
                cds_counter.clear()

                data[old_gene_name]["CDSs"] = as_events.copy()
                as_events.clear()
                i = 0

            # extract start
            start = int(line.split("\t")[3])

            # in case of multiple alternative sequencing events, add last CDSs collection to as_events
            if strand == "+":
                if gene_id == old_gene_id and start < old_start:
                    as_events[i] = CDSs.copy()
                    as_events[str(i) + "a"] = extracted_sequences.copy()
                    as_events[str(i) + "b"] = aa_sequences.copy()
                    extracted_sequences.clear()
                    aa_sequences.clear()
                    CDSs.clear()
                    j = 0
                    i += 1
            else:
                if gene_id == old_gene_id and start > old_start:
                    as_events[i] = CDSs.copy()
                    as_events[str(i) + "a"] = extracted_sequences.copy()
                    as_events[str(i) + "b"] = aa_sequences.copy()
                    extracted_sequences.clear()
                    aa_sequences.clear()
                    CDSs.clear()
                    j = 0
                    i += 1

            # extract information
            stop = int(line.split("\t")[4])
            frame = int(line.split("\t")[7])

            CDSs[j] = {"start": start, "stop": stop, "frame": frame}

            # filter out chromosome names I cannot map
            if data[gene_name]["chromosome"] != "-":
                if data[gene_name]["chromosome"] != "X|Y":
                    # extract sequence out of genome fasta file
                    for record in sequences:
                        if record.id == chromosome_to_head_mapping[data[gene_name]["chromosome"]]:
                            region_seq = record.seq[start-1:stop]
                            extracted_sequences[j] = str(region_seq)

                            # convert extracted sequence to protein sequence
                            aa_sequences[j] = str(convert_nuc2aa(region_seq, frame))
                            break

            # remember old gene id and increment CDS counter
            old_gene_id = gene_id
            old_gene_name = gene_name
            old_start = start
            j += 1



# convert to pandas dataframe
df = pd.DataFrame.from_dict(data, orient="index")
df.replace(to_replace="", value="NA")
df.index.name = "gene_name"
df.to_csv("output/uniprot_genebank_homo_sapiens.tsv", sep="\t")
