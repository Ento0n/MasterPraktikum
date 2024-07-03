from typing import List

import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import argparse


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


# get gene names, cath supaerfamily, pdb ID and protein sequence from UniProt
def read_uniprot_tsv(path: str, data_fct: dict):
    with open(path) as f:
        for line in tqdm(f):
            if line.startswith("Entry"):
                continue

            uniprot_gene_name = line.split("\t")[4].split(" ")[0].upper()
            cath_superfamily = line.split("\t")[7]
            pdb_id = line.split("\t")[8].strip()
            uniprot_sequence = line.split("\t")[9].strip()
            org = line.split("\t")[2].split("_")[1]

            if uniprot_gene_name.endswith(";"):
                uniprot_gene_name = uniprot_gene_name[:-1]

            data_fct[uniprot_gene_name] = dict(cath_superfamily=cath_superfamily, pdb_id=pdb_id,
                                               uniprot_sequence=uniprot_sequence, organism=org)


def read_genome_fasta(fna_file_path):
    sequences_fct = list(SeqIO.parse(fna_file_path, "fasta"))

    return sequences_fct


# extract gene ID number from gff file
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


# extract gene name from gff file
def extract_gene_name_gff(line: str) -> str:
    gene_name = None

    attributes = line.split("\t")[8]

    for item in attributes.split(";"):
        if item.split("=")[0] == "gene":
            gene_name = item.split("=")[1].upper()  # all gene names are handled as upper case
            break

    return gene_name


# translate nucleotide sequence to protein sequence
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


# convert a dna string of the "+" strand to "-" strand
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


def cds2seq(start: int, stop: int, sequences_fct, sequence_region: str, strand: str):
    seq = ""

    # extract sequence out of genome fasta file
    for record in sequences_fct:
        if record.id == sequence_region:
            region_seq = record.seq[start - 1:stop]

            # Turn around sequence for reverse strand
            if strand == "+":
                seq = str(region_seq)
            else:
                seq = reverse_complement(str(region_seq))

            # record id found end loop
            break

    return seq


def collect_prot_seqs(as_events: dict) -> str:
    prot_seqs = []

    for i in range(int(len(as_events) / 2)):
        full_seq = ""
        for seq in as_events[str(i) + "a"].values():
            seq = str(seq)  # not necessary needed, but pycharm gives warning otherwise
            full_seq = full_seq + seq

        # check whether there is even 1 CDS entry
        if full_seq != "":
            pro_seq = convert_nuc2aa(full_seq, int(as_events[i][0]["frame"]))

            if pro_seq.endswith("*") or pro_seq.endswith("U") or pro_seq.endswith("-"):
                pro_seq = pro_seq[:-1]

            prot_seqs.append(pro_seq)

        else:
            prot_seqs.append("-")

    if len(prot_seqs) > 1:
        prot_seq = ",".join(prot_seqs)
    elif len(prot_seqs) == 1:
        prot_seq = prot_seqs[0]
    else:
        prot_seq = "-"

    return prot_seq


def read_gff(path: str, data_fct: dict, sequences_fct):
    with open(path) as f:
        cdss = {}
        cds_chain = False
        extracted_sequences = {}
        as_events = {}
        old_gene_name = ""
        i = 0
        j = 0
        for line in tqdm(f):
            if line.startswith("#"):
                continue

            if line.split("\t")[2] == "pseudogene":
                gene_name = extract_gene_name_gff(line)

                if gene_name not in data_fct.keys():
                    continue

                data_fct[gene_name]["pseudogene"] = True

                cds_chain = False

            elif line.split("\t")[2] == "CDS":
                gene_name = extract_gene_name_gff(line)

                if gene_name not in data_fct.keys():
                    continue

                # in case an entry is already present, skip this one...
                if "CDSs" in data_fct[gene_name].keys():
                    continue

                data_fct[gene_name]["pseudogene"] = False
                strand = line.split("\t")[6]
                data_fct[gene_name]["strand"] = strand
                sequence_region = line.split("\t")[0]
                data_fct[gene_name]["gff_sequence_region"] = sequence_region

                # clear the record of CDSs in case new gene is handled
                if gene_name != old_gene_name:
                    # add last collected CDSs
                    as_events[i] = cdss.copy()
                    as_events[str(i) + "a"] = extracted_sequences.copy()
                    extracted_sequences.clear()
                    cdss.clear()
                    j = 0

                    # add counter for coding sequences
                    cds_counter = []
                    for i in range(int(len(as_events) / 2)):  # divide by 2 because of 0a and so on
                        cds_counter.append(len(as_events[i]))
                    data_fct[old_gene_name]["cds_counter"] = cds_counter.copy()
                    cds_counter.clear()

                    # Add protein sequences as new column of df
                    data_fct[old_gene_name]["protein_sequences"] = collect_prot_seqs(as_events)

                    # add all collected alternative splicing events to the collected data
                    data_fct[old_gene_name]["CDSs"] = as_events.copy()
                    as_events.clear()
                    i = 0

                # extract information
                start = int(line.split("\t")[3])
                stop = int(line.split("\t")[4])
                frame = int(line.split("\t")[7])

                # in case of multiple alternative sequencing events, add last CDSs collection to as_events
                if not cds_chain and gene_name == old_gene_name:
                    as_events[i] = cdss.copy()
                    as_events[str(i) + "a"] = extracted_sequences.copy()
                    extracted_sequences.clear()
                    cdss.clear()
                    j = 0
                    i += 1

                cdss[j] = {"start": start, "stop": stop, "frame": frame}

                # extract sequence out of genome fasta file
                extracted_sequences[j] = cds2seq(start, stop, sequences_fct, sequence_region, strand)

                # remember old gene name and increment CDS counter, also indicated that cds has been processed
                cds_chain = True
                old_gene_name = gene_name
                j += 1

            else:
                cds_chain = False


# convert to pandas dataframe and save as .tsv-file
def print_data(path: str, data_fct: dict):
    df = pd.DataFrame.from_dict(data_fct, orient="index")
    df.replace(to_replace="", value="NA")
    df.index.name = "gene_name"
    df.to_csv(path, sep="\t")


def get_organism_paths(org) -> (str, str, str):
    possible_organisms: list[str] = ["human", "mouse"]

    if org == "human":
        uniprot_path = "data/uniprot_gene3d_human.tsv"
        gff_file_path = "data/genomes/homo_sapiens/GCF_000001405.40_GRCh38.p14_genomic.gff"
        fna_file_path = "data/genomes/homo_sapiens/GCF_000001405.40_GRCh38.p14_genomic.fna"
        out_file_path = "output/uniprot_genbank_homo_sapiens.tsv"
    elif org == "mouse":
        uniprot_path = "data/uniprot_gene3d_mouse.tsv"
        gff_file_path = "data/genomes/mus_musculus/GCF_000001635.27_GRCm39_genomic.gff"
        fna_file_path = "data/genomes/mus_musculus/GCF_000001635.27_GRCm39_genomic.fna"
        out_file_path = "output/uniprot_genbank_mus_musculus.tsv"
    else:
        raise Exception(f"given organism cannot be processed, possible organisms: {possible_organisms}")

    return uniprot_path, gff_file_path, fna_file_path, out_file_path


if __name__ == "__main__":
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--organism",
        required=True,
        type=str,
        help=(
            "Organism f.e. human, mouse... for which the data should be collected."
        )
    )
    args = parser.parse_args()

    # get organism to work with
    organism = args.organism

    # get paths for given organism
    uni_path, gff_path, fna_path, out_path = get_organism_paths(organism)

    data = dict()
    read_uniprot_tsv(uni_path, data)

    sequences = read_genome_fasta(fna_path)
    read_gff(gff_path, data, sequences)

    print_data(out_path, data)



