import pandas as pd
import numpy as np
import argparse


def get_organism_info(org: str):
    possible_organisms: list[str] = ["human", "mouse"]

    if org == "human":
        df = pd.read_csv("output/uniprot_genbank_homo_sapiens.tsv", sep="\t", index_col=0)
        out_path = "output/uniprot_genbank_homo_sapiens.tsv"
    elif org == "mouse":
        df = pd.read_csv("output/uniprot_genbank_mus_musculus.tsv", sep="\t", index_col=0)
        out_path = "output/uniprot_genbank_mus_musculus.tsv"
    else:
        raise Exception(f"given organism cannot be processed, possible organisms: {possible_organisms}")

    return df, out_path


def check_correctly_translated(df: pd.DataFrame, v: bool):
    wrong_sequences = set()
    counter = 0
    gene_names_correct = set()
    first_three_bases_correct = set()

    # boolean value whether catching for each gene whether gene is correctly translated or not
    correctly_translated = dict()
    correct_index = dict()

    for i, row in df.iterrows():
        correct_flag = False
        first_3_flag = False
        uniprot_seq = row["uniprot_sequence"]

        if row["protein_sequences"] is np.nan:
            correctly_translated[i] = False
            correct_index[i] = "-"
            continue

        strand = row["strand"]

        counter += 1

        for j, sequence in enumerate(row["protein_sequences"].split(",")):
            if sequence == uniprot_seq:
                correct_flag = True
                correct_index[i] = j
            elif sequence[:3] == uniprot_seq[:3] and strand == "+":
                first_3_flag = True

        if correct_flag:
            correctly_translated[i] = True
            gene_names_correct.add(i)
        elif first_3_flag:
            correctly_translated[i] = False
            correct_index[i] = "-"
            first_three_bases_correct.add(i)
            wrong_sequences.add(i)
        else:
            correctly_translated[i] = False
            correct_index[i] = "-"
            wrong_sequences.add(i)

    if v:
        print(f"sequences in general: {counter}")
        print(f"correct sequences: {len(gene_names_correct)}")
        # print(f"gene names correct sequences: {gene_names_correct}")
        print(f"wrong sequences: {len(wrong_sequences)}")
        print(f"wrong sequences: {list(wrong_sequences)[:10]}")
        print(f"only first 3 match, rest uncertain: {first_three_bases_correct}")
        print(f"length of first 3 match: {len(first_three_bases_correct)}")

    return correctly_translated, correct_index


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
    parser.add_argument(
        "-v",
        "--verbose",
        required=False,
        action="store_true",
        help=(
            "Enables printing # of correct and wrong sequences etc."
        )
    )
    args = parser.parse_args()

    # get organism to work with
    organism = args.organism
    verbose = args.verbose

    # get df
    data_frame, output_path = get_organism_info(organism)

    correctly_trans, correct_ind = check_correctly_translated(data_frame, verbose)

    # add column to df and save it as .tsv-file
    data_frame["correctly_translated"] = correctly_trans
    data_frame["correct_index"] = correct_ind
    data_frame.to_csv(output_path, sep="\t")
