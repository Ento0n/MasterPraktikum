import pygenomeviz
import argparse
import pandas as pd


def extract_collected_info(orgs: [str], superf: str):
    # convert organism path to df
    dfs = []
    for org in orgs:
        dfs.append(pd.read_csv(org, sep="\t", index_col=0))

    # extract entries that have the corresponding superfamily
    selections = []
    for df in dfs:
        selections.append(df[df["cath_superfamily"].apply(lambda x: superf in x.split(";"))])

    # filter out which genes are given in all selections
    gene_names = set(selections[0].index.tolist())
    for selection in selections:
        selection_gene_names = set(selection.index.tolist())
        gene_names = gene_names & selection_gene_names







def get_organism_paths(orgs: [str]):
    possible_organisms: list[str] = ["human", "mouse"]

    paths = []
    for org in orgs:
        if org == "human":
            paths.append("output/uniprot_genbank_homo_sapiens.tsv")

        if org == "mouse":
            paths.append("output/uniprot_genbank_mus_musculus.tsv")

        if org not in possible_organisms:
            raise Exception(f"given organism cannot be processed, possible organisms: {possible_organisms}")

    return paths


if __name__ == "__main__":
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-os",
        "--organisms",
        required=True,
        type=str,
        help=(
            "Organism f.e. human, mouse... for which the data should be collected."
        )
    )
    parser.add_argument(
        "-sf",
        "--superfamily",
        required=True,
        type=str,
        help=(
            "superfamily for which the synteny plot should be created"
        )
    )
    args = parser.parse_args()

    # get organisms and superfamily to work with
    organisms = args.organisms
    superfamily = args.superfamily

    # convert organisms into list
    organisms = organisms.split(",")

    # get paths for wanted organisms
    org_paths = get_organism_paths(organisms)

    # extract needed infos from .tsv-file
    extract_collected_info(org_paths, superfamily)



