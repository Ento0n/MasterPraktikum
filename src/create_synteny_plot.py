import sys
import ast

from pygenomeviz import GenomeViz
import argparse
import pandas as pd


def coordinates2segments(starts: list, stops: list):
    starts.sort()
    stops.sort()

    # buffer around the starts and stop for nicer displaying
    buffer = 5000

    segments = list()
    next_segment_start = starts[0] - buffer
    for i, stop in enumerate(stops):
        # check whether all are processed
        if i+1 == len(stops):
            # add last segment
            segments.append((next_segment_start, stop + buffer))

            # end loop, all found
            break

        # get second coordinate
        start_next = starts[i+1]

        # check whether distance is higher than 60000, otherwise split into segments, avg length gene 62.000
        if start_next - stop > 60000:
            # create new segment
            new_segment = (next_segment_start, stop + buffer)
            segments.append(new_segment)

            # save next segment start for next new segment
            next_segment_start = start_next - buffer

    if not segments:
        segments.append((starts[0], stops[0]))

    return segments


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

    # get gene start and stop and collect everything sorted by sequence region
    sequence_region_list = dict()
    for selection in selections:
        # go through intersected genes
        for gene in gene_names:
            try:
                correct_index = int(selection.at[gene, "correct_index"])
            except ValueError:
                print(f"gene {gene} of organism {selection.at[gene, 'organism']} could not be correctly"
                      f"translated and is therefore not represented in the synteny plot!")
                continue

            start = sys.maxsize
            stop = - sys.maxsize
            # go through CDSs of the correct as event
            tmp_dict = ast.literal_eval(selection.at[gene, "CDSs"])
            for cds in tmp_dict[correct_index].values():
                # check start, and stop whether it is smaller or bigger than already extracted start and stop
                if cds["start"] < start:
                    start = cds["start"]
                if cds["stop"] > stop:
                    stop = cds["stop"]

            # sort stuff by gff sequence region
            sequence_region = selection.at[gene, "gff_sequence_region"]

            if sequence_region in sequence_region_list.keys():
                # check whether start is larger or stop is smaller
                if sequence_region_list[sequence_region]["start"] > start:
                    sequence_region_list[sequence_region]["start"] = start
                if sequence_region_list[sequence_region]["stop"] < stop:
                    sequence_region_list[sequence_region]["stop"] = stop

                # add gene to features
                sequence_region_list[sequence_region]["features"].append(dict(
                    name=gene + f"_{selection.at[gene, 'organism']}", strand=selection.at[gene, "strand"],
                    start=start, stop=stop
                ))
            else:
                sequence_region_list[sequence_region] = dict()
                sequence_region_list[sequence_region]["start"] = start
                sequence_region_list[sequence_region]["stop"] = stop
                sequence_region_list[sequence_region]["features"] = [dict(
                    name=gene + f"_{selection.at[gene, 'organism']}", strand=selection.at[gene, "strand"],
                    start=start, stop=stop
                )]

        # divide starts and stops of the genes into segments -> otherwise region is way too big
        for sequence_region in sequence_region_list:
            starts = list()
            stops = list()
            for i, _ in enumerate(sequence_region_list[sequence_region]["features"]):
                starts.append(sequence_region_list[sequence_region]["features"][i]["start"])
                stops.append(sequence_region_list[sequence_region]["features"][i]["stop"])

            # process starts and stops into segments
            segments = coordinates2segments(starts, stops)
            sequence_region_list[sequence_region]["segments"] = segments

        # sort features by coordinates
        for attributes in sequence_region_list.values():
            attributes["features"] = sorted(attributes["features"], key=lambda x: x["start"])

    return sequence_region_list


def create_plot(sequence_region_list: dict):
    gv = GenomeViz(track_align_type="left")

    for sequence_region, attributes in sequence_region_list.items():
        track = gv.add_feature_track(sequence_region, segments=attributes["segments"])

        # change seperator of the segments
        track.set_segment_sep(symbol="//")  # is shifted around like crazy...only when track_align_type = center!

        # add sub label for segments
        for segment in track.segments:
            segment.add_sublabel(size=6)

        for feature in attributes["features"]:
            if feature["strand"] == "+":
                strand = +1
            else:
                strand = -1

            # go through segments and add to according segment
            for segment in track.segments:
                if segment.start <= feature["start"] <= segment.start + segment.size:
                    segment.add_feature(feature["start"], feature["stop"], strand, label=feature["name"])
                    break

    gv.savefig("test_plot.png")


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
    seq_reg_list = extract_collected_info(org_paths, superfamily)

    # create the synteny plot
    create_plot(seq_reg_list)



