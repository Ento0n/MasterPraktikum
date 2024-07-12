import sys
import ast

import pygenomeviz.exception
from pygenomeviz import GenomeViz
import argparse
import pandas as pd


"python src/create_synteny_plot.py -sf 2.40.10.10 -os human,mouse,bovin,pig,chimpanzee,chicken -srs NC_000019.10,NC_072418.2,NC_037345.1,NC_010448.4,NC_000073.7,NC052532.1,NC_052532.1 -o output/KLK/KLK.png -kos CTRB2,HP,LOC100848132,LOC100621820,MASP2,CELA2A,OVCH2,PRSS53,PRSS8,PRSS36,PRSS23,F7,F10,PROZ,TMPRSS3,TMPRSS2,PRSS54,GZMM,PRSS57,AZU1,PRTN3,ELANE,CFD,TMPRSS9,TMPRSS15,LOC418421,TMPRSS7 -cdss"
"python src/create_synteny_plot.py -os human,chimpanzee,bovin,pig,mouse,chicken,anole,platypus -srs NC_000008.11,NC_072405.2,NC_037341.1,NC_010446.5,NC_000081.7,NC_052533.1,NC_041731.1,NC_085844.1 -sf 2.10.60.10 -o output/LY6/LY6.png"


def coordinates2segments(starts: list, stops: list):
    starts.sort()
    stops.sort()

    # buffer around the starts and stop for nicer displaying
    buffer = 2000

    segments = list()
    next_segment_start = starts[0] - buffer
    gene2exon_index = 0
    index_list = list()
    for i, stop in enumerate(stops):
        # check whether all are processed
        if i+1 == len(stops):
            # add last segment
            segments.append((next_segment_start, stop + buffer))

            # add index of last gene
            index_list.append(gene2exon_index)

            # end loop, all found
            break

        # get second coordinate
        start_next = starts[i+1]

        # add index to index list before check for new segment, always looking at next start
        index_list.append(gene2exon_index)

        # check whether distance is higher than 60000, otherwise split into segments, avg length gene 62.000
        if start_next - stop > 40000:
            # create new segment
            new_segment = (next_segment_start, stop + buffer)
            segments.append(new_segment)

            # save next segment start for next new segment
            next_segment_start = start_next - buffer

            # new segment, increment counter
            gene2exon_index += 1

    if not segments:
        segments.append((starts[0], stops[0]))

    return segments, index_list


def extract_collected_info(orgs: [str], superf: str, ko_genes: list, wanted_genes: list, v: bool):
    # convert organism path to df
    dfs = []
    for org in orgs:
        dfs.append(pd.read_csv(org, sep="\t", index_col=0))

    # extract entries that have the corresponding superfamily
    selections = []
    for df in dfs:
        selections.append(df[df["cath_superfamily"].apply(lambda x: superf in x.split(";"))])

    # get gene start and stop and collect everything sorted by sequence region
    sequence_region_list = dict()
    for selection in selections:
        # go through intersected genes
        for gene in selection.index:
            # skip the unwanted genes
            if ko_genes:
                if gene in ko_genes:
                    continue

            if wanted_genes:
                if gene not in wanted_genes:
                    continue

            try:
                correct_index = int(selection.at[gene, "correct_index"])
                correct = True

            except ValueError:
                correct = False

                if v:
                    print(f"gene {gene} of organism {selection.at[gene, 'organism']} could not be correctly"
                          f"translated, therefore simmply the first alternative splicing event is displayed!")
                try:
                    cds_counter = ast.literal_eval(selection.at[gene, "cds_counter"])
                    correct_index = cds_counter.index(max(cds_counter))
                except ValueError:
                    correct_index = 0

            # skip the entries where everythig is missing
            if pd.isna(selection.at[gene, "CDSs"]):
                continue

            start = sys.maxsize
            stop = - sys.maxsize
            cdss = list()
            # go through CDSs of the correct as event
            tmp_dict = ast.literal_eval(selection.at[gene, "CDSs"])
            for cds in tmp_dict[correct_index].values():
                # check start, and stop whether it is smaller or bigger than already extracted start and stop
                if cds["start"] < start:
                    start = cds["start"]
                if cds["stop"] > stop:
                    stop = cds["stop"]

                # add start and stop to CDS list
                cdss.append((cds["start"], cds["stop"]))

            # sort stuff by gff sequence region
            sequence_region = selection.at[gene, "gff_sequence_region"]

            # check whether sequence region (chromosome) is already added or not
            if sequence_region in sequence_region_list.keys():
                # check whether start is larger or stop is smaller
                if sequence_region_list[sequence_region]["start"] > start:
                    sequence_region_list[sequence_region]["start"] = start
                if sequence_region_list[sequence_region]["stop"] < stop:
                    sequence_region_list[sequence_region]["stop"] = stop

                # add gene to features
                sequence_region_list[sequence_region]["features"].append(dict(
                    name=gene + f"_{selection.at[gene, 'organism']}", strand=selection.at[gene, "strand"],
                    start=start, stop=stop, cdss=cdss, correct=correct
                ))
            else:
                sequence_region_list[sequence_region] = dict()
                sequence_region_list[sequence_region]["start"] = start
                sequence_region_list[sequence_region]["stop"] = stop
                sequence_region_list[sequence_region]["features"] = [dict(
                    name=gene + f"_{selection.at[gene, 'organism']}", strand=selection.at[gene, "strand"],
                    start=start, stop=stop, cdss=cdss, correct=correct
                )]

        # sort features by coordinates
        for attributes in sequence_region_list.values():
            attributes["features"] = sorted(attributes["features"], key=lambda x: x["start"])

        # divide starts and stops of the genes into segments -> otherwise region is way too big
        for sequence_region in sequence_region_list:
            starts = list()
            stops = list()
            for i, _ in enumerate(sequence_region_list[sequence_region]["features"]):
                starts.append(sequence_region_list[sequence_region]["features"][i]["start"])
                stops.append(sequence_region_list[sequence_region]["features"][i]["stop"])

            # process starts and stops into segments
            segments, index_list = coordinates2segments(starts, stops)
            sequence_region_list[sequence_region]["segments"] = segments

            # add index to data
            for i, index in enumerate(index_list):
                sequence_region_list[sequence_region]["features"][i]["index"] = index

    return sequence_region_list


def create_track(gv, attributes: dict, sequence_region: str, cdss: bool):
    track = gv.add_feature_track(sequence_region, segments=attributes["segments"])

    # change seperator of the segments
    track.set_segment_sep(symbol="//")  # is shifted around like crazy...only when track_align_type = center!

    # add sub label for segments
    for i, segment in enumerate(track.segments):
        if i % 2 == 0:
            y_margin = 0.2
        else:
            y_margin = 0.7

        segment.add_sublabel(size=7, ymargin=y_margin)

    # add features to track
    for feature in attributes["features"]:
        if feature["strand"] == "+":
            strand = +1
        else:
            strand = -1

        # go through segments and add to according segment
        for segment in track.segments:
            if segment.start <= feature["start"] <= segment.start + segment.size:
                if feature["correct"]:
                    color = "green"
                else:
                    color = "red"

                if cdss:
                    segment.add_exon_feature(feature["cdss"], strand, intron_patch_kws=dict(ec=color),
                                             label=feature["name"] + f"_{len(feature['cdss'])}" +
                                             f"_{feature['correct']}")
                else:
                    segment.add_feature(feature["start"], feature["stop"], strand, label=feature["name"])

                break


def add_link(gv, feature: dict, sequence_region: str, tmp_sequence_region: str, gene: str, tmp_attributes: dict,
             v: bool):
    # skip the sequence region handled
    if tmp_sequence_region == sequence_region:
        return

    for tmp_feature in tmp_attributes["features"]:
        tmp_gene = tmp_feature["name"].split("_")[0]

        # if same gene is found add link
        if tmp_gene == gene:
            try:
                # turn around 1 start and stop if direction of genes is different
                if feature["strand"] == "+" and tmp_feature["strand"] == "-" or feature["strand"] == "-" and \
                        tmp_feature["strand"] == "+":
                    gv.add_link((sequence_region, "seg" + str(feature["index"] + 1),
                                 feature["stop"], feature["start"]),
                                (tmp_sequence_region, "seg" + str(tmp_feature["index"] + 1),
                                 tmp_feature["start"], tmp_feature["stop"]), curve=True, inverted_color="skyblue",
                                color="lightgrey")
                else:
                    gv.add_link((sequence_region, "seg" + str(feature["index"] + 1),
                                 feature["start"], feature["stop"]),
                                (tmp_sequence_region, "seg" + str(tmp_feature["index"] + 1),
                                 tmp_feature["start"], tmp_feature["stop"]), curve=True, inverted_color="skyblue",
                                color="lightgrey")
            except pygenomeviz.exception.LinkTrackNotFoundError:
                if v:
                    print(
                        f"Link between {sequence_region} and {tmp_sequence_region} "
                        f"for {gene} not possible, tracks not adjacent")

            # set flag for outer loops and break inner loop
            break


def create_plot(sequence_region_list: dict, wanted_sequence_regions, out_path: str, cdss: bool, v: bool):
    gv = GenomeViz(track_align_type="left")

    if wanted_sequence_regions is None:
        for sequence_region, attributes in sequence_region_list.items():
            create_track(gv, attributes, sequence_region, cdss)

    else:
        for sequence_region in wanted_sequence_regions:
            # skip not present ones
            if sequence_region not in sequence_region_list.keys():
                continue

            attributes = sequence_region_list[sequence_region]

            create_track(gv, attributes, sequence_region, cdss)

    # add links to the plot
    if wanted_sequence_regions is None:
        for sequence_region, attributes in sequence_region_list.items():
            for feature in attributes["features"]:
                gene = feature["name"].split("_")[0]

                # go through other sequence regions
                for tmp_sequence_region, tmp_attributes in sequence_region_list.items():
                    add_link(gv, feature, sequence_region, tmp_sequence_region, gene, tmp_attributes, v)

    else:
        for sequence_region in wanted_sequence_regions:
            # skip not present ones
            if sequence_region not in sequence_region_list.keys():
                continue

            attributes = sequence_region_list[sequence_region]

            for feature in attributes["features"]:
                gene = feature["name"].split("_")[0]

                # set up flag to leave loops when match is found
                flag = False

                # go through other sequence regions
                for tmp_sequence_region in wanted_sequence_regions:
                    # skip not present ones
                    if tmp_sequence_region not in sequence_region_list.keys():
                        continue

                    tmp_attributes = sequence_region_list[tmp_sequence_region]

                    add_link(gv, feature, sequence_region, tmp_sequence_region, gene, tmp_attributes, flag)

    gv.savefig(out_path)


def get_organism_paths(orgs: [str]):
    possible_organisms: list[str] = ["human", "mouse", "chicken", "bovin", "zebrafish", "pig", "chimpanzee", "platypus",
                                     "anole"]

    paths = []
    for org in orgs:
        if org == "human":
            paths.append("output/uniprot_genbank_homo_sapiens.tsv")

        if org == "mouse":
            paths.append("output/uniprot_genbank_mus_musculus.tsv")

        if org == "chicken":
            paths.append("output/uniprot_genbank_gallus_gallus.tsv")

        if org == "bovin":
            paths.append("output/uniprot_genbank_bos_taurus.tsv")

        if org == "zebrafish":
            paths.append("output/uniprot_genbank_danio_rerio.tsv")

        if org == "pig":
            paths.append("output/uniprot_genbank_sus_scrofa.tsv")

        if org == "chimpanzee":
            paths.append("output/uniprot_genbank_pan_troglodytes.tsv")

        if org == "platypus":
            paths.append("output/uniprot_genbank_platypus.tsv")

        if org == "anole":
            paths.append("output/uniprot_genbank_anole.tsv")

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
    parser.add_argument(
        "-srs",
        "--sequence_regions",
        required=False,
        type=str,
        help=(
            "sequence regions to be displayed in the synteny plot"
        )
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=str,
        help=(
            "output path for the synteny plot"
        )
    )
    parser.add_argument(
        "-cdss",
        "--coding_sequences",
        required=False,
        action="store_true",
        help=(
            "Enables displaying of cdss, genes are displayed when disabled"
        )
    )
    parser.add_argument(
        "-kos",
        "--ko_genes",
        required=False,
        type=str,
        help=(
            "genes that shouldn't be displayed in the synteny plot"
        )
    )
    parser.add_argument(
        "--genes",
        required=False,
        type=str
    )
    parser.add_argument(
        "-v",
        "--verbose",
        required=False,
        action="store_true",
        help=(
            "Enables printing of warnings"
        )
    )
    args = parser.parse_args()

    # get organisms and superfamily to work with
    organisms = args.organisms
    superfamily = args.superfamily
    sequence_regions = args.sequence_regions
    output_path = args.output
    cdss_wanted = args.coding_sequences
    knock_out_genes = args.ko_genes
    genes = args.genes
    verbose = args.verbose

    # convert organisms and sequence regions into list
    organisms = organisms.split(",")
    if sequence_regions:
        sequence_regions = sequence_regions.split(",")
    if knock_out_genes:
        knock_out_genes = knock_out_genes.split(",")
    if genes:
        genes = genes.split(",")

    # get paths for wanted organisms
    org_paths = get_organism_paths(organisms)

    # extract needed infos from .tsv-file
    seq_reg_list = extract_collected_info(org_paths, superfamily, knock_out_genes, genes, verbose)

    # create the synteny plot
    create_plot(seq_reg_list, sequence_regions, output_path, cdss_wanted, verbose)



