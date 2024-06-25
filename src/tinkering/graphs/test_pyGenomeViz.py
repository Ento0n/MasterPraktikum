from pygenomeviz import GenomeViz
from pygenomeviz.parser import Gff

gff = Gff("data/sample_data/UCHL1_human.gff")

gv = GenomeViz()

track = gv.add_feature_track("UCHL1 on human chromosome 4", (41256928 - 3000, 41268455 + 1000))
cds_features = gff.extract_features(feature_type="CDS")
track.add_features(cds_features)
track.add_sublabel()

gff = Gff("data/sample_data/UCHL1_mouse.gff")

track = gv.add_feature_track("UCHL1 on mouse chromosome 5", (66833464 - 3000, 66844577 + 1000))
cds_features = gff.extract_features(feature_type="CDS")
track.add_features(cds_features)
track.add_sublabel()

gff = Gff("data/sample_data/BAP1_human.gff")

track = gv.add_feature_track("BAP1 on human chromosome 3", (52401008 - 3000, 52410008 + 1000))
cds_features = gff.extract_features(feature_type="CDS")
track.add_features(cds_features)
track.add_sublabel()

gff = Gff("data/sample_data/BAP1_mouse.gff")

track = gv.add_feature_track("BAP1 on mouse chromosome 14", (30973446 - 3000, 30981886 + 1000))
cds_features = gff.extract_features(feature_type="CDS")
track.add_features(cds_features)
track.add_sublabel()

gv.add_link(("UCHL1 on human chromosome 4", 41256977, 41268073), ("UCHL1 on mouse chromosome 5", 66833633, 66844271))
gv.add_link(("BAP1 on human chromosome 3", 52409878, 52402288), ("BAP1 on mouse chromosome 14", 30973575, 30980796))

gv.savefig("gff_sample_features.png")
