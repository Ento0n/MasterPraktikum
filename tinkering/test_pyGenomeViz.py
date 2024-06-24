from pygenomeviz import GenomeViz
from pygenomeviz.parser import Gff

gff = Gff("data/sample.gff")

gv = GenomeViz()

# segment adapted to range of UCHL1 gene
track = gv.add_feature_track(name=gff.name, segments=(41256928, 41268455))

# plot cds features
cds_features = gff.extract_features(feature_type="CDS")
track.add_features(cds_features)
track.add_sublabel()

gv.savefig("gff_sample_features.png")
