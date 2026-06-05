library(phylo2vec)
library(ape)

ultra1k <- read.tree("ultra1k.nwk")

p2vecu1k <- from_newick("phylo2vecu1k.nwk", with_mapping = TRUE)
