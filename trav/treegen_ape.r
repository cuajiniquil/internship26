library("ape")
library("ggplot2")

set.seed(273)

trsetu1 = lapply(1:1000, function(i) rcoal(1000))
#trsetu10 = lapply(1:1000, function(i) rcoal(10000))
trseth1 = lapply(1:1000, function(i) rtree(1000))
#trseth10 = lapply(1:1000, function(i) rtree(10000))

class(trsetu1) <- "multiPhylo"
class(trseth1) <- "multiPhylo"

write.tree(trsetu1, file = "ultra1k.nwk")
#write.tree(trsetu10, file = "ultra10k.nwk")
write.tree(trseth1, file = "hetero1k.nwk")
#write.tree(trseth10, file = "hetero10k.nwk")