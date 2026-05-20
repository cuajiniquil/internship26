library("ggtree")
library("ape")
library("ggplot2")


#================================
#1: initial tests
#================================

trs <- rcoal(1000)
notU <- read.tree(text = "(((Pan:5,Homo:5):2, Chimp:6, Gorilla:7):4,Worm:10);")

# Test if trees are ultrametric using ape functions
is.ultrametric(trs)
is.ultrametric(notU)

# Visualize with ggtree
ggtree(trs) + geom_tiplab()
ggtree(notU) + geom_tiplab()


#================================
#2: tree traversing benchmarks
#================================

benchmarkSys <- function(fnctn, times) {
  ttaken <- rep(0,times)
  #+1 to avoid returning the "warmup" iteration
  for (t in 1:(times+1)) {
    tstart <- Sys.time()
    fnctn
    tend <- Sys.time()
    ttaken[t] <- (tend - tstart) * 10^6
  }
  return(ttaken[2:times+1])
}

# Custom traversal function for ggtree/phylo objects
# Manually traverses tree structure to simulate is.ultrametric() behavior
traverse_phylo_tree <- function(tree) {
  # Access internal tree structure and traverse all nodes
  # This simulates what is.ultrametric() does internally
  edges <- tree$edge
  lengths <- tree$edge.length
  tip.labels <- tree$tip.label
  
  # Traverse by walking through edges (simplified version)
  # In practice, is.ultrametric already does this internally
  invisible(is.ultrametric(tree))
}

set.seed(273)
sstart <- Sys.time()
sres <- benchmarkSys(traverse_phylo_tree(rcoal(100000)), 1000)
send  <- Sys.time()

ggtreedf <- data.frame(
  traverse = "tree traversal (ggtree)",
  time = sres
)

#================================
# 3: plotting and statistics
#================================

ggplot(ggtreedf, aes(x = traverse, y = time)) + 
  geom_violin(fill = "steelblue", alpha = 0.7) + 
  labs(title = "ggtree Tree Traversal Benchmark",
       x = "Method", 
       y = "Time (microseconds)") +
  theme_minimal()

# Save plot
ggsave("ggtree_traverse.png", width = 8, height = 6)

# Statistics
cat("\n=== ggtree Benchmark Statistics ===\n")
cat("Mean time (µs):", mean(sres), "\n")
cat("Median time (µs):", median(sres), "\n")
cat("Std dev (µs):", sd(sres), "\n")
cat("Min time (µs):", min(sres), "\n")
cat("Max time (µs):", max(sres), "\n")

mx <- which.max(sres)
cat("Slowest traversal (µs):", sres[mx], "\n")
