#!/usr/bin/env Rscript
library("castor")
library("ggplot2")


#================================
# 1: Constants
#================================

FIXED_PHI0 <- 1.0
FIXED_R <- 0.023882220114924


#================================
# 2: Helper functions
#================================

benchmarkSys <- function(trset, fn, times) {
  ttaken <- rep(0, times)
  n <- length(trset)
  
  fn(trset[[1]])  # warmup
  
  for (t in 1:times) {
    idx <- ((t - 1) %% n) + 1
    tstart <- Sys.time()
    fn(trset[[idx]])
    tend <- Sys.time()
    ttaken[t] <- as.numeric(tend - tstart, units = "secs") * 1e6
  }
  return(ttaken)
}


load_newick_trees <- function(file_path) {
  lines <- readLines(file_path)
  lines <- lines[nzchar(lines)]
  lapply(lines, function(line) read_tree(string = line))
}


#================================
# 3: Coalescent likelihood
#================================

extract_coalescent_intervals <- function(tree) {
  n_tips  <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode
  
  # native root-to-tip distances via castor's level-order traversal
  dists <- castor::get_all_distances_to_root(tree, as_edge_count = FALSE)
  dists <- dists[seq_len(n_nodes)]
  
  present <- max(dists[seq_len(n_tips)])
  ages <- present - dists
  
  events <- data.frame(
    time = ages,
    type = c(rep("sample", n_tips), rep("coalescent", tree$Nnode)),
    name = c(tree$tip.label, as.character((n_tips + 1):n_nodes)),
    stringsAsFactors = FALSE
  )
  events <- events[order(events$time), ]
  rownames(events) <- NULL
  
  n_active <- 0
  intervals <- list()
  for (i in seq_len(nrow(events) - 1)) {
    n_active <- n_active + ifelse(events$type[i] == "sample", 1L, -1L)
    intervals[[i]] <- list(
      a         = events$time[i],
      b         = events$time[i + 1],
      k         = n_active,
      ends_coal = events$type[i + 1] == "coalescent"
    )
  }
  
  intervals
}


coalescent_loglikelihood <- function(intervals,
                                     phi0 = FIXED_PHI0,
                                     r    = FIXED_R) {
  log_phi0 <- log(phi0)
  logL     <- 0.0
  
  for (iv in intervals) {
    k  <- iv$k
    a  <- iv$a
    b  <- iv$b
    c2 <- k * (k - 1) / 2
    
    if (c2 <= 0) next
    
    # survival term
    if (r == 0) {
      logL <- logL - c2 * (b - a) / phi0
    } else {
      logL <- logL - (c2 / (phi0 * r)) * (exp(r * b) - exp(r * a))
    }
    
    # event term
    if (iv$ends_coal) {
      logL <- logL + log(c2) - log_phi0 + r * b
    }
  }
  
  logL
}


loglik_castor <- function(tree,
                          phi0 = FIXED_PHI0,
                          r    = FIXED_R) {
  intervals <- extract_coalescent_intervals(tree)
  coalescent_loglikelihood(intervals, phi0, r)
}


#================================
# 4: Sanity checks
#================================

sntree <- read_tree(string = "((t1:0.9383623928,t2:0.5619198801):0.531734891,t3:0.286919398);")
cat("Sanity check logL:", loglik_castor(sntree, 10, 0.1), "\n")
cat("Expected:          -4.457106\n")


#================================
# 5: Load tree sets
#================================

cat("\nLoading trees...\n")
ultra1k  <- load_newick_trees("ultra1k.nwk")
hetero1k <- load_newick_trees("hetero1k.nwk")
cat("Loaded", length(ultra1k), "ultrametric trees and", length(hetero1k), "heterochronic trees\n")

cat("\nFirst 5 ultrametric log-likelihoods:\n")
for (i in 1:5) {
  cat("  Ultra#", i, ": ", round(loglik_castor(ultra1k[[i]]), 4), "\n", sep = "")
}

cat("\nFirst 5 heterochronic log-likelihoods:\n")
for (i in 1:5) {
  cat("  Hetero#", i, ": ", round(loglik_castor(hetero1k[[i]]), 4), "\n", sep = "")
}


#================================
# 6: Benchmarks
#================================

cat("\nStarting benchmarking (ultrametric)...\n")
ultra_results <- list(
  interval_extraction = benchmarkSys(ultra1k, extract_coalescent_intervals, 1000),
  log_likelihood      = benchmarkSys(ultra1k, loglik_castor, 1000)
)

cat("\nStarting benchmarking (heterochronic)...\n")
hetero_results <- list(
  interval_extraction = benchmarkSys(hetero1k, extract_coalescent_intervals, 1000),
  log_likelihood      = benchmarkSys(hetero1k, loglik_castor, 1000)
)


#================================
# 7: Statistics
#================================

print_stats <- function(tree_label, results) {
  cat("\n=== castor Benchmark Statistics Coalescent likelihood (", tree_label, ") ===\n", sep = "")
  for (step in names(results)) {
    step_times <- results[[step]]
    cat("\n", step, ":\n", sep = "")
    cat("  Mean time (µs):   ", round(mean(step_times), 4), "\n")
    cat("  Median time (µs): ", round(median(step_times), 4), "\n")
    cat("  Std dev (µs):     ", round(sd(step_times), 4), "\n")
    cat("  Min time (µs):    ", round(min(step_times), 4), "\n")
    cat("  Max time (µs):    ", round(max(step_times), 4), "\n")
  }
}

print_stats("ultrametric", ultra_results)
print_stats("heterochronic", hetero_results)


#================================
# 8: Plotting
#================================

save_plot <- function(results, tree_label, filename) {
  df <- data.frame(
    step = rep(c("interval extraction", "log likelihood"), each = 1000),
    time = c(results$interval_extraction, results$log_likelihood)
  )
  
  p <- ggplot(df, aes(x = step, y = time)) +
    geom_violin(fill = "steelblue", alpha = 0.7) +
    labs(x = "Method", y = "Time (microseconds)",
         title = paste("castor Coalescent Likelihood -", tree_label)) +
    theme_minimal()
  
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 100)
  cat("Plot saved to '", filename, "'\n", sep = "")
}

save_plot(ultra_results, "Ultrametric", "castor_coal_likelihood_ultra.png")
save_plot(hetero_results, "Heterochronic", "castor_coal_likelihood_hetero.png")
