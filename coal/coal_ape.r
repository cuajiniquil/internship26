library("ape")
library("ggplot2")

#args <- commandArgs(trailingOnly = TRUE)
#setwd(args[1])

benchmarkSys <- function(trset, fn, times) {
  ttaken <- rep(0, times)
  fn(trset[[1]])

  for (t in 1:times) {
    tree <- trset[[((t - 1) %% length(trset)) + 1]]
    tstart <- Sys.time()
    fn(tree)
    tend <- Sys.time()
    ttaken[t] <- (tend - tstart) * 10^6
  }

  ttaken
}

load_newick_trees <- function(file_path) {
  lines <- readLines(file_path)
  lines <- lines[nzchar(lines)]
  lapply(lines, function(line) read.tree(text = line))
}

children_map <- function(tree) {
  split(tree$edge[, 2], tree$edge[, 1])
}

edge_length_map <- function(tree) {
  stats::setNames(tree$edge.length, tree$edge[, 2])
}

traverse_preorder <- function(tree) {
  root <- Ntip(tree) + 1L
  children <- children_map(tree)
  stack <- c(root)
  count <- 0L

  while (length(stack) > 0) {
    node <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    count <- count + 1L
    kids <- children[[as.character(node)]]
    if (!is.null(kids)) {
      stack <- c(stack, rev(kids))
    }
  }

  count
}

traverse_postorder <- function(tree) {
  root <- Ntip(tree) + 1L
  children <- children_map(tree)
  stack <- list(list(node = root, visited = FALSE))
  count <- 0L

  while (length(stack) > 0) {
    frame <- stack[[length(stack)]]
    stack <- stack[-length(stack)]

    if (frame$visited) {
      count <- count + 1L
    } else {
      stack[[length(stack) + 1L]] <- list(node = frame$node, visited = TRUE)
      kids <- children[[as.character(frame$node)]]
      if (!is.null(kids)) {
        for (child in rev(kids)) {
          stack[[length(stack) + 1L]] <- list(node = child, visited = FALSE)
        }
      }
    }
  }

  count
}

traverse_levelorder <- function(tree) {
  root <- Ntip(tree) + 1L
  children <- children_map(tree)
  queue <- c(root)
  count <- 0L

  while (length(queue) > 0) {
    node <- queue[[1]]
    queue <- queue[-1]
    count <- count + 1L
    kids <- children[[as.character(node)]]
    if (!is.null(kids)) {
      queue <- c(queue, kids)
    }
  }

  count
}

is_ultrametric_manual <- function(tree) {
  root <- Ntip(tree) + 1L
  children <- children_map(tree)
  lengths <- edge_length_map(tree)
  distances <- numeric(0)

  visit <- function(node, dist) {
    kids <- children[[as.character(node)]]
    if (is.null(kids)) {
      distances <<- c(distances, dist)
      return()
    }

    for (child in kids) {
      branch_len <- lengths[[as.character(child)]]
      if (is.null(branch_len)) {
        branch_len <- 0
      }
      visit(child, dist + branch_len)
    }
  }

  visit(root, 0)
  if (length(distances) <= 1) {
    return(TRUE)
  }
  max(distances) - min(distances) < 1e-6
}

# -------------------------------------------------------
# Coalescent log-likelihood under exponential growth
# N(t) = phi * exp(-r*t),  t = backward time from present
# Tutorial eq. (8)
# -------------------------------------------------------
coal_exp_loglik <- function(tree, phi = 10, r = 0.1) {
  Ntips  <- length(tree$tip.label)
  Nnodes <- tree$Nnode
  
  # castor native: vector[1..Ntips] = tip distances,
  #               vector[(Ntips+1)..(Ntips+Nnodes)] = internal-node distances
  all_dists  <- get_all_distances_to_root(tree)
  tip_dists  <- all_dists[seq_len(Ntips)]
  node_dists <- all_dists[(Ntips + 1L):(Ntips + Nnodes)]
  
  max_rtt   <- max(tip_dists)
  tip_ages  <- max_rtt - tip_dists    # sampling events (backward time)
  node_ages <- max_rtt - node_dists   # coalescent events (backward time)
  
  # Merge and sort all events
  evt_time <- c(tip_ages, node_ages)
  evt_type <- c(rep("S", Ntips), rep("C", Nnodes))   # S=sample, C=coal
  ord      <- order(evt_time)
  evt_time <- evt_time[ord]
  evt_type <- evt_type[ord]
  n_events <- length(evt_time)
  
  # Closed-form integral of 1/N(t) over [a, b]
  int_invN <- function(a, b) {
    if (abs(r) < 1e-12) (b - a) / phi
    else (exp(r * b) - exp(r * a)) / (phi * r)
  }
  
  n_active <- 0L
  ll <- 0.0
  
  for (i in seq_len(n_events)) {
    if (evt_type[[i]] == "S") n_active <- n_active + 1L
    else                       n_active <- n_active - 1L
    
    if (i < n_events) {
      k   <- n_active
      kc2 <- k * (k - 1L) / 2L
      if (kc2 == 0L) next          # k < 2: no coalescence possible
      
      a <- evt_time[[i]]
      b <- evt_time[[i + 1L]]
      
      ll <- ll - kc2 * int_invN(a, b)              # survival term
      if (evt_type[[i + 1L]] == "C")               # rate term at coalescent
        ll <- ll + log(kc2) - log(phi) + r * b
    }
  }
  ll
}

ultra1k <- load_newick_trees("ultra1k.nwk")
hetero1k <- load_newick_trees("hetero1k.nwk")

benchmark_tree_set <- function(tree_set) {
  coal_fn <- function(tree) coal_exp_loglik(tree, phi = 10, r = 0.1)
  list(
    preorder       = benchmarkSys(tree_set, traverse_preorder,     1000),
    postorder      = benchmarkSys(tree_set, traverse_postorder,    1000),
    levelorder     = benchmarkSys(tree_set, traverse_levelorder,   1000),
    is_ultrametric = benchmarkSys(tree_set, is_ultrametric_manual, 1000),
    coal_loglik    = benchmarkSys(tree_set, coal_fn,               1000)
  )
}

print_stats <- function(tree_label, results) {
  cat(sprintf("\n=== castor Benchmark Statistics (%s) ===\n", tree_label))
  for (method_name in c("preorder", "postorder", "levelorder",
                        "is_ultrametric", "coal_loglik")) {
    cat(sprintf("\n%s:\n", method_name))
    print(summary(results[[method_name]]))
  }
}

save_plot <- function(results, title, filename) {
  methods <- c("preorder", "postorder", "levelorder",
               "is_ultrametric", "coal_loglik")
  plotdf <- data.frame(
    method = rep(methods, each = 1000),
    time   = c(results$preorder, results$postorder, results$levelorder,
               results$is_ultrametric, results$coal_loglik)
  )
  p <- ggplot(plotdf, aes(x = method, y = time)) +
    geom_violin() +
    labs(x = "Method", y = "Time (microseconds)", title = title)
  ggsave(filename, plot = p, width = 10, height = 6)
}

benchmark_tree_set(ultra1k)