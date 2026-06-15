library("castor")
library("ggplot2")

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])

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
  lapply(lines, function(line) read_tree(string = line))
}

children_map <- function(tree) {
  split(tree$edge[, 2], tree$edge[, 1])
}

edge_length_map <- function(tree) {
  stats::setNames(tree$edge.length, tree$edge[, 2])
}

traverse_preorder <- function(tree) {
  root <- tree$root
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
  root <- tree$root
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
  root <- tree$root
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
  root <- tree$root
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

ultra1k <- load_newick_trees("ultra1k.nwk")
hetero1k <- load_newick_trees("hetero1k.nwk")

benchmark_tree_set <- function(tree_set) {
  list(
    preorder = benchmarkSys(tree_set, traverse_preorder, 1000),
    postorder = benchmarkSys(tree_set, traverse_postorder, 1000),
    levelorder = benchmarkSys(tree_set, traverse_levelorder, 1000),
    is_ultrametric = benchmarkSys(tree_set, is_ultrametric_manual, 1000)
  )
}

print_stats <- function(tree_label, results) {
  cat(sprintf("\n=== castor Benchmark Statistics (%s) ===\n", tree_label))
  for (method_name in c("preorder", "postorder", "levelorder", "is_ultrametric")) {
    method_times <- results[[method_name]]
    cat(sprintf("\n%s:\n", method_name))
    print(summary(method_times))
  }
}

save_plot <- function(results, title, filename) {
  plotdf <- data.frame(
    method = rep(c("preorder", "postorder", "levelorder", "is_ultrametric"), each = 1000),
    time = c(results$preorder, results$postorder, results$levelorder, results$is_ultrametric)
  )

  p <- ggplot(plotdf, aes(x = method, y = time)) +
    geom_violin() +
    labs(x = "Method", y = "Time (microseconds)", title = title)

  ggsave(filename, plot = p, width = 8, height = 6)
}

ultra_results <- benchmark_tree_set(ultra1k)
hetero_results <- benchmark_tree_set(hetero1k)

print_stats("Ultrametric", ultra_results)
print_stats("Heterochronic", hetero_results)
save_plot(ultra_results, "castor Tree Traversal Methods Comparison - Ultrametric", "castor_traverse_ultra.png")
save_plot(hetero_results, "castor Tree Traversal Methods Comparison - Heterochronic", "castor_traverse_hetero.png")