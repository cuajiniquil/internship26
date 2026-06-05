#!/usr/bin/env python3
"""
Tree traversal benchmark for cogent3 package.
Generates ultrametric trees and measures traversal time.
"""

import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from cogent3 import make_tree
from cogent3.evolve.models import get_model


#================================
#1: initial tests
#================================

print("Generating initial test trees with cogent3...")

# Create an ultrametric tree using cogent3 newick format
# This is a simple balanced binary tree where all leaves are equidistant from root
ultra_newick = "((A:1,B:1):1,(C:1,D:1):1);"
trs = make_tree(ultra_newick)

# Create a non-ultrametric tree for comparison
non_ultra_newick = "(A:1,B:5);"
notU = make_tree(non_ultra_newick)

def is_ultrametric(tree):
    """
    Check if all leaf nodes are equidistant from root.
    """
    if tree.is_tip():
        return True
    
    def get_leaf_distances(node, dist=0.0):
        if node.is_tip():
            return [dist]
        distances = []
        for child in node.children:
            child_dist = dist + (child.length if child.length else 0)
            distances.extend(get_leaf_distances(child, child_dist))
        return distances
    
    distances = get_leaf_distances(tree)
    if len(distances) <= 1:
        return True
    # Use approximate equality for floating point comparison
    return np.allclose(distances, distances[0])


print("Ultrametric test tree is ultrametric:", is_ultrametric(trs))
print("Non-ultrametric test tree is ultrametric:", is_ultrametric(notU))

# Print tree info
print("\nTest tree structure:")
print(trs)
print("\nTest tree tip names:", [tip.name for tip in trs.tips()])


#================================
#2: tree traversing benchmarks
#================================

def benchmark_sys(tree_set, fn, times):
    """Benchmark function execution time over a list of trees."""
    ttaken = []

    fn(tree_set[0])

    for t in range(times):
        tree = tree_set[t % len(tree_set)]
        tstart = time.perf_counter()
        fn(tree)
        tend = time.perf_counter()
        ttaken.append((tend - tstart) * 1e6)

    return ttaken


def create_balanced_ultrametric_tree(num_taxa):
    """
    Create a balanced ultrametric tree using cogent3.
    Recursively builds a binary tree with equal branch lengths.
    """
    def build_balanced(taxa_list, branch_len=1.0):
        if len(taxa_list) == 1:
            return f"{taxa_list[0]}"

        mid = len(taxa_list) // 2
        left_tree = build_balanced(taxa_list[:mid], branch_len)
        right_tree = build_balanced(taxa_list[mid:], branch_len)

        # Combine with equal branch lengths to maintain ultrametricity
        combined = f"({left_tree}:{branch_len},{right_tree}:{branch_len})"
        return combined

    taxa_names = [f"taxon_{i}" for i in range(num_taxa)]
    newick_str = build_balanced(taxa_names) + ";"
    return make_tree(newick_str)


def traverse_tree(tree):
    """Traverse all nodes in tree."""
    count = 0
    for node in tree.preorder():
        count += 1
    return count


def traverse_tree_postorder(tree):
    """Traverse all nodes in tree (postorder)."""
    count = 0
    for node in tree.postorder():
        count += 1
    return count


def traverse_tree_levelorder(tree):
    """Traverse all nodes in tree (level order)."""
    count = 0
    for node in tree.levelorder():
        count += 1
    return count


print("\nExtracting trees from nwk files...")

def load_newick_trees(file_path):
    trees = []
    with open(file_path) as handle:
        for line in handle:
            line = line.strip()
            if line:
                trees.append(make_tree(line))
    return trees


ultra1k = load_newick_trees("ultra1k.nwk")
hetero1k = load_newick_trees("hetero1k.nwk")

print("Starting benchmarks...")
sstart = time.perf_counter()

ultra_results = {
    "preorder": benchmark_sys(ultra1k, traverse_tree, 1000),
    "postorder": benchmark_sys(ultra1k, traverse_tree_postorder, 1000),
    "levelorder": benchmark_sys(ultra1k, traverse_tree_levelorder, 1000),
    "is_ultrametric": benchmark_sys(ultra1k, is_ultrametric, 1000),
}

hetero_results = {
    "preorder": benchmark_sys(hetero1k, traverse_tree, 1000),
    "postorder": benchmark_sys(hetero1k, traverse_tree_postorder, 1000),
    "levelorder": benchmark_sys(hetero1k, traverse_tree_levelorder, 1000),
    "is_ultrametric": benchmark_sys(hetero1k, is_ultrametric, 1000),
}

send = time.perf_counter()

print(f"Benchmark completed in {send - sstart:.2f} seconds")


#================================
# 3: plotting and statistics
#================================


def print_stats(tree_label, results):
    print(f"\n=== cogent3 Benchmark Statistics ({tree_label}) ===")
    for method_name in ['preorder', 'postorder', 'levelorder', 'is_ultrametric']:
        method_times = results[method_name]
        print(f"\n{method_name}:")
        print(f"  Mean time (µs): {np.mean(method_times):.4f}")
        print(f"  Median time (µs): {np.median(method_times):.4f}")
        print(f"  Std dev (µs): {np.std(method_times):.4f}")
        print(f"  Min time (µs): {np.min(method_times):.4f}")
        print(f"  Max time (µs): {np.max(method_times):.4f}")


def save_plot(results, title, filename):
    fig, ax = plt.subplots(figsize=(10, 6))
    data_by_method = [results['preorder'], results['postorder'], results['levelorder'], results['is_ultrametric']]
    positions = [1, 2, 3, 4]
    ax.violinplot(data_by_method, positions=positions, widths=0.7)
    ax.set_ylabel('Time (microseconds)')
    ax.set_title(title)
    ax.set_xticks(positions)
    ax.set_xticklabels(['preorder', 'postorder', 'levelorder', 'is_ultrametric'])
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=100)
    print(f"\nPlot saved to '{filename}'")
    plt.close()


print_stats('ultrametric', ultra_results)
print_stats('heterochronic', hetero_results)
save_plot(ultra_results, 'cogent3 Tree Traversal Methods Comparison - Ultrametric', 'cogent3_traverse_ultra.png')
save_plot(hetero_results, 'cogent3 Tree Traversal Methods Comparison - Heterochronic', 'cogent3_traverse_hetero.png')
