#!/usr/bin/env python3
"""
Tree traversal benchmark for biopython Phylo package.
Generates ultrametric trees and measures traversal time.
"""

import time
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree as PhyloTree, Clade
from io import StringIO


#================================
#1: initial tests
#================================

print("Generating initial test trees with biopython...")

# Create an ultrametric tree using Newick format
ultra_newick = "((A:1,B:1):1,(C:1,D:1):1);"
handle = StringIO(ultra_newick)
trs = Phylo.read(handle, "newick")

# Create a non-ultrametric tree for comparison
non_ultra_newick = "((A:5,B:5):2,C:7);"
handle = StringIO(non_ultra_newick)
notU = Phylo.read(handle, "newick")

def is_ultrametric(tree):
    """
    Check if all leaf nodes are equidistant from root.
    biopython doesn't have native is_ultrametric, so we implement it.
    """
    def get_leaf_distances(clade, dist=0.0):
        if clade.is_terminal():
            return [dist]
        distances = []
        for child in clade.clades:
            child_dist = dist + (child.branch_length if child.branch_length else 0)
            distances.extend(get_leaf_distances(child, child_dist))
        return distances

    distances = get_leaf_distances(tree.root)
    if len(distances) <= 1:
        return True
    return np.allclose(distances, distances[0])


print("Ultrametric test tree is ultrametric:", is_ultrametric(trs))
print("Non-ultrametric test tree is ultrametric:", is_ultrametric(notU))

# Print tree info
print("\nTest tree structure:")
Phylo.draw_ascii(trs)
print("Test tree leaf count:", len(trs.get_terminals()))


#================================
#2: tree traversing benchmarks
#================================

def benchmark_sys(fn, times):
    """
    Benchmark function execution time.
    fn: callable to benchmark
    times: number of iterations (warmup + actual)
    Returns: list of execution times in microseconds (excluding warmup)
    """
    ttaken = []

    # +1 to account for warmup iteration
    for t in range(times + 1):
        tstart = time.perf_counter()
        fn()
        tend = time.perf_counter()

        # Convert to microseconds
        elapsed_us = (tend - tstart) * 1e6
        ttaken.append(elapsed_us)

    # Return all except warmup iteration
    return ttaken[1:]


def create_balanced_ultrametric_tree(num_taxa):
    """
    Create a balanced ultrametric tree using biopython.
    Recursively builds a binary tree with equal branch lengths.
    """
    def build_balanced(taxa_list, branch_len=1.0):
        if len(taxa_list) == 1:
            return Clade(name=taxa_list[0], branch_length=0)

        mid = len(taxa_list) // 2
        left_clade = build_balanced(taxa_list[:mid], branch_len)
        right_clade = build_balanced(taxa_list[mid:], branch_len)

        left_clade.branch_length = branch_len
        right_clade.branch_length = branch_len
        parent = Clade(clades=[left_clade, right_clade])
        return parent

    taxa_names = [f"taxon_{i}" for i in range(num_taxa)]
    root = build_balanced(taxa_names)
    root.branch_length = 0
    return PhyloTree(root=root)


def traverse_tree_preorder(tree):
    """Traverse all nodes in tree (preorder)."""
    count = 0
    for clade in tree.find_clades(order="preorder"):
        count += 1
    return count


def traverse_tree_postorder(tree):
    """Traverse all nodes in tree (postorder)."""
    count = 0
    for clade in tree.find_clades(order="postorder"):
        count += 1
    return count


def traverse_tree_levelorder(tree):
    """Traverse all nodes in tree (level order)."""
    count = 0
    for clade in tree.find_clades(order="level"):
        count += 1
    return count

# Create test tree
np.random.seed(273)
print("\nGenerating large ultrametric tree for benchmark...")
test_tree = create_balanced_ultrametric_tree(10000)

print("Starting benchmarks...")
sstart = time.perf_counter()
sres_preorder = benchmark_sys(lambda: traverse_tree_preorder(test_tree), 1000)
sres_postorder = benchmark_sys(lambda: traverse_tree_postorder(test_tree), 1000)
sres_levelorder = benchmark_sys(lambda: traverse_tree_levelorder(test_tree), 1000)
sres_ultrametric = benchmark_sys(lambda: is_ultrametric(test_tree), 1000)
send = time.perf_counter()

print(f"Benchmark completed in {send - sstart:.2f} seconds")

# Combine all results
all_times = sres_preorder + sres_postorder + sres_levelorder + sres_ultrametric
all_methods = (
    ['preorder'] * len(sres_preorder) +
    ['postorder'] * len(sres_postorder) +
    ['levelorder'] * len(sres_levelorder) +
    ['is_ultrametric'] * len(sres_ultrametric)
)

biopython_df = {
    'method': all_methods,
    'time': all_times
}


#================================
# 3: plotting and statistics
#================================

print("\n=== biopython Benchmark Statistics ===")
for method_name, method_times in [
    ('preorder', sres_preorder),
    ('postorder', sres_postorder),
    ('levelorder', sres_levelorder),
    ('is_ultrametric', sres_ultrametric)
]:
    print(f"\n{method_name}:")
    print(f"  Mean time (µs): {np.mean(method_times):.4f}")
    print(f"  Median time (µs): {np.median(method_times):.4f}")
    print(f"  Std dev (µs): {np.std(method_times):.4f}")
    print(f"  Min time (µs): {np.min(method_times):.4f}")
    print(f"  Max time (µs): {np.max(method_times):.4f}")

# Create violin plot comparing all methods
fig, ax = plt.subplots(figsize=(10, 6))
data_by_method = [sres_preorder, sres_postorder, sres_levelorder, sres_ultrametric]
positions = [1, 2, 3, 4]
ax.violinplot(data_by_method, positions=positions, widths=0.7)
ax.set_ylabel('Time (microseconds)')
ax.set_title('biopython Tree Traversal Methods Comparison')
ax.set_xticks(positions)
ax.set_xticklabels(['preorder', 'postorder', 'levelorder', 'is_ultrametric'])
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('biopython_traverse.png', dpi=100)
print("\nPlot saved to 'biopython_traverse.png'")
plt.close()
