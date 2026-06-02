#!/usr/bin/env python3
"""
Tree traversal benchmark for ete3 package.
Generates ultrametric trees and measures traversal time.
"""

import time
import numpy as np
import matplotlib.pyplot as plt
from ete3 import Tree


#================================
#1: initial tests
#================================

print("Generating initial test trees with ete3...")

# Create an ultrametric tree (Newick format)
ultra_newick = "((A:1,B:1):1,(C:1,D:1):1);"
trs = Tree(ultra_newick)

# Create a non-ultrametric tree for comparison
non_ultra_newick = "(A:1,B:5);"
notU = Tree(non_ultra_newick)

def is_ultrametric_slow(tree):
    """
    Check if all leaf nodes are equidistant from root.
    """
    distances = []
    for leaf in tree.get_leaves():
        dist = tree.get_distance(leaf)
        distances.append(dist)
    if len(distances) <= 1:
        return True
    return np.allclose(distances, distances[0])

#hgih search time fix:
def is_ultrametric(tree):
    ref = None
    # stack contient (noeud, distance_accumulée_depuis_racine)
    stack = [(tree, tree.dist)]
    while stack:
        node, dist = stack.pop()
        if node.is_leaf():
            if ref is None:
                ref = dist
            elif not np.isclose(dist, ref):
                return False  # early exit possible sur arbre non-ultrametrique
        else:
            for child in node.children:
                stack.append((child, dist + child.dist))
    return True
    
print("Ultrametric test tree is ultrametric:", is_ultrametric(trs))
print("Non-ultrametric test tree is ultrametric:", is_ultrametric(notU))

# Print tree info
print("\nTest tree structure:")
print(trs.get_ascii(show_internal=True))
print("Test tree leaf count:", len(trs.get_leaves()))


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
    Create a balanced ultrametric tree using ete3.
    Recursively builds a binary tree with equal branch lengths.
    """
    def build_balanced(taxa_list, branch_len=1.0):
        if len(taxa_list) == 1:
            return f"{taxa_list[0]}"

        mid = len(taxa_list) // 2
        left_tree = build_balanced(taxa_list[:mid], branch_len)
        right_tree = build_balanced(taxa_list[mid:], branch_len)

        combined = f"({left_tree}:{branch_len},{right_tree}:{branch_len})"
        return combined

    taxa_names = [f"taxon_{i}" for i in range(num_taxa)]
    newick_str = build_balanced(taxa_names) + ";"
    return Tree(newick_str)


def traverse_tree_preorder(tree):
    """Traverse all nodes in tree (preorder)."""
    count = 0
    for node in tree.traverse("preorder"):
        count += 1
    return count


def traverse_tree_postorder(tree):
    """Traverse all nodes in tree (postorder)."""
    count = 0
    for node in tree.traverse("postorder"):
        count += 1
    return count


def traverse_tree_levelorder(tree):
    """Traverse all nodes in tree (levelorder)."""
    count = 0
    for node in tree.traverse("levelorder"):
        count += 1
    return count


# Create test tree
print("\nExtracting trees from nwk files...")
# Trees genereated in R:
ultra1k = open("ultra1k.nwk").readlines()
parsedu1k = Tree(ultra1k[0])

print("Starting benchmarks...")
sstart = time.perf_counter()
sres_preorder = benchmark_sys(lambda: traverse_tree_preorder(parsedu1k), 1000)
sres_postorder = benchmark_sys(lambda: traverse_tree_postorder(parsedu1k), 1000)
sres_levelorder = benchmark_sys(lambda: traverse_tree_levelorder(parsedu1k), 1000)
sres_ultrametric = benchmark_sys(lambda: is_ultrametric(parsedu1k), 1000)
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

ete3_df = {
    'method': all_methods,
    'time': all_times
}


#================================
# 3: plotting and statistics
#================================

print("\n=== ete3 Benchmark Statistics ===")
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
ax.set_title('ete3 Tree Traversal Methods Comparison')
ax.set_xticks(positions)
ax.set_xticklabels(['preorder', 'postorder', 'levelorder', 'is_ultrametric'])
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('ete3_traverse.png', dpi=100)
print("\nPlot saved to 'ete3_traverse.png'")
plt.close()
