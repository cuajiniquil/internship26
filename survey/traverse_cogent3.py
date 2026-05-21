#!/usr/bin/env python3
"""
Tree traversal benchmark for cogent3 package.
Generates ultrametric trees and measures traversal time.
"""

import time
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
non_ultra_newick = "((A:5,B:5):2,C:7);"
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


# Create test tree with smaller size for cogent3 (tree construction can be slow)
np.random.seed(273)
print("\nGenerating large ultrametric tree for benchmark...")

# Note: For cogent3, using a smaller tree size as recursive tree building is slower
test_tree = create_balanced_ultrametric_tree(1000)

print("Starting benchmarks...")
sstart = time.perf_counter()
sres_preorder = benchmark_sys(lambda: traverse_tree(test_tree), 1000)
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

cogent3_df = {
    'method': all_methods,
    'time': all_times
}


#================================
# 3: plotting and statistics
#================================

print("\n=== cogent3 Benchmark Statistics ===")
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
ax.set_title('cogent3 Tree Traversal Methods Comparison')
ax.set_xticks(positions)
ax.set_xticklabels(['preorder', 'postorder', 'levelorder', 'is_ultrametric'])
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('cogent3_traverse.png', dpi=100)
print("\nPlot saved to 'cogent3_traverse.png'")
plt.close()
