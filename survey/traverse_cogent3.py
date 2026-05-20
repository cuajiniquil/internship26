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
            child_dist = dist + (child.branch_length if child.branch_length else 0)
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
            return make_tree(taxa_list[0])
        
        mid = len(taxa_list) // 2
        left_tree = build_balanced(taxa_list[:mid], branch_len)
        right_tree = build_balanced(taxa_list[mid:], branch_len)
        
        # Combine with equal branch lengths to maintain ultrametricity
        left_str = str(left_tree) if hasattr(left_tree, '__str__') else left_tree
        right_str = str(right_tree) if hasattr(right_tree, '__str__') else right_tree
        
        combined = f"({left_str}:{branch_len},{right_str}:{branch_len});"
        return make_tree(combined)
    
    taxa_names = [f"taxon_{i}" for i in range(num_taxa)]
    return build_balanced(taxa_names)


def traverse_tree(tree):
    """Traverse all nodes in tree."""
    count = 0
    for node in tree.traverse(include_root=True):
        count += 1
    return count


# Create test tree with smaller size for cogent3 (tree construction can be slow)
np.random.seed(273)
print("\nGenerating large ultrametric tree for benchmark...")

# Note: For cogent3, using a smaller tree size as recursive tree building is slower
test_tree = create_balanced_ultrametric_tree(1000)

print("Starting benchmarks...")
sstart = time.perf_counter()
sres = benchmark_sys(lambda: traverse_tree(test_tree), 1000)
send = time.perf_counter()

print(f"Benchmark completed in {send - sstart:.2f} seconds")

cogent3_df = {
    'traverse': ['tree traversal (cogent3)'] * len(sres),
    'time': sres
}


#================================
# 3: plotting and statistics
#================================

print("\n=== cogent3 Benchmark Statistics ===")
print(f"Mean time (µs): {np.mean(sres):.4f}")
print(f"Median time (µs): {np.median(sres):.4f}")
print(f"Std dev (µs): {np.std(sres):.4f}")
print(f"Min time (µs): {np.min(sres):.4f}")
print(f"Max time (µs): {np.max(sres):.4f}")

mx_idx = np.argmax(sres)
print(f"Slowest traversal (µs): {sres[mx_idx]:.4f}")

# Create violin plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.violinplot(sres, positions=[1], widths=0.5)
ax.set_ylabel('Time (microseconds)')
ax.set_title('cogent3 Tree Traversal Benchmark')
ax.set_xticks([1])
ax.set_xticklabels(['tree traversal (cogent3)'])
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('cogent3_traverse.png', dpi=100)
print("\nPlot saved to 'cogent3_traverse.png'")
plt.close()
