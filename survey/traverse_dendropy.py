#!/usr/bin/env python3
"""
Tree traversal benchmark for dendropy package.
Generates ultrametric trees and measures traversal time.
"""

import dendropy
import time
import numpy as np
import matplotlib.pyplot as plt
from dendropy.simulate import treesim


#================================
#1: initial tests
#================================

# Generate an ultrametric tree (coalescent)
print("Generating initial test trees...")
trs = treesim.pure_kingman_tree(
    taxon_namespace=dendropy.TaxonNamespace(
        ["taxon_{}".format(i) for i in range(1000)]
    ),
    num_lineages=1000
)

# Check if tree is ultrametric
def is_ultrametric(tree):
    """Check if all leaf nodes are equidistant from root."""
    if not tree.seed_node:
        return True
    
    def get_leaf_distances(node, dist=0.0):
        if node.is_leaf():
            return [dist]
        distances = []
        for child in node.child_node_iter():
            child_dist = dist + (child.edge_length if child.edge_length else 0)
            distances.extend(get_leaf_distances(child, child_dist))
        return distances
    
    distances = get_leaf_distances(tree.seed_node)
    if len(distances) <= 1:
        return True
    return len(set(distances)) == 1


# Create a non-ultrametric tree for comparison
non_ultra_str = "((A:5,B:5):2,C:7);"
non_ultra = dendropy.Tree.get(data=non_ultra_str, schema="newick")

print("Ultrametric test tree is ultrametric:", is_ultrametric(trs))
print("Non-ultrametric test tree is ultrametric:", is_ultrametric(non_ultra))

# Print tree structure
print("\nTest tree root:", trs.seed_node)
print("Test tree number of leaves:", len(trs.leaf_nodes()))


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


def create_ultrametric_tree(num_taxa):
    """Create an ultrametric tree using dendropy."""
    taxon_ns = dendropy.TaxonNamespace(
        ["taxon_{}".format(i) for i in range(num_taxa)]
    )
    return treesim.pure_kingman_tree(
        taxon_namespace=taxon_ns,
        num_lineages=num_taxa
    )


def traverse_tree(tree):
    """Traverse all nodes in tree."""
    for node in tree.preorder_node_iter():
        pass


# Create test tree
np.random.seed(273)
test_tree = create_ultrametric_tree(100000)

print("\nStarting benchmarks...")
sstart = time.perf_counter()
sres = benchmark_sys(lambda: traverse_tree(test_tree), 1000)
send = time.perf_counter()

print(f"Benchmark completed in {send - sstart:.2f} seconds")

dendropy_df = {
    'traverse': ['tree traversal (dendropy)'] * len(sres),
    'time': sres
}


#================================
# 3: plotting and statistics
#================================

print("\n=== dendropy Benchmark Statistics ===")
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
ax.set_title('dendropy Tree Traversal Benchmark')
ax.set_xticks([1])
ax.set_xticklabels(['tree traversal (dendropy)'])
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('dendropy_traverse.png', dpi=100)
print("\nPlot saved to 'dendropy_traverse.png'")
plt.close()
