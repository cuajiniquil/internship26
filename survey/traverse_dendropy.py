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
    )
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
    return np.allclose(distances, distances[0])


# Create a non-ultrametric tree for comparison
non_ultra_str = "(A:1,B:5);"
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
        taxon_namespace=taxon_ns
    )


def traverse_tree(tree):
    """Traverse all nodes in tree (preorder)."""
    for node in tree.preorder_node_iter():
        pass


def traverse_tree_postorder(tree):
    """Traverse all nodes in tree (postorder)."""
    for node in tree.postorder_node_iter():
        pass


def traverse_tree_levelorder(tree):
    """Traverse all nodes in tree (level order)."""
    for node in tree.levelorder_node_iter():
        pass


# Create test tree
print("\nExtracting trees from nwk files...")
# Trees generated in R:
ultra1k = open("ultra1k.nwk").readlines()
parsedu1k = dendropy.Tree.get(data=ultra1k[0], schema="newick")
test_tree = parsedu1k

print("\nStarting benchmarks...")
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

dendropy_df = {
    'method': all_methods,
    'time': all_times
}


#================================
# 3: plotting and statistics
#================================

print("\n=== dendropy Benchmark Statistics ===")
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
ax.set_title('dendropy Tree Traversal Methods Comparison')
ax.set_xticks(positions)
ax.set_xticklabels(['preorder', 'postorder', 'levelorder', 'is_ultrametric'])
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('dendropy_traverse.png', dpi=100)
print("\nPlot saved to 'dendropy_traverse.png'")
plt.close()
