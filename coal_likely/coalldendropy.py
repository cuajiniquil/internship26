#!/usr/bin/env python3
"""
Tree traversal benchmark for dendropy package.
Generates ultrametric trees and measures traversal time.
"""

import dendropy
import time
import math
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


# Create test trees
print("\nExtracting trees from nwk files...")

def load_newick_trees(file_path):
    trees = []
    with open(file_path) as handle:
        for line in handle:
            line = line.strip()
            if line:
                trees.append(dendropy.Tree.get(data=line, schema="newick"))
    return trees


ultra1k = load_newick_trees("ultra1k.nwk")
hetero1k = load_newick_trees("hetero1k.nwk")

print("\nStarting benchmarks...")
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
# 3: Coalescent Likelihood
#================================

FIXED_PHI0 = 1.0
FIXED_R    = 0.023882220114924

def extract_coalescent_intervals(tree):
    tree.calc_node_root_distances(return_leaf_distances_only=False)
    present = max(nd.root_distance for nd in tree.leaf_node_iter())

    events = []
    for nd in tree.leaf_node_iter():
        events.append({
            "time": present - nd.root_distance,
            "type": "sample",
            "name": nd.taxon.label if nd.taxon else None
        })
    for nd in tree.postorder_internal_node_iter():
        events.append({
            "time": present - nd.root_distance,
            "type": "coalescent",
            "name": None
        })

    events.sort(key=lambda e: e["time"])

    intervals = []
    n_active = 0
    for i, ev in enumerate(events[:-1]):
        n_active += 1 if ev["type"] == "sample" else -1
        intervals.append({
            "a": ev["time"],
            "b": events[i + 1]["time"],
            "k": n_active,
            "ends_coal": events[i + 1]["type"] == "coalescent"
        })

    return intervals


def coalescent_loglikelihood(intervals, phi0=FIXED_PHI0, r=FIXED_R):
    """
    Compute coalescent log-likelihood from intervals.

    intervals: iterable of objects with attributes t0, t1, k
    """
    logL = 0.0
    log_phi0 = math.log(phi0)

    for iv in intervals:
        k, a, b, ends = iv["k"], iv["a"], iv["b"], iv["ends_coal"]
        c2 = k * (k - 1) / 2.0
        if c2 <= 0:
            continue

        # survival
        if r == 0:
            logL -= c2 * (b - a) / phi0
        else:
            logL -= (c2 / (phi0 * r)) * (
                math.exp(r * b) - math.exp(r * a)
            )

        # event
        if ends:
            logL += math.log(c2) - log_phi0 + r * b
            
    return logL


def loglik_dendropy(tree, phi0=FIXED_PHI0, r=FIXED_R):
    """
    Wraper: 
    Extract intervals using the existing dendropy extractor and compute
    the coalescent log-likelihood.
    """
    intervals = extract_coalescent_intervals(tree)
    return coalescent_loglikelihood(intervals, phi0, r)


coalikely_results = {
    "intv_extract": benchmark_sys(hetero1k, loglik_dendropy, 1000),
}

liketest = loglik_dendropy(hetero1k[0])

#================================
# 4: plotting and statistics
#================================


def print_stats(tree_label, results):
    print(f"\n=== dendropy Benchmark Statistics ({tree_label}) ===")
    for method_name in ['preorder', 'postorder', 'levelorder', 'is_ultrametric']:
        method_times = results[method_name]
        print(f"\n{method_name}:")
        print(f"  Mean time (µs): {np.mean(method_times):.4f}")
        print(f"  Median time (µs): {np.median(method_times):.4f}")
        print(f"  Std dev (µs): {np.std(method_times):.4f}")
        print(f"  Min time (µs): {np.min(method_times):.4f}")
        print(f"  Max time (µs): {np.max(method_times):.4f}")

def print_step(tree_label, results):
    print(f"\n=== cogent3 Benchmark Statistics ({tree_label}) ===")
    for step in ["intv_extract"]:
        step_times = results[step]
        print(f"\n{step}:")
        print(f"  Mean time (µs): {np.mean(step_times):.4f}")
        print(f"  Median time (µs): {np.median(step_times):.4f}")
        print(f"  Std dev (µs): {np.std(step_times):.4f}")
        print(f"  Min time (µs): {np.min(step_times):.4f}")
        print(f"  Max time (µs): {np.max(step_times):.4f}")

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
print_step('interval extraction', coalikely_results)
times = np.asarray(coalikely_results["intv_extract"])

mx = int(np.argmax(times))
print(f"Max time: {times[mx]:.4f} µs at index {mx}")

slow_tree = hetero1k[mx]
#print(slow_tree)

print("Coalescent likelihood test:")
for i in range (0,5):
    print(loglik_dendropy(hetero1k[i]))

sntree = dendropy.Tree.get(
        data="((t1:0.9383623928,t2:0.5619198801):0.531734891,t3:0.286919398);",
        schema="newick")
print("Sanity check")
print(loglik_dendropy(sntree, 10, 0.1))

save_plot(ultra_results, 'dendropy Tree Traversal Methods Comparison - Ultrametric', 'dendropy_traverse_ultra.png')
save_plot(hetero_results, 'dendropy Tree Traversal Methods Comparison - Heterochronic', 'dendropy_traverse_hetero.png')
