#!/usr/bin/env python3
"""
Coalescent likelihood benchmark for toytree package.
"""

import time
import math
import numpy as np
import matplotlib.pyplot as plt
import toytree
from ete3 import Tree


FIXED_PHI0 = 1.0
FIXED_R = 0.023882220114924


def benchmark_sys(tree_set, fn, times):
    """Benchmark function execution time over a list of trees."""
    ttaken = []

    fn(tree_set[0])  # warmup

    n = len(tree_set)
    for t in range(times):
        tree = tree_set[t % n]
        tstart = time.perf_counter()
        fn(tree)
        tend = time.perf_counter()
        ttaken.append((tend - tstart) * 1e6)

    return ttaken


def extract_coalescent_intervals(tree):
    """
    Extract coalescent intervals from a toytree.
    """
    # compute root-to-node depths
    depths = {}
    for node in tree.traverse("preorder"):
        if node.is_root():
            depths[node] = 0.0
        else:
            depths[node] = depths[node.up] + (node.dist if node.dist else 0.0)

    # present time = max root-to-tip distance
    tips = [n for n in tree.traverse() if n.is_leaf()]
    present = max(depths[t] for t in tips)

    # build events from native preorder traversal
    events = []
    for node in tree.traverse("preorder"):
        events.append({
            "time": present - depths[node],
            "type": "sample" if node.is_leaf() else "coalescent",
            "name": node.name
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


def loglik_toytree(tree, phi0=FIXED_PHI0, r=FIXED_R):
    """
    Wrapper: Extract intervals and compute the coalescent log-likelihood.
    """
    intervals = extract_coalescent_intervals(tree)
    return coalescent_loglikelihood(intervals, phi0, r)


def load_trees(path):
    """
    Load Newick trees into toytree, using ete3 as a robust intermediate parser.
    """
    trees = []
    with open(path) as fh:
        for line_num, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                trees.append(toytree.tree(line))
            except Exception:
                # fallback: parse with ete3 and re-export normalized Newick
                try:
                    ete_tree = Tree(line, format=1)
                    clean_newick = ete_tree.write(format=1)
                    trees.append(toytree.tree(clean_newick))
                except Exception as e:
                    print(f"Skipping bad line {line_num}: {e}")
                    print(f"  Preview: {line[:200]}")
    return trees


def print_stats(tree_label, results):
    print(f"\n=== toytree Benchmark Statistics Coalescent likelihood ({tree_label}) ===")
    for step, step_times in results.items():
        print(f"\n{step}:")
        print(f"  Mean time (µs): {np.mean(step_times):.4f}")
        print(f"  Median time (µs): {np.median(step_times):.4f}")
        print(f"  Std dev (µs): {np.std(step_times):.4f}")
        print(f"  Min time (µs): {np.min(step_times):.4f}")
        print(f"  Max time (µs): {np.max(step_times):.4f}")


def save_plot(results, title, filename):
    fig, ax = plt.subplots(figsize=(10, 6))
    data_by_method = [results["interval extraction"], results["log likelihood"]]
    positions = [1, 2]
    ax.violinplot(data_by_method, positions=positions, widths=0.7)
    ax.set_ylabel("Time (microseconds)")
    ax.set_title(title)
    ax.set_xticks(positions)
    ax.set_xticklabels(["interval extraction", "log likelihood"])
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=100)
    print(f"\nPlot saved to '{filename}'")
    plt.close()


if __name__ == "__main__":
    print("Loading trees...")
    ultra1k = load_trees("ultra1k.nwk")
    hetero1k = load_trees("hetero1k.nwk")

    sntree = toytree.tree("((t1:0.9383623928,t2:0.5619198801):0.531734891,t3:0.286919398);")

    print("Starting sanity checks...")
    for i in range(5):
        print(f"Coal Likelihood (Ultra#{i+1}):  {loglik_toytree(ultra1k[i]):.2f}")
    for i in range(5):
        print(f"Coal Likelihood (Hetero#{i+1}): {loglik_toytree(hetero1k[i]):.2f}")
    print(f"Test tree logL: {loglik_toytree(sntree, 10, 0.1):.6f}")
    print("Expected value: -4.457106")

    print("\nStarting benchmarking (heterochronic)...")
    coal_results = {
        "interval extraction": benchmark_sys(hetero1k, extract_coalescent_intervals, 1000),
        "log likelihood": benchmark_sys(hetero1k, loglik_toytree, 1000),
    }

    print_stats("heterochronic", coal_results)
    save_plot(coal_results, "toytree Coalescent likelihood", "toytree_coal_lik.png")
