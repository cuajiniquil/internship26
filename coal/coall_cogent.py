import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import math
from cogent3 import make_tree

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
    
FIXED_PHI0 = 1.0
FIXED_R    = 0.023882220114924

def node_depths(tree):
    depths = {}
    def walk(node, depth):
        depths[node] = depth
        for child in node.children:
            edge = child.length if child.length is not None else 0.0
            walk(child, depth + edge)
    walk(tree, 0.0)
    return depths


def extract_coalescent_intervals(tree):
    depths = {tree: 0.0}

    for node in tree.preorder():
        for child in node.children:
            depths[child] = depths[node] + (child.length or 0.0)

    present = max(depths[t] for t in tree.tips())

    events = []
    for node in tree.preorder():
        events.append({
            "time": present - depths[node],
            "type": "sample" if node.is_tip() else "coalescent",
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


def loglik_cogent3(tree, phi0=FIXED_PHI0, r=FIXED_R):
    """
    Wraper: Extract intervals using the existing cogent3 extractor and compute
    the coalescent log-likelihood.
    """
    intervals = extract_coalescent_intervals(tree)
    return coalescent_loglikelihood(intervals, phi0, r)


def load_trees(path):
    with open(path) as fh:
        return [make_tree(line.strip()) for line in fh if line.strip()]

def print_stats(tree_label, results):
    print(f"\n=== cogent3 Benchmark Statistics Coalescent likelihood ({tree_label}) ===")
    for step, step_times in results.items():
        print(f"\n{step}:")
        print(f"  Mean time (µs): {np.mean(step_times):.4f}")
        print(f"  Median time (µs): {np.median(step_times):.4f}")
        print(f"  Std dev (µs): {np.std(step_times):.4f}")
        print(f"  Min time (µs): {np.min(step_times):.4f}")
        print(f"  Max time (µs): {np.max(step_times):.4f}")

 
       
def save_plot(results, title, filename):
    fig, ax = plt.subplots(figsize=(10, 6))
    data_by_method = [results['interval extraction'], results['log likelihood']]
    positions = [1, 2]
    ax.violinplot(data_by_method, positions=positions, widths=0.7)
    ax.set_ylabel('Time (microseconds)')
    ax.set_title(title)
    ax.set_xticks(positions)
    ax.set_xticklabels(['interval extraction', 'log likelihood'])
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=100)
    print(f"\nPlot saved to '{filename}'")
    plt.close()

if __name__ == "__main__":
    print("Loading trees...")
    ultra1k  = load_trees("ultra1k.nwk")
    hetero1k = load_trees("hetero1k.nwk")
    sntree = make_tree("((t1:0.9383623928,t2:0.5619198801):0.531734891,t3:0.286919398);")
    print("Starting sanity checks...")
    for i in range(5):
        print(f"Coal Likelihood (Ultra#{i+1}):  {loglik_cogent3(ultra1k[i]):.2f}")
    for i in range(5):
        print(f"Coal Likelihood (Hetero#{i+1}): {loglik_cogent3(hetero1k[i]):.2f}")
    print(f"Test tree logL: {loglik_cogent3(sntree,10,0.1)}")
    print("Expected value : -4.457106")
    print("Starting benchmarking...")
    coal_results = {
    "interval extraction": benchmark_sys(hetero1k, extract_coalescent_intervals, 1000),
    "log likelihood": benchmark_sys(hetero1k, loglik_cogent3, 1000),
    }
    print_stats('interval extraction', coal_results)
    save_plot(coal_results, 'cogent3 Coalescent likelihood', 'cogent3_coal_lik.png')

