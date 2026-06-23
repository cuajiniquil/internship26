import math
from cogent3 import make_tree

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


def extract_coalescent_intervals_ultrametric(tree):
    """
    Reference-style extraction: treat all tips as sampled at t=0.
    Only coalescence events change k.
    """
    depths = node_depths(tree)
    present = max(depths[t] for t in tree.tips())

    # Coalescence events only: internal nodes sorted by time before present.
    coal_events = sorted(
        ({"time": present - depths[node], "node": node}
         for node in tree.preorder() if not node.is_tip()),
        key=lambda e: e["time"],
    )

    n_tips = len(list(tree.tips()))
    k = n_tips
    logL_terms = []

    # Interval from present (t=0) to first coalescence
    if not coal_events:
        return []

    intervals = []
    prev_time = 0.0
    for ev in coal_events:
        intervals.append({
            "a": prev_time,
            "b": ev["time"],
            "k": k,
            "ends_coal": True,
        })
        k -= 1
        prev_time = ev["time"]

    return intervals


def coalescent_loglikelihood(intervals):
    logL = 0.0
    log_phi0 = math.log(FIXED_PHI0)
    for iv in intervals:
        k, a, b = iv["k"], iv["a"], iv["b"]
        c2 = k * (k - 1) / 2.0
        if c2:
            logL += math.log(c2) - log_phi0 - FIXED_R * b
        logL -= c2 * (math.exp(-FIXED_R * a) - math.exp(-FIXED_R * b)) / FIXED_PHI0
    return logL


def loglik_cogent3(tree):
    return coalescent_loglikelihood(extract_coalescent_intervals_ultrametric(tree))


def load_trees(path):
    with open(path) as fh:
        return [make_tree(line.strip()) for line in fh if line.strip()]


if __name__ == "__main__":
    ultra1k  = load_trees("ultra1k.nwk")
    hetero1k = load_trees("hetero1k.nwk")

    for i in range(5):
        print(f"Coal Likelihood (Ultra#{i+1}):  {loglik_cogent3(ultra1k[i]):.2f}")
    for i in range(5):
        print(f"Coal Likelihood (Hetero#{i+1}): {loglik_cogent3(hetero1k[i]):.2f}")
