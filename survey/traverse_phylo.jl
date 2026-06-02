#!/usr/bin/env julia
"""
Tree traversal benchmark for Phylo.jl package.
Generates ultrametric trees and measures traversal time.
"""

using Phylo
using StatsBase
using PyPlot
import Random

# ================================
# 1: initial tests
# ================================

println("Generating initial test trees with Phylo.jl...")

# Create an ultrametric tree from Newick format
ultra_newick = "((A:1,B:1):1,(C:1,D:1):1);"
trs = parsenewick(ultra_newick)

# Create a non-ultrametric tree for comparison
non_ultra_newick = "(A:1,B:5);"
notU = parsenewick(non_ultra_newick)

function is_ultrametric(tree::AbstractTree)::Bool
    """
    Check if all leaf nodes are equidistant from root.
    """
    function get_leaf_distances(node, dist=0.0)
        if isleaf(tree, node)
            return [dist]
        end
        distances = Float64[]
        for child in getchildren(tree, node)
            branch = getinbound(tree, child)
            child_dist = dist + getlength(tree, branch)
            append!(distances, get_leaf_distances(child, child_dist))
        end
        return distances
    end

    distances = get_leaf_distances(getroot(tree))
    if length(distances) <= 1
        return true
    end
    # Check if all distances are approximately equal
    return maximum(distances) - minimum(distances) < 1e-10
end


println("Ultrametric test tree is ultrametric: $(is_ultrametric(trs))")
println("Non-ultrametric test tree is ultrametric: $(is_ultrametric(notU))")

# Print tree info
println("\nTest tree leaf count: $(nleaves(trs))")


# ================================
# 2: tree traversing benchmarks
# ================================

function benchmark_sys(fn::Function, times::Int)::Vector{Float64}
    """
    Benchmark function execution time.
    Returns: vector of execution times in microseconds (excluding warmup)
    """
    ttaken = Float64[]

    for t in 1:(times + 1)
        tstart = time_ns()
        fn()
        tend = time_ns()

        elapsed_us = (tend - tstart) / 1000.0
        push!(ttaken, elapsed_us)
    end

    return ttaken[2:end]
end


function create_balanced_ultrametric_tree(num_taxa::Int)
    return rand(Ultrametric(num_taxa))
end


function traverse_tree_preorder(tree::AbstractTree, stack::Vector)::Int
    empty!(stack)
    push!(stack, getroot(tree))
    count = 0
    while !isempty(stack)
        node = pop!(stack)
        count += 1
        for child in getchildren(tree, node)
            push!(stack, child)
        end
    end
    return count
end


function traverse_tree_postorder(tree::AbstractTree, stack::Vector, output::Vector)::Int
    empty!(stack)
    empty!(output)
    push!(stack, getroot(tree))
    while !isempty(stack)
        node = pop!(stack)
        push!(output, node)
        for child in getchildren(tree, node)
            push!(stack, child)
        end
    end
    for node in Iterators.reverse(output)
        # visit
    end
    return length(output)
end


function traverse_tree_levelorder(tree::AbstractTree, queue::Vector)::Int
    empty!(queue)
    push!(queue, getroot(tree))
    count = 0
    while !isempty(queue)
        node = popfirst!(queue)
        count += 1
        for child in getchildren(tree, node)
            push!(queue, child)
        end
    end
    return count
end


# Create test tree
println("\nExtracting trees from nwk files...")
# Trees generated in R:
ultra1k = readlines("ultra1k.nwk")
parsedu1k = parsenewick(ultra1k[1])
test_tree = parsedu1k

# Pre-allocate buffers — reused across all benchmark runs, zero heap allocation per call
preorder_stack   = []
postorder_stack  = []
postorder_output = []
levelorder_queue = []

preorder_fn    = () -> traverse_tree_preorder(test_tree, preorder_stack)
postorder_fn   = () -> traverse_tree_postorder(test_tree, postorder_stack, postorder_output)
levelorder_fn  = () -> traverse_tree_levelorder(test_tree, levelorder_queue)
ultrametric_fn = () -> is_ultrametric(test_tree)

println("Starting benchmarks...")
sstart = time_ns()
sres_preorder    = benchmark_sys(preorder_fn, 1000)
sres_postorder   = benchmark_sys(postorder_fn, 1000)
sres_levelorder  = benchmark_sys(levelorder_fn, 1000)
sres_ultrametric = benchmark_sys(ultrametric_fn, 1000)
send = time_ns()

benchmark_time_s = (send - sstart) / 1e9
println("Benchmark completed in $(round(benchmark_time_s, digits=2)) seconds")

# ================================
# 3: plotting and statistics
# ================================

println("\n=== Phylo.jl Benchmark Statistics ===")
for (method_name, method_times) in [
    ("preorder", sres_preorder),
    ("postorder", sres_postorder),
    ("levelorder", sres_levelorder),
    ("is_ultrametric", sres_ultrametric)
]
    println("\n$method_name:")
    println("  Mean time (µs): $(round(mean(method_times), digits=4))")
    println("  Median time (µs): $(round(median(method_times), digits=4))")
    println("  Std dev (µs): $(round(std(method_times), digits=4))")
    println("  Min time (µs): $(round(minimum(method_times), digits=4))")
    println("  Max time (µs): $(round(maximum(method_times), digits=4))")
end

# Create violin plot comparing all methods
fig, ax = subplots(figsize=(10, 6))
data_by_method = [sres_preorder, sres_postorder, sres_levelorder, sres_ultrametric]
positions = [1, 2, 3, 4]
ax.violinplot(data_by_method, positions=positions, widths=0.7)
ax.set_ylabel("Time (microseconds)")
ax.set_title("Phylo.jl Tree Traversal Methods Comparison")
ax.set_xticks(positions)
ax.set_xticklabels(["preorder", "postorder", "levelorder", "is_ultrametric"])
ax.grid(true, axis="y", alpha=0.3)

tight_layout()
savefig("phylo_jl_traverse.png", dpi=100)
println("\nPlot saved to 'phylo_jl_traverse.png'")
close()