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


#replace with @benchmark 
function benchmark_sys(tree_set, fn::Function, times::Int)::Vector{Float64}
    """Benchmark function execution time over a list of trees."""
    ttaken = Float64[]

    fn(tree_set[1])

    for t in 1:times
        tree = tree_set[mod1(t, length(tree_set))]
        tstart = time_ns()
        fn(tree)
        tend = time_ns()

        elapsed_us = (tend - tstart) / 1000.0
        push!(ttaken, elapsed_us)
    end

    return ttaken
end


function create_balanced_ultrametric_tree(num_taxa::Int)
    return rand(Ultrametric(num_taxa))
end


function traverse_tree_preorder_man(tree::AbstractTree)::Int
    count = 0
    stack = [getroot(tree)]
    while !isempty(stack)
        node = pop!(stack)
        count += 1
        for child in getchildren(tree, node)
            push!(stack, child)
        end
    end
    return count
end


function traverse_tree_postorder_man(tree::AbstractTree)::Int
    """Traverse all nodes in tree (postorder)."""
    count = 0
    stack = [getroot(tree)]
    visited = Set()

    while !isempty(stack)
        node = last(stack)
        if node in visited
            pop!(stack)
            count += 1
        else
            push!(visited, node)
            for child in getchildren(tree, node)
                push!(stack, child)
            end
        end
    end
    return count
end


function traverse_tree_levelorder_man(tree::AbstractTree)::Int
    """Traverse all nodes in tree (levelorder)."""
    count = 0
    queue = [getroot(tree)]

    while !isempty(queue)
        node = popfirst!(queue)
        count += 1
        for child in getchildren(tree, node)
            push!(queue, child)
        end
    end
    return count
end


function traverse_tree_preorder(tree::AbstractTree)::Int
    count = 0
    for _ in traversal(tree, preorder)
        count += 1
    end
    return count
end


function traverse_tree_postorder(tree::AbstractTree)::Int
    count = 0
    for _ in traversal(tree, postorder)
        count += 1
    end
    return count
end


function traverse_tree_levelorder(tree::AbstractTree)::Int
    count = 0
    for _ in traversal(tree, breadthfirst)
        count += 1
    end
    return count
end

#load the function
#read trees (generated in ape)
#do loops 

# Create test trees
println("\nExtracting trees from nwk files...")

function load_newick_trees(file_path::String)
    trees = Any[]
    for line in readlines(file_path)
        stripped = strip(line)
        if !isempty(stripped)
            push!(trees, parsenewick(stripped))
        end
    end
    return trees
end

ultra1k = load_newick_trees("ultra1k.nwk")
hetero1k = load_newick_trees("hetero1k.nwk")

preorder_fn = traverse_tree_preorder
postorder_fn = traverse_tree_postorder
levelorder_fn = traverse_tree_levelorder
ultrametric_fn = is_ultrametric

println("Starting benchmarks...")
sstart = time_ns()

ultra_results = Dict(
    "preorder" => benchmark_sys(ultra1k, preorder_fn, 1000),
    "postorder" => benchmark_sys(ultra1k, postorder_fn, 1000),
    "levelorder" => benchmark_sys(ultra1k, levelorder_fn, 1000),
    "is_ultrametric" => benchmark_sys(ultra1k, ultrametric_fn, 1000),
)

hetero_results = Dict(
    "preorder" => benchmark_sys(hetero1k, preorder_fn, 1000),
    "postorder" => benchmark_sys(hetero1k, postorder_fn, 1000),
    "levelorder" => benchmark_sys(hetero1k, levelorder_fn, 1000),
    "is_ultrametric" => benchmark_sys(hetero1k, ultrametric_fn, 1000),
)

send = time_ns()

benchmark_time_s = (send - sstart) / 1e9
println("Benchmark completed in $(round(benchmark_time_s, digits=2)) seconds")

# ================================
# 3: plotting and statistics
# ================================

function print_stats(tree_label, results)
    println("\n=== Phylo.jl Benchmark Statistics ($tree_label) ===")
    for method_name in ["preorder", "postorder", "levelorder", "is_ultrametric"]
        method_times = results[method_name]
        println("\n$method_name:")
        println("  Mean time (µs): $(round(mean(method_times), digits=4))")
        println("  Median time (µs): $(round(median(method_times), digits=4))")
        println("  Std dev (µs): $(round(std(method_times), digits=4))")
        println("  Min time (µs): $(round(minimum(method_times), digits=4))")
        println("  Max time (µs): $(round(maximum(method_times), digits=4))")
    end
end

function save_plot(results, title::String, filename::String)
    fig, ax = subplots(figsize=(10, 6))
    data_by_method = [results["preorder"], results["postorder"], results["levelorder"], results["is_ultrametric"]]
    positions = [1, 2, 3, 4]
    ax.violinplot(data_by_method, positions=positions, widths=0.7)
    ax.set_ylabel("Time (microseconds)")
    ax.set_title(title)
    ax.set_xticks(positions)
    ax.set_xticklabels(["preorder", "postorder", "levelorder", "is_ultrametric"])
    ax.grid(true, axis="y", alpha=0.3)

    tight_layout()
    savefig(filename, dpi=100)
    println("\nPlot saved to '$filename'")
    close()
end

print_stats("Ultrametric", ultra_results)
print_stats("Heterochronic", hetero_results)
save_plot(ultra_results, "Phylo.jl Tree Traversal Methods Comparison - Ultrametric", "phylo_jl_traverse_ultra.png")
save_plot(hetero_results, "Phylo.jl Tree Traversal Methods Comparison - Heterochronic", "phylo_jl_traverse_hetero.png")
