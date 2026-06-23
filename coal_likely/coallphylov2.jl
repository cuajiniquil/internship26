#!/usr/bin/env julia
"""
Phylo.jl coalescent likelihood (serial, rate-growing-into-the-past convention).
Uses native Phylo.jl traversal functions.
"""

using Phylo

const FIXED_PHI0 = 1.0
const FIXED_R    = 0.023882220114924

# ---------------------------------
# Tree loading
# ---------------------------------

function load_newick_trees(file_path::String)
    trees = Any[]
    for line in readlines(file_path)
        stripped = strip(line)
        if !isempty(stripped)
            push!(trees, parsenewick(String(stripped)))
        end
    end
    return trees
end

# ---------------------------------
# Depth computation (manual stack-based preorder)
# ---------------------------------

function compute_node_depths(tree::AbstractTree)
    depths = Dict{Any, Float64}()
    root   = getroot(tree)
    depths[root] = 0.0

    stack = [root]
    while !isempty(stack)
        node = pop!(stack)
        for child in getchildren(tree, node)
            branch = getinbound(tree, child)
            bl     = branch === nothing ? 0.0 : getlength(tree, branch)
            depths[child] = depths[node] + bl
            push!(stack, child)
        end
    end

    return depths
end

# ---------------------------------
# Interval extraction (native preorder traversal for events)
# ---------------------------------

function extract_coalescent_intervals(tree::AbstractTree)
    depths = compute_node_depths(tree)

    # collect leaves using native preorder traversal
    leaves  = [n for n in traversal(tree, preorder) if isleaf(tree, n)]
    present = maximum(depths[l] for l in leaves)

    # Build events from native preorder traversal
    events = []
    for node in traversal(tree, preorder)
        time   = present - depths[node]
        is_tip = isleaf(tree, node)
        push!(events, (time, is_tip, node))
    end

    sort!(events, by = e -> e[1])

    intervals = []
    k = 0
    for i in 1:length(events)-1
        time_i, is_tip_i, _ = events[i]
        time_next, is_tip_next, _ = events[i+1]

        k += is_tip_i ? 1 : -1

        push!(intervals, (
            a         = time_i,
            b         = time_next,
            k         = k,
            ends_coal = !is_tip_next,
        ))
    end

    return intervals
end

# ---------------------------------
# Coalescent log-likelihood
# ---------------------------------

function coalescent_loglikelihood(intervals;
                                  phi0::Float64 = FIXED_PHI0,
                                  r::Float64    = FIXED_R)::Float64
    logL = 0.0
    log_phi0 = log(phi0)

    for iv in intervals
        k, a, b, ends = iv.k, iv.a, iv.b, iv.ends_coal
        c2 = k * (k - 1) / 2.0
        c2 <= 0.0 && continue

        # Survival term: λ(t) = (1/φ0) * exp(r*t)
        if r == 0.0
            logL -= c2 * (b - a) / phi0
        else
            logL -= (c2 / (phi0 * r)) * (exp(r * b) - exp(r * a))
        end

        # Event term only if interval ends in a coalescence
        if ends
            logL += log(c2) - log_phi0 + r * b
        end
    end

    return logL
end

function loglik_phylo(tree::AbstractTree;
                      phi0::Float64 = FIXED_PHI0,
                      r::Float64    = FIXED_R)::Float64
    intervals = extract_coalescent_intervals(tree)
    return coalescent_loglikelihood(intervals; phi0=phi0, r=r)
end

# ---------------------------------
# Main
# ---------------------------------

function main()
    sntree = parsenewick("((t1:0.9383623928,t2:0.5619198801):0.531734891,t3:0.286919398);")
    sntree_ll = loglik_phylo(sntree; phi0=10.0, r=0.1)
    println("sntree logL (phi0=10, r=0.1): ", round(sntree_ll, digits=6))
    println("Expected:                     -4.457106")

    println("\nLoading tree sets...")
    ultra1k  = load_newick_trees("ultra1k.nwk")
    hetero1k = load_newick_trees("hetero1k.nwk")

    println("\nUltra first 5:")
    for i in 1:min(5, length(ultra1k))
        ll = loglik_phylo(ultra1k[i])
        println("  Ultra#$i: ", round(ll, digits=4))
    end

    println("\nHetero first 5:")
    for i in 1:min(5, length(hetero1k))
        ll = loglik_phylo(hetero1k[i])
        println("  Hetero#$i: ", round(ll, digits=4))
    end

    println("\nsntree intervals:")
    for iv in extract_coalescent_intervals(sntree)
        println("  a=$(round(iv.a, digits=6))  b=$(round(iv.b, digits=6))  ",
                "k=$(iv.k)  ends_coal=$(iv.ends_coal)")
    end
end

main()
