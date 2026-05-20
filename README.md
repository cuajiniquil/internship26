# internship26

There are two angles to your project that I like and which I would like you to explore over the next few weeks. The first, is a survey of methods to manipulate and visualise trees programmatically. If you get a chance, have a look at ggtree, ETE3, dendropy, phylo.jl, Biopython, and there are possibly others. Ideally we would like some sort of summary of what they do and what they don't. Then we could write some small tests where we simulate very large trees and obtain the time each takes to read it, traverse it, obtain annotations, etc...

The goal is to benchmark how long it takes to run a script similar to this NELSI one: 

library(NELSI)
tr <- rtree(1000)
tr$tip.label <- paste(tr$tip.label, "_", sample(c("A", "B)")
find.monophyletic(tr, "_A")

They are first tested with a small tree (ape.tree), and then tested with a large tree; the benchmarking should be separated by function (reading, traversing,...) so ideally those are all done separetely to make the benchmarking easier

Tools tested in this survey:
    -ggtree
    -ape
    -ete3
    -dendropy 
    -biopython
    -toytree
    -cogent3
    -phylo.jl 
    -gotree (go - ligne cmd developpe ici) 
    -bio++ (c++)
    
CURRENT GOAL:    
use rcoal(10000) to generate ULTRAMETRIC coalescent trees
then benchmark is.ultrametric 


find equivalent scripts/functions with  the rest of the tools 
-tr$edge -> to find the nodes? but what were testing id voal
