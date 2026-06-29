using Phylo

# Check the file
filepath = "ultra1k.nwk"
lines = readlines(filepath)
println("Number of lines: ", length(lines))
println("First line length: ", length(lines[1]))
println("First line preview: ", lines[1][1:min(200, length(lines[1]))])

# Try to parse the first line
try
    tree = parsenewick(lines[1])
    println("Parse OK: nleaves = ", nleaves(tree))
catch e
    println("Parse FAILED: ", e)
end

# Try with strip
try
    tree = parsenewick(strip(lines[1]))
    println("Stripped parse OK: nleaves = ", nleaves(tree))
catch e
    println("Stripped parse FAILED: ", e)
end
