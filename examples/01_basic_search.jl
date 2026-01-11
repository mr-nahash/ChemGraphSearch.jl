# examples/01_basic_search.jl
# Run with: julia --project=../.. 01_basic_search.jl

using ChemGraphSearch

println("=== Basic substructure search example ===\n")

# Small built-in dataset for demo
smiles = [
    "c1ccccc1",           # benzene
    "CCO",                # ethanol
    "O=c1ccccc1",         # benzaldehyde
    "c1ccncc1",           # pyridine
    "C1CCCCC1"            # cyclohexane
]

ids = ["mol$i" for i in 1:length(smiles)]

println("Building small demo index...")
idx = build_index(smiles, ids; verbose=true)

println("\nSearching for benzene (c1ccccc1)...")
results = search(idx, "c1ccccc1"; verbose=true)

println("\nResults:")
for (id, _) in results
    println("  • $id")
end

println("\nDone! You can run this file directly:")
println("    julia --project=../.. 01_basic_search.jl")