# Show atom mapping
using ChemGraphSearch

smiles = ["c1ccccc1OC", "c1ccccc1"]  # anisole + benzene
ids = ["anisole", "benzene"]

idx = build_index(smiles, ids)

println("Searching for methoxybenzene pattern: c1ccccc1OC")
results = search(idx, "c1ccccc1OC"; return_mappings=true)

for (id, mapping) in results
    println("\nMatch found in: $id")
    println("Mapping (query atom → target atom):")
    for (q, t) in enumerate(mapping)
        q == 0 && continue
        println("  q$q → t$t")
    end
end