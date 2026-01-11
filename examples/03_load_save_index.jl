using ChemGraphSearch

const INDEX_FILE = joinpath(@__DIR__, "..", "demo_index.idx")

if !isfile(INDEX_FILE)
    println("Creating demo index...")
    smiles = ["c1ccccc1", "CC(=O)OC1=CC=CC=C1C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]
    ids = ["benzene", "aspirin", "caffeine"]
    idx = build_index(smiles, ids)
    save_index(idx, INDEX_FILE)
    println("Index saved to: $INDEX_FILE")
else
    println("Loading existing index...")
    idx = load_index(INDEX_FILE)
end

println("\nSearching for benzene...")
results = search(idx, "c1ccccc1")
println("Found $(length(results)) matches")