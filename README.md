# ChemGraphSearch.jl
ChemGraphSearch.jl is a small, pure-Julia substructure search engine for SMILES-based molecule libraries. It builds a local index using fast path-based fingerprints to prune candidates, then verifies matches with a VF2-style subgraph isomorphism search (with optional atom mappings). It also includes automatic Kekulé ↔ aromatic normalization so equivalent representations match reliably. The project is intentionally lightweight and hackable—ideal for teaching, prototyping, and Julia-native workflows on small-to-medium collections.

No external dependencies · MIT licensed · Fully hackable Julia codebase

---

## Why ChemGraphSearch?

ChemGraphSearch is built for situations where you want:

* Full **control** over chemical matching logic
* **Transparent**, inspectable algorithms
* No heavyweight dependencies (RDKit, OpenBabel, Java, C++)
* Tight integration with **Julia workflows**
* A **teachable**, research-friendly codebase

It is especially useful for:

* Medicinal chemistry prototyping
* Scaffold exploration
* Teaching chemical graph algorithms
* Small proprietary compound collections
* Algorithm development & experimentation

---
## What happens when you search?

1. ChemGraphSearch compiles each SMILES into a graph (atoms/bonds).
2. It computes a 2048-bit fingerprint to quickly reject most molecules.
3. It runs a VF2-style exact subgraph match only on the remaining candidates.
4. Optionally, it returns an atom mapping (query atom → target atom indices).

## Current Status (January 2026)

⚠️ **Experimental / research / prototyping grade**

ChemGraphSearch is **not** a drop-in replacement for RDKit, but it is:

* Stable enough for real work on small–medium libraries
* Actively developed
* Well-suited for learning and experimentation
* Easy to modify and extend

---

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/mr-nahash/ChemGraphSearch.jl.git")

using ChemGraphSearch
```

### Quick smoke test

```julia
using ChemGraphSearch
println("Fingerprint size: ", ChemGraphSearch.FP_BITS)
```

Expected output:

```
Fingerprint size: 2048
```

---

## Quick Demo – Copy & Paste (< 60 seconds)

This demo shows:

* Index construction
* Exact vs pharma-friendly search
* Kekulé ↔ aromatic equivalence
* Saving & reloading an index

```julia
using ChemGraphSearch

smiles = [
    "c1ccccc1", "c1ccncc1", "O=c1ccccc1", "C1=CC=CC=C1",
    "C1CCCCC1", "CCO", "c1ccc(cc1)O", "c1ccc2ccccc2c1"
]
ids = ["benzene","pyridine","benzaldehyde","benzene_kekule","cyclohexane","ethanol","phenol","naphthalene"]

println("Building index (once)...")
idx = build_index(smiles, ids; verbose=false)
println("✓ indexed ", length(idx.mols), " molecules\n")

query = "c1ccccc1"
println("Query = ", query, "\n")

function show_hits(title, hits; mappings=false)
    println(title)
    isempty(hits) && return println("  (no matches)\n")
    for h in hits
        if mappings && hasproperty(h, :mapping) && h.mapping !== nothing
            println("  • ", h.id, " | mapping = ", h.mapping)
        else
            println("  • ", h.id)
        end
    end
    println()
end

hits_exact = search(idx, query; mode=ChemGraphSearch.EXACT, verbose=false)
show_hits("EXACT (strict element matching):", hits_exact)

hits_gen = search(idx, query; mode=ChemGraphSearch.GENERALIZED, verbose=false)
show_hits("GENERALIZED (aromatic C can match aromatic N):", hits_gen)

hits_map = search(idx, query; mode=ChemGraphSearch.GENERALIZED, return_mappings=true, verbose=false)
show_hits("GENERALIZED + atom mappings:", hits_map; mappings=true)

save_index(idx, "demo_index.idx")
idx2 = load_index("demo_index.idx")
hits_reload = search(idx2, query; mode=ChemGraphSearch.GENERALIZED, verbose=false)
show_hits("After reload (GENERALIZED):", hits_reload)

println("Done.")
```
## Core Usage (Cheat Sheet)

```julia
# Build index (slowest step)
idx = build_index(smiles, ids)

# Save for later
save_index(idx, "my_collection.idx")

# Reload instantly
idx = load_index("my_collection.idx")

# Search
search(idx, "c1ccccc1")
search(idx, "c1ccncc1C(=O)O"; return_mappings=true)
```

---

## Testing

```bash
julia --project=test -e 'using Pkg; Pkg.test("ChemGraphSearch")'
```

---

## Realistic Expectations (January 2026)

### What works well

* Pure Julia, no dependencies
* Fast queries after indexing (often 1–200 ms)
* Deterministic atom mapping
* Kekulé/aromatic equivalence
* Ideal for teaching & prototyping
* Fully inspectable internals

### Known limitations

* Indexing is slow (≈ 0.5–5 s per molecule)
* High memory usage for very large libraries
* Limited element set (organic subset)
* No stereochemistry
* SMILES only (no SMARTS)
* No similarity search (yet)

---

## Contributing

Contributions are very welcome, especially:

* Edge-case SMILES & bug reports
* Performance improvements
* Macrocycles, tautomers, charged systems
* Documentation & teaching material
* Algorithmic ideas

MIT License — free for academic, personal, and commercial use.

**If you have weird molecules, we want them.**
Open an issue — pathological SMILES are our favorite 🧪
