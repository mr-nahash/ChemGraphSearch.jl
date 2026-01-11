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
<div align="center">
  <img width="400" height="500"
       alt="ChemSearchDiagram"
       src="https://github.com/user-attachments/assets/9fc981cf-f79c-4ca8-aaa7-c30897ca072f" />
</div>

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

---

## Quick Demo – Copy & Paste (< 60 seconds)

This demo shows:

* Index construction
* Exact vs pharma-friendly search
* Kekulé ↔ aromatic equivalence
* Saving & reloading an index

```julia
using ChemGraphSearch

# ============================================================
# What you’ll learn:
#   1) Build a tiny index from SMILES
#   2) Search a substructure (benzene ring)
#   3) EXACT vs GENERALIZED matching
#   4) Optional atom mappings (query atom → target atom)
#   5) Save & reload an index
# ============================================================

println("\n== 1) Build a tiny demo index ==")

smiles = [
    "c1ccccc1",            # benzene
    "c1ccncc1",            # pyridine (one N in aromatic ring)
    "O=c1ccccc1",          # benzaldehyde
    "C1=CC=CC=C1",         # benzene (Kekulé form)
    "C1CCCCC1",            # cyclohexane (non-aromatic)
    "CCO",                 # ethanol
    "c1ccc(cc1)O",         # phenol
    "c1ccc2ccccc2c1"       # naphthalene
]

ids = [
    "benzene",
    "pyridine",
    "benzaldehyde",
    "benzene_kekule",
    "cyclohexane",
    "ethanol",
    "phenol",
    "naphthalene"
]

# TIP: verbose=true prints internal steps; keep false for a clean demo.
idx = build_index(smiles, ids; verbose=false)
println("✓ Built index with ", length(idx.mols), " molecules")

println("\n== 2) Define a query substructure ==")
query = "c1ccccc1"  # benzene ring in aromatic form
println("Query SMILES: ", query)

# ------------------------------------------------------------
# Helper: print results in a beginner-friendly way
# (works with SearchHit objects; mappings optional)
# ------------------------------------------------------------
function show_hits(title, hits; show_mapping=false)
    println("\n" * title)
    if isempty(hits)
        println("  (no matches)")
        return
    end
    for h in hits
        # SearchHit usually has fields: id, mapping (mapping may be nothing)
        if show_mapping && hasproperty(h, :mapping) && h.mapping !== nothing
            println("  • ", h.id, " | mapping = ", h.mapping)
        else
            println("  • ", h.id)
        end
    end
end

println("\n== 3) EXACT mode (strict chemistry) ==")
println("In EXACT mode, aromatic carbon must match aromatic carbon (C ≠ N).")
hits_exact = search(idx, query; mode=ChemGraphSearch.EXACT, verbose=false)
show_hits("Matches (EXACT):", hits_exact)

println("\n== 4) GENERALIZED mode (scaffold-like) ==")
println("In GENERALIZED mode, aromatic C in the *query* may match aromatic C or aromatic N.")
println("That means a benzene query can also match pyridine.")
hits_gen = search(idx, query; mode=ChemGraphSearch.GENERALIZED, verbose=false)
show_hits("Matches (GENERALIZED):", hits_gen)

println("\n== 5) OPTIONAL: return atom mappings ==")
println("Atom mapping shows which target atoms correspond to each query atom index.")
hits_map = search(idx, query; mode=ChemGraphSearch.GENERALIZED, return_mappings=true, verbose=false)
show_hits("Matches (GENERALIZED + mappings):", hits_map; show_mapping=true)

println("\n== 6) Save & reload the index ==")
save_index(idx, "demo_index.idx")
idx2 = load_index("demo_index.idx")

hits_reload = search(idx2, query; mode=ChemGraphSearch.GENERALIZED, verbose=false)
show_hits("Matches after reload (GENERALIZED):", hits_reload)

println("\nDone ✅")
```

## Core Usage (Cheat Sheet)

```julia
## Core Usage (Cheat Sheet)

using ChemGraphSearch

# ------------------------------------------------------------
# 1) Build an index (slowest step — do this once)
# ------------------------------------------------------------
idx = build_index(smiles, ids; verbose=false)     # smiles::Vector{String}, ids::Vector{String}

# ------------------------------------------------------------
# 2) Save / reload (recommended)
# ------------------------------------------------------------
save_index(idx, "my_collection.idx")              # saved under ~/.chemgraphsearch/ unless absolute path
idx = load_index("my_collection.idx")             # reload fast

# ------------------------------------------------------------
# 3) Search (two matching modes)
# ------------------------------------------------------------

# EXACT (default-like): strict element matching
hits = search(idx, "c1ccccc1"; mode=ChemGraphSearch.EXACT, verbose=false)

# GENERALIZED: scaffold-like (aromatic query C can match aromatic C or N)
hits = search(idx, "c1ccccc1"; mode=ChemGraphSearch.GENERALIZED, verbose=false)

# ------------------------------------------------------------
# 4) Search + atom mappings (query atom → target atom index)
# ------------------------------------------------------------
hits_map = search(idx, "c1ccncc1C(=O)O";
                  mode=ChemGraphSearch.EXACT,
                  return_mappings=true,
                  verbose=false)

# hits / hits_map are vectors of SearchHit objects.
# Each hit has: h.id  and (optionally) h.mapping
for h in hits_map
    println(h.id, "  mapping=", h.mapping)
end
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
