# ChemGraphSearch.jl

**Pure Julia · Lightweight chemical substructure search engine**
*Designed for small-to-medium molecular libraries (≈ up to 100k–300k compounds)*

ChemGraphSearch.jl is a **Julia-native cheminformatics engine** for **substructure search on molecular graphs**.

It combines:

* **Path-based fingerprints (2048 bits)** for fast candidate pruning
* **VF2-style subgraph isomorphism** for exact chemical matching
* **Automatic Kekulé ↔ aromatic normalization**
* **Pharma-friendly “generalized” matching modes**

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

# ------------------------------------------------------------
# 1) Small demo dataset
# ------------------------------------------------------------
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

println("Building demo index...")
idx = build_index(smiles, ids; verbose=false)
println("Index contains ", length(idx.mols), " molecules\n")

# ------------------------------------------------------------
# Helper: normalize results (SearchHit OR tuple) into (id, mapping)
# ------------------------------------------------------------
function _hit_id_mapping(hit)
    # New API: SearchHit(id=..., mapping=...)
    if hasproperty(hit, :id)
        id = getproperty(hit, :id)
        mapping = hasproperty(hit, :mapping) ? getproperty(hit, :mapping) : nothing
        return id, mapping
    end

    # Old API: (id, mapping)
    if hit isa Tuple && length(hit) == 2
        return hit[1], hit[2]
    end

    # Fallback
    return string(hit), nothing
end

function hit_ids(results)
    [first(_hit_id_mapping(h)) for h in results]
end

function show_results(title, results; show_mapping=false)
    println(title)
    if isempty(results)
        println("  (no matches)\n")
        return
    end

    for h in results
        id, mapping = _hit_id_mapping(h)
        if show_mapping && mapping !== nothing
            println("  • ", id, "  | mapping = ", mapping)
        else
            println("  • ", id)
        end
    end
    println()
end

# ------------------------------------------------------------
# 2) Query: benzene ring
# ------------------------------------------------------------
query = "c1ccccc1"

# ------------------------------------------------------------
# 3) EXACT mode (strict chemistry)
#    Aromatic carbon must stay carbon
# ------------------------------------------------------------
res_exact = search(idx, query; mode=ChemGraphSearch.EXACT, verbose=false)
show_results("EXACT mode (strict element matching):", res_exact)

# ------------------------------------------------------------
# 4) GENERALIZED mode (pharma-friendly)
#    Aromatic C can match aromatic C or aromatic N
# ------------------------------------------------------------
res_gen = search(idx, query; mode=ChemGraphSearch.GENERALIZED, verbose=false)
show_results("GENERALIZED mode (scaffold-like matching):", res_gen)

extra = setdiff(hit_ids(res_gen), hit_ids(res_exact))
println("New hits in GENERALIZED: ",
        isempty(extra) ? "(none)" : join(extra, ", "))
println()

# ------------------------------------------------------------
# 5) Atom mappings (optional)
# ------------------------------------------------------------
res_map = search(idx, query;
    mode=ChemGraphSearch.GENERALIZED,
    return_mappings=true,
    verbose=false
)
show_results("GENERALIZED + atom mappings:", res_map; show_mapping=true)

# ------------------------------------------------------------
# 6) Save & reload index
# ------------------------------------------------------------
save_index(idx, "demo_index.idx")
idx2 = load_index("demo_index.idx")

res_reload = search(idx2, query; mode=ChemGraphSearch.GENERALIZED, verbose=false)
show_results("GENERALIZED after reload:", res_reload)

println("Done.")
---

## EXACT vs GENERALIZED — Why it Matters

| Mode            | Meaning                              | Typical Use                                     |
| --------------- | ------------------------------------ | ----------------------------------------------- |
| **EXACT**       | Element-strict substructure matching | Patents, filtering, formal substructure queries |
| **GENERALIZED** | Pharma-style scaffold matching       | Hit expansion, SAR, scaffold hopping            |

Example:

* Benzene **does not** match pyridine in EXACT mode
* Benzene **does** match pyridine in GENERALIZED mode

This mirrors how medicinal chemists actually reason about scaffolds.

---

## Ready-to-Run Examples

All examples live in `examples/`:

```bash
cd examples
julia --project=../.. 01_basic_search.jl
```

| File                       | What it demonstrates                       | Best for                  |
| -------------------------- | ------------------------------------------ | ------------------------- |
| `01_basic_search.jl`       | Simple substructure queries                | First steps               |
| `02_mappings.jl`           | Atom-to-atom mapping                       | Visualization & SAR       |
| `03_realistic_workflow.jl` | Read `.smi` → index → save → load → search | Real-world usage          |
| `04_kekule_aromatic.jl`    | Kekulé ↔ aromatic normalization            | Understanding equivalence |
| `05_read_smi_file.jl`      | Loading ChEMBL-style `.smi` files          | Your own datasets         |

Each example is **self-contained** and uses the correct project environment.

---

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
