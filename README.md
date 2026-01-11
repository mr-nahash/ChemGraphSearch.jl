Below is an **enhanced, clearer, newcomer-friendly version** of your README section.
I’ve kept your technical accuracy intact, but improved:

* **Narrative flow**
* **First-time reader friendliness**
* **Clear mental model (what / why / when)**
* **Reduced intimidation**
* **More explicit pharma relevance**

You can copy-paste this directly over your current README.

---

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
println("Index contains ", length(idx.mols), " molecules")

# ------------------------------------------------------------
# Helper: pretty printing
# ------------------------------------------------------------
function show_results(title, results)
    println("\n" * title)
    isempty(results) && return println("  (no matches)")
    for (name, mapping) in results
        println(mapping === nothing ? "  • $name" :
                "  • $name | mapping = $mapping")
    end
end

# ------------------------------------------------------------
# 2) Query: benzene ring
# ------------------------------------------------------------
query = "c1ccccc1"

# ------------------------------------------------------------
# 3) EXACT mode (strict chemistry)
#    Aromatic carbon must stay carbon
# ------------------------------------------------------------
res_exact = search(idx, query; mode=EXACT, verbose=false)
show_results("EXACT search (strict element matching):", res_exact)

# ------------------------------------------------------------
# 4) GENERALIZED mode (pharma-friendly)
#    Aromatic C can match aromatic C or N
# ------------------------------------------------------------
res_gen = search(idx, query; mode=GENERALIZED, verbose=false)
show_results("GENERALIZED search (scaffold-like matching):", res_gen)

# ------------------------------------------------------------
# 5) Atom mappings (optional)
# ------------------------------------------------------------
res_map = search(idx, query;
                 mode=GENERALIZED,
                 return_mappings=true,
                 verbose=false)

show_results("GENERALIZED search with atom mappings:", res_map)

# ------------------------------------------------------------
# 6) Save & reload index
# ------------------------------------------------------------
save_index(idx, "demo_index.idx")
idx2 = load_index("demo_index.idx")

res2 = search(idx2, query; mode=GENERALIZED, verbose=false)
show_results("GENERALIZED search after reload:", res2)

println("\nDone.")
```

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
