# ChemGraphSearch.jl

**Pure Julia • Lightweight chemical substructure search engine**  
(for small-to-medium sized chemical libraries – up to ~100–300k compounds)

Path-based fingerprints (2048 bits) + VF2-style subgraph isomorphism matching  
Automatic Kekulé ↔ aromatic normalization  
No external dependencies · MIT licensed

**Current status (January 2026):**  
Experimental / research / prototyping grade tool  
Great for teaching, small proprietary collections, Julia-native workflows

## Installation

```julia
using Pkg
 Pkg.add(url="https://github.com/mr-nahash/ChemGraphSearch.jl.git")

using ChemGraphSearch
```

### Quick smoke test after installation

```julia
using ChemGraphSearch
println("Fingerprint size: ", ChemGraphSearch.FP_BITS)   # should print 2048
```

## Quick Demo – Copy & Paste in < 60 seconds

```julia
using ChemGraphSearch

# Tiny built-in example dataset
smiles = [
    "c1ccccc1",           # benzene
    "CCO",                # ethanol
    "O=c1ccccc1",         # benzaldehyde
    "c1ccncc1",           # pyridine
    "C1CCCCC1"            # cyclohexane
]
ids = ["benzene", "ethanol", "benzald", "pyridine", "cyclohexane"]

println("Building small demo index...")
idx = build_index(smiles, ids)

println("\nSearching for any benzene ring...")
results = search(idx, "c1ccccc1")

println("\nFound matches:")
foreach(r -> println("  • ", r[1]), results)
# Expected output includes: benzene, benzald, pyridine
```

Save for next time (much faster):

```julia
save_index(idx, "demo_index.idx")
idx = load_index("demo_index.idx")   # ← reload in < 1 second
```

## Ready-to-Run Examples

All examples are in the `examples/` folder — run them like this:

```bash
cd examples
julia --project=../.. 01_basic_search.jl
```

| File                        | What it shows                                      | Best for                          |
|-----------------------------|----------------------------------------------------|-----------------------------------|
| `01_basic_search.jl`        | Multiple simple queries on small dataset           | First steps                       |
| `02_mappings.jl`            | Atom-by-atom mapping between query & target        | Visualization / SAR analysis      |
| `03_realistic_workflow.jl`  | Read .smi → build → save → load → search           | Real-world typical usage          |
| `04_kekule_aromatic.jl`     | Kekulé form matches aromatic query automatically   | Understanding normalization       |
| `05_read_smi_file.jl`       | Loading real ChEMBL-style .smi files               | Working with your own data        |

Every example is self-contained and uses the correct project environment.

## Core Usage in One Table

```julia
# Build once (slowest step)
idx = build_index(smiles_vector, ids_vector; verbose=true)

# Persistence (highly recommended!)
save_index(idx, "my_collection.idx")
idx = load_index("my_collection.idx")          # fast!

# Search
search(idx, "c1ccccc1")                        # → list of matching ids
search(idx, "c1ccncc1C(=O)O", return_mappings=true)  # → ids + atom mappings
```

## Run the Tests

```bash
# From project root
julia --project=test -e 'using Pkg; Pkg.test("ChemGraphSearch")'

# or directly:
julia --project=test test/runtests.jl
```

## Realistic Expectations – January 2026

### What currently works very well

- Completely pure Julia → easy install, modify, deploy
- Fast queries after indexing (usually 1–200 ms)
- Atom mapping included by default
- Kekulé/aromatic equivalence handling
- No cloud/external services needed
- Great for teaching, prototyping, small–medium local collections

### Current main limitations

- Index creation is slow (~0.5–5 seconds per molecule)
- Memory usage grows quickly (>5–10 GB for 200k+ molecules)
- Only basic organic elements supported (B,C,N,O,F,P,S,halogens)
- No stereochemistry handling
- No SMARTS (only plain SMILES queries)
- No similarity/Tanimoto search (yet)

## Contributing

We especially welcome:

- Bug reports with tricky/problematic SMILES
- Performance improvement ideas
- Additional test cases (macrocycles, unusual charges, tautomers…)
- Better documentation & more examples

MIT License — free for academic, personal and commercial use.

**Feedback & weird molecules are very welcome!**  
Open an issue — we love collecting challenging SMILES cases 🧪
