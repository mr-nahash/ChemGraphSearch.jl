
# ChemGraphSearch.jl

**Pure-Julia lightweight chemical substructure search engine**  
(for small to medium-sized chemical libraries – up to ~200–500k compounds depending on hardware)

Current status (Jan 2026): experimental / research-grade tool  
Good for: proof-of-concept, teaching, small proprietary collections, rapid prototyping  
Not (yet) suitable for: production systems with millions of compounds

## Quick Start (3–5 minutes)

```julia
using Pkg
Pkg.add(url="https://github.com/yourusername/ChemGraphSearch.jl.git")

using ChemGraphSearch

# ────────────────────────────── One-time index creation ──────────────────────────────

# Option A: from .smi file (ChEMBL-style: ID<tab>SMILES or just SMILES)
ids, smiles_list = read_smi_file("compounds_50k.smi")

index = build_index(smiles_list, ids; verbose=true)

save_index(index, "chem_index_50k.idx")   # strongly recommended!

# Option B: already have index file from previous run
index = load_index("chem_index_50k.idx")

# ────────────────────────────── Search examples ──────────────────────────────

# Basic substructure search
results = search(index, "c1ccccc1")           # benzene
results = search(index, "c1ccccc1C(=O)O")     # benzoic acid
results = search(index, "c1ccncc1")           # pyridine

println("Found $(length(results)) matches:")
foreach(r -> println("  • ", r[1]), results)

# Get atom mappings (slower but very useful for visualization/debugging)
results_with_map = search(index, "c1ccccc1NC(=O)C", return_mappings=true)

for (id, mapping) in results_with_map
    println("Match in $id → query atom → target atom:")
    for (q,t) in enumerate(mapping)
        println("  q$q → t$t")
    end
end
```

## Installation

```julia
using Pkg

# Most common ways:
Pkg.add(url="https://github.com/yourusername/ChemGraphSearch.jl.git")
# or for development:
# Pkg.develop(path="/path/to/local/ChemGraphSearch")
```

**Requirements**

- Julia 1.9 – 1.11 (1.10+ recommended)
- No binary dependencies
- RAM: ~1–4 GB per 100k average-sized organic molecules (rough estimate)

## Main API Overview

| Function                          | Purpose                                          | Performance note                     | Returns                             |
|-----------------------------------|--------------------------------------------------|--------------------------------------|-------------------------------------|
| `read_smi_file(path)`             | Read standard .smi file                          | fast                                 | `(ids::Vector{String}, smiles)`     |
| `compile_mol(id, smiles)`         | Parse & prepare single molecule                  | ~5–50 ms/molecule                    | `Molecule` struct                   |
| `build_index(smiles, ids)`        | Create full search index                         | slowest step (~0.5–5 s/molecule)     | `Index`                             |
| `save_index(index, path)`         | Serialize index to disk                          | fast                                 | path (String)                       |
| `load_index(path)`                | Restore index from disk                          | very fast                            | `Index`                             |
| `search(index, query_smiles; …)`  | Main substructure search                         | 1–500 ms/query after FP prefilter    | `Vector{Tuple{String, Union{Nothing,Vector{Int}}}}` |

**Search flags**

```julia
search(index, smiles;
       return_mappings = false,   # atom→atom mapping (much slower)
       verbose         = true     # detailed step-by-step logging
)
```

## Documentation

### 1. Molecule Representation (`Molecule` struct)

```julia
struct Molecule
    id::String
    natoms::Int
    atoms::Vector{Atom}
    edge_src_offsets::Vector{Int32}   # CSR-like format
    edge_dst::Vector{Int32}
    edge_btype::Vector{UInt8}         # 1=single,2=double,3=triple,4=aromatic
    edge_ring::BitVector              # is this edge part of any ring?
    degree::Vector{UInt8}
    valence::Vector{UInt8}            # sum of bond orders (aromatic=1)
    ringatom::BitVector
    neigh_hash::Vector{UInt32}        # local neighborhood hash (for pruning)
    fp::Vector{UInt64}                # 2048-bit path fingerprint
end
```

### 2. Supported Chemistry Features (2026 status)

| Feature                          | Supported? | Comments / Limitations                                                                 |
|----------------------------------|------------|----------------------------------------------------------------------------------------|
| Elements                         | Partial    | B,C,N,O,F,P,S,Cl,Br,I only                                                             |
| Charges                          | Yes        | simple +1/−1, +2/−2, etc. via bracket notation                                        |
| Explicit hydrogens               | Yes        | `[CH3]`, `[NH3+]`, etc.                                                                |
| Aromaticity                      | Basic      | lowercase symbols + Kekulé → aromatic heuristic for 6-membered rings                  |
| Ring perception                  | Yes        | Any size, but mostly used for 5/6 membered + aromatic detection                       |
| Bond types                       | Yes        | single, double, triple, aromatic                                                       |
| Stereochemistry                  | No         | @ / @@ ignored                                                                         |
| Isotopes                         | No         | ignored                                                                                |
| Coordinate bonds / metals        | No         | not designed for organometallics / inorganics                                          |
| SMARTS                           | No         | only plain SMILES queries for now                                                      |

### 3. How the search algorithm works (simplified)

1. **Fingerprint pre-filter**  
   2048-bit hashed path fingerprint (length ≤ 7)  
   Very fast rejection (~95–99% of non-matches eliminated quickly)

2. **VF2-style subgraph isomorphism** (with many heuristics)

   - Atom compatibility: atomic number, charge, aromaticity, degree, valence, ring membership
   - Neighborhood hash pruning
   - Dynamic ordering: query atoms with fewest candidates first
   - Feasibility checking during backtrack (bond & ring consistency)
   - Early pruning using already mapped neighbors

3. **Result**  
   Either simple list of matching IDs or list of `(id, mapping)` where mapping is Vector{Int}  
   (query atom index → target atom index)

### 4. Performance Expectations (rough numbers – 2026 laptops)

| Database size | Index build time | Memory (after build) | Query time (typical) | Query time (hard cases) |
|---------------|------------------|----------------------|----------------------|--------------------------|
| 10,000        | 1–5 min          | ~300–800 MB          | 1–30 ms              | 100–500 ms               |
| 100,000       | 15–90 min        | 2–8 GB               | 5–100 ms             | 0.5–5 s                  |
| 500,000       | 2–10 hours       | 10–40 GB             | 10–500 ms            | 5–60 s                   |


### 5. Strengths & Current Major Limitations (January 2026)

#### Strengths – what this library does really well right now

| Advantage                                      | Why it matters                                                                 | Practical benefit                                      |
|------------------------------------------------|--------------------------------------------------------------------------------|--------------------------------------------------------|
| **Completely pure Julia**                      | No Python/RDKit dependency, no Conda, no binaries                             | Easy deployment, reproducible environments, works on clusters without hassle |
| **Very fast loading & querying once indexed**  | Load once (~seconds), then most queries < 100 ms                              | Excellent for interactive exploration & repeated searches |
| **Good fingerprint pre-filter**                | 2048-bit path fingerprints reject most non-matches extremely quickly          | Saves huge amounts of time on large-ish libraries      |
| **Transparent & hackable code**                | ~1000 lines, relatively clean, no deep magic                                  | Great for teaching, prototyping, and extending        |
| **Atom mapping out of the box**                | Returns query → target atom correspondence                                    | Very useful for visualization, SAR analysis, debugging |
| **Kekulé → aromatic normalization**            | `C1=CC=CC=C1` matches `c1ccccc1`                                              | Much better behavior on real-world datasets            |
| **No external services / cloud required**      | Everything runs locally                                                       | Works offline, in air-gapped environments, privacy-sensitive projects |
| **MIT license**                                | Commercial use allowed without restrictions                                   | Suitable for startups & proprietary screening          |

#### Current Major Limitations & Known Pain Points

| Limitation / Pain Point                        | Severity     | Workaround / Comment                                                                 |
|------------------------------------------------|--------------|--------------------------------------------------------------------------------------|
| **Very slow index creation**                   | High         | Takes 0.5–5 seconds per molecule → hours/days for big libraries                     |
| **High memory usage**                          | High         | ~2–8 GB for 100k molecules, can become problematic >300k                            |
| **Limited element support**                    | Medium–High  | Only B,C,N,O,F,P,S,halogens — many metals, Si, Se, etc. missing                     |
| **No stereochemistry**                         | High         | @ / @@ completely ignored — chiral matches are not distinguished                   |
| **Weak aromaticity model**                     | Medium       | Only simple Kekulé→aromatic + lowercase → may fail on complicated fused systems    |
| **No SMARTS / query features**                 | High         | Only exact SMILES substructure — no ring-size, atom lists, recursive queries, etc. |
| **No multi-threading**                         | Medium–High  | Both index build and search are single-threaded                                     |
| **Scalability ceiling**                        | Medium–High  | Practical limit ~200–600k compounds depending on hardware & molecule size          |
| **No similarity / Tanimoto search**            | Medium       | Only yes/no substructure — no ranking by similarity                                 |
| **Macrocycles & unusual valences can hurt**    | Medium       | Ring detection & valence perception may give wrong results in edge cases           |

### Quick Summary Table (how people usually compare tools)

| Criterion                     | ChemGraphSearch.jl | RDKit (Python) | Open Babel | Conclusion for this library |
|-------------------------------|--------------------|----------------|------------|------------------------------|
| Speed of single query         | ★★★★☆             | ★★★★★          | ★★★★☆      | Very good after indexing     |
| Index creation speed          | ★☆☆☆☆             | ★★★★☆          | ★★★☆☆      | Major weakness               |
| Memory efficiency             | ★★☆☆☆             | ★★★★☆          | ★★★★☆      | Quite hungry                 |
| Ease of deployment            | ★★★★★             | ★★☆☆☆          | ★★★☆☆      | Excellent advantage          |
| Customizability / hacking     | ★★★★★             | ★★★☆☆          | ★★☆☆☆      | Very strong                  |
| Feature completeness          | ★★☆☆☆             | ★★★★★          | ★★★★☆      | Still very early stage       |
| Best for                      | Prototyping, teaching, small–medium local collections, Julia-native workflows | Production screening, very large libraries, full-featured cheminformatics | Quick & dirty conversion & search | — |

**Bottom line (Jan 2026)**  
If you need something **fast to start**, **easy to understand/modify**, **completely Julia-native**, and you work with **< 200–300k compounds** → this library can be surprisingly pleasant to use.

If you need **million-scale databases**, **stereo**, **SMARTS**, **very fast indexing**, or **broad element coverage** → you should currently prefer RDKit (via PythonCall.jl if you must stay in Julia).

Hope this balanced view helps people decide whether to try it or not! 🧪

### 6. Development Roadmap – 2026/2027 ideas (not ordered)

- [ ] Multi-threaded index construction
- [ ] Tanimoto similarity search (on same fingerprints)
- [ ] Better aromaticity perception (Hückel / more ring types)
- [ ] Support more elements (Si, Se, etc.)
- [ ] Basic SMARTS subset (ring size, connectivity queries…)
- [ ] Query caching / batch searching
- [ ] Integration with Julia visualization packages (Makie, GraphMakie?)
- [ ] Memory-efficient storage (maybe memory-mapped index?)

## Contributing

We especially welcome:

- Bug reports with problematic SMILES
- Performance improvement ideas/PRs
- Additional test cases (especially macrocycles, tautomers, charged species, weird valences)
- Better documentation / examples

## License

MIT License

Use freely in academic, personal and commercial projects.

---

**Questions / ideas / horror SMILES?**  
Open an issue — we love tricky molecules! 🧪
```

This version tries to balance three goals:

1. Quick-start remains very easy for newcomers
2. Detailed technical documentation for people who want to understand/modify/extend the code
3. Honest communication of current limitations (very important for credibility in cheminformatics)

Feel free to adjust the level of optimism/pessimism and the timeline according to your actual plans. 😄
