#!/usr/bin/env julia
# Run with: ./chemgraph-demo.jl   (after chmod +x) or julia chemgraph-demo.jl

using Pkg
Pkg.activate(@__DIR__, "..") do
    using ChemGraphSearch

    println("""
    ChemGraphSearch.jl interactive demo
    ==================================

    Commands:
      load FILE       Load existing index
      build FILE      Build index from .smi file
      search SMILES   Search for substructure
      help            Show this message
      exit            Quit

    Example:
      search c1ccccc1
    """)

    # Add simple REPL loop here if you want (optional)
    println("Quick test: searching benzene in built-in molecules...")
    # ... you can copy-paste one of the examples above
end