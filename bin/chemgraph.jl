#!/usr/bin/env julia
using ChemGraphSearch

ids, smiles = read_smi_file(ARGS[1])
idx = build_index(smiles, ids)
save_index(idx, ARGS[2])
