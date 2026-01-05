module ChemGraphSearch

export Atom, Molecule, Index,
       compile_mol, build_index, search,
       save_index, load_index, read_smi_file


# (paste almost all your code here)

# IMPORTANT: remove `using Test` from the module file.
# IMPORTANT: do NOT keep the CLI `if abspath(PROGRAM_FILE) == @__FILE__` inside the module.

using Serialization

# -----------------------------
# Constants
# -----------------------------
const FP_BITS = 2048
const FP_WORDS = FP_BITS ÷ 64
const MAX_PATH_LEN = 7

const ATOM_SYMBOL_TO_Z = Dict(
    "B" => 5, "C" => 6, "N" => 7, "O" => 8, "F" => 9,
    "P" => 15, "S" => 16, "Cl" => 17, "Br" => 35, "I" => 53
)

const VALENCE_TABLE = Dict(
    5 => [3],        # B
    6 => [4],        # C
    7 => [3, 5],     # N
    8 => [2],        # O
    9 => [1],        # F
    15 => [3, 5],    # P
    16 => [2, 4, 6], # S
    17 => [1],       # Cl
    35 => [1],       # Br
    53 => [1]        # I
)

const BOND_SINGLE   = UInt8(1)
const BOND_DOUBLE   = UInt8(2)
const BOND_TRIPLE   = UInt8(3)
const BOND_AROMATIC = UInt8(4)

# -----------------------------
# Structs
# -----------------------------
struct Atom
    z::UInt8
    charge::Int8
    aromatic::Bool
    h_count::UInt8  # explicit or implicit
end

struct Molecule
    id::String
    natoms::Int
    atoms::Vector{Atom}
    edge_src_offsets::Vector{Int32}  # length natoms+1
    edge_dst::Vector{Int32}
    edge_btype::Vector{UInt8}
    edge_ring::BitVector
    degree::Vector{UInt8}
    valence::Vector{UInt8}
    ringatom::BitVector
    neigh_hash::Vector{UInt32}
    fp::Vector{UInt64}  # FP_WORDS long
end

const Query = Molecule

struct Index
    mols::Vector{Molecule}
end

# -----------------------------
# Hash
# -----------------------------
function simple_hash(x)::UInt32
    x = UInt32(x)
    x = ((x >> 16) ⊻ x) * 0x45d9f3b
    x = ((x >> 16) ⊻ x) * 0x45d9f3b
    x = (x >> 16) ⊻ x
    return x
end

# -----------------------------
# SMILES Parser
# -----------------------------
mutable struct ParserState
    s::String
    pos::Int
    atoms::Vector{Atom}
    adj::Vector{Vector{Tuple{Int, UInt8}}}  # (dst, btype)
    ring_map::Dict{Int, Int}
    branch_stack::Vector{Int}
    prev_atom::Int
    curr_bond::UInt8
    ParserState(s::String) = new(s, 1, Atom[], Vector{Vector{Tuple{Int,UInt8}}}(),
                                Dict{Int,Int}(), Int[], 0, BOND_SINGLE)
end

function parse_smiles(smiles::String)::Tuple{Vector{Atom}, Vector{Vector{Tuple{Int, UInt8}}}}
    state = ParserState(smiles)
    while state.pos <= lastindex(state.s)
        parse_atom_or_branch_or_ring!(state)
    end
    if !isempty(state.branch_stack)
        error("Unclosed branch")
    end
    if !isempty(state.ring_map)
        error("Unclosed rings: $(collect(keys(state.ring_map)))")
    end
    return state.atoms, state.adj
end

function skip_whitespace!(state::ParserState)
    while state.pos <= lastindex(state.s) && isspace(state.s[state.pos])
        state.pos += 1
    end
end

function peek_char(state::ParserState)::Char
    state.pos > lastindex(state.s) ? '\0' : state.s[state.pos]
end

function consume!(state::ParserState)::Char
    c = peek_char(state)
    if c != '\0'
        state.pos += 1
    end
    return c
end

function parse_atom_or_branch_or_ring!(state::ParserState)
    skip_whitespace!(state)
    c = peek_char(state)
    if c == '\0'
        return
    elseif c == '('
        consume!(state)
        push!(state.branch_stack, state.prev_atom)
        return
    elseif c == ')'
        consume!(state)
        isempty(state.branch_stack) && error("Unmatched ) at pos $(state.pos)")
        state.prev_atom = pop!(state.branch_stack)
        return
    elseif c == '.'
        consume!(state)
        state.prev_atom = 0
        state.curr_bond = 0x00
        return
    elseif c == '['
        parse_bracket_atom!(state)
        return
    elseif isdigit(c) || c == '%'
        parse_ring!(state)
        return
    elseif c in ['-', '=', '#', ':']
        parse_bond!(state)
        return
    else
        parse_atom!(state)
        return
    end
end

function parse_bond!(state::ParserState)
    c = consume!(state)
    state.curr_bond = if c == '-'
        BOND_SINGLE
    elseif c == '='
        BOND_DOUBLE
    elseif c == '#'
        BOND_TRIPLE
    elseif c == ':'
        BOND_AROMATIC
    else
        error("Invalid bond $c")
    end
end

function parse_digit!(state::ParserState)::Int
    c = consume!(state)
    isdigit(c) || error("Expected digit, got $c at pos $(state.pos)")
    return Int(c - '0')
end

function parse_ring_number!(state::ParserState)::Int
    if peek_char(state) == '%'
        consume!(state)
        tens = parse_digit!(state)
        ones = parse_digit!(state)
        return 10 * tens + ones
    else
        return parse_digit!(state)
    end
end

function add_bond!(state::ParserState, a1::Int, a2::Int, btype::UInt8)
    push!(state.adj[a1], (a2, btype))
    push!(state.adj[a2], (a1, btype))
end

function parse_ring!(state::ParserState)
    num = parse_ring_number!(state)
    if haskey(state.ring_map, num)
        open_atom = state.ring_map[num]
        add_bond!(state, open_atom, state.prev_atom, state.curr_bond == 0x00 ? BOND_SINGLE : state.curr_bond)
        delete!(state.ring_map, num)
    else
        state.ring_map[num] = state.prev_atom
    end
    state.curr_bond = BOND_SINGLE
end

function parse_atom!(state::ParserState)
    first = consume!(state)
    aromatic = islowercase(first)

    # normalize symbol
    symbol = uppercase(string(first))

    # two-letter only for standard organic subset (Cl, Br)
    if isuppercase(first) && islowercase(peek_char(state))
        candidate = symbol * string(peek_char(state))
        if get(ATOM_SYMBOL_TO_Z, candidate, 0) != 0
            consume!(state)
            symbol = candidate
        end
    end

    z = get(ATOM_SYMBOL_TO_Z, symbol, 0)
    z == 0 && error("Unknown atom $symbol at pos $(state.pos)")

    atom = Atom(UInt8(z), Int8(0), aromatic, UInt8(0))
    push!(state.atoms, atom)
    push!(state.adj, Vector{Tuple{Int,UInt8}}())

    atom_id = length(state.atoms)
    if state.prev_atom > 0 && state.curr_bond > 0x00
        add_bond!(state, state.prev_atom, atom_id, state.curr_bond)
    end
    state.prev_atom = atom_id
    state.curr_bond = BOND_SINGLE
end

function parse_bracket_atom!(state::ParserState)
    consume!(state)  # '['

    # Very simplified bracket parsing:
    # [nH], [O-], [NH3+], [C@@H], etc.
    first = consume!(state)
    first == '\0' && error("Unclosed bracket at end of string")

    aromatic = islowercase(first)
    symbol = uppercase(string(first))

    if isuppercase(first) && islowercase(peek_char(state))
        candidate = symbol * string(peek_char(state))
        if get(ATOM_SYMBOL_TO_Z, candidate, 0) != 0
            consume!(state)
            symbol = candidate
        end
    end

    z = get(ATOM_SYMBOL_TO_Z, symbol, 0)
    z == 0 && error("Unknown bracket atom $symbol")

    h_count::Int = 0
    charge::Int = 0

    while true
        c = peek_char(state)
        c == '\0' && error("Unclosed bracket atom")
        if c == ']'
            consume!(state)
            break
        elseif c == 'H'
            consume!(state)
            if isdigit(peek_char(state))
                h_count = parse_digit!(state)
            else
                h_count = 1
            end
        elseif c == '+' || c == '-'
            consume!(state)
            sign = (c == '+') ? 1 : -1
            if isdigit(peek_char(state))
                charge = sign * parse_digit!(state)
            else
                charge = sign
            end
        else
            # skip chiral, isotope, etc.
            consume!(state)
        end
    end

    atom = Atom(UInt8(z), Int8(charge), aromatic, UInt8(h_count))
    push!(state.atoms, atom)
    push!(state.adj, Vector{Tuple{Int,UInt8}}())

    atom_id = length(state.atoms)
    if state.prev_atom > 0 && state.curr_bond > 0x00
        add_bond!(state, state.prev_atom, atom_id, state.curr_bond)
    end
    state.prev_atom = atom_id
    state.curr_bond = BOND_SINGLE
end

# -----------------------------
# Perception
# -----------------------------
function aromatize_bonds!(atoms::Vector{Atom}, adj::Vector{Vector{Tuple{Int, UInt8}}})
    n = length(atoms)
    for i in 1:n
        for k in eachindex(adj[i])
            (j, b) = adj[i][k]
            if atoms[i].aromatic && atoms[j].aromatic
                # If not explicitly double/triple/aromatic, treat as aromatic
                if b == BOND_SINGLE
                    adj[i][k] = (j, BOND_AROMATIC)
                end
            end
        end
    end
end

function perceive!(atoms::Vector{Atom}, adj::Vector{Vector{Tuple{Int, UInt8}}})
    aromatize_bonds!(atoms, adj)

    n = length(atoms)
    for i in 1:n
        # only assign implicit H if not explicitly specified in bracket
        if atoms[i].h_count == 0
            valences = get(VALENCE_TABLE, atoms[i].z, [0])
            # aromatic bonds count as 1 for valence here
            curr_val = 0
            for (_, btype) in adj[i]
                curr_val += (btype == BOND_AROMATIC) ? 1 : Int(btype)
            end
            min_h = max(0, minimum(valences) - curr_val)
            atoms[i] = Atom(atoms[i].z, atoms[i].charge, atoms[i].aromatic, UInt8(min_h))
        end
    end
    return nothing
end

# -----------------------------
# Ring detection (simple DFS marking)
# -----------------------------
sort_tuple(t::Tuple{Int,Int}) = (t[1] < t[2]) ? t : (t[2], t[1])

function detect_rings(n::Int, adj::Vector{Vector{Int}}, edge_indices::Vector{Tuple{Int,Int}})
    ringatom = falses(n)
    edge_ring = falses(length(edge_indices))
    visited = falses(n)
    parent = fill(-1, n)
    for start in 1:n
        if !visited[start]
            dfs_ring!(start, visited, parent, ringatom, edge_ring, adj, edge_indices)
        end
    end
    return ringatom, edge_ring
end

function dfs_ring!(u, visited, parent, ringatom, edge_ring, adj, edge_indices)
    visited[u] = true
    for v in adj[u]
        if !visited[v]
            parent[v] = u
            dfs_ring!(v, visited, parent, ringatom, edge_ring, adj, edge_indices)
        elseif v != parent[u]
            mark_cycle!(u, v, parent, ringatom, edge_ring, edge_indices)
        end
    end
end

function mark_cycle!(u, v, parent, ringatom, edge_ring, edge_indices)
    path = Int[]
    curr = u
    while curr != v
        push!(path, curr)
        curr = parent[curr]
        if curr == -1
            return
        end
    end
    push!(path, v)

    for a in path
        ringatom[a] = true
    end

    for i in 1:(length(path)-1)
        e = sort_tuple((path[i], path[i+1]))
        edge_idx = findfirst(x -> x == e, edge_indices)
        if edge_idx !== nothing
            edge_ring[edge_idx] = true
        end
    end
end

# -----------------------------
# Fingerprint generation
# -----------------------------
function set_bit!(fp::Vector{UInt64}, bit::Int)
    word = (bit - 1) ÷ 64 + 1
    pos = (bit - 1) % 64
    fp[word] |= (UInt64(1) << pos)
end

function hash_path(path)::Int
    h = UInt32(0)
    for p in path
        h = simple_hash(h ⊻ UInt32(p))
    end
    return (Int(h % UInt32(FP_BITS)) + 1)
end

function dfs_paths(curr::Int, path::Vector{Int}, visited_edges::Vector{Tuple{Int,Int}},
                   adj, atoms, fp::Vector{UInt64})
    if length(path) > MAX_PATH_LEN
        return
    end
    set_bit!(fp, hash_path(path))
    for (neigh, btype) in adj[curr]
        edge = sort_tuple((curr, neigh))
        if edge in visited_edges
            continue
        end
        new_path = copy(path)
        push!(new_path, Int(btype))
        push!(new_path, Int(atoms[neigh].z) * 100 + (atoms[neigh].aromatic ? 1 : 0))
        new_visited = copy(visited_edges)
        push!(new_visited, edge)
        dfs_paths(neigh, new_path, new_visited, adj, atoms, fp)
    end
end

function generate_fp(atoms, adj)::Vector{UInt64}
    fp = zeros(UInt64, FP_WORDS)
    n = length(atoms)
    for start in 1:n
        seed = Int(atoms[start].z) * 100 + (atoms[start].aromatic ? 1 : 0)
        dfs_paths(start, [seed], Tuple{Int,Int}[], adj, atoms, fp)
    end
    return fp
end

function fp_subset(qfp::Vector{UInt64}, tfp::Vector{UInt64})::Bool
    for i in 1:FP_WORDS
        if (qfp[i] & tfp[i]) != qfp[i]
            return false
        end
    end
    return true
end

# -----------------------------
# Compile molecule
# -----------------------------
function compile_mol(id::String, smiles::String)::Molecule
    atoms, adj = parse_smiles(smiles)
    perceive!(atoms, adj)
    n = length(atoms)

    # Build CSR-like half-edge storage: store each undirected edge once under min(a,b)
    edge_dst = Int32[]
    edge_btype = UInt8[]
    edge_ring = BitVector()
    edge_src_offsets = zeros(Int32, n+1)
    edge_indices = Tuple{Int,Int}[]

    for src in 1:n
        edge_src_offsets[src] = Int32(length(edge_dst) + 1)
        for (dst, btype) in sort(adj[src]; by = x -> x[1])
            if src < dst
                push!(edge_dst, Int32(dst))
                push!(edge_btype, btype)
                push!(edge_ring, false) # placeholder, replaced after ring detection
                push!(edge_indices, (src, dst))
            end
        end
    end
    edge_src_offsets[n+1] = Int32(length(edge_dst) + 1)

    # Rings
    adj_list = [ [d for (d, _) in adj[i]] for i in 1:n ]
    ringatom, edge_ring2 = detect_rings(n, adj_list, edge_indices)
    edge_ring = edge_ring2
    
    # Degree / valence
    degree = UInt8[length(adj[i]) for i in 1:n]
    valence = UInt8[
        sum(((b == BOND_AROMATIC) ? 1 : Int(b)) for (_, b) in adj[i]; init=0)
        for i in 1:n
    ]
    # Neighborhood hash
    neigh_hash = zeros(UInt32, n)
    for i in 1:n
        neighs = sort([(atoms[j].z, b, atoms[j].aromatic) for (j, b) in adj[i]])
        h = UInt32(0)
        for (z, b, a) in neighs
            val = UInt32(UInt32(z) * 100 + UInt32(b) * 10 + (a ? 1 : 0))
            h ⊻= simple_hash(val)
        end
        neigh_hash[i] = h
    end

    # Fingerprint
    fp = generate_fp(atoms, adj)

    return Molecule(id, n, atoms, edge_src_offsets, edge_dst, edge_btype,
                    edge_ring, degree, valence, ringatom, neigh_hash, fp)
end

# -----------------------------
# Matcher (VF2-ish)
# -----------------------------
function compatible(qatom::Atom, tatom::Atom,
                    qdeg::UInt8, tdeg::UInt8,
                    qval::UInt8, tval::UInt8,
                    qring::Bool, tring::Bool,
                    qhash::UInt32, thash::UInt32)
    return (qatom.z == tatom.z) &&
           (qatom.charge == tatom.charge) &&
           (qatom.aromatic == tatom.aromatic) &&
           (qdeg <= tdeg) &&
           (qval <= tval) &&
           (!qring || tring)
           # && (qhash == thash)  # optionally enable for extra pruning
end

function get_bond_type(mol::Molecule, a1::Int, a2::Int)::UInt8
    a = min(a1, a2)
    b = max(a1, a2)
    start = Int(mol.edge_src_offsets[a])
    stop  = Int(mol.edge_src_offsets[a+1]) - 1
    for i in start:stop
        if mol.edge_dst[i] == b
            return mol.edge_btype[i]
        end
    end
    return 0x00
end

function edge_in_ring(mol::Molecule, a1::Int, a2::Int)::Bool
    a = min(a1, a2)
    b = max(a1, a2)
    start = Int(mol.edge_src_offsets[a])
    stop  = Int(mol.edge_src_offsets[a+1]) - 1
    for i in start:stop
        if mol.edge_dst[i] == b
            return mol.edge_ring[i]
        end
    end
    return false
end

function get_neighbors(mol::Molecule, src::Int)::Vector{Tuple{Int, UInt8}}
    neigh = Tuple{Int,UInt8}[]
    # edges stored under min endpoint, so gather:
    # 1) edges where src is min endpoint (direct in src segment)
    start = Int(mol.edge_src_offsets[src])
    stop  = Int(mol.edge_src_offsets[src+1]) - 1
    for i in start:stop
        push!(neigh, (Int(mol.edge_dst[i]), mol.edge_btype[i]))
    end
    # 2) edges where src is max endpoint (scan mins < src)
    for a in 1:(src-1)
        s = Int(mol.edge_src_offsets[a])
        e = Int(mol.edge_src_offsets[a+1]) - 1
        for i in s:e
            if Int(mol.edge_dst[i]) == src
                push!(neigh, (a, mol.edge_btype[i]))
            end
        end
    end
    return neigh
end

function check_feasibility(query::Molecule, target::Molecule,
                           mapping::Vector{Int}, used::BitVector,
                           q::Int, t::Int,
                           depth::Int, order::Vector{Int}, candidates::Vector{Vector{Int}})
    for (qn, qb) in get_neighbors(query, q)
        # If neighbor already mapped, enforce bond & ring consistency
        mapped = false
        for d in 1:(depth-1)
            qq = order[d]
            if qq == qn
                tn = mapping[qq]
                tb = get_bond_type(target, tn, t)
                if tb != qb
                    return false
                end
                if edge_in_ring(query, q, qn) && !edge_in_ring(target, t, tn)
                    return false
                end
                mapped = true
                break
            end
        end

        # If not mapped yet, ensure there exists at least one viable future candidate
        if !mapped
            has_cand = false
            for tc in candidates[qn]
                if !used[tc] && get_bond_type(target, tc, t) == qb
                    if !edge_in_ring(query, q, qn) || edge_in_ring(target, t, tc)
                        has_cand = true
                        break
                    end
                end
            end
            if !has_cand
                return false
            end
        end
    end
    return true
end

function vf2_backtrack!(query::Molecule, target::Molecule,
                        order::Vector{Int}, depth::Int,
                        mapping::Vector{Int}, used::BitVector,
                        candidates::Vector{Vector{Int}})::Bool
    if depth > length(order)
        return true
    end
    q = order[depth]
    for t in candidates[q]
        if !used[t] && check_feasibility(query, target, mapping, used, q, t, depth, order, candidates)
            mapping[q] = t
            used[t] = true
            if vf2_backtrack!(query, target, order, depth + 1, mapping, used, candidates)
                return true
            end
            used[t] = false
            mapping[q] = 0
        end
    end
    return false
end

function substructure_match(query::Query, target::Molecule; return_mapping::Bool=false)::Union{Bool, Vector{Int}}
    # Candidate prefilter per query atom
    candidates = [Int[] for _ in 1:query.natoms]
    for q in 1:query.natoms
        for t in 1:target.natoms
            if compatible(query.atoms[q], target.atoms[t],
                          query.degree[q], target.degree[t],
                          query.valence[q], target.valence[t],
                          query.ringatom[q], target.ringatom[t],
                          query.neigh_hash[q], target.neigh_hash[t])
                push!(candidates[q], t)
            end
        end
    end

    # Order query atoms: smallest candidate set first, then higher degree, then ring atoms
    order = sort(1:query.natoms; by = q -> (length(candidates[q]), -Int(query.degree[q]), query.ringatom[q] ? 0 : 1))

    mapping = fill(0, query.natoms)
    used = falses(target.natoms)

    if vf2_backtrack!(query, target, order, 1, mapping, used, candidates)
        return return_mapping ? mapping : true
    end
    return false
end

# -----------------------------
# Index + Search
# -----------------------------
function build_index(smiles_list::Vector{String}, ids::Vector{String})::Index
    mols = Molecule[]
    for (sm, id) in zip(smiles_list, ids)
        try
            mol = compile_mol(id, sm)
            push!(mols, mol)
        catch e
            @warn "Failed to compile $id: $e"
        end
    end
    return Index(mols)
end

function search(index::Index, query_smiles::String; return_mappings::Bool=false)::Vector{Tuple{String, Union{Nothing, Vector{Int}}}}
    query = compile_mol("query", query_smiles)
    results = Tuple{String, Union{Nothing, Vector{Int}}}[]

    # FP prefilter
    candidates = Molecule[]
    for mol in index.mols
        if fp_subset(query.fp, mol.fp)
            push!(candidates, mol)
        end
    end

    # Verify with exact match
    for mol in candidates
        match = substructure_match(query, mol; return_mapping=return_mappings)
        if match !== false
            push!(results, (mol.id, return_mappings ? match : nothing))
        end
    end
    return results
end

# -----------------------------
# Persistence
# -----------------------------

# Pick a safe default directory for data files
function default_data_dir()::String
    dir = joinpath(homedir(), ".chemgraphsearch")
    isdir(dir) || mkpath(dir)
    return dir
end

# If user passes a relative path like "test.idx", write it under default_data_dir()
function normalize_path(path::AbstractString)::String
    p = String(path)
    return isabspath(p) ? p : joinpath(default_data_dir(), p)
end

function save_index(index::Index, path::AbstractString)
    p = normalize_path(path)
    open(p, "w") do io
        serialize(io, index)
    end
    return p
end

function load_index(path::AbstractString)::Index
    p = normalize_path(path)
    open(p) do io
        return deserialize(io)
    end
end

# -----------------------------
# Utility: read SMILES file
# -----------------------------
function read_smi_file(path::String)::Tuple{Vector{String}, Vector{String}}
    ids = String[]
    smiles = String[]
    for line in readlines(path)
        parts = split(strip(line))
        isempty(parts) && continue
        if length(parts) == 1
            push!(smiles, parts[1])
            push!(ids, "mol$(length(ids)+1)")
        else
            push!(ids, parts[1])
            push!(smiles, parts[2])
        end
    end
    return ids, smiles
end
end # module