module ChemGraphSearch

export Atom, Molecule, Index, MatchMode, SearchHit,
       compile_mol, build_index, search,
       save_index, load_index, read_smi_file,
       hit_ids, show_hits, EXACT, GENERALIZED

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

# Match mode:
# - EXACT: element-typed substructure (default; benzene won't match pyridine)
# - GENERALIZED: pharma-friendly scaffold relaxations (benzene can match pyridine)
@enum MatchMode begin
    EXACT
    GENERALIZED
end

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
    fp_exact::Vector{UInt64}  # FP_WORDS long
    fp_gen::Vector{UInt64}    # FP_WORDS long (generalized)
end

const Query = Molecule

struct Index
    mols::Vector{Molecule}
end

"""
A single search match.

- `id` is the target molecule id.
- `mapping` is either `nothing` or a vector mapping query atom index -> target atom index.
"""
struct SearchHit
    id::String
    mapping::Union{Nothing, Vector{Int}}
end

"Return just the molecule ids from a vector of hits."
hit_ids(hits::Vector{SearchHit}) = [h.id for h in hits]

"""
Pretty-print search results in a user-friendly way.
"""
function show_hits(hits::Vector{SearchHit}; title::AbstractString="Matches")
    println(title, " (", length(hits), "):")
    isempty(hits) && return println("  (no matches)")
    for h in hits
        if h.mapping === nothing
            println("  • ", h.id)
        else
            println("  • ", h.id, " | mapping = ", h.mapping)
        end
    end
    return nothing
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

function parse_smiles(smiles::String; verbose::Bool=true)::Tuple{Vector{Atom}, Vector{Vector{Tuple{Int, UInt8}}}}
    if verbose
        println("Starting SMILES parsing for: $smiles")
    end
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
    if verbose
        println("SMILES parsing completed: $(length(state.atoms)) atoms found")
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
function aromatize_bonds!(atoms::Vector{Atom}, adj::Vector{Vector{Tuple{Int, UInt8}}}; verbose::Bool=true)
    if verbose
        println("Aromatizing bonds for aromatic atoms")
    end
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

function perceive!(atoms::Vector{Atom}, adj::Vector{Vector{Tuple{Int, UInt8}}}; verbose::Bool=true)
    aromatize_bonds!(atoms, adj; verbose=verbose)

    if verbose
        println("Perceiving implicit hydrogens based on valence")
    end

    n = length(atoms)
    for i in 1:n
        # Only assign implicit H if not explicitly specified in bracket
        if atoms[i].h_count == 0
            expected_valence = 0

            if atoms[i].aromatic && atoms[i].z == UInt8(6)
                # Aromatic carbon: treat as valence 3
                expected_valence = 3
            elseif atoms[i].aromatic && atoms[i].z == UInt8(7)
                # Aromatic nitrogen:
                # - "n" (pyridine-like) has no H unless explicitly [nH]
                # - we still use valence 3 for counting purposes
                expected_valence = 3
            else
                valences = get(VALENCE_TABLE, atoms[i].z, [0])
                expected_valence = minimum(valences)
            end

            curr_val = sum(((btype == BOND_AROMATIC) ? 1 : Int(btype)) for (_, btype) in adj[i]; init=0)
            min_h = max(0, expected_valence - curr_val)

            atoms[i] = Atom(atoms[i].z, atoms[i].charge, atoms[i].aromatic, UInt8(min_h))
        end
    end

    return nothing
end


# -----------------------------
# Ring detection (simple DFS marking)
# -----------------------------
sort_tuple(t::Tuple{Int,Int}) = (t[1] < t[2]) ? t : (t[2], t[1])

function detect_rings(n::Int, adj::Vector{Vector{Int}}, edge_indices::Vector{Tuple{Int,Int}}; verbose::Bool=true)
    if verbose
        println("Detecting rings using DFS")
    end
    ringatom = falses(n)
    edge_ring = falses(length(edge_indices))
    visited = falses(n)
    parent = fill(-1, n)
    for start in 1:n
        if !visited[start]
            dfs_ring!(start, visited, parent, ringatom, edge_ring, adj, edge_indices)
        end
    end
    if verbose
        println("Ring detection complete: $(sum(ringatom)) ring atoms, $(sum(edge_ring)) ring edges")
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

# Mode-dependent atom code for fingerprints
function fp_atom_code(a::Atom, mode::MatchMode)::Int
    if mode == GENERALIZED && a.aromatic && (a.z == UInt8(6) || a.z == UInt8(7))
        # collapse aromatic C/N into a single code so generalized queries don't get pruned
        return 600 + 1
    end
    return Int(a.z) * 100 + (a.aromatic ? 1 : 0)
end

function dfs_paths(curr::Int, path::Vector{Int}, visited_edges::Vector{Tuple{Int,Int}},
                   adj, atoms, fp::Vector{UInt64}, mode::MatchMode)
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
        push!(new_path, fp_atom_code(atoms[neigh], mode))
        new_visited = copy(visited_edges)
        push!(new_visited, edge)
        dfs_paths(neigh, new_path, new_visited, adj, atoms, fp, mode)
    end
end

function generate_fp(atoms, adj; verbose::Bool=true, mode::MatchMode=EXACT)::Vector{UInt64}
    if verbose
        println("Generating fingerprint using path hashing (up to length $MAX_PATH_LEN) [mode=$(mode)]")
    end
    fp = zeros(UInt64, FP_WORDS)
    n = length(atoms)
    for start in 1:n
        seed = fp_atom_code(atoms[start], mode)
        dfs_paths(start, [seed], Tuple{Int,Int}[], adj, atoms, fp, mode)
    end
    if verbose
        println("Fingerprint generation complete")
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
# Kekulé -> Aromatic normalization (chemical equivalence)
# -----------------------------

function _adj_bond_type(adj::Vector{Vector{Tuple{Int,UInt8}}}, a::Int, b::Int)::UInt8
    for (nbr, bt) in adj[a]
        if nbr == b
            return bt
        end
    end
    return 0x00
end

function _adj_set_bond!(adj::Vector{Vector{Tuple{Int,UInt8}}}, a::Int, b::Int, bt::UInt8)
    for k in eachindex(adj[a])
        if adj[a][k][1] == b
            adj[a][k] = (b, bt)
            break
        end
    end
    for k in eachindex(adj[b])
        if adj[b][k][1] == a
            adj[b][k] = (a, bt)
            break
        end
    end
    return nothing
end

function _canon6(cyc::NTuple{6,Int})::NTuple{6,Int}
    v = collect(cyc)

    function rotmin(w)
        mins = findall(==(minimum(w)), w)
        candidates = NTuple{6,Int}[]
        for idx in mins
            r = [w[idx:end]; w[1:idx-1]]
            push!(candidates, (r[1],r[2],r[3],r[4],r[5],r[6]))
        end
        return minimum(candidates)
    end

    a = rotmin(v)
    b = rotmin(reverse(v))
    return min(a, b)
end

function _find_6cycles(adj::Vector{Vector{Tuple{Int,UInt8}}}, n::Int; verbose::Bool=true)
    if verbose
        println("Finding 6-member cycles for aromatic normalization")
    end
    cycles = Set{NTuple{6,Int}}()

    function dfs(start::Int, curr::Int, path::Vector{Int})
        if length(path) == 6
            if any(t -> t[1] == start, adj[curr])
                cyc = (path[1], path[2], path[3], path[4], path[5], path[6])
                push!(cycles, _canon6(cyc))
            end
            return
        end

        for (nbr, _) in adj[curr]
            if nbr == start
                continue
            end
            if nbr in path
                continue
            end
            dfs(start, nbr, [path; nbr])
        end
    end

    for s in 1:n
        dfs(s, s, [s])
    end

    if verbose
        println("Found $(length(cycles)) unique 6-cycles")
    end
    return collect(cycles)
end

function _is_alternating_1_2(bonds::Vector{UInt8})::Bool
    @assert length(bonds) == 6
    ok12 = true
    ok21 = true
    for i in 1:6
        expected12 = isodd(i) ? BOND_SINGLE : BOND_DOUBLE
        expected21 = isodd(i) ? BOND_DOUBLE : BOND_SINGLE
        ok12 &= (bonds[i] == expected12)
        ok21 &= (bonds[i] == expected21)
    end
    return ok12 || ok21
end

"""
    normalize_kekule_aromatic!(atoms, adj)

Heuristic chemical equivalence:
- Detect 6-member rings with alternating single/double bonds
- Convert ring bonds to aromatic and mark ring atoms aromatic

This makes `C1=CC=CC=C1` behave like `c1ccccc1` for matching/search.
"""
function normalize_kekule_aromatic!(atoms::Vector{Atom}, adj::Vector{Vector{Tuple{Int,UInt8}}}; verbose::Bool=true)
    n = length(atoms)
    n < 6 && return nothing

    cycles = _find_6cycles(adj, n; verbose=verbose)

    normalized_count = 0
    for cyc in cycles
        nodes = collect(cyc)

        if any(i -> atoms[i].aromatic, nodes)
            continue
        end

        if any(i -> !(atoms[i].z in (UInt8(6), UInt8(7))), nodes)
            continue
        end

        bonds = UInt8[]
        ok = true
        for i in 1:6
            a = nodes[i]
            b = nodes[i == 6 ? 1 : i+1]
            bt = _adj_bond_type(adj, a, b)
            if bt != BOND_SINGLE && bt != BOND_DOUBLE
                ok = false
                break
            end
            push!(bonds, bt)
        end
        ok || continue

        _is_alternating_1_2(bonds) || continue

        for i in nodes
            atoms[i] = Atom(atoms[i].z, atoms[i].charge, true, atoms[i].h_count)
        end
        for i in 1:6
            a = nodes[i]
            b = nodes[i == 6 ? 1 : i+1]
            _adj_set_bond!(adj, a, b, BOND_AROMATIC)
        end
        normalized_count += 1
    end

    if verbose
        println("Normalized $normalized_count Kekulé rings to aromatic")
    end
    return nothing
end

# -----------------------------
# Compile molecule
# -----------------------------
function compile_mol(id::String, smiles::String; verbose::Bool=true)::Molecule
    if verbose
        println("Compiling molecule '$id' from SMILES: $smiles")
    end
    atoms, adj = parse_smiles(smiles; verbose=verbose)
    if verbose
        println("Parsed $(length(atoms)) atoms")
    end
    perceive!(atoms, adj; verbose=verbose)
    normalize_kekule_aromatic!(atoms, adj; verbose=verbose)
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
    if verbose
        println("Built CSR edge storage with $(length(edge_dst)) undirected edges")
    end

    # Rings
    adj_list = [ [d for (d, _) in adj[i]] for i in 1:n ]
    ringatom, edge_ring2 = detect_rings(n, adj_list, edge_indices; verbose=verbose)
    edge_ring = edge_ring2

    # Degree / valence
    degree = UInt8[length(adj[i]) for i in 1:n]
    valence = UInt8[
        sum(((b == BOND_AROMATIC) ? 1 : Int(b)) for (_, b) in adj[i]; init=0)
        for i in 1:n
    ]

    # Neighborhood hash (kept for future use; not used in matching right now)
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
    if verbose
        println("Computed degrees, valences, ring flags, and neighbor hashes")
    end

    # Fingerprints (exact + generalized)
    fp_exact = generate_fp(atoms, adj; verbose=verbose, mode=EXACT)
    fp_gen   = generate_fp(atoms, adj; verbose=verbose, mode=GENERALIZED)

    if verbose
        println("Molecule compilation complete for '$id'")
    end
    return Molecule(id, n, atoms, edge_src_offsets, edge_dst, edge_btype,
                    edge_ring, degree, valence, ringatom, neigh_hash,
                    fp_exact, fp_gen)
end

# -----------------------------
# Matcher (VF2-ish)
# -----------------------------

# Atom compatibility under a mode (pharma-friendly defaults)
function atom_ok(q::Atom, t::Atom, mode::MatchMode)::Bool
    if mode == EXACT
        return q.z == t.z
    end

    # GENERALIZED:
    # - Query aromatic carbon matches aromatic carbon or aromatic nitrogen (benzene ↔ pyridine)
    if q.aromatic && q.z == UInt8(6)
        return t.aromatic && (t.z == UInt8(6) || t.z == UInt8(7))
    end

    # Otherwise conservative: require element equality
    return q.z == t.z
end

function compatible(qatom::Atom, tatom::Atom,
                    qdeg::UInt8, tdeg::UInt8,
                    qval::UInt8, tval::UInt8,
                    qring::Bool, tring::Bool,
                    qhash::UInt32, thash::UInt32,
                    mode::MatchMode)::Bool

    return atom_ok(qatom, tatom, mode) &&
           (qatom.charge == tatom.charge) &&
           (qatom.aromatic == tatom.aromatic) &&
           (qdeg <= tdeg) &&
           (qval <= tval) &&
           (!qring || tring)
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

function substructure_match(query::Query, target::Molecule;
                            return_mapping::Bool=false,
                            verbose::Bool=true,
                            mode::MatchMode=EXACT)::Union{Bool, Vector{Int}}
    if verbose
        println("Performing substructure matching [mode=$(mode)] for query ($(query.natoms) atoms) against target '$(target.id)' ($(target.natoms) atoms)")
    end

    # Candidate prefilter per query atom
    candidates = [Int[] for _ in 1:query.natoms]
    for q in 1:query.natoms
        for t in 1:target.natoms
            if compatible(query.atoms[q], target.atoms[t],
                          query.degree[q], target.degree[t],
                          query.valence[q], target.valence[t],
                          query.ringatom[q], target.ringatom[t],
                          query.neigh_hash[q], target.neigh_hash[t],
                          mode)
                push!(candidates[q], t)
            end
        end
    end
    if verbose
        println("Prefiltered candidates per query atom: $([length(c) for c in candidates])")
    end

    # Order query atoms: smallest candidate set first, then higher degree, then ring atoms
    order = sort(1:query.natoms; by = q -> (length(candidates[q]), -Int(query.degree[q]), query.ringatom[q] ? 0 : 1))
    if verbose
        println("Query atom matching order: $order")
    end

    mapping = fill(0, query.natoms)
    used = falses(target.natoms)

    if vf2_backtrack!(query, target, order, 1, mapping, used, candidates)
        if verbose
            println("Substructure match found")
        end
        return return_mapping ? mapping : true
    end
    if verbose
        println("No substructure match found")
    end
    return false
end

# -----------------------------
# Index + Search
# -----------------------------
function build_index(smiles_list::Vector{String}, ids::Vector{String}; verbose::Bool=true, explain::Bool=false)::Index
    n = length(smiles_list)
    explain && println("ChemGraphSearch.build_index: compiling $n molecules")

    mols = Molecule[]
    failed = 0

    for (i, (sm, id)) in enumerate(zip(smiles_list, ids))
        if explain
            println("  [$i/$n] $id")
        elseif verbose
            println("[$i/$n] Compiling molecule '$id': $sm")
        end

        try
            mol = compile_mol(id, sm; verbose=verbose)
            push!(mols, mol)
        catch e
            failed += 1
            @warn "Failed to compile $id: $e"
        end
    end

    if explain
        println("  ok:     ", length(mols))
        println("  failed: ", failed)
    elseif verbose
        println("Index built with $(length(mols)) molecules")
    end

    return Index(mols)
end

function search(index::Index, query_smiles::String;
                return_mappings::Bool=false,
                verbose::Bool=true,
                mode::MatchMode=EXACT,
                explain::Bool=false)::Vector{SearchHit}

    # explain = "structured transparency" without noisy internals
    if explain
        println("ChemGraphSearch.search")
        println("  query:  ", query_smiles)
        println("  mode:   ", mode)
        println("  index:  ", length(index.mols), " molecules")
        println("  maps?:  ", return_mappings)
    elseif verbose
        println("Starting search for substructure from SMILES: $query_smiles [mode=$(mode)]")
        println("Index contains $(length(index.mols)) molecules")
    end

    query = compile_mol("query", query_smiles; verbose=verbose)

    # FP prefilter (mode-aware)
    qfp = (mode == EXACT) ? query.fp_exact : query.fp_gen

    fp_pass = Molecule[]
    for mol in index.mols
        tfp = (mode == EXACT) ? mol.fp_exact : mol.fp_gen
        if fp_subset(qfp, tfp)
            push!(fp_pass, mol)
        end
    end

    if explain
        println("  fp_pass: ", length(fp_pass), " candidates after fingerprint filter")
    elseif verbose
        println("Fingerprint prefilter reduced to $(length(fp_pass)) candidates")
    end

    # Verify with substructure matching
    hits = SearchHit[]
    checked = 0
    for mol in fp_pass
        checked += 1
        m = substructure_match(query, mol;
                               return_mapping=return_mappings,
                               verbose=verbose,
                               mode=mode)
        if m !== false
            push!(hits, SearchHit(mol.id, return_mappings ? m : nothing))
        end
    end

    if explain
        println("  verified: ", checked, " checked")
        println("  hits:     ", length(hits))
    elseif verbose
        println("Search complete: $(length(hits)) matches found")
    end

    return hits
end


# -----------------------------
# Persistence
# -----------------------------

function default_data_dir()::String
    dir = joinpath(homedir(), ".chemgraphsearch")
    isdir(dir) || mkpath(dir)
    return dir
end

function normalize_path(path::AbstractString)::String
    p = String(path)
    return isabspath(p) ? p : joinpath(default_data_dir(), p)
end

function save_index(index::Index, path::AbstractString; verbose::Bool=true)
    p = normalize_path(path)

    # Keep output short & useful
    if verbose
        println("ChemGraphSearch.save_index")
        println("  file:   ", p)
        println("  mols:   ", length(index.mols))
        println("  note:   index includes fp_exact + fp_gen (rebuild older saved indexes)")
    end

    open(p, "w") do io
        serialize(io, index)
    end

    verbose && println("  status: OK")
    return p
end

function load_index(path::AbstractString; verbose::Bool=true)::Index
    p = normalize_path(path)

    if verbose
        println("ChemGraphSearch.load_index")
        println("  file:   ", p)
    end

    idx = open(p) do io
        deserialize(io)
    end

    if verbose
        println("  mols:   ", length(idx.mols))
        println("  status: OK")
    end

    return idx
end
# -----------------------------
# Utility: read SMILES file
# -----------------------------
function read_smi_file(path::String; verbose::Bool=true)::Tuple{Vector{String}, Vector{String}}
    if verbose
        println("Reading SMILES file: $path")
    end
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
    if verbose
        println("Read $(length(ids)) molecules from file")
    end
    return ids, smiles
end

# -----------------------------
# User-friendly display (REPL / notebooks)
# -----------------------------

# Index: short summary
Base.show(io::IO, ::MIME"text/plain", idx::Index) =
    print(io, "ChemGraphSearch.Index(", length(idx.mols), " molecules)")

# Molecule: concise human-readable summary
Base.show(io::IO, ::MIME"text/plain", mol::Molecule) = begin
    arom = count(a -> a.aromatic, mol.atoms)
    print(io,
        "ChemGraphSearch.Molecule(\"",
        mol.id,
        "\", ",
        mol.natoms,
        " atoms, ",
        arom,
        " aromatic)"
    )
end

# Search result
Base.show(io::IO, ::MIME"text/plain", h::SearchHit) = begin
    if h.mapping === nothing
        print(io, "SearchHit(", h.id, ")")
    else
        print(io, "SearchHit(", h.id, ", mapping=", h.mapping, ")")
    end
end

end # module

