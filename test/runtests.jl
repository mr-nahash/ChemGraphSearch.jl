using Test
using ChemGraphSearch

@testset "SMILES Parser" begin
    atoms, adj = parse_smiles("C1CCCCC1")
    @test length(atoms) == 6

    atoms2, adj2 = parse_smiles("c1ccccc1")
    @test length(atoms2) == 6
    @test all(a -> a.aromatic, atoms2)

    atoms3, adj3 = parse_smiles("[nH]1cccc1")
    @test length(atoms3) == 5
    @test atoms3[1].aromatic == true
end

@testset "Aromatic matching" begin
    idx = build_index(["c1ccccc1", "CCO", "O=c1ccccc1"], ["benz", "eth", "benzald"])
    res = search(idx, "c1ccccc1")
    @test length(res) >= 1
end

@testset "Aromatic/kekule equivalence" begin
    idxA = build_index(["c1ccccc1", "C1=CC=CC=C1"], ["arom", "kekule"])
    res = search(idxA, "c1ccccc1")
    ids = sort(first.(res))
    @test ids == ["arom", "kekule"]
end
