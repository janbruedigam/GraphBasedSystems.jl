using GraphBasedSystems
using GraphBasedSystems: srand, ranges
using LinearAlgebra
using OrderedCollections

include("adjacency_matrix.jl")


for i=1:5
    system = System{Float64}(A, rand(0:3, N), symmetric=true)
    initialize!(system)
    F = full_matrix(system)
    while !isposdef(F)
        for i=1:N
            system.matrix_entries[i,i].value += I 
        end
        F = full_matrix(system)
    end

    F = full_matrix(system)
    f = full_vector(system)
    llt_solve!(system)
    @test maximum(abs.(full_vector(system)-F\f)) < 1e-3
end

for i=1:5
    system = System{Float64}(A, rand(0:3, N), symmetric=true)
    initialize!(system)
    F = full_matrix(system)
    while !isposdef(F)
        for i=1:N
            system.matrix_entries[i,i].value += I 
        end
        F = full_matrix(system)
    end
    system.actives = srand(Bool,N)

    range = vcat(sort!(OrderedDict(ranges(system; actives=system.actives))).vals...)

    F = full_matrix(system)[range, range]
    f = full_vector(system)[range]
    llt_solve!(system)
    @test maximum(abs.(full_vector(system)[range]-F\f)) < 1e-3
end

for i=1:5
    system = System{Float64}(A, rand(0:3, N), symmetric=true)
    initialize!(system)
    Bmat = deepcopy(system.matrix_entries)
    F2 = full_matrix(Bmat,false,system.dims,system.dims)
    initialize!(system)
    F1 = full_matrix(system)
    while !isposdef(F1)
        for i=1:N
            system.matrix_entries[i,i].value += I 
        end
        F1 = full_matrix(system)
    end
    F1 = full_matrix(system)

    Cmat = llt_matrix_solve!(system, Bmat)
    @test maximum(abs.(full_matrix(Cmat,false,system.dims,system.dims)-F1\F2)) < 1e-3
end

system = System{Float64}(A, rand(0:3, N), symmetric=true)

display(system)