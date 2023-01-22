using GraphBasedSystems
using LinearAlgebra

include("adjacency_matrix.jl")


for i=1:5
    system = System{Float64}(A, rand(0:3, size(A)[1]))
    initialize!(system)

    F = full_matrix(system)
    f = full_vector(system)
    lu_solve!(system)
    @test maximum(abs.(full_vector(system)-F\f)) < 1e-3
end

for i=1:5
    system = System{Float64}(A, rand(0:3, size(A)[1]))
    initialize!(system)
    Bmat = deepcopy(system.matrix_entries)
    F2 = full_matrix(Bmat,false,system.dims,system.dims)
    initialize!(system)
    F1 = full_matrix(system)

    Cmat = lu_matrix_solve!(system, Bmat)
    @test maximum(abs.(full_matrix(Cmat,false,system.dims,system.dims)-F1\F2)) < 1e-3
end

for i=1:5
    system = System{Float64}(A, rand(0:3, size(A)[1]))
    initialize!(system)
    Bmat = deepcopy(system.matrix_entries)
    F2 = full_matrix(Bmat,false,system.dims,system.dims)
    initialize!(system)
    F1 = full_matrix(system)

    Cmat = system \ Bmat
    @test maximum(abs.(full_matrix(Cmat,false,system.dims,system.dims)-F1\F2)) < 1e-3
end

system = System{Float64}(A, rand(0:3, size(A)[1]))

display(system)