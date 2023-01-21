using GraphBasedSystems
using LinearAlgebra

include("adjacency_matrix.jl")


system = System{Float64}(A, rand(0:3, size(A)[1]))

for i=1:10
    initialize!(system)

    F = full_matrix(system)
    f = full_vector(system)
    lu_solve!(system)
    @test maximum(abs.(full_vector(system)-F\f)) < 1e-3
end

display(system)