using GraphBasedSystems
using LinearAlgebra

include("adjacency_matrix.jl")


system = System{Float64}(A, rand(0:3,size(A)[1]), symmetric=true)

for i=1:10
    randomize!(system)
    F = full_matrix(system)
    while !isposdef(F)
        for i=1:size(A)[1]
            system.matrix_entries[i,i].value += I 
        end
        F = full_matrix(system)
    end

    F = full_matrix(system)
    f = full_vector(system)
    llt_solve!(system)
    @test maximum(abs.(full_vector(system)-F\f)) < 1e-3
end

