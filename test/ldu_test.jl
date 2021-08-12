using GraphBasedSystems
using LinearAlgebra
GBS = GraphBasedSystems

Z = zeros(Int64,10,10)
A = [
    0 1 0 1 1 0 1 0 1 0
    1 0 1 0 0 0 0 0 0 0
    0 1 0 0 0 0 0 0 0 0
    1 0 0 0 1 0 0 0 0 0
    1 0 0 1 0 1 0 0 0 0
    0 0 0 0 1 0 0 0 0 0
    1 0 0 0 0 0 0 1 0 0
    0 0 0 0 0 0 1 0 1 1
    1 0 0 0 0 0 0 1 0 0
    0 0 0 0 0 0 0 1 0 0
]

# Graph 1 is disconnected, 2-3-4-5 are connected, 6 is disconnected

A = [
    A Z Z Z Z Z
    Z A Z Z Z Z
    Z Z A Z Z Z
    Z Z Z A Z Z
    Z Z Z Z A Z
    Z Z Z Z Z A
]

A[13,21] = A[21,13] = 1
A[16,31] = A[31,16] = 1
A[20,41] = A[41,20] = 1


system = System{Float64}(A, [rand(1:3) for i=1:size(A)[1]])

for i=1:10
    for entry in collect(values(system.matrix_entries))
        GBS.randomize!(entry)
    end
    for entry in collect(values(system.vector_entries))
        GBS.randomize!(entry)
    end

    B = full_matrix(system)
    b = full_vector(system)
    ldu_solve!(system)
    @test maximum(abs.(full_vector(system)-B\b)) < 1e-3
end

