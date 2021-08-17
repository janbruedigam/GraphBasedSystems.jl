using GraphBasedSystems
using LinearAlgebra
GBS = GraphBasedSystems

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

B = [
    0 1 0 0 0 1 0 0 0
    1 0 0 0 0 1 0 0 0
    0 0 0 0 0 1 1 0 1
    0 0 0 0 0 0 1 1 0
    0 0 0 0 0 0 0 1 1
    1 1 1 0 0 0 1 0 0
    0 0 1 1 0 1 0 1 1
    0 0 0 1 1 0 1 0 1
    0 0 1 0 1 0 1 1 0
]

C = [
    0 1 1 0 0 0
    1 0 0 1 0 0
    1 0 0 1 1 0
    0 1 1 0 0 1
    0 0 1 0 0 1
    0 0 0 1 1 0
]

ZAA = zeros(Int64,10,10)
ZAB = zeros(Int64,10,9)
ZBA = ZAB'
ZAC = zeros(Int64,10,6)
ZCA = ZAC'
ZBC = zeros(Int64,9,6)
ZCB = ZBC'

# Graph 1 is disconnected, 2-3-4-5 are connected, 6 is disconnected

A = [
    A   ZAA ZAB ZAC ZAA ZAA
    ZAA A   ZAB ZAC ZAA ZAA
    ZBA ZBA B   ZBC ZBA ZBA
    ZCA ZCA ZCB C   ZCA ZCA
    ZAA ZAA ZAB ZAC A   ZAA
    ZAA ZAA ZAB ZAC ZAA A
]

A[13,21] = A[21,13] = 1
A[16,30] = A[30,16] = 1
A[20,36] = A[36,20] = 1


system = System{Float64}(A, rand(0:3,size(A)[1]))

for i=1:10
    for entry in system.matrix_entries.nzval
        GBS.randomize!(entry)
    end
    for entry in system.vector_entries
        GBS.randomize!(entry)
    end

    F = full_matrix(system)
    f = full_vector(system)
    ldu_solve!(system)
    @test maximum(abs.(full_vector(system)-F\f)) < 1e-3
end

