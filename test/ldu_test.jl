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

D = [
    0 1 1 0 1 0
    1 0 0 1 0 1
    1 0 0 1 0 0
    0 1 1 0 0 0
    1 0 0 0 0 1
    0 1 0 0 1 0]

ZAA = zeros(Int64,10,10)
ZAB = zeros(Int64,10,9)
ZAC = zeros(Int64,10,6)
ZAD = zeros(Int64,10,6)
ZBB = zeros(Int64,9,9)
ZBC = zeros(Int64,9,6)
ZBD = zeros(Int64,9,6)
ZCC = zeros(Int64,6,6)
ZCD = zeros(Int64,6,6)

ZBA = ZAB'
ZCA = ZAC'
ZDA = ZAD'
ZCB = ZBC'
ZDB = ZBD'
ZDC = ZCD'


# Graph 1 is disconnected, 2-3-4-5 are connected, 6 is disconnected, 7 is disconnected

A = [
    A   ZAA ZAB ZAC ZAA ZAA ZAD
    ZAA A   ZAB ZAC ZAA ZAA ZAD
    ZBA ZBA B   ZBC ZBA ZBA ZBD
    ZCA ZCA ZCB C   ZCA ZCA ZCD
    ZAA ZAA ZAB ZAC A   ZAA ZAD
    ZAA ZAA ZAB ZAC ZAA A   ZAD
    ZDA ZDA ZDB ZDC ZDA ZDA D
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

