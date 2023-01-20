using GraphBasedSystems
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


system = System{Float64}(A, ones(Int,size(A)[1])*6)

function randomize_posdef!(system)
    randomize!(system)
    F = full_matrix(system)
    while !isposdef(F)
        for i=1:size(A)[1]
            system.matrix_entries[i,i].value += I 
        end
        F = full_matrix(system)
    end
end

SUITE["ldu"] = @benchmarkable ldu_solve!($system) samples=2 setup=(randomize!($system))
SUITE["lu"] = @benchmarkable lu_solve!($system) samples=2 setup=(randomize!($system))
SUITE["ldlt"] = @benchmarkable ldlt_solve!($system) samples=2 setup=(randomize!($system))
SUITE["llt"] = @benchmarkable llt_solve!($system) samples=2 setup=(randomize_posdef!($system))