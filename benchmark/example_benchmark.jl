using GraphBasedSystems
using LinearAlgebra

include("adjacency_matrix.jl")


function initialize!_posdef!(system::System{N}) where N 
    initialize!(system,rand)
    for i=1:N
        system.matrix_entries[i,i].value += 100*I 
    end
    # display(isposdef(full_matrix(system)))
end


system = System{Float64}(A, ones(Int,size(A)[1])*3)
systemldlt = System{Float64}(A, ones(Int,size(A)[1])*3, symmetric=true)
systemllt = System{Float64}(A, ones(Int,size(A)[1])*3, symmetric=true)

SUITE["sparse_ldu"] = @benchmarkable ldu_solve!($system) setup=(initialize!($system))
SUITE["sparse_lu"] = @benchmarkable lu_solve!($system) setup=(initialize!($system))
SUITE["dense_lu"] = @benchmarkable lu(F)\f setup=(initialize!($system);F=full_matrix($system);f=full_vector($system))
SUITE["sparse_ldlt"] = @benchmarkable ldlt_solve!($systemldlt) setup=(initialize!($systemldlt))
SUITE["dense_ldlt"] = @benchmarkable bunchkaufman(F)\f setup=(initialize!($systemldlt);F=full_matrix($systemldlt);f=full_vector($systemldlt))
SUITE["sparse_llt"] = @benchmarkable llt_solve!($systemllt) setup=(initialize!_posdef!($systemllt))
SUITE["dense_llt"] = @benchmarkable cholesky(F)\f setup=(initialize!_posdef!($systemllt);F=full_matrix($systemllt);f=full_vector($systemllt))

SUITE["sparse_add"] = @benchmarkable +($system,$system) setup=(initialize!($system))
SUITE["dense_add"] = @benchmarkable (+(F,F);+(f+f)) setup=(initialize!($system);F=full_matrix($system);f=full_vector($system))
SUITE["sparse_mul"] = @benchmarkable *($system,$system) setup=(initialize!($system))
SUITE["dense_mul"] = @benchmarkable *(F,F) setup=(initialize!($system);F=full_matrix($system))
SUITE["sparse_solve"] = @benchmarkable \($system,Bmat) setup=(initialize!($system);Bmat=deepcopy(system.matrix_entries);initialize!($system))
SUITE["dense_solve"] = @benchmarkable \(F1,F2) setup=(initialize!($system);F2=full_matrix($system);initialize!($system);F1=full_matrix($system))
