module GraphBasedSystems

using LinearAlgebra
using SparseArrays
using SparseArrays: widelength
using SparseArrays.HigherOrderFns: _sumnnzs, _allocres, _map_zeropres!, _densestructure!
using StaticArrays
using Graphs


export System,
    Entry,

    full_matrix,
    full_vector,
    initialize!,
    reset_inverse_diagonals!,
    
    children,
    connections,
    parents,

    ldu_solve!,
    ldu_matrix_solve!,
    ldu_factorization!,
    ldu_backsubstitution!,
    lu_solve!,
    lu_matrix_solve!,
    lu_factorization!,
    lu_backsubstitution!,
    ldlt_solve!,
    ldlt_matrix_solve!,
    ldlt_factorization!,
    ldlt_backsubstitution!,
    llt_solve!,
    llt_matrix_solve!,
    llt_factorization!,
    llt_backsubstitution!


include(joinpath("util", "custom_static.jl"))

include(joinpath("system", "entry.jl"))
include(joinpath("system", "system.jl"))
include(joinpath("system", "graph_functions.jl"))

include(joinpath("system", "interface.jl"))
include(joinpath("system", "dense.jl"))

include(joinpath("solvers", "matrix.jl"))
include(joinpath("solvers", "lu.jl"))
include(joinpath("solvers", "llt.jl"))
include(joinpath("solvers", "ldlt.jl"))
include(joinpath("solvers", "ldu.jl"))

end
