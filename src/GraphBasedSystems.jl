module GraphBasedSystems

using SparseArrays
using StaticArrays
using Graphs


export System,
    full_matrix,
    full_vector,
    children,
    connections,
    parents,

    ldu_solve!


include(joinpath("util", "custom_static.jl"))

include(joinpath("system", "entry.jl"))
include(joinpath("system", "system.jl"))
include(joinpath("system", "setup_functions.jl"))

include(joinpath("solvers", "lu.jl"))
include(joinpath("solvers", "llt.jl"))
include(joinpath("solvers", "ldlt.jl"))
include(joinpath("solvers", "ldu.jl"))

end
