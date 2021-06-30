module GraphBasedSystems

using StaticArrays
using LightGraphs


export System
    ldu_factorization!,
    ldu_backsubstitution!,
    ldu_solve!

include("custom_static.jl")

include("system.jl")
include("graph_functions.jl")

include("ldu.jl")

end
