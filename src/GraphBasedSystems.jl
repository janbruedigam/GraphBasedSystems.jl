module GraphBasedSystems

using StaticArrays
using LightGraphs


export System,
    full_matrix,
    full_vector,
    children,
    connections,
    parents,

    ldu_solve!


include("custom_static.jl")

include("entry.jl")
include("system.jl")
include("setup_functions.jl")

# include("lu.jl")
# include("llt.jl")
# include("ldlt.jl")
include("ldu.jl")

end
