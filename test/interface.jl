using GraphBasedSystems
using LinearAlgebra

include("adjacency_matrix.jl")

system = System{Float64}(A, rand(0:3, size(A)[1]))

children(system, 1)
connections(system, 1)
parents(system, 1)
children(system, 10)
connections(system, 10)
parents(system, 10)
@test true

@test GraphBasedSystems.ranges(system) == GraphBasedSystems.ranges(system.dims)

initialize!(system, rand)
initialize!(system, randn)
initialize!(system, ones)
initialize!(system, zeros)

@test true

reset_inverse_diagonals!(system)

@test true