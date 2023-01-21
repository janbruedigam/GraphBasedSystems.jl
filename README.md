# GraphBasedSystems
[![Build Status](https://github.com/janbruedigam/GraphBasedSystems.jl/workflows/CI/badge.svg)](https://github.com/janbruedigam/GraphBasedSystems.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/janbruedigam/GraphBasedSystems.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/janbruedigam/GraphBasedSystems.jl)

A package for efficiently solving linear systems of equations based on graph properties. Any matrix representing an underlying real-world system can be represented by a graph, for example a mechanical system, chemical molecules, neural networks, ... The code of the package is currently used for the robotics simulator [Dojo.jl](https://github.com/dojo-sim/Dojo.jl)

By providing the adjacency matrix of a graph-based system, the `GraphBasedSystems` package can automatically exploit the existing sparsity when solving the linear system to speed up calculations. Currently, the LDU, LU, LDLt, and LLt decomposition/backsubstitution are implemented. LLt is currently slow due to memory allocations.

```julia
using GraphBasedSystems

graph_matrix = [  # The adjacency matrix for the underlying graph
  0 1 1 0
  1 0 1 0
  1 1 0 1
  0 0 1 0
]

dimensions = [2; 3; 0; 1] # The dimension of row/column

system = System{Float64}(graph_matrix, dimensions; symmetric=false) # The resulting linear system. Set symmetric=true for symmetric systems

initialize!(system, randn) # initialize all system entries randomly
system.matrix_entries[1,2].value = rand(2,3) # Directly set the value of a matrix entry
system.vector_entries[4].value = rand(1) # Directly set the value of a vector entry

A = full_matrix(system) # Inspect the matrix
b = full_vector(system) # Inspect the vector

ldu_solve!(system) # Solve the system inplace

x = full_vector(system) # Inspect the result
x - A\b # Compare to classical implementation

```
