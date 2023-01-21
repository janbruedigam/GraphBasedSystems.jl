@inline children(system::System, v) = outneighbors(system.dfs_graph, v) # all direct children of v
@inline connections(system::System, v) = neighbors(system.graph, v)     # all connected nodes of v
@inline parents(system::System, v) = inneighbors(system.dfs_graph, v)   # same elements as system.parents[v], but potentially different order

function ranges(system::System{N}) where N
    dims = system.dims
    range_dict = Dict(1=>1:dims[1])
    for i=2:N
        range_dict[i] = last(range_dict[i-1])+1:sum(dims[1:i])
    end

    return range_dict
end
function ranges(dims::AbstractVector)
    range_dict = Dict(1=>1:dims[1])
    for i=2:size(dims)[1]
        range_dict[i] = last(range_dict[i-1])+1:sum(dims[1:i])
    end

    return range_dict
end

function initialize!(system::System, init_function = randn)
    for entry in system.matrix_entries.nzval
        initialize!(entry, init_function)
    end
    for entry in system.vector_entries
        initialize!(entry, init_function)
    end
end

function initialize!(system::System{N,<:Symmetric}, init_function = randn) where N
    matrix_entries = system.matrix_entries

    for entry in matrix_entries.nzval
        initialize!(entry, init_function)
    end
    for i=1:N
        matrix_entries[i,i].value += matrix_entries[i,i].value'
    end

    for entry in system.vector_entries
        initialize!(entry, init_function)
    end
end

function reset_inverse_diagonals!(system::System)
    for entry in system.diagonal_inverses
        entry.isinverted = false
    end
end

adjacency(system::System) = adjacency(system.matrix_entries)
function adjacency(M::SparseMatrixCSC{Entry, Int64})
    N = size(M)[1]
    adj = zeros(Int64,N,N)

    for i=1:N
        for j=1:N
            i == j && continue
            Mij = M[i,j]
            Mij isa Entry{Nothing} ? adj[i,j] = 0 : adj[i,j] = 1
        end
    end

    return adj
end