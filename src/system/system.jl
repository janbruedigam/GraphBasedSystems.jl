abstract type Symmetric end                     # symmetric = true
abstract type BlockSymmetric <: Symmetric end   # symmetric = true
abstract type SparsitySymmetric end             # symmetric = false, Any matrix is sparsity symmetric since zero entries can be assumed nonzero

struct System{N,T}
    matrix_entries::SparseMatrixCSC{GraphBasedSystems.Entry, Int64}
    vector_entries::Vector{Entry}
    diagonal_inverses::Vector{Entry}
    acyclic_children::Vector{Vector{Int64}} # Contains direct children that are not part of a cycle
    cyclic_children::Vector{Vector{Int64}}  # Contains direct and indirect children that are part of a cycle
    parents::Vector{Vector{Int64}}          # Contains direct and cycle-opening parents
    dfs_list::SVector{N,Int64}
    graph::SimpleGraph{Int64}
    dfs_graph::SimpleDiGraph{Int64}

    function System{T}(A, dims; force_static = false, symmetric = false) where T
        N = length(dims)
        if symmetric
            all(dims .== 1) ? S = Symmetric : S = BlockSymmetric
        else
            S = SparsitySymmetric
        end
        static = force_static || all(dims.<=10)

        full_graph = Graph(A)

        matrix_entries = spzeros(Entry,N,N)

        for (i,dimi) in enumerate(dims)
                matrix_entries[i,i] = Entry{T}(dimi, dimi, static = static)
        end

        vector_entries = [Entry{T}(dim, static = static) for dim in dims]
        diagonal_inverses = [Entry{T}(dim, dim, static = static) for dim in dims]

        graphs, roots = split_adjacency(A)
        dfs_list = Int64[]
        acyclic_children = [Int64[] for i=1:N]
        cycles = [Vector{Int64}[] for i=1:N]
        parents = [Int64[] for i=1:N]
        edgelist = Graphs.SimpleEdge{Int64}[]

        for (i,graph) in enumerate(graphs)
            root = roots[i]
            dfs_graph = dfs_tree(graph, root)
            sub_dfs_list, cycle_closures = dfs(graph, root)
            append!(dfs_list, sub_dfs_list)

            cycle_dfs_graph = copy(dfs_graph)
            cycle_dfs_graph_reverse = copy(dfs_graph)
            for cycle_closure in cycle_closures
                add_edge!(cycle_dfs_graph, cycle_closure...)
                add_edge!(cycle_dfs_graph_reverse, reverse(cycle_closure)...)
            end
            append!(edgelist, collect(edges(cycle_dfs_graph_reverse)))

            for v in sub_dfs_list
                acyclic_children[v] = neighbors(dfs_graph, v)
                parents[v] = neighbors(reverse(dfs_graph), v)
            end   

            for cyclic_members in simplecycles(cycle_dfs_graph)
                v, cyclic_children = cycle_parent_children(cyclic_members, parents)
                cyclic_children = reverse(cyclic_children)
                lump_cycles!(cycles[v], cyclic_children)

                acyclic_children[v] = setdiff(acyclic_children[v], cyclic_children)
                for c in cyclic_children
                    matrix_entries[v,c] = Entry{T}(dims[v], dims[c], static = static)
                    !symmetric && (matrix_entries[c,v] = Entry{T}(dims[c], dims[v], static = static))

                    v ∉ parents[c] && push!(parents[c],v)
                end
            end  
            for v in sub_dfs_list
                for c in acyclic_children[v]
                    matrix_entries[v,c] = Entry{T}(dims[v], dims[c], static = static)
                    !symmetric && (matrix_entries[c,v] = Entry{T}(dims[c], dims[v], static = static))
                end
            end
        end

        full_dfs_graph = SimpleDiGraph(edgelist)
        cyclic_children = [unique(vcat(cycles[i]...)) for i=1:N]
        cyclic_children = [intersect(dfs_list,cyclic_children[i]) for i=1:N]

        new{N,S}(matrix_entries, vector_entries, diagonal_inverses, acyclic_children, cyclic_children, parents, dfs_list, full_graph, full_dfs_graph)
    end

    System(A, dims; force_static = false, symmetric = false) = System{Float64}(A, dims; force_static, symmetric)
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, system::System{N,S}) where {N,S}
    if S <: BlockSymmetric
        print(io, "BlockSymmetric ")
    elseif S <: Symmetric
        print(io, "Symmetric ")
    elseif S <: SparsitySymmetric
        print(io, "SparsitySymmetric ")
    end
    println(io, "System with "*string(N)*" nodes.")
    SparseArrays._show_with_braille_patterns(io, system.matrix_entries)
end


@inline children(system, v) = outneighbors(system.dfs_graph, v)
@inline connections(system, v) = neighbors(system.graph, v)
@inline parents(system, v) = inneighbors(system.dfs_graph, v)

# There probably exists a smarter way of getting the dense matrix from the spares one 
function full_matrix(system::System{N}) where N
    dims = [length(system.vector_entries[i].value) for i=1:N]

    range = [1:dims[1]]

    for (i,dim) in enumerate(collect(Iterators.rest(dims, 2)))
        push!(range,sum(dims[1:i])+1:sum(dims[1:i])+dim)
    end
    A = zeros(sum(dims),sum(dims))

    for (i,row) in enumerate(system.matrix_entries.rowval)
        col = findfirst(x->i<x,system.matrix_entries.colptr)-1
            A[range[row],range[col]] = system.matrix_entries[row,col].value
    end
    return A
end

# function full_matrix(system::System{N,<:Symmetric}) where N
#     dims = [length(system.vector_entries[i].value) for i=1:N]

#     range = [1:dims[1]]

#     for (i,dim) in enumerate(collect(Iterators.rest(dims, 2)))
#         push!(range,sum(dims[1:i])+1:sum(dims[1:i])+dim)
#     end
#     A = zeros(sum(dims),sum(dims))

#     for (i,row) in enumerate(system.matrix_entries.rowval)
#         col = findfirst(x->i<x,system.matrix_entries.colptr)-1
#         A[range[row],range[col]] = system.matrix_entries[row,col].value
#         if col != row
#             A[range[col],range[row]] = system.matrix_entries[row,col].value'
#         end
#     end
#     return A
# end

full_vector(system) = vcat(getfield.(system.vector_entries,:value)...)

reordered_matrix(system) = full_matrix(system)[system.dfs_list,system.dfs_list]
reordered_vector(system) = full_vector(system)[system.dfs_list]