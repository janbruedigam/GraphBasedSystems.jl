mutable struct Entry{ET}
    value::ET

    function Entry{T}(dims...; static = true) where T
        static ? value = srand(T,dims...) : value = rand(T,dims...)
        new{typeof(value)}(value)
    end

    Entry(dims...; static = true) = Entry{Float64}(dims...; static = static)
end

struct System{N}
    matrixentries::Dict{Tuple{Int64, Int64}, Entry}
    vectorentries::Vector{Entry}
    graph::SimpleGraph{Int64}
    dfs_graph::SimpleDiGraph{Int64}
    reverse_dfs_graph::SimpleDiGraph{Int64}
    dfs_list::SVector{N,Int64}
    reverse_dfs_list::SVector{N,Int64}

    function System{T}(A, dims; force_static = false) where T
        N = length(dims)
        static = force_static || all(dims.<=10)

        graph = Graph(A)
        @assert nv(graph) == N

        matrixentries = Dict{Tuple{Int64, Int64}, Entry}()

        for (i,dimi) in enumerate(dims)
            for (j,dimj) in enumerate(dims)
                if i==j
                    matrixentries[(i,j)] = Entry{T}(dimi, dimj, static = static)
                elseif j ∈ all_neighbors(graph,i)
                    matrixentries[(i,j)] = Entry{T}(dimi, dimj, static = static)
                    matrixentries[(j,i)] = Entry{T}(dimj, dimi, static = static)
                end
            end
        end

        vectorentries = [Entry{T}(dim, static = static) for dim in dims]

        dfs_graph = dfs_tree(graph,1)
        reverse_dfs_graph = reverse(dfs_graph)
        dfs_list = dfs(graph,1)
        reverse_dfs_list = reverse(dfs_list)

        new{N}(matrixentries,vectorentries,graph,dfs_graph,reverse_dfs_graph,dfs_list,reverse_dfs_list)
    end

    System(A, dims; force_static = false) = System{Float64}(A, dims; force_static = force_static)
end


function full_matrix(system)
    n = nv(system.graph)
    dims = [length(system.vectorentries[i].value) for i=1:n]

    range = [1:dims[1]]

    for (i,dim) in enumerate(collect(Iterators.rest(dims, 2)))
        push!(range,sum(dims[1:i])+1:sum(dims[1:i])+dim)
    end
    A = zeros(sum(dims),sum(dims))
    for key in keys(system.matrixentries)
        A[range[key[1]],range[key[2]]] = system.matrixentries[key].value
    end
    return A
end

full_vector(system) = vcat(getfield.(system.vectorentries,:value)...)