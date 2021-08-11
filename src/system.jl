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
    acyclic_children::Vector{Vector{Int64}}         # Contains direct children that are not part of a cycle
    cycles::Vector{Vector{Vector{Int64}}}   # Contains direct and indirect children that are part of a cycle (sorted in cycles)
    parents::Vector{Vector{Int64}}          # Contains direct and cycle-opening parent (2 elemtents at most)
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
                if i == j
                    matrixentries[(i,j)] = Entry{T}(dimi, dimj, static = static)
                elseif j ∈ all_neighbors(graph,i)
                    matrixentries[(i,j)] = Entry{T}(dimi, dimj, static = static)
                    matrixentries[(j,i)] = Entry{T}(dimj, dimi, static = static)
                end
            end
        end

        vectorentries = [Entry{T}(dim, static = static) for dim in dims]

        dfs_graph = dfs_tree(graph,1)
        dfs_list, cycle_closures = dfs(graph,1)
        reverse_dfs_list = reverse(dfs_list)
        cycle_dfs_graph = copy(dfs_graph)
        for cycle_closure in cycle_closures
            add_edge!(cycle_dfs_graph, cycle_closure...)
        end

        acyclic_children = [Int64[] for i=1:N]
        parents = [Int64[] for i=1:N]
        cycles = [Vector{Int64}[] for i=1:N]

        for v = 1:N
            acyclic_children[v] = neighbors(dfs_graph, v)
            parents[v] = neighbors(reverse(dfs_graph), v)
        end     
        for cyclemembers in simplecycles(cycle_dfs_graph)
            v = cyclemembers[1]
            cyclechildren = reverse(cyclemembers[2:end])
            push!(cycles[v], cyclechildren)

            acyclic_children[v] = setdiff(acyclic_children[v],cyclechildren)
            for c in cyclechildren
                matrixentries[(v,c)] = Entry{T}(dims[v], dims[c], static = static)
                matrixentries[(c,v)] = Entry{T}(dims[c], dims[v], static = static)

                v != parents[c][1] && push!(parents[c],v)
            end
        end

        new{N}(matrixentries,vectorentries,acyclic_children,cycles,parents,dfs_list,reverse_dfs_list)
    end

    System(A, dims; force_static = false) = System{Float64}(A, dims; force_static = force_static)
end


function full_matrix(system::System{N}) where N
    dims = [length(system.vectorentries[i].value) for i=1:N]

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