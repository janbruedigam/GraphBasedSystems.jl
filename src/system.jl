struct System{N}
    matrix_entries::Dict{Tuple{Int64, Int64}, Entry}
    vector_entries::Dict{Int64, Entry}
    acyclic_children::Dict{Int64,Vector{Int64}} # Contains direct children that are not part of a cycle
    cycles::Dict{Int64,Vector{Vector{Int64}}}   # Contains direct and indirect children that are part of a cycle (sorted in cycles)
    parents::Dict{Int64,Vector{Int64}}          # Contains direct and cycle-opening parent (2 elemtents at most)
    dfs_list::SVector{N,Int64}
    reverse_dfs_list::SVector{N,Int64}

    function System{T}(A, dims; ids = collect(1:length(dims)), force_static = false) where T
        N = length(dims)
        static = force_static || all(dims.<=10)

        full_graph = reindex_graph(Graph(A), ids)

        matrix_entries = Dict{Tuple{Int64, Int64}, Entry}()

        for (_i,dimi) in enumerate(dims)
            i = ids[_i]
            for (_j,dimj) in enumerate(dims)
                j = ids[_j]
                if i == j
                    matrix_entries[(i,j)] = Entry{T}(dimi, dimj, static = static)
                elseif j ∈ all_neighbors(full_graph,i)
                    matrix_entries[(i,j)] = Entry{T}(dimi, dimj, static = static)
                    matrix_entries[(j,i)] = Entry{T}(dimj, dimi, static = static)
                end
            end
        end

        vector_entries = Dict{Int64, Entry}()

        for (i,dim) in enumerate(dims)
            vector_entries[ids[i]] = Entry{T}(dim, static = static)
        end

        graphs, roots = split_adjacency(A)
        dfs_list = Int64[]
        acyclic_children = Dict{Int64,Vector{Int64}}()
        cycles = Dict([ids[i] => Vector{Int64}[] for i=1:N]...)
        parents = Dict{Int64,Vector{Int64}}()
        dims = Dict([ids[i] => dims[i] for i=1:N]...)

        for (i,graph) in enumerate(graphs)
            graph = reindex_graph(graph, ids)
            root = ids[roots[i]]
            dfs_graph = dfs_tree(graph, root)
            sub_dfs_list, cycle_closures = dfs(graph, root)
            append!(dfs_list, sub_dfs_list)

            cycle_dfs_graph = copy(dfs_graph)
            for cycle_closure in cycle_closures
                add_edge!(cycle_dfs_graph, cycle_closure...)
            end

            for v in sub_dfs_list
                acyclic_children[v] = neighbors(dfs_graph, v)
                parents[v] = neighbors(reverse(dfs_graph), v)
            end   

            for cyclic_members in simplecycles(cycle_dfs_graph)
                v, cyclic_children = cycle_parent_children(cyclic_members, parents)
                cyclic_children = reverse(cyclic_children)
                push!(cycles[v], cyclic_children)

                acyclic_children[v] = setdiff(acyclic_children[v], cyclic_children)
                for c in cyclic_children
                    matrix_entries[(v,c)] = Entry{T}(dims[v], dims[c], static = static)
                    matrix_entries[(c,v)] = Entry{T}(dims[c], dims[v], static = static)

                    v != parents[c][1] && push!(parents[c],v)
                end
            end            
        end

        reverse_dfs_list = reverse(dfs_list)

        new{N}(matrix_entries, vector_entries, acyclic_children, cycles, parents, dfs_list, reverse_dfs_list)
    end

    System(A, dims; force_static = false) = System{Float64}(A, dims; force_static = force_static)
end


function split_adjacency(A)
    graphs = SimpleGraph{Int64}[]

    subinds = connected_components(Graph(A))
    for subset in subinds
        subA = zeros(Int64,size(A)...)
        subA[subset,subset] = A[subset,subset]
        push!(graphs, Graph(subA))
    end

    roots = [subinds[i][1] for i=1:length(subinds)]

    return graphs, roots
end

function reindex_graph(g, ids)
    originaledges = collect(edges(g))
    newedges = LightGraphs.SimpleEdge{Int64}[]
    for edge in originaledges
        push!(newedges, LightGraphs.SimpleEdge(ids[edge.src], ids[edge.dst]))
    end

    return SimpleGraphFromIterator(newedges)
end

function cycle_parent_children(cyclic_members, parents)
    parent = -1
    for member in cyclic_members
        !isempty(parents[member]) && parents[member][1] ∈ cyclic_members && continue
        parent = member
        break
    end

    return parent, setdiff(cyclic_members, parent)
end


function full_matrix(system::System{N}) where N
    ids = convert(Vector,sort(system.dfs_list)) # SA for some reason gives wrong 'rest' in Iterators
    dims = Dict([id => length(system.vector_entries[id].value) for id in ids]...)

    range = Dict(ids[1] => 1:dims[ids[1]])
    
    for (i,id) in enumerate(collect(Iterators.rest(ids, 2)))
        lastind = last(range[ids[i]])
        range[id] = lastind+1:lastind+dims[id]
    end
    A = zeros(sum(values(dims)),sum(values(dims)))
    for key in keys(system.matrix_entries)
        A[range[key[1]],range[key[2]]] = system.matrix_entries[key].value
    end
    return A
end

function full_vector(system::System{N}) where N
    ids = convert(Vector,sort(system.dfs_list)) # SA for some reason gives wrong 'rest' in Iterators
    dims = Dict([id => length(system.vector_entries[id].value) for id in ids]...)

    range = Dict(ids[1] => 1:dims[ids[1]])
    
    for (i,id) in enumerate(collect(Iterators.rest(ids, 2)))
        lastind = last(range[ids[i]])
        range[id] = lastind+1:lastind+dims[id]
    end
    b = zeros(sum(values(dims)))
    for key in keys(system.vector_entries)
        b[range[key[1]]] = system.vector_entries[key].value
    end
    return b
end
