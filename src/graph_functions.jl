function dfs(graph,v)
    n = nv(graph)
    visited = zeros(Bool, n)
    list = Int64[]
    cycleclosures = Vector{Int64}[]
    _dfs(graph,v,0,list,cycleclosures,visited)

    cycleclosures = cycleclosures[sortperm(sort.(cycleclosures))[1:2:length(cycleclosures)]] # removes double entries of cycle connections and keeps the first found pair

    return list, cycleclosures
end

function _dfs(graph,v,p,list,cycleclosures,visited)
    if !visited[v]
        visited[v] = true
        for node in all_neighbors(graph,v)
            if node != p && visited[node]
                push!(cycleclosures,[v;node])
            end
            _dfs(graph,node,v,list,cycleclosures,visited)
        end
        push!(list,v)
    end

    return
end
