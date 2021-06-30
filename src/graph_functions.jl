function dfs(graph,v)
    n = nv(graph)
    visited = zeros(Bool, n)
    list = Int64[]
    dfs2(graph,v,list,visited)

    return list
end

function dfs2(graph,v,list,visited)
    if !visited[v]
        visited[v] = true
        for node in all_neighbors(graph,v)
            dfs2(graph,node,list,visited)
        end
        push!(list,v)
    end

    return
end

children(system, v) = neighbors(system.dfs_graph, v)
parents(system, v) = neighbors(system.reverse_dfs_graph, v)