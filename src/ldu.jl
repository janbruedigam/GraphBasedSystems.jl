function ldu_factorization_l!(offdiagonal, diagonal)
    offdiagonal.value = offdiagonal.value/diagonal.value
    return
end
function ldu_factorization_u!(offdiagonal, diagonal)
    offdiagonal.value = diagonal.value\offdiagonal.value
    return
end
function ldu_factorization_d!(offdiagonal_l, offdiagonal_u, diagonal_c, diagonal_v)
    diagonal_v.value -= offdiagonal_l.value*diagonal_c.value*offdiagonal_u.value
    return
end

function ldu_factorization!(system)
    matrixentries = system.matrixentries

    for v in system.dfs_list
        for c in children(system,v) 
            ldu_factorization_l!(matrixentries[(v,c)], matrixentries[(c,c)])
            ldu_factorization_u!(matrixentries[(c,v)], matrixentries[(c,c)])
            ldu_factorization_d!(matrixentries[(v,c)], matrixentries[(c,v)], matrixentries[(c,c)], matrixentries[(v,v)])
        end
    end
    return
end

function ldu_backsubstitution_l!(vector_v, vector_c, offdiagonal)
    vector_v.value -= offdiagonal.value*vector_c.value
    return
end
function ldu_backsubstitution_u!(vector_v, vector_p, offdiagonal)
    vector_v.value -= offdiagonal.value*vector_p.value
    return
end
function ldu_backsubstitution_d!(vector, diagonal)
    vector.value = diagonal.value\vector.value
    return
end

function ldu_backsubstitution!(system)
    matrixentries = system.matrixentries
    vectorentries = system.vectorentries

    for v in system.dfs_list
        for c in children(system,v) # children
            ldu_backsubstitution_l!(vectorentries[v], vectorentries[c], matrixentries[(v,c)])
        end
    end
    for v in system.reverse_dfs_list
        ldu_backsubstitution_d!(vectorentries[v], matrixentries[(v,v)])
        for p in parents(system,v) # parent
            ldu_backsubstitution_u!(vectorentries[v], vectorentries[p], matrixentries[(v,p)])
        end
    end
end

function ldu_solve!(system)
    ldu_factorization!(system)
    ldu_backsubstitution!(system)
    return
end  
