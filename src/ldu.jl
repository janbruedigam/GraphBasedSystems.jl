function ldu_factorization_l!(offdiagonal, diagonal)
    offdiagonal.value = offdiagonal.value/diagonal.value
    return
end
function ldu_factorization_u!(offdiagonal, diagonal)
    offdiagonal.value = diagonal.value\offdiagonal.value
    return
end
function ldu_factorization_lud!(entry_v, offdiagonal_lu, diagonal_c, offdiagonal_ul)
    entry_v.value -= offdiagonal_lu.value*diagonal_c.value*offdiagonal_ul.value
    return
end

function ldu_factorization!(system)
    matrix_entries = system.matrix_entries
    acyclic_children = system.acyclic_children
    cycles = system.cycles

    for v in system.dfs_list
        for cyclic_children in cycles[v]
            for c in cyclic_children
                for cc in cyclic_children
                    cc == c && break 
                    cc ∉ acyclic_children[c] && continue
                    ldu_factorization_lud!(matrix_entries[(v,c)], matrix_entries[(v,cc)], matrix_entries[(cc,cc)], matrix_entries[(cc,c)])
                    ldu_factorization_lud!(matrix_entries[(c,v)], matrix_entries[(c,cc)], matrix_entries[(cc,cc)], matrix_entries[(cc,v)])
                end
                ldu_factorization_l!(matrix_entries[(v,c)], matrix_entries[(c,c)])
                ldu_factorization_u!(matrix_entries[(c,v)], matrix_entries[(c,c)])
                ldu_factorization_lud!(matrix_entries[(v,v)], matrix_entries[(v,c)], matrix_entries[(c,c)], matrix_entries[(c,v)])
            end
        end
        for c in acyclic_children[v]
            ldu_factorization_l!(matrix_entries[(v,c)], matrix_entries[(c,c)])
            ldu_factorization_u!(matrix_entries[(c,v)], matrix_entries[(c,c)])
            ldu_factorization_lud!(matrix_entries[(v,v)], matrix_entries[(v,c)], matrix_entries[(c,c)], matrix_entries[(c,v)])
        end
    end
    return
end

function ldu_backsubstitution_l!(vector_v, offdiagonal, vector_c)
    vector_v.value -= offdiagonal.value*vector_c.value
    return
end
function ldu_backsubstitution_u!(vector_v, offdiagonal, vector_p)
    vector_v.value -= offdiagonal.value*vector_p.value
    return
end
function ldu_backsubstitution_d!(vector, diagonal)
    vector.value = diagonal.value\vector.value
    return
end

function ldu_backsubstitution!(system)
    matrix_entries = system.matrix_entries
    vector_entries = system.vector_entries
    acyclic_children = system.acyclic_children
    cycles = system.cycles
    parents = system.parents

    for v in system.dfs_list
        for cyclic_children in cycles[v]
            for c in cyclic_children
                ldu_backsubstitution_l!(vector_entries[v], matrix_entries[(v,c)], vector_entries[c])
            end
        end
        for c in acyclic_children[v]
            ldu_backsubstitution_l!(vector_entries[v], matrix_entries[(v,c)], vector_entries[c])
        end
    end
    for v in system.reverse_dfs_list
        ldu_backsubstitution_d!(vector_entries[v], matrix_entries[(v,v)])
        for p in parents[v]
            ldu_backsubstitution_u!(vector_entries[v], matrix_entries[(v,p)], vector_entries[p])
        end
    end
end

function ldu_solve!(system)
    ldu_factorization!(system)
    ldu_backsubstitution!(system)
    return
end  
