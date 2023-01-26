# LDU factorization for unsymmetric systems

function ldu_factorization_acyclic!(diagonal_v, offdiagonal_l, diagonal_c, offdiagonal_u, diagonal_inverse_c)
    if diagonal_inverse_c.isinverted
        invdiagonal_c = diagonal_inverse_c.value
    else
        invdiagonal_c = inv(diagonal_c.value)
        diagonal_inverse_c.value = invdiagonal_c
        diagonal_inverse_c.isinverted = true
    end
    offdiagonal_l.value = offdiagonal_l.value * invdiagonal_c
    offdiagonal_u.value = invdiagonal_c * offdiagonal_u.value
    diagonal_v.value -= offdiagonal_l.value * diagonal_c.value * offdiagonal_u.value
    return
end
function ldu_factorization_cyclic!(entry_lu, offdiagonal_lu, diagonal_c, offdiagonal_ul)
    entry_lu.value -= offdiagonal_lu.value * diagonal_c.value * offdiagonal_ul.value
    return
end

function ldu_factorization!(system::System)
    matrix_entries = system.matrix_entries
    diagonal_inverses = system.diagonal_inverses
    acyclic_children = system.acyclic_children
    cyclic_children = system.cyclic_children
    actives = system.actives

    reset_inverse_diagonals!(system)

    for v in system.dfs_list
        !actives[v] && continue
        for c in acyclic_children[v]
            !actives[c] && continue
            ldu_factorization_acyclic!(matrix_entries[v,v], matrix_entries[v,c], matrix_entries[c,c], matrix_entries[c,v], diagonal_inverses[c])
        end
        for c in cyclic_children[v]
            !actives[c] && continue
            for cc in cyclic_children[v]
                !actives[cc] && continue
                cc == c && break 
                (cc ∉ acyclic_children[c] && cc ∉ cyclic_children[c]) && continue
                ldu_factorization_cyclic!(matrix_entries[v,c], matrix_entries[v,cc], matrix_entries[cc,cc], matrix_entries[cc,c])
                ldu_factorization_cyclic!(matrix_entries[c,v], matrix_entries[c,cc], matrix_entries[cc,cc], matrix_entries[cc,v])
            end
            ldu_factorization_acyclic!(matrix_entries[v,v], matrix_entries[v,c], matrix_entries[c,c], matrix_entries[c,v], diagonal_inverses[c])
        end
    end
    return
end

function ldu_backsubstitution_l!(vector_v, offdiagonal, vector_c)
    vector_v.value -= offdiagonal.value * vector_c.value
    return
end
function ldu_backsubstitution_u!(vector_v, offdiagonal, vector_p)
    vector_v.value -= offdiagonal.value * vector_p.value
    return
end
function ldu_backsubstitution_d!(vector, diagonal, diagonal_inverse)
    if diagonal_inverse.isinverted
        vector.value = diagonal_inverse.value * vector.value
    else
        vector.value = diagonal.value \ vector.value
    end
    return
end

function ldu_backsubstitution!(system::System)
    matrix_entries = system.matrix_entries
    diagonal_inverses = system.diagonal_inverses
    vector_entries = system.vector_entries
    acyclic_children = system.acyclic_children
    cyclic_children = system.cyclic_children
    parents = system.parents
    dfs_list = system.dfs_list
    actives = system.actives

    for v in dfs_list
        !actives[v] && continue
        for c in cyclic_children[v]
            !actives[c] && continue
            ldu_backsubstitution_l!(vector_entries[v], matrix_entries[v,c], vector_entries[c])
        end
        for c in acyclic_children[v]
            !actives[c] && continue
            ldu_backsubstitution_l!(vector_entries[v], matrix_entries[v,c], vector_entries[c])
        end
    end
    for v in reverse(dfs_list)
        !actives[v] && continue
        ldu_backsubstitution_d!(vector_entries[v], matrix_entries[v,v], diagonal_inverses[v])
        for p in parents[v]
            !actives[p] && continue
            ldu_backsubstitution_u!(vector_entries[v], matrix_entries[v,p], vector_entries[p])
        end
    end
end

function ldu_solve!(system::System)
    ldu_factorization!(system)
    ldu_backsubstitution!(system)
    return
end  

function ldu_matrix_solve!(system::System, matrix::SparseMatrixCSC{Entry, Int64}; keep_vector = true)
    keep_vector && (vector_entries = deepcopy(system.vector_entries))
    ldu_factorization!(system)
    C = matrix_backsubsitution!(system, matrix, ldu_backsubstitution!)
    keep_vector && (system.vector_entries .= vector_entries)

    return C
end