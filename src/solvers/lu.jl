# LU factorization for unsymmetric systems

function lu_factorization_acyclic!(diagonal_v, offdiagonal_l, diagonal_c, offdiagonal_u, diagonal_inverse_c)
    if diagonal_inverse_c.isinverted
        invdiagonal_c = diagonal_inverse_c.value
    else
        invdiagonal_c = inv(diagonal_c.value)
        diagonal_inverse_c.value = invdiagonal_c
        diagonal_inverse_c.isinverted = true
    end
    offdiagonal_l.value = offdiagonal_l.value * invdiagonal_c
    diagonal_v.value -= offdiagonal_l.value*offdiagonal_u.value
    return
end
function lu_factorization_cyclic!(entry_lu, offdiagonal_lu, offdiagonal_ul)
    entry_lu.value -= offdiagonal_lu.value*offdiagonal_ul.value
    return
end

function lu_factorization!(system::System)
    matrix_entries = system.matrix_entries
    diagonal_inverses = system.diagonal_inverses
    acyclic_children = system.acyclic_children
    cyclic_children = system.cyclic_children

    for v in system.dfs_list
        for c in acyclic_children[v]
            lu_factorization_acyclic!(matrix_entries[v,v], matrix_entries[v,c], matrix_entries[c,c], matrix_entries[c,v], diagonal_inverses[c])
        end
        for c in cyclic_children[v]
            for cc in cyclic_children[v]
                cc == c && break 
                (cc ∉ acyclic_children[c] && cc ∉ cyclic_children[c]) && continue
                lu_factorization_cyclic!(matrix_entries[v,c], matrix_entries[v,cc], matrix_entries[cc,c])
                lu_factorization_cyclic!(matrix_entries[c,v], matrix_entries[c,cc], matrix_entries[cc,v])
            end
            lu_factorization_acyclic!(matrix_entries[v,v], matrix_entries[v,c], matrix_entries[c,c], matrix_entries[c,v], diagonal_inverses[c])
        end
    end
    return
end

function lu_backsubstitution_l!(vector_v, offdiagonal, vector_c)
    vector_v.value -= offdiagonal.value * vector_c.value
    return
end
function lu_backsubstitution_u!(vector_v, offdiagonal, vector_p)
    vector_v.value -= offdiagonal.value * vector_p.value
    return
end
function lu_backsubstitution_d!(vector, diagonal, diagonal_inverse)
    if diagonal_inverse.isinverted
        vector.value = diagonal_inverse.value * vector.value
    else
        vector.value = diagonal.value \ vector.value
    end
    return
end

function lu_backsubstitution!(system::System)
    matrix_entries = system.matrix_entries
    diagonal_inverses = system.diagonal_inverses
    vector_entries = system.vector_entries
    acyclic_children = system.acyclic_children
    cyclic_children = system.cyclic_children
    parents = system.parents
    dfs_list = system.dfs_list

    for v in dfs_list
        for c in cyclic_children[v]
            lu_backsubstitution_l!(vector_entries[v], matrix_entries[v,c], vector_entries[c])
        end
        for c in acyclic_children[v]
            lu_backsubstitution_l!(vector_entries[v], matrix_entries[v,c], vector_entries[c])
        end
    end
    for v in reverse(dfs_list)
        for p in parents[v]
            lu_backsubstitution_u!(vector_entries[v], matrix_entries[v,p], vector_entries[p])
        end
        lu_backsubstitution_d!(vector_entries[v], matrix_entries[v,v], diagonal_inverses[v])
    end
end

function lu_solve!(system::System)
    reset_inverse_diagonals!(system)
    lu_factorization!(system)
    lu_backsubstitution!(system)
    return
end
