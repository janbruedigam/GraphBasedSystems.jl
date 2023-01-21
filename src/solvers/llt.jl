# LLᵀ (Cholesky) factorization for symmetric systems

function llt_factorization_acyclic!(diagonal_v, offdiagonal_l, diagonal_c, diagonal_inverse_c)
    if diagonal_inverse_c.isinverted
        invdiagonal_c = diagonal_inverse_c.value
    else
        invdiagonal_c = inv(diagonal_c.value)
        diagonal_inverse_c.value = invdiagonal_c
        diagonal_inverse_c.isinverted = true
    end
    offdiagonal_l.value = offdiagonal_l.value * invdiagonal_c
    diagonal_v.value -= offdiagonal_l.value * offdiagonal_l.value'
    return
end
function llt_factorization_cyclic!(entry_lu, offdiagonal_lu, offdiagonal_ul)
    entry_lu.value -= offdiagonal_lu.value * offdiagonal_ul.value'
    return
end
# TODO currently allocates memory
function llt_diagonal_root!(diagonal)
    diagonal.value = sqrt(diagonal.value)
end


function llt_factorization!(system::System)
    matrix_entries = system.matrix_entries
    diagonal_inverses = system.diagonal_inverses
    acyclic_children = system.acyclic_children
    cyclic_children = system.cyclic_children

    reset_inverse_diagonals!(system)

    for v in system.dfs_list
        for c in acyclic_children[v]
            llt_factorization_acyclic!(matrix_entries[v,v], matrix_entries[v,c], matrix_entries[c,c], diagonal_inverses[c])
        end
        for c in cyclic_children[v]
            for cc in cyclic_children[v]
                cc == c && break 
                (cc ∉ acyclic_children[c] && cc ∉ cyclic_children[c]) && continue
                llt_factorization_cyclic!(matrix_entries[v,c], matrix_entries[v,cc], matrix_entries[c,cc])
            end
            llt_factorization_acyclic!(matrix_entries[v,v], matrix_entries[v,c], matrix_entries[c,c], diagonal_inverses[c])
        end
        llt_diagonal_root!(matrix_entries[v,v])
    end
    return
end

function llt_backsubstitution_l!(vector_v, offdiagonal, vector_c)
    vector_v.value -= offdiagonal.value * vector_c.value
    return
end
function llt_backsubstitution_lt!(vector_v, offdiagonal, vector_p)
    vector_v.value -= offdiagonal.value' * vector_p.value
    return
end
function llt_backsubstitution_d!(vector, diagonal, diagonal_inverse)
    if diagonal_inverse.isinverted
        vector.value = diagonal_inverse.value * vector.value
    else
        invdiagonal = inv(diagonal.value)
        diagonal_inverse.value = invdiagonal
        diagonal_inverse.isinverted = true
        vector.value = diagonal_inverse.value * vector.value
    end
    return
end

function llt_backsubstitution!(system::System)
    matrix_entries = system.matrix_entries
    diagonal_inverses = system.diagonal_inverses
    vector_entries = system.vector_entries
    acyclic_children = system.acyclic_children
    cyclic_children = system.cyclic_children
    parents = system.parents
    dfs_list = system.dfs_list

    for v in dfs_list
        for c in cyclic_children[v]
            llt_backsubstitution_l!(vector_entries[v], matrix_entries[v,c], vector_entries[c])
        end
        for c in acyclic_children[v]
            llt_backsubstitution_l!(vector_entries[v], matrix_entries[v,c], vector_entries[c])
        end
        llt_backsubstitution_d!(vector_entries[v], matrix_entries[v,v], diagonal_inverses[v])
    end
    for v in reverse(dfs_list)
        for p in parents[v]
            llt_backsubstitution_lt!(vector_entries[v], matrix_entries[p,v], vector_entries[p])
        end
        llt_backsubstitution_d!(vector_entries[v], matrix_entries[v,v], diagonal_inverses[v])
    end
end

function llt_solve!(system::System)
    llt_factorization!(system)
    llt_backsubstitution!(system)
    return
end  

function llt_matrix_solve!(A::System, B::SparseMatrixCSC{Entry, Int64}; keep_vector = true)
    keep_vector && (vector_entries = deepcopy(A.vector_entries))
    llt_factorization!(A)
    C = matrix_backsubsitution!(A, B, llt_backsubstitution!)
    keep_vector && (A.vector_entries .= vector_entries)

    return C
end