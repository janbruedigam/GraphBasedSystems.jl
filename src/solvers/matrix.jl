# Assumes quadratic matrices of same size
function matrix_backsubsitution!(systemA::System{N}, B::SparseMatrixCSC{Entry, Int64}, backsubsitution!) where N
    vector_entries = systemA.vector_entries
    dims = systemA.dims

    maxnnzC = Int(widelength(systemA.matrix_entries))
    C = _allocres((N,N), Int64,  Entry, maxnnzC)
    # _densestructure!(C)

    for i = 1:N
        for j = 1:N
            Bji = B[j,i]
            Bji isa Entry{Nothing} ? vector_entries[j] = Entry(dims[j], dims[i]) : vector_entries[j] = Bji
        end
        backsubsitution!(systemA)
        for j = 1:N
            C[j,i] = vector_entries[j]
        end
    end

    return C
end