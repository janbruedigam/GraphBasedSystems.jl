full_matrix(system::System) = full_matrix(system.matrix_entries, issymmetric(system), system.dims, system.dims)
# There probably exists a smarter way of getting the dense matrix from the spares one 
function full_matrix(matrix_entries::SparseMatrixCSC, symmetric::Bool, dimensions_rows, dimensions_cols)
    range_dict_rows = ranges(dimensions_rows)
    range_dict_cols = ranges(dimensions_cols)
    A = zeros(sum(dimensions_rows), sum(dimensions_cols))

    for (i,row) in enumerate(matrix_entries.rowval)
        col = findfirst(x -> i<x, matrix_entries.colptr)-1
        A[range_dict_rows[row],range_dict_cols[col]] = matrix_entries[row,col].value
        if symmetric && col != row
            A[range_dict_cols[col],range_dict_rows[row]] = matrix_entries[row,col].value'
        end
    end
    return A
end

full_vector(system::System) = full_vector(system.vector_entries, system.dims)
full_vector(vector_entries::AbstractVector, dimensions_rows) = vcat([getfield(vector_entries[i],:value) for i=1:size(dimensions_rows)[1]]...)