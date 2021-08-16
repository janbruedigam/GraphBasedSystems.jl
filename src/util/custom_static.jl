szeros(::Type{T}, N) where T = @SVector zeros(T, N)
szeros(N)= @SVector zeros(N)
szeros(::Type{T}, N1, N2) where T = @SMatrix zeros(T, N1, N2)
szeros(N1, N2)= @SMatrix zeros(N1, N2)

sones(::Type{T}, N) where T = @SVector ones(T, N)
sones(N)= @SVector ones(N)
sones(::Type{T}, N1, N2) where T = @SMatrix ones(T, N1, N2)
sones(N1, N2)= @SMatrix ones(N1, N2)

srand(::Type{T}, N) where T = @SVector rand(T, N)
srand(N)= @SVector rand(N)
srand(::Type{T}, N1, N2) where T = @SMatrix rand(T, N1, N2)
srand(N1, N2)= @SMatrix rand(N1, N2)


# To fix StaticArray bug
Base.:*(::SMatrix{N,0,T,0},::SVector{0}) where {T,N} = szeros(T,N)
Base.:*(::SMatrix{N,0,T,0},::SMatrix{0,M,T,0}) where {T,N,M} = szeros(T,N,M)
Base.inv(::SMatrix{0,0,T,0}) where T = szeros(T,0,0)
