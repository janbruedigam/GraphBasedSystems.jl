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
zerodimslash(A,B) = A/B
zerodimslash(::SMatrix{0,N,T,0},::SMatrix{N,N,T,N2}) where {T,N,N2} = SMatrix{0,N,T,0}()
zerodimslash(::SMatrix{N,0,T,0},::SMatrix{0,0,T,0}) where {T,N} = SMatrix{N,0,T,0}()
zerodimslash(::SMatrix{0,0,T,0},::SMatrix{0,0,T,0}) where {T} = SMatrix{0,0,T,0}()
zerodimbackslash(A,B) = A\B
zerodimbackslash(::SMatrix{N,N,T,N2},::SMatrix{N,0,T,0}) where {T,N,N2} = SMatrix{N,0,T,0}()
zerodimbackslash(::SMatrix{0,0,T,0},::SMatrix{0,N,T,0}) where {T,N} = SMatrix{0,N,T,0}()
zerodimbackslash(::SMatrix{0,0,T,0},::SMatrix{0,0,T,0}) where {T} = SMatrix{0,0,T,0}()
Base.:*(::SMatrix{N,0,T,0},::SVector{0}) where {T,N} = szeros(T,N)
Base.:*(::SMatrix{N,0,T,0},::SMatrix{0,M,T,0}) where {T,N,M} = szeros(T,N,M)