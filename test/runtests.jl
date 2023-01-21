using Test
using SafeTestsets

@safetestset "LDU Tests" begin
    include("ldu_test.jl")
end

@safetestset "LU Tests" begin
    include("lu_test.jl")
end

@safetestset "LDLᵀ Tests" begin
    include("ldlt_test.jl")
end

@safetestset "LLᵀ Tests" begin
    include("llt_test.jl")
end

@safetestset "Algebra Tests" begin
    include("algebra.jl")
end

@safetestset "Interface Tests" begin
    include("interface.jl")
end