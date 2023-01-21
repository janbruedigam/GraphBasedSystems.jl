using GraphBasedSystems
using LinearAlgebra

include("adjacency_matrix.jl")

for i=1:5
    r1 = randn()
    r2 = randn()
    system = System{Float64}(A, rand(0:3, size(A)[1]))
    initialize!(system)

    F = full_matrix(system)
    f = full_vector(system)

    system = r1*system+r2*system
    F = r1*F+r2*F
    f = r1*f+r2*f
    @test maximum(abs.(full_vector(system)-f)) < 1e-3
    @test maximum(abs.(full_matrix(system)-F)) < 1e-3

    system = r1*system-r2*system
    F = r1*F-r2*F
    f = r1*f-r2*f
    @test maximum(abs.(full_vector(system)-f)) < 1e-3
    @test maximum(abs.(full_matrix(system)-F)) < 1e-3

    system = system*r1
    F = F*r1
    f = f*r1
    @test maximum(abs.(full_vector(system)-f)) < 1e-3
    @test maximum(abs.(full_matrix(system)-F)) < 1e-3

    system = r1*system
    F = r1*F
    f = r1*f
    @test maximum(abs.(full_vector(system)-f)) < 1e-3
    @test maximum(abs.(full_matrix(system)-F)) < 1e-3

    system = system/r1
    F = F/r1
    f = f/r1
    @test maximum(abs.(full_vector(system)-f)) < 1e-3
    @test maximum(abs.(full_matrix(system)-F)) < 1e-3

    system = r1\system
    F = r1\F
    f = r1\f
    @test maximum(abs.(full_vector(system)-f)) < 1e-3
    @test maximum(abs.(full_matrix(system)-F)) < 1e-3

    system = r1*system*r2*system
    F = r1*F*r2*F
    @test maximum(abs.(full_matrix(system)-F)) < 1e-3
end