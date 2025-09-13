using Test
using StaticArrays
using LinearAlgebra
using MoM3D

const EPS = 1e-10

@testset "GreenFunctions basic tests" begin
    # simple far-field check: for r large, green_3d ~ exp(-ikR)/(4πR)
    r_obs = SVector{3,Float64}(1.0, 0.0, 0.0)
    r_src = SVector{3,Float64}(0.0, 0.0, 0.0)
    k = 2.0
    g = green_3d(r_obs, r_src, k)
    R = norm(r_obs - r_src)
    @test isapprox(g, cis(-k*R) / (4π*R); atol=1e-12, rtol=1e-12)

    # gradient should be roughly along R_vec direction (not exported by MoM3D)
    grad = MoM3D.GreenFunctions.green_3d_gradient(r_obs, r_src, k)
    # ensure gradient is a 3-vector of complex numbers
    @test length(grad) == 3
    @test typeof(grad[1]) <: Complex
    # directional check: gradient should be parallel (or anti-parallel) to R_vec
    R_vec = r_obs - r_src
    proj = dot(real.(grad), R_vec) / (norm(R_vec) * norm(real.(grad)) + EPS)
    @test abs(proj) > 0.1

    # singular extraction: when points coincide, singular part ~ 1/(4πR)
    r_same = SVector{3,Float64}(0.0, 0.0, 0.0)
    se = MoM3D.GreenFunctions.green_3d_singular_extraction(r_same, r_same, k)
    # singular should be present and be a Float64
    @test typeof(se.singular) <: Float64

    # non-singular pair should split into singular + regular that sums to full green
    r2 = SVector{3,Float64}(0.1, 0.2, -0.1)
    parts = MoM3D.GreenFunctions.green_3d_singular_extraction(r2, r_src, k)
    full = green_3d(r2, r_src, k)
    reconstructed = parts.singular + parts.regular
    @test isapprox(reconstructed, full; atol=1e-12, rtol=1e-12)
end
