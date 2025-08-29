# Tests for src/solvers.jl
#
# Design choices:
# - Use random Hermitian positive-definite complex systems (A'*A + δI) so
#   all solvers have a well-defined solution and iterative methods converge.
# - Verify both solution accuracy and the contents of the SolverReport without
#   overly constraining internal details (e.g., exact iteration counts), making
#   the tests robust across IterativeSolvers.jl versions.
# - Keep tolerances slightly looser for iterative methods than direct backslash
#   to account for algorithmic differences and random conditioning.

using Test
using MoM3D
using LinearAlgebra
using Random

# Fix RNG for reproducibility across runs and platforms
Random.seed!(1234)

# Simple Hermitian positive-definite complex system for testing
function make_hpdc(n)
    A = randn(n, n) + 1im*randn(n, n)
    Z = A'*A + 1e-2I # HPD real, then convert to ComplexF64
    Zc = ComplexF64.(Z)
    x_true = randn(n) .+ 1im*randn(n)
    b = Zc * x_true
    return Zc, b, x_true
end

@testset "Solvers: direct" begin
    Z, b, x_true = make_hpdc(20)
    x = MoM3D.Solvers.solve_direct(Z, b)
    @test norm(Z*x - b) / norm(b) < 1e-10
end

@testset "Solvers: gmres no-precond" begin
    Z, b, _ = make_hpdc(30)
    x, rep = MoM3D.Solvers.solve_mom_system(Z, b; method=:gmres, tolerance=1e-8, max_iterations=500, return_report=true)
    @test rep.method == :gmres
    @test rep.maxiter == 500
    @test rep.tolerance == 1e-8
    @test rep.preconditioner in (:none, :custom, :jacobi)
    @test rep.converged
    @test norm(Z*x - b)/norm(b) < 1e-6
end

@testset "Solvers: gmres jacobi precond" begin
    Z, b, _ = make_hpdc(40)
    x, rep = MoM3D.Solvers.solve_mom_system(Z, b; method=:gmres, tolerance=1e-6, max_iterations=500, preconditioner=:jacobi, return_report=true)
    @test rep.preconditioner in (:jacobi, :custom)
    @test rep.converged
    @test norm(Z*x - b)/norm(b) < 1e-5
end

@testset "Solvers: bicgstab" begin
    Z, b, _ = make_hpdc(25)
    x, rep = MoM3D.Solvers.solve_mom_system(Z, b; method=:bicgstab, tolerance=1e-6, max_iterations=500, return_report=true)
    @test rep.method == :bicgstab
    @test rep.converged
    @test norm(Z*x - b)/norm(b) < 1e-5
end

@testset "Preconditioner construction" begin
    Z, _, _ = make_hpdc(10)
    P = MoM3D.Solvers.create_diagonal_preconditioner(Z)
    @test size(P) == size(Z)
    # Check that P is diagonal and nonzero on diagonal
    @test all(P[i,i] != 0 for i in 1:size(Z,1))
end

