module Solvers

using LinearAlgebra
using IterativeSolvers
using SparseArrays

export solve_mom_system, solve_direct, solve_iterative

function solve_mom_system(Z::Matrix{ComplexF64}, b::Vector{ComplexF64}; 
                         method::Symbol=:direct, tolerance::Float64=1e-6, max_iterations::Int=1000)
    if method == :direct
        return solve_direct(Z, b)
    elseif method == :gmres
        return solve_iterative(Z, b, :gmres, tolerance, max_iterations)
    elseif method == :bicgstab
        return solve_iterative(Z, b, :bicgstab, tolerance, max_iterations)
    else
        error("Unknown solver method: $method")
    end
end

function solve_direct(Z::Matrix{ComplexF64}, b::Vector{ComplexF64})
    try
        return Z \ b
    catch e
        @warn "Direct solver failed, trying with regularization"
        n = size(Z, 1)
        regularization = 1e-12 * I(n)
        return (Z + regularization) \ b
    end
end

function solve_iterative(Z::Matrix{ComplexF64}, b::Vector{ComplexF64}, 
                        method::Symbol, tolerance::Float64, max_iterations::Int)
    
    if method == :gmres
        result = gmres(Z, b, rtol=tolerance, maxiter=max_iterations, verbose=true)
    elseif method == :bicgstab
        result = bicgstabl(Z, b, rtol=tolerance, maxiter=max_iterations, verbose=true)
    else
        error("Unknown iterative method: $method")
    end
    
    return result
end

function condition_number_estimate(Z::Matrix{ComplexF64})
    return cond(Z)
end

function create_diagonal_preconditioner(Z::Matrix{ComplexF64})
    diagonal_elements = diag(Z)
    return Diagonal(1.0 ./ diagonal_elements)
end

end
