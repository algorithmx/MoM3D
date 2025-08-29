"""
Solvers

Numerical linear solvers for MoM3D.

Provides direct and iterative Krylov solvers (GMRES, BiCGSTAB(l)) via IterativeSolvers.jl,
with optional preconditioning and structured reporting. Functions are written to be
backward-compatible: by default they return only the solution; set `return_report=true`
for a `SolverReport` with convergence information.
"""

module Solvers

using LinearAlgebra
using IterativeSolvers
using SparseArrays
"""
    struct SolverReport

Compact metadata describing a linear solve:
- method::Symbol           : :direct | :gmres | :bicgstab
- converged::Bool          : whether the solver signaled convergence
- iterations::Int          : number of iterations (best-effort, parsed from history)
- tolerance::Float64       : requested relative tolerance
- maxiter::Int             : maximum iterations requested
- residual_norm::Float64   : ||Z*x - b||_2 on return
- preconditioner::Symbol   : :none | :jacobi | :custom
- time_sec::Float64        : wall-clock seconds for the call
"""


export solve_mom_system, solve_direct, solve_iterative, create_diagonal_preconditioner

# Lightweight report to make solver behavior observable without breaking API
struct SolverReport
    method::Symbol
    converged::Bool
    iterations::Int
    tolerance::Float64
    maxiter::Int
    residual_norm::Float64
    preconditioner::Symbol
    time_sec::Float64
end

# Public entry point. Backward-compatible: returns just the solution by default.
"""
    solve_mom_system(Z, b; method=:direct, tolerance=1e-6, max_iterations=1000,
                      preconditioner=:none, Pl=nothing, Pr=nothing,
                      restart=0, initial_guess=nothing, return_report=false)

Solve the MoM linear system Z*x = b using either a direct solver or an iterative
Krylov method from IterativeSolvers.jl.

Arguments:
- Z::AbstractMatrix{ComplexF64}
- b::AbstractVector{ComplexF64}

Keywords:
- method::Symbol = :direct | :gmres | :bicgstab
- tolerance::Float64 = relative tolerance used by the iterative solver
- max_iterations::Int = iteration cap
- preconditioner::Symbol = :none | :jacobi (ignored if Pl/Pr provided)
- Pl, Pr = optional left/right preconditioners compatible with IterativeSolvers
- restart::Int = GMRES restart parameter (0 for no restart)
- initial_guess = optional initial guess vector
- return_report::Bool = when true, also return a SolverReport

Returns:
- x or (x, report::SolverReport)
"""

function solve_mom_system(Z::AbstractMatrix{ComplexF64}, b::AbstractVector{ComplexF64};
                          method::Symbol = :direct,
                          tolerance::Float64 = 1e-6,
                          max_iterations::Int = 1000,
                          # preconditioning
                          preconditioner::Symbol = :none,  # :none | :jacobi
                          Pl = nothing, Pr = nothing,
                          # GMRES options
                          restart::Int = 0,
                          initial_guess::Union{Nothing,AbstractVector{ComplexF64}} = nothing,
                          return_report::Bool = false)
    if method == :direct
        x, rep = solve_direct(Z, b; return_report=true)
    elseif method == :gmres || method == :bicgstab
        x, rep = solve_iterative(Z, b, method, tolerance, max_iterations;
                                 preconditioner=preconditioner, Pl=Pl, Pr=Pr,
                                 restart=restart, initial_guess=initial_guess,
                                 return_report=true)
    else
        error("Unknown solver method: $method")
    end
    return return_report ? (x, rep) : x
end

"""
    solve_direct(Z, b; return_report=false)

Solve Z*x = b using a direct backslash solve. If the solve fails (e.g., near-singular),
a tiny Tikhonov regularization λI is added and the solve is retried.

Returns x or (x, SolverReport) when return_report=true.
"""

function solve_direct(Z::AbstractMatrix{ComplexF64}, b::AbstractVector{ComplexF64};
                      return_report::Bool=false)
    t0 = time()
    try
        x = Z \ b
        t = time() - t0
        rep = SolverReport(:direct, true, 1, 0.0, 1, norm(Z * x - b), :none, t)
        return return_report ? (x, rep) : x
    catch e
        @warn "Direct solver failed, trying with Tikhonov regularization" exception=e
        n = size(Z, 1)
        λ = 1e-12
        x = (Z + λ * I(n)) \ b
        t = time() - t0
        rep = SolverReport(:direct, true, 1, 0.0, 1, norm(Z * x - b), :none, t)
        return return_report ? (x, rep) : x
    end
end

"""
    solve_iterative(Z, b, method, tolerance, max_iterations; preconditioner=:none,
                    Pl=nothing, Pr=nothing, restart=0, initial_guess=nothing,
                    return_report=false)

Internal helper powering iterative solves (GMRES, BiCGSTAB(l)) with optional
preconditioning. See `solve_mom_system` for the public API.

Returns x or (x, SolverReport) when return_report=true.
"""

# Iterative methods with optional preconditioning and logging
function solve_iterative(Z::AbstractMatrix{ComplexF64}, b::AbstractVector{ComplexF64},
                         method::Symbol, tolerance::Float64, max_iterations::Int;
                         preconditioner::Symbol = :none, Pl = nothing, Pr = nothing,
                         restart::Int = 0,
                         initial_guess::Union{Nothing,AbstractVector{ComplexF64}} = nothing,
                         return_report::Bool=false)
    # Build preconditioner if requested and not provided
    Pl_eff = Pl
    Pr_eff = Pr
    pc_used = :none
    if Pl === nothing && Pr === nothing
        if preconditioner == :jacobi
            Pl_eff = create_diagonal_preconditioner(Z)
            pc_used = :jacobi
        elseif preconditioner == :none
            pc_used = :none
        else
            @warn "Unknown preconditioner symbol $preconditioner; proceeding without preconditioning"
            pc_used = :none
        end
    else
        pc_used = :custom
    end

    x = initial_guess === nothing ? zeros(ComplexF64, size(Z, 2)) : copy(initial_guess)
    t0 = time()

    # Run selected method with logging enabled (use in-place variants to control initial guess)
    if method == :gmres
        restart_val = max(restart, 0)
        # Try reltol (older IterativeSolvers) then fall back to rtol
        try
            hist = gmres!(x, Z, b; Pl=Pl_eff, Pr=Pr_eff, reltol=tolerance,
                          maxiter=max_iterations, restart=restart_val, log=true)
        catch
            try
                hist = gmres!(x, Z, b; Pl=Pl_eff, reltol=tolerance,
                              maxiter=max_iterations, restart=restart_val, log=true)
            catch
                # Final fallback for newer API: rtol
                try
                    hist = gmres!(x, Z, b; Pl=Pl_eff, Pr=Pr_eff, rtol=tolerance,
                                  maxiter=max_iterations, restart=restart_val, log=true)
                catch
                    hist = gmres!(x, Z, b; Pl=Pl_eff, rtol=tolerance,
                                  maxiter=max_iterations, restart=restart_val, log=true)
                end
            end
        end
        converged, iters = _extract_convergence(hist)
    elseif method == :bicgstab
        # Older IterativeSolvers expect reltol and max_mv_products, not rtol/maxiter
        max_mv_products = max(2 * max_iterations, max_iterations)
        try
            hist = bicgstabl!(x, Z, b; Pl=Pl_eff, reltol=tolerance,
                              max_mv_products=max_mv_products, log=true)
        catch
            # Fallback without Pl if keyword unsupported
            hist = bicgstabl!(x, Z, b; reltol=tolerance,
                              max_mv_products=max_mv_products, log=true)
        end
        converged, iters = _extract_convergence(hist)
    else
        error("Unknown iterative method: $method")
    end

    t = time() - t0
    resn = norm(Z * x - b)
    rep = SolverReport(method, converged, iters, tolerance, max_iterations, resn, pc_used, t)
    return return_report ? (x, rep) : x
end

function condition_number_estimate(Z::AbstractMatrix{ComplexF64})
    return cond(Matrix(Z))
end

# Jacobi (diagonal) preconditioner: P = diag(Z)
# Note: we build P, not P^{-1}; solvers use P \ v internally.
function create_diagonal_preconditioner(Z::AbstractMatrix{ComplexF64})
    n = min(size(Z,1), size(Z,2))
    d = [Z[i,i] for i in 1:n]
    # Guard against zeros on the diagonal
    repaired = similar(d)
    @inbounds for i in eachindex(d)
        repaired[i] = abs(d[i]) < eps(real(1.0)) ? one(ComplexF64) : d[i]
    end
    return Diagonal(repaired)
end

# --- helpers ---
# Try to robustly extract convergence info from IterativeSolvers history
function _extract_convergence(hist)
    # Attempt a few common patterns across IterativeSolvers versions
    iters = 0
    try
        if hasproperty(hist, :iters)
            iters = getproperty(hist, :iters)
        elseif hasproperty(hist, :numiter)
            iters = getproperty(hist, :numiter)
        else
            # Fallback: parse from string like "Converged after N iterations"
            m = match(r"(\d+)", string(hist))
            iters = m === nothing ? 0 : parse(Int, m.captures[1])
        end
    catch
        iters = 0
    end
    # Convergence boolean: search for "Converged" token or check status
    converged = occursin("Converged", string(hist))
    return converged, iters
end

end
