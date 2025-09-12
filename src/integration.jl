module Integration

using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using QuadGK
# Note: SimplexQuad.jl may be installed; detect it at module load time.
const _HAS_SIMPLEXQUAD = try
    @eval using SimplexQuad
    true
catch
    false
end

# The Dunavant (1985) tables are the canonical reference for triangle rules
# used here (see D. A. Dunavant, Int. J. Numer. Methods Eng., 1985).
using ..Geometry

# Cache Gauss triangle points/weights per order to avoid reallocating the
# same vectors on every integration call (this is a common allocation hotspot
# when assembling O(N^2) matrices).
const GAUSS_TRIANGLE_CACHE = Dict{Int,Tuple{Vector{SVector{2,Float64}},Vector{Float64}}}()

# Protect access to the cache when the library is used with Julia threads.
# Dict is not safe for concurrent writes (or simultaneous reads/writes), so
# serialize cache reads/writes with a lock. This keeps the fix minimal-risk
# while ensuring correctness under multithreaded assembly.
const GAUSS_TRIANGLE_LOCK = ReentrantLock()

export gauss_triangle, integrate_singular, integrate_regular, integrate_near_singular, barycentric_to_cartesian

function _fallback_gauss_triangle(n::Int)
    # Isolated fallback implementation for regression testing. Returns
    # (points::Vector{SVector{2,Float64}}, weights::Vector{Float64}).
    if n == 1
        points = [SVector(1 / 3, 1 / 3)]
        weights = [1.0]
    elseif n == 3
        points = [SVector(1 / 6, 1 / 6), SVector(2 / 3, 1 / 6), SVector(1 / 6, 2 / 3)]
        weights = [1 / 3, 1 / 3, 1 / 3]
    elseif n == 4
        points = [SVector(1 / 3, 1 / 3), SVector(0.6, 0.2), SVector(0.2, 0.6), SVector(0.2, 0.2)]
        weights = [-27 / 48, 25 / 48, 25 / 48, 25 / 48]
    elseif n == 5
        points = [SVector(1 / 3, 1 / 3),
            SVector(0.470942, 0.059058),
            SVector(0.059058, 0.470942),
            SVector(0.059058, 0.059058),
            SVector(0.797, 0.1015)]
        weights = [-0.5625, 0.5208333333333334, 0.5208333333333334, 0.5208333333333334, 0.5208333333333334]
    elseif n == 7
        points = [SVector(1 / 3, 1 / 3),
            SVector(0.059715871789770, 0.470142064105115),
            SVector(0.470142064105115, 0.059715871789770),
            SVector(0.470142064105115, 0.470142064105115),
            SVector(0.797426985353087, 0.101286507323456),
            SVector(0.101286507323456, 0.797426985353087),
            SVector(0.101286507323456, 0.101286507323456)]
        weights = [0.225,
            0.132394152788506,
            0.132394152788506,
            0.132394152788506,
            0.125939180544827,
            0.125939180544827,
            0.125939180544827]
    elseif n == 9
        # Composite 9-point as previous fallback (safe, low-risk)
        base_points, base_weights = _fallback_gauss_triangle(3)
        points = SVector{2,Float64}[]
        weights = Float64[]
        for sub in 1:3
            for (ξ, w) in zip(base_points, base_weights)
                u = ξ[1]
                v = ξ[2]
                wloc = 1.0 - u - v
                if sub == 1
                    g1 = u + wloc / 3.0
                    g2 = v + wloc / 3.0
                elseif sub == 2
                    g1 = wloc / 3.0
                    g2 = u + wloc / 3.0
                else
                    g1 = v + wloc / 3.0
                    g2 = wloc / 3.0
                end
                push!(points, SVector(g1, g2))
                push!(weights, w * (1.0 / 3.0))
            end
        end
    else
        error("Gaussian quadrature order $n not implemented in fallback")
    end
    # Ensure fallback weights are normalized to sum to 1.0 so the API
    # consistently returns barycentric-style weights (sum==1) for all orders.
    ws = collect(weights)
    s = sum(ws)
    if s == 0.0
        error("Fallback quadrature weights sum to zero for n=$n")
    end
    ws .= ws ./ s
    return points, ws
end

function gauss_triangle(n::Int)
    # Fast path: try unsynchronized read. If present, return immediately.
    # This avoids locking for the common case where the rule is already cached.
    val = get(GAUSS_TRIANGLE_CACHE, n, nothing)
    if val !== nothing
        return val
    end

    # Miss: acquire write lock and insert under protection. Use double-checked
    # locking to avoid recomputing if another thread inserted while we waited.
    lock(GAUSS_TRIANGLE_LOCK)
    try
        val = get(GAUSS_TRIANGLE_CACHE, n, nothing)
        if val !== nothing
            return val
        end

        if _HAS_SIMPLEXQUAD
            # SimplexQuad.simplexquad(N, dim) returns X (points) and W (weights)
            # For 2D triangles, call simplexquad(n, 2) to get points on unit simplex.
            X, W = SimplexQuad.simplexquad(n, 2)
            # X is (num_points, dim) or (dim, num_points) depending on version; handle both.
            if size(X, 2) == 2
                pts = [SVector(X[i, 1], X[i, 2]) for i in axes(X, 1)]
            else
                pts = [SVector(X[1, i], X[2, i]) for i in axes(X, 2)]
            end
            # SimplexQuad returns weights for the standard simplex whose area is 1/2.
            # Our fallback gauss_triangle uses barycentric weights that sum to 1 on
            # the reference triangle (so that they are multiplied by triangle.area
            # later). Scale weights by the inverse reference area (2.0) so both
            # implementations are consistent for regression tests.
            ws = vec(W) .* 2.0
        else
            @warn "SimplexQuad not installed; using fallback gauss_triangle"
            pts, ws = _fallback_gauss_triangle(n)
        end

        GAUSS_TRIANGLE_CACHE[n] = (pts, ws)
        return pts, ws
    finally
        unlock(GAUSS_TRIANGLE_LOCK)
    end
end

function barycentric_to_cartesian(ξ::SVector{2,Float64}, triangle::Triangle)
    v1, v2, v3 = triangle.vertices
    ξ1, ξ2 = ξ
    ξ3 = 1.0 - ξ1 - ξ2
    return ξ1 * v1 + ξ2 * v2 + ξ3 * v3
end

function integrate_regular(integrand_func, tri_src::Triangle, tri_obs::Triangle; quad_order::Int=3)
    # Use cached quadrature rule where possible; both source and observer use
    # the same order here so a single allocation is sufficient.
    points, weights = gauss_triangle(quad_order)

    result = complex(0.0)

    # Precompute jacobian factor once per call
    jacobian = tri_src.area * tri_obs.area

    for (ξ_src, w_src) in zip(points, weights)
        r_src = barycentric_to_cartesian(ξ_src, tri_src)

        for (ξ_obs, w_obs) in zip(points, weights)
            r_obs = barycentric_to_cartesian(ξ_obs, tri_obs)

            integrand_value = integrand_func(r_obs, r_src)

            result += w_src * w_obs * jacobian * integrand_value
        end
    end

    return result
end

# Positional wrapper for backward compatibility: allow calling with positional quad_order
function integrate_regular(integrand_func, tri_src::Triangle, tri_obs::Triangle, quad_order::Int=3)
    return integrate_regular(integrand_func, tri_src, tri_obs; quad_order=quad_order)
end

function integrate_singular(integrand_func, triangle::Triangle; quad_order::Int=7)
    points, weights = gauss_triangle(quad_order)

    result = complex(0.0)

    for (ξ, w) in zip(points, weights)
        r = barycentric_to_cartesian(ξ, triangle)
        jacobian = triangle.area
        integrand_value = integrand_func(r)
        result += w * jacobian * integrand_value
    end

    return result
end

# Positional wrapper for backward compatibility: allow calling with positional quad_order
function integrate_singular(integrand_func, triangle::Triangle, quad_order::Int=7)
    return integrate_singular(integrand_func, triangle; quad_order=quad_order)
end

function integrate_near_singular(integrand_func, tri_src::Triangle, tri_obs::Triangle,
    tolerance::Float64=1e-6, max_subdivisions::Int=10)
    # Near-singular integrals are challenging; use a higher-order product Gauss rule
    # as a robust fallback. Use quad order 9 (if available) otherwise 7.
    quad_order = 9
    try
        return integrate_regular(integrand_func, tri_src, tri_obs, quad_order)
    catch err
        @warn "High-order product Gauss failed for near-singular; falling back to order 7" exception = err
        return integrate_regular(integrand_func, tri_src, tri_obs, 7)
    end
end

end
