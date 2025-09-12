module IntegralEquations
using ProgressMeter
using Printf
using LinearAlgebra
using StaticArrays
using ..Geometry
using ..BasisFunctions
using ..GreenFunctions: green_3d, green_3d_singular_extraction, wavenumber, green_3d_fast, FOUR_PI, μ0, ε0
using ..Integration

export assemble_efie_matrix, assemble_mfie_matrix, assemble_cfie_matrix
export compute_excitation_vector

# Module-level lock for safely building/storing per-mesh integration caches
const INTEGRATION_CACHE_LOCK = ReentrantLock()

function assemble_efie_matrix(mesh::Mesh3D, frequency::Float64; progress::Bool=false, parallel::Bool=false, chunk_size::Int=0)
    k = wavenumber(frequency)
    n_edges = mesh.num_edges
    Z = zeros(ComplexF64, n_edges, n_edges)
    p = progress ? ProgressMeter.Progress(n_edges; desc="Assembling EFIE") : nothing

    # Precompute RWG basis and warm common quadrature rules on the main thread
    # to avoid races or duplicated construction when threading.
    rwgs = get_rwgs(mesh)

    # If mesh contains a persistent integration cache, use it; otherwise build
    # the caches and persist them on the mesh. Double-checked locking ensures
    # only one thread builds and stores the cache if multiple assembly calls
    # happen concurrently.
    if mesh.integration_cache !== nothing
        caches = mesh.integration_cache
    else
        lock(INTEGRATION_CACHE_LOCK)
        try
            if mesh.integration_cache !== nothing
                caches = mesh.integration_cache
            else
                orders = (3, 5, 7, 9)
                gauss_rules = Dict{Int,Tuple{Vector{SVector{2,Float64}},Vector{Float64}}}()
                for ord in orders
                    gauss_rules[ord] = gauss_triangle(ord)
                end

                n_tri = length(mesh.triangles)
                # Map triangle -> list of RWG indices referencing it
                triangle_rwgs = [Int[] for _ in 1:n_tri]
                for (i, rwg) in enumerate(rwgs)
                    push!(triangle_rwgs[rwg.triangle_plus], i)
                    push!(triangle_rwgs[rwg.triangle_minus], i)
                end

                # tri_cart_pts[ord] = Vector of cartesian points per triangle
                tri_cart_pts = Dict{Int,Vector{Vector{SVector{3,Float64}}}}()
                # tri_rwg_vals[ord] = Vector per triangle of Dict{rwg_idx => Vector{SVector{3,Float64}}}
                tri_rwg_vals = Dict{Int,Vector{Dict{Int,Vector{SVector{3,Float64}}}}}()

                for ord in orders
                    bary_pts, weights = gauss_rules[ord]
                    cart_pts_per_tri = Vector{Vector{SVector{3,Float64}}}(undef, n_tri)
                    rwgvals_per_tri = Vector{Dict{Int,Vector{SVector{3,Float64}}}}(undef, n_tri)

                    Threads.@threads for tid in 1:n_tri
                        tri = mesh.triangles[tid]
                        # compute cartesian quadrature points for this triangle and order
                        cps = Vector{SVector{3,Float64}}(undef, length(bary_pts))
                        for (j, ξ) in enumerate(bary_pts)
                            cps[j] = barycentric_to_cartesian(ξ, tri)
                        end
                        cart_pts_per_tri[tid] = cps

                        # precompute RWG evaluations for RWGs that reference this triangle
                        dict = Dict{Int,Vector{SVector{3,Float64}}}()
                        for rwg_idx in triangle_rwgs[tid]
                            vals = Vector{SVector{3,Float64}}(undef, length(cps))
                            for (j, pt) in enumerate(cps)
                                vals[j] = evaluate_rwg(rwgs[rwg_idx], pt, mesh)
                            end
                            dict[rwg_idx] = vals
                        end
                        rwgvals_per_tri[tid] = dict
                    end

                    tri_cart_pts[ord] = cart_pts_per_tri
                    tri_rwg_vals[ord] = rwgvals_per_tri
                end

                caches = Dict{Symbol,Any}(
                    :gauss_rules => gauss_rules,
                    :tri_cart_pts => tri_cart_pts,
                    :tri_rwg_vals => tri_rwg_vals
                )

                mesh.integration_cache = caches
            end
        finally
            unlock(INTEGRATION_CACHE_LOCK)
        end
    end

    if parallel && Threads.nthreads() > 1
        # Keep outer loop single-threaded so ProgressMeter remains correct.
        # Parallelize the inner loop over `n` so each row is assembled in parallel
        # across threads while the main thread advances the progress meter.
        # Determine block size: if chunk_size <= 0, auto compute based on n_edges
        if chunk_size <= 0
            # Aim for ~8-16 blocks per thread to balance load and overhead
            nthreads = Threads.nthreads()
            nblocks = max(16, 8 * nthreads)
            bs = max(1, div(n_edges + nblocks - 1, nblocks))
        else
            bs = max(1, chunk_size)
        end

        for m in 1:n_edges
            rwg_m = rwgs[m]

            nblocks = div(n_edges + bs - 1, bs)
            Threads.@threads for bi in 1:nblocks
                i1 = (bi - 1) * bs + 1
                i2 = min(n_edges, bi * bs)
                for n in i1:i2
                    rwg_n = rwgs[n]
                    Z[m, n] = compute_efie_matrix_element(rwg_m, rwg_n, mesh, k)
                end
            end

            if p !== nothing
                ProgressMeter.next!(p)
            end
        end
    else
        for m in 1:n_edges
            rwg_m = rwgs[m]
            for n in 1:n_edges
                rwg_n = rwgs[n]
                Z[m, n] = compute_efie_matrix_element(rwg_m, rwg_n, mesh, k; caches=caches)
            end
            if p !== nothing
                ProgressMeter.next!(p)
            end
        end
    end

    return Z
end

function compute_efie_matrix_element(rwg_m::RWGFunction, rwg_n::RWGFunction, mesh::Mesh3D, k::Float64; caches=nothing)
    η = sqrt(μ0 / ε0)

    tri_m_plus_idx = rwg_m.triangle_plus
    tri_m_minus_idx = rwg_m.triangle_minus
    tri_n_plus_idx = rwg_n.triangle_plus
    tri_n_minus_idx = rwg_n.triangle_minus

    tri_m_plus = mesh.triangles[tri_m_plus_idx]
    tri_m_minus = mesh.triangles[tri_m_minus_idx]
    tri_n_plus = mesh.triangles[tri_n_plus_idx]
    tri_n_minus = mesh.triangles[tri_n_minus_idx]

    result = complex(0.0)

    triangle_pairs = [
        (tri_m_plus, tri_n_plus, 1.0, 1.0),
        (tri_m_plus, tri_n_minus, 1.0, -1.0),
        (tri_m_minus, tri_n_plus, -1.0, 1.0),
        (tri_m_minus, tri_n_minus, -1.0, -1.0)
    ]

    for (tri_obs, tri_src, sign_m, sign_n) in triangle_pairs
        # Determine triangle indices for cache access
        tri_obs_idx = tri_obs === tri_m_plus ? tri_m_plus_idx : (tri_obs === tri_m_minus ? tri_m_minus_idx : (tri_obs === tri_n_plus ? tri_n_plus_idx : tri_n_minus_idx))
        tri_src_idx = tri_src === tri_m_plus ? tri_m_plus_idx : (tri_src === tri_m_minus ? tri_m_minus_idx : (tri_src === tri_n_plus ? tri_n_plus_idx : tri_n_minus_idx))

        if triangles_overlap(tri_obs, tri_src)
            element = compute_singular_efie_element(tri_obs_idx, tri_src_idx, tri_obs, tri_src, rwg_m, rwg_n, mesh, k; caches=caches)
        elseif triangles_are_adjacent(tri_obs, tri_src)
            element = compute_near_singular_efie_element(tri_obs_idx, tri_src_idx, tri_obs, tri_src, rwg_m, rwg_n, mesh, k; caches=caches)
        else
            element = compute_regular_efie_element(tri_obs_idx, tri_src_idx, tri_obs, tri_src, rwg_m, rwg_n, mesh, k; caches=caches)
        end

        result += sign_m * sign_n * element
    end

    return 1im * η * k * result
end

function compute_regular_efie_element(tri_obs_idx::Int, tri_src_idx::Int, tri_obs::Triangle, tri_src::Triangle,
    rwg_m::RWGFunction, rwg_n::RWGFunction, mesh::Mesh3D, k::Float64; caches=nothing)

    quad_order = 3
    if caches !== nothing && haskey(caches, :tri_cart_pts) && haskey(caches, :tri_rwg_vals)
        # Use precomputed points and RWG values when available
        tri_cart = caches[:tri_cart_pts][quad_order][tri_src_idx]
        obs_cart = caches[:tri_cart_pts][quad_order][tri_obs_idx]
        tri_rwg = caches[:tri_rwg_vals][quad_order]

        result = complex(0.0)
        jacobian = tri_src.area * tri_obs.area
        weights = caches[:gauss_rules][quad_order][2]

        # Hoist constant/divergence evaluations and cache lookups out of inner loops
        div_m = evaluate_rwg_divergence(rwg_m)
        div_n = evaluate_rwg_divergence(rwg_n)

        src_rwg_dict = tri_rwg[tri_src_idx]
        obs_rwg_dict = tri_rwg[tri_obs_idx]

        src_vals = haskey(src_rwg_dict, rwg_n.edge_index) ? src_rwg_dict[rwg_n.edge_index] : nothing
        obs_vals = haskey(obs_rwg_dict, rwg_m.edge_index) ? obs_rwg_dict[rwg_m.edge_index] : nothing

        # If either side is not cached, precompute the missing side once
        if src_vals === nothing
            src_vals = Vector{SVector{3,Float64}}(undef, length(tri_cart))
            @inbounds for i in eachindex(tri_cart)
                src_vals[i] = evaluate_rwg(rwg_n, tri_cart[i], mesh)
            end
        end
        if obs_vals === nothing
            obs_vals = Vector{SVector{3,Float64}}(undef, length(obs_cart))
            @inbounds for j in eachindex(obs_cart)
                obs_vals[j] = evaluate_rwg(rwg_m, obs_cart[j], mesh)
            end
        end

        # Precompute k^2 to avoid repeated division
        k_squared = k^2
        div_product = (div_m * div_n) / k_squared
        
        @inbounds for i in eachindex(tri_cart)
            f_n = src_vals[i]
            wi = weights[i]
            r_src = tri_cart[i]
            for j in eachindex(obs_cart)
                f_m = obs_vals[j]
                wj = weights[j]
                r_obs = obs_cart[j]
                
                # Use fast inlined Green's function
                g = green_3d_fast(r_obs, r_src, k)
                
                # Use StaticArrays optimized dot product
                dot_product = dot(f_m, f_n)
                
                # Combine terms efficiently
                weight_jac = wi * wj * jacobian
                vector_contrib = dot_product * g
                scalar_contrib = (div_product * g)
                
                result += weight_jac * (vector_contrib + scalar_contrib)
            end
        end
        return result
    else
        function integrand(r_obs, r_src)
            f_m = evaluate_rwg(rwg_m, r_obs, mesh)
            f_n = evaluate_rwg(rwg_n, r_src, mesh)
            g = green_3d(r_obs, r_src, k)

            vector_term = dot(f_m, f_n) * g

            div_m = evaluate_rwg_divergence(rwg_m)
            div_n = evaluate_rwg_divergence(rwg_n)
            scalar_term = (div_m * div_n * g) / (k^2)

            return vector_term + scalar_term
        end

        return integrate_regular(integrand, tri_src, tri_obs, quad_order)
    end
end

function compute_singular_efie_element(tri_obs_idx::Int, tri_src_idx::Int, tri_obs::Triangle, tri_src::Triangle,
    rwg_m::RWGFunction, rwg_n::RWGFunction, mesh::Mesh3D, k::Float64; caches=nothing)

    # Use order 7 for singular extraction and 5 for the regular remainder
    # If caches available, reuse precomputed cartesian points and RWG vals.
    if caches !== nothing && haskey(caches, :tri_cart_pts)
        # singular part uses barycentric points on obs triangle
        ord_sing = 7
        ord_reg = 5
        obs_pts = caches[:tri_cart_pts][ord_sing][tri_obs_idx]
        reg_src_pts = caches[:tri_cart_pts][ord_reg][tri_src_idx]
        reg_obs_pts = caches[:tri_cart_pts][ord_reg][tri_obs_idx]

        singular_part = complex(0.0)
        # Precompute divergences and cached RWG arrays (or compute once)
        div_m = evaluate_rwg_divergence(rwg_m)
        div_n = evaluate_rwg_divergence(rwg_n)

        sing_rwg_dict = caches[:tri_rwg_vals][ord_sing][tri_obs_idx]
        sing_m_vals = haskey(sing_rwg_dict, rwg_m.edge_index) ? sing_rwg_dict[rwg_m.edge_index] : nothing
        sing_n_vals = haskey(sing_rwg_dict, rwg_n.edge_index) ? sing_rwg_dict[rwg_n.edge_index] : nothing

        if sing_m_vals === nothing || sing_n_vals === nothing
            sing_m_vals = Vector{SVector{3,Float64}}(undef, length(obs_pts))
            sing_n_vals = Vector{SVector{3,Float64}}(undef, length(obs_pts))
            @inbounds for i in eachindex(obs_pts)
                sing_m_vals[i] = evaluate_rwg(rwg_m, obs_pts[i], mesh)
                sing_n_vals[i] = evaluate_rwg(rwg_n, obs_pts[i], mesh)
            end
        end

        # Precompute constant factor
        area_over_4pi = tri_obs.area / FOUR_PI
        weights_sing = caches[:gauss_rules][ord_sing][2]
        
        for i in eachindex(obs_pts)
            # Manual dot product for performance
            dot_product = sing_m_vals[i][1] * sing_n_vals[i][1] + sing_m_vals[i][2] * sing_n_vals[i][2] + sing_m_vals[i][3] * sing_n_vals[i][3]
            singular_part += weights_sing[i] * area_over_4pi * dot_product
        end

        regular_part = complex(0.0)
        jacobian = tri_src.area * tri_obs.area
        weights = caches[:gauss_rules][ord_reg][2]
        tri_rwg_reg = caches[:tri_rwg_vals][ord_reg]

        src_dict_reg = tri_rwg_reg[tri_src_idx]
        obs_dict_reg = tri_rwg_reg[tri_obs_idx]

        src_vals_reg = haskey(src_dict_reg, rwg_n.edge_index) ? src_dict_reg[rwg_n.edge_index] : nothing
        obs_vals_reg = haskey(obs_dict_reg, rwg_m.edge_index) ? obs_dict_reg[rwg_m.edge_index] : nothing

        if src_vals_reg === nothing
            src_vals_reg = Vector{SVector{3,Float64}}(undef, length(reg_src_pts))
            @inbounds for j in eachindex(reg_src_pts)
                src_vals_reg[j] = evaluate_rwg(rwg_n, reg_src_pts[j], mesh)
            end
        end
        if obs_vals_reg === nothing
            obs_vals_reg = Vector{SVector{3,Float64}}(undef, length(reg_obs_pts))
            @inbounds for i in eachindex(reg_obs_pts)
                obs_vals_reg[i] = evaluate_rwg(rwg_m, reg_obs_pts[i], mesh)
            end
        end

        @inbounds for i in eachindex(reg_obs_pts)
            f_m = obs_vals_reg[i]
            wi = weights[i]
            r_obs = reg_obs_pts[i]
            for j in eachindex(reg_src_pts)
                f_n = src_vals_reg[j]
                wj = weights[j]
                r_src = reg_src_pts[j]
                g_parts = green_3d_singular_extraction(r_obs, r_src, k)
                vector_term = dot(f_m, f_n) * g_parts.regular
                scalar_term = (div_m * div_n * g_parts.regular) / (k^2)
                regular_part += wi * wj * jacobian * (vector_term + scalar_term)
            end
        end

        return singular_part + regular_part
    else
        function integrand_singular(r)
            f_m = evaluate_rwg(rwg_m, r, mesh)
            f_n = evaluate_rwg(rwg_n, r, mesh)
            return dot(f_m, f_n) / FOUR_PI
        end

        function integrand_regular(r_obs, r_src)
            f_m = evaluate_rwg(rwg_m, r_obs, mesh)
            f_n = evaluate_rwg(rwg_n, r_src, mesh)
            g_parts = green_3d_singular_extraction(r_obs, r_src, k)

            vector_term = dot(f_m, f_n) * g_parts.regular

            div_m = evaluate_rwg_divergence(rwg_m)
            div_n = evaluate_rwg_divergence(rwg_n)
            scalar_term = (div_m * div_n * g_parts.regular) / (k^2)

            return vector_term + scalar_term
        end

        singular_part = integrate_singular(integrand_singular, tri_obs, 7)
        regular_part = integrate_regular(integrand_regular, tri_src, tri_obs, 5)

        return singular_part + regular_part
    end
end

function compute_near_singular_efie_element(tri_obs_idx::Int, tri_src_idx::Int, tri_obs::Triangle, tri_src::Triangle,
    rwg_m::RWGFunction, rwg_n::RWGFunction, mesh::Mesh3D, k::Float64; caches=nothing)

    quad_order = 9
    if caches !== nothing && haskey(caches, :tri_cart_pts)
        # Try high-order product Gauss using precomputed points
        points_obs = caches[:tri_cart_pts][quad_order][tri_obs_idx]
        points_src = caches[:tri_cart_pts][quad_order][tri_src_idx]
        weights = caches[:gauss_rules][quad_order][2]
        tri_rwg = caches[:tri_rwg_vals][quad_order]
        jacobian = tri_src.area * tri_obs.area

        result = complex(0.0)

        # hoist divergences
        div_m = evaluate_rwg_divergence(rwg_m)
        div_n = evaluate_rwg_divergence(rwg_n)

        src_dict = tri_rwg[tri_src_idx]
        obs_dict = tri_rwg[tri_obs_idx]

        src_vals = haskey(src_dict, rwg_n.edge_index) ? src_dict[rwg_n.edge_index] : nothing
        obs_vals = haskey(obs_dict, rwg_m.edge_index) ? obs_dict[rwg_m.edge_index] : nothing

        # Precompute missing RWG evaluations
        if src_vals === nothing
            src_vals = Vector{SVector{3,Float64}}(undef, length(points_src))
            @inbounds for i in eachindex(points_src)
                src_vals[i] = evaluate_rwg(rwg_n, points_src[i], mesh)
            end
        end
        if obs_vals === nothing
            obs_vals = Vector{SVector{3,Float64}}(undef, length(points_obs))
            @inbounds for j in eachindex(points_obs)
                obs_vals[j] = evaluate_rwg(rwg_m, points_obs[j], mesh)
            end
        end

        # Precompute k^2 to avoid repeated division
        k_squared = k^2
        div_product = (div_m * div_n) / k_squared
        
        @inbounds for i in eachindex(points_src)
            f_n = src_vals[i]
            wi = weights[i]
            r_src = points_src[i]
            for j in eachindex(points_obs)
                f_m = obs_vals[j]
                wj = weights[j]
                r_obs = points_obs[j]
                
                # Use fast inlined Green's function
                g = green_3d_fast(r_obs, r_src, k)
                
                # Use StaticArrays optimized dot product
                dot_product = dot(f_m, f_n)
                
                # Combine terms efficiently
                weight_jac = wi * wj * jacobian
                vector_contrib = dot_product * g
                scalar_contrib = (div_product * g)
                
                result += weight_jac * (vector_contrib + scalar_contrib)
            end
        end
        return result
    else
        function integrand(r_obs, r_src)
            f_m = evaluate_rwg(rwg_m, r_obs, mesh)
            f_n = evaluate_rwg(rwg_n, r_src, mesh)
            g = green_3d(r_obs, r_src, k)

            vector_term = dot(f_m, f_n) * g

            div_m = evaluate_rwg_divergence(rwg_m)
            div_n = evaluate_rwg_divergence(rwg_n)
            scalar_term = (div_m * div_n * g) / (k^2)

            return vector_term + scalar_term
        end

        return integrate_near_singular(integrand, tri_src, tri_obs, 1e-8, 15)
    end
end

function triangles_overlap(tri1::Triangle, tri2::Triangle)
    return tri1 === tri2
end

function triangles_are_adjacent(tri1::Triangle, tri2::Triangle)
    shared_vertices = 0
    for v1 in tri1.vertices
        for v2 in tri2.vertices
            if norm(v1 - v2) < 1e-12
                shared_vertices += 1
            end
        end
    end
    return shared_vertices >= 2
end

function compute_excitation_vector(mesh::Mesh3D, incident_field_func, frequency::Float64)
    n_edges = mesh.num_edges
    b = zeros(ComplexF64, n_edges)
    # Use cached RWG basis stored on the mesh
    rwgs = get_rwgs(mesh)

    for m in 1:n_edges
        rwg_m = rwgs[m]
        b[m] = compute_excitation_element(rwg_m, mesh, incident_field_func)
    end

    return b
end

function compute_excitation_element(rwg_m::RWGFunction, mesh::Mesh3D, incident_field_func)
    tri_plus = mesh.triangles[rwg_m.triangle_plus]
    tri_minus = mesh.triangles[rwg_m.triangle_minus]

    function integrand_plus(r)
        f_m = evaluate_rwg(rwg_m, r, mesh)
        e_inc = incident_field_func(r)
        return dot(f_m, e_inc)
    end

    function integrand_minus(r)
        f_m = evaluate_rwg(rwg_m, r, mesh)
        e_inc = incident_field_func(r)
        return dot(f_m, e_inc)
    end

    integral_plus = integrate_singular(integrand_plus, tri_plus, 3)
    integral_minus = integrate_singular(integrand_minus, tri_minus, 3)

    return integral_plus + integral_minus
end

end
