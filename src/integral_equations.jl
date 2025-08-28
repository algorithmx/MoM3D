module IntegralEquations

using LinearAlgebra
using StaticArrays
using ..Geometry
using ..BasisFunctions
using ..GreenFunctions
using ..Integration

export assemble_efie_matrix, assemble_mfie_matrix, assemble_cfie_matrix
export compute_excitation_vector

function assemble_efie_matrix(mesh::Mesh3D, frequency::Float64)
    k = wavenumber(frequency)
    n_edges = mesh.num_edges
    Z = zeros(ComplexF64, n_edges, n_edges)
    
    for m in 1:n_edges
        rwg_m = RWGFunction(mesh.edges[m], mesh)
        
        for n in 1:n_edges
            rwg_n = RWGFunction(mesh.edges[n], mesh)
            Z[m, n] = compute_efie_matrix_element(rwg_m, rwg_n, mesh, k)
        end
    end
    
    return Z
end

function compute_efie_matrix_element(rwg_m::RWGFunction, rwg_n::RWGFunction, mesh::Mesh3D, k::Float64)
    η = sqrt(μ0 / ε0)
    
    tri_m_plus = mesh.triangles[rwg_m.triangle_plus]
    tri_m_minus = mesh.triangles[rwg_m.triangle_minus]
    tri_n_plus = mesh.triangles[rwg_n.triangle_plus]
    tri_n_minus = mesh.triangles[rwg_n.triangle_minus]
    
    result = complex(0.0)
    
    triangle_pairs = [
        (tri_m_plus, tri_n_plus, 1.0, 1.0),
        (tri_m_plus, tri_n_minus, 1.0, -1.0),
        (tri_m_minus, tri_n_plus, -1.0, 1.0),
        (tri_m_minus, tri_n_minus, -1.0, -1.0)
    ]
    
    for (tri_obs, tri_src, sign_m, sign_n) in triangle_pairs
        if triangles_overlap(tri_obs, tri_src)
            element = compute_singular_efie_element(tri_obs, tri_src, rwg_m, rwg_n, k)
        elseif triangles_are_adjacent(tri_obs, tri_src)
            element = compute_near_singular_efie_element(tri_obs, tri_src, rwg_m, rwg_n, k)
        else
            element = compute_regular_efie_element(tri_obs, tri_src, rwg_m, rwg_n, k)
        end
        
        result += sign_m * sign_n * element
    end
    
    return 1im * η * k * result
end

function compute_regular_efie_element(tri_obs::Triangle, tri_src::Triangle, 
                                     rwg_m::RWGFunction, rwg_n::RWGFunction, k::Float64)
    
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
    
    return integrate_regular(integrand, tri_src, tri_obs, 3)
end

function compute_singular_efie_element(tri_obs::Triangle, tri_src::Triangle,
                                      rwg_m::RWGFunction, rwg_n::RWGFunction, k::Float64)
    
    function integrand_singular(r)
        f_m = evaluate_rwg(rwg_m, r, mesh)
        f_n = evaluate_rwg(rwg_n, r, mesh)
        return dot(f_m, f_n) / (4π)
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

function compute_near_singular_efie_element(tri_obs::Triangle, tri_src::Triangle,
                                           rwg_m::RWGFunction, rwg_n::RWGFunction, k::Float64)
    
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
    
    for m in 1:n_edges
        rwg_m = RWGFunction(mesh.edges[m], mesh)
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
