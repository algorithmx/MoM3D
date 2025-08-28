module PostProcess

using LinearAlgebra
using StaticArrays
using ..Geometry
using ..BasisFunctions
using ..GreenFunctions

export compute_near_field, compute_far_field, compute_rcs, compute_current_density

function compute_near_field(current_coeffs::Vector{ComplexF64}, mesh::Mesh3D, 
                           observation_points::Vector{SVector{3, Float64}}, frequency::Float64)
    k = wavenumber(frequency)
    η = sqrt(μ0 / ε0)
    
    E_field = [SVector(complex(0.0), complex(0.0), complex(0.0)) for _ in observation_points]
    H_field = [SVector(complex(0.0), complex(0.0), complex(0.0)) for _ in observation_points]
    
    for (obs_idx, r_obs) in enumerate(observation_points)
        for edge_idx in 1:mesh.num_edges
            rwg = RWGFunction(mesh.edges[edge_idx], mesh)
            I_n = current_coeffs[edge_idx]
            
            E_contribution = compute_electric_field_contribution(rwg, I_n, r_obs, mesh, k, η)
            H_contribution = compute_magnetic_field_contribution(rwg, I_n, r_obs, mesh, k)
            
            E_field[obs_idx] += E_contribution
            H_field[obs_idx] += H_contribution
        end
    end
    
    return E_field, H_field
end

function compute_electric_field_contribution(rwg::RWGFunction, current_coeff::ComplexF64,
                                           r_obs::SVector{3, Float64}, mesh::Mesh3D, k::Float64, η::Float64)
    
    tri_plus = mesh.triangles[rwg.triangle_plus]
    tri_minus = mesh.triangles[rwg.triangle_minus]
    
    function integrand_plus(r_src)
        f_n = evaluate_rwg(rwg, r_src, mesh)
        g = green_3d(r_obs, r_src, k)
        grad_g = green_3d_gradient(r_obs, r_src, k)
        
        vector_term = 1im * η * k * f_n * g
        scalar_term = (1im * η / k) * evaluate_rwg_divergence(rwg) * grad_g
        
        return vector_term + scalar_term
    end
    
    function integrand_minus(r_src)
        f_n = evaluate_rwg(rwg, r_src, mesh)
        g = green_3d(r_obs, r_src, k)
        grad_g = green_3d_gradient(r_obs, r_src, k)
        
        vector_term = -1im * η * k * f_n * g
        scalar_term = -(1im * η / k) * evaluate_rwg_divergence(rwg) * grad_g
        
        return vector_term + scalar_term
    end
    
    integral_plus = integrate_regular(integrand_plus, tri_plus, Triangle(r_obs, r_obs, r_obs), 3)
    integral_minus = integrate_regular(integrand_minus, tri_minus, Triangle(r_obs, r_obs, r_obs), 3)
    
    return current_coeff * (integral_plus + integral_minus)
end

function compute_magnetic_field_contribution(rwg::RWGFunction, current_coeff::ComplexF64,
                                           r_obs::SVector{3, Float64}, mesh::Mesh3D, k::Float64)
    
    tri_plus = mesh.triangles[rwg.triangle_plus]
    tri_minus = mesh.triangles[rwg.triangle_minus]
    
    function integrand_plus(r_src)
        f_n = evaluate_rwg(rwg, r_src, mesh)
        grad_g = green_3d_gradient(r_obs, r_src, k)
        return cross(f_n, grad_g)
    end
    
    function integrand_minus(r_src)
        f_n = evaluate_rwg(rwg, r_src, mesh)
        grad_g = green_3d_gradient(r_obs, r_src, k)
        return -cross(f_n, grad_g)
    end
    
    integral_plus = integrate_regular(integrand_plus, tri_plus, Triangle(r_obs, r_obs, r_obs), 3)
    integral_minus = integrate_regular(integrand_minus, tri_minus, Triangle(r_obs, r_obs, r_obs), 3)
    
    return current_coeff * (integral_plus + integral_minus)
end

function compute_far_field(current_coeffs::Vector{ComplexF64}, mesh::Mesh3D,
                          theta_angles::Vector{Float64}, phi_angles::Vector{Float64}, frequency::Float64)
    k = wavenumber(frequency)
    η = sqrt(μ0 / ε0)
    
    E_theta = zeros(ComplexF64, length(theta_angles), length(phi_angles))
    E_phi = zeros(ComplexF64, length(theta_angles), length(phi_angles))
    
    for (i, theta) in enumerate(theta_angles)
        for (j, phi) in enumerate(phi_angles)
            r_hat = SVector(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))
            theta_hat = SVector(cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta))
            phi_hat = SVector(-sin(phi), cos(phi), 0.0)
            
            E_far = SVector(complex(0.0), complex(0.0), complex(0.0))
            
            for edge_idx in 1:mesh.num_edges
                rwg = RWGFunction(mesh.edges[edge_idx], mesh)
                I_n = current_coeffs[edge_idx]
                
                E_contribution = compute_far_field_contribution(rwg, I_n, r_hat, mesh, k, η)
                E_far += E_contribution
            end
            
            E_theta[i, j] = dot(E_far, theta_hat)
            E_phi[i, j] = dot(E_far, phi_hat)
        end
    end
    
    return E_theta, E_phi
end

function compute_far_field_contribution(rwg::RWGFunction, current_coeff::ComplexF64,
                                       r_hat::SVector{3, Float64}, mesh::Mesh3D, k::Float64, η::Float64)
    
    tri_plus = mesh.triangles[rwg.triangle_plus]
    tri_minus = mesh.triangles[rwg.triangle_minus]
    
    function integrand_plus(r_src)
        f_n = evaluate_rwg(rwg, r_src, mesh)
        phase = exp(1im * k * dot(r_hat, r_src))
        return f_n * phase
    end
    
    function integrand_minus(r_src)
        f_n = evaluate_rwg(rwg, r_src, mesh)
        phase = exp(1im * k * dot(r_hat, r_src))
        return -f_n * phase
    end
    
    integral_plus = integrate_regular(integrand_plus, tri_plus, tri_plus, 3)
    integral_minus = integrate_regular(integrand_minus, tri_minus, tri_minus, 3)
    
    total_integral = integral_plus + integral_minus
    
    return current_coeff * (-1im * η * k / (4π)) * (total_integral - dot(total_integral, r_hat) * r_hat)
end

function compute_rcs(E_theta::Matrix{ComplexF64}, E_phi::Matrix{ComplexF64})
    return 4π * (abs2.(E_theta) + abs2.(E_phi))
end

function compute_current_density(current_coeffs::Vector{ComplexF64}, mesh::Mesh3D,
                                evaluation_points::Vector{SVector{3, Float64}})
    
    current_density = [SVector(complex(0.0), complex(0.0), complex(0.0)) for _ in evaluation_points]
    
    for (point_idx, r) in enumerate(evaluation_points)
        for edge_idx in 1:mesh.num_edges
            rwg = RWGFunction(mesh.edges[edge_idx], mesh)
            I_n = current_coeffs[edge_idx]
            
            f_n = evaluate_rwg(rwg, r, mesh)
            current_density[point_idx] += I_n * f_n
        end
    end
    
    return current_density
end

end
