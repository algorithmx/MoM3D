module Integration

using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using QuadGK
using ..Geometry

export gauss_triangle, integrate_singular, integrate_regular

function gauss_triangle(n::Int)
    if n == 1
        points = [SVector(1/3, 1/3)]
        weights = [1.0]
    elseif n == 3
        points = [SVector(1/6, 1/6), SVector(2/3, 1/6), SVector(1/6, 2/3)]
        weights = [1/3, 1/3, 1/3]
    elseif n == 4
        points = [SVector(1/3, 1/3), SVector(0.6, 0.2), SVector(0.2, 0.6), SVector(0.2, 0.2)]
        weights = [-27/48, 25/48, 25/48, 25/48]
    else
        error("Gaussian quadrature order $n not implemented")
    end
    
    return points, weights
end

function barycentric_to_cartesian(ξ::SVector{2, Float64}, triangle::Triangle)
    v1, v2, v3 = triangle.vertices
    ξ1, ξ2 = ξ
    ξ3 = 1.0 - ξ1 - ξ2
    return ξ1 * v1 + ξ2 * v2 + ξ3 * v3
end

function integrate_regular(integrand_func, tri_src::Triangle, tri_obs::Triangle, quad_order::Int=3)
    points_src, weights_src = gauss_triangle(quad_order)
    points_obs, weights_obs = gauss_triangle(quad_order)
    
    result = complex(0.0)
    
    for (i, (ξ_src, w_src)) in enumerate(zip(points_src, weights_src))
        r_src = barycentric_to_cartesian(ξ_src, tri_src)
        
        for (j, (ξ_obs, w_obs)) in enumerate(zip(points_obs, weights_obs))
            r_obs = barycentric_to_cartesian(ξ_obs, tri_obs)
            
            jacobian = tri_src.area * tri_obs.area
            integrand_value = integrand_func(r_obs, r_src)
            
            result += w_src * w_obs * jacobian * integrand_value
        end
    end
    
    return result
end

function integrate_singular(integrand_func, triangle::Triangle, quad_order::Int=7)
    points, weights = gauss_triangle(quad_order)
    
    result = complex(0.0)
    
    for (ξ, w) in zip(points, weights)
        r = barycentric_to_cartesian(ξ, triangle)
        jacobian = triangle.area
        integrand_value = integrand_func(r, r)
        result += w * jacobian * integrand_value
    end
    
    return result
end

function integrate_near_singular(integrand_func, tri_src::Triangle, tri_obs::Triangle, 
                                tolerance::Float64=1e-6, max_subdivisions::Int=10)
    function adaptive_integrand(ξ_flat)
        n_points = length(ξ_flat) ÷ 4
        result = complex(0.0)
        
        for i in 1:n_points
            idx = 4*(i-1)
            ξ_src = SVector(ξ_flat[idx+1], ξ_flat[idx+2])
            ξ_obs = SVector(ξ_flat[idx+3], ξ_flat[idx+4])
            
            r_src = barycentric_to_cartesian(ξ_src, tri_src)
            r_obs = barycentric_to_cartesian(ξ_obs, tri_obs)
            
            result += integrand_func(r_obs, r_src)
        end
        
        return result * tri_src.area * tri_obs.area
    end
    
    integral, _ = quadgk(adaptive_integrand, [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0], 
                        rtol=tolerance, maxevals=max_subdivisions*1000)
    
    return integral
end

end
