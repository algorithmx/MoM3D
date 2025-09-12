module GreenFunctions

using LinearAlgebra
using StaticArrays
using SpecialFunctions

export c0, μ0, ε0, green_3d, green_3d_gradient, green_3d_singular_extraction, wavenumber, green_3d_fast

const c0 = 2.99792458e8
const μ0 = 4π * 1e-7
const ε0 = 1.0 / (μ0 * c0^2)
const FOUR_PI = 4π

function wavenumber(frequency::Float64)
    return 2π * frequency / c0
end

function green_3d(r_obs::SVector{3, Float64}, r_src::SVector{3, Float64}, k::Float64)
    # Compute R^2 first to avoid sqrt when possible
    dx = r_obs[1] - r_src[1]
    dy = r_obs[2] - r_src[2] 
    dz = r_obs[3] - r_src[3]
    R_squared = dx*dx + dy*dy + dz*dz
    
    if R_squared < 1e-24  # 1e-12^2
        return complex(0.0)
    end
    
    R = sqrt(R_squared)
    kr = k * R
    
    # Use cis for exp(-im*kr) = cos(kr) - im*sin(kr)
    return cis(-kr) / (FOUR_PI * R)
end

# Fast inlined version for hot loop computations
@inline function green_3d_fast(r_obs::SVector{3, Float64}, r_src::SVector{3, Float64}, k::Float64)
    # Compute R^2 first to avoid sqrt when possible
    dx = r_obs[1] - r_src[1]
    dy = r_obs[2] - r_src[2] 
    dz = r_obs[3] - r_src[3]
    R_squared = dx*dx + dy*dy + dz*dz
    
    if R_squared < 1e-24  # 1e-12^2
        return complex(0.0)
    end
    
    R = sqrt(R_squared)
    kr = k * R
    
    # Use cis for exp(-im*kr) = cos(kr) - im*sin(kr)
    return cis(-kr) / (FOUR_PI * R)
end

function green_3d_gradient(r_obs::SVector{3, Float64}, r_src::SVector{3, Float64}, k::Float64)
    R_vec = r_obs - r_src
    R = norm(R_vec)
    
    if R < 1e-12
        return SVector(complex(0.0), complex(0.0), complex(0.0))
    end
    
    R_hat = R_vec / R
    g = green_3d(r_obs, r_src, k)
    factor = (-1im * k - 1.0 / R) * g / R
    
    return factor * R_hat
end

function green_3d_singular_extraction(r_obs::SVector{3, Float64}, r_src::SVector{3, Float64}, k::Float64)
    R = norm(r_obs - r_src)
    if R < 1e-12
        return (singular=1.0/(FOUR_PI*R), regular=complex(0.0))
    end
    
    singular_part = 1.0 / (FOUR_PI * R)
    regular_part = (exp(-1im * k * R) - 1.0) / (FOUR_PI * R)
    
    return (singular=singular_part, regular=regular_part)
end

end
