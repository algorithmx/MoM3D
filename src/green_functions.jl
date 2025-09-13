module GreenFunctions

using LinearAlgebra
using StaticArrays
using SpecialFunctions

export c0, μ0, ε0, green_3d, green_3d_gradient, green_3d_singular_extraction, wavenumber, green_3d

const c0 = 2.99792458e8
const μ0 = 4π * 1e-7
const ε0 = 1.0 / (μ0 * c0^2)
const FOUR_PI = 4π
const GREEN_FUNCTION_CUTOFF = 1e-12
const GREEN_FUNCTION_CUTOFF2 = 1e-24

function wavenumber(frequency::Float64)
    return 2π * frequency / c0
end

@inline function green_3d(R_vec::SVector{3,Float64}, k::Float64)
    R_squared = dot(R_vec, R_vec)
    if R_squared < GREEN_FUNCTION_CUTOFF2
        return complex(0.0)
    end
    R = sqrt(R_squared)
    # Use cis for exp(-im*kr) = cos(kr) - im*sin(kr)
    return cis(-k * R) / (FOUR_PI * R)
end


@inline green_3d(r_obs::SVector{3, Float64}, r_src::SVector{3, Float64}, k::Float64) = green_3d(r_obs - r_src, k)

function green_3d_gradient(r_obs::SVector{3, Float64}, r_src::SVector{3, Float64}, k::Float64)
    R_vec = r_obs - r_src
    R_squared = dot(R_vec, R_vec)
    if R_squared < GREEN_FUNCTION_CUTOFF2
        return SVector(complex(0.0), complex(0.0), complex(0.0))
    end
    factor = (-1im * k - 1.0 / sqrt(R_squared)) * green_3d(R_vec, k) / R_squared
    return factor * R_vec
end

function green_3d_singular_extraction(r_obs::SVector{3, Float64}, r_src::SVector{3, Float64}, k::Float64)
    R = norm(r_obs - r_src)
    if R < GREEN_FUNCTION_CUTOFF
        return (singular=1.0/(FOUR_PI*R), regular=complex(0.0))
    end
    
    singular_part = 1.0 / (FOUR_PI * R)
    regular_part = (exp(-1im * k * R) - 1.0) / (FOUR_PI * R)
    
    return (singular=singular_part, regular=regular_part)
end

end

