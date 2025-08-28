module GreenFunctions

using LinearAlgebra
using StaticArrays
using SpecialFunctions

export green_3d, green_3d_gradient, wavenumber

const c0 = 2.99792458e8
const μ0 = 4π * 1e-7
const ε0 = 1.0 / (μ0 * c0^2)

function wavenumber(frequency::Float64)
    return 2π * frequency / c0
end

function green_3d(r_obs::SVector{3, Float64}, r_src::SVector{3, Float64}, k::Float64)
    R = norm(r_obs - r_src)
    if R < 1e-12
        return complex(0.0)
    end
    return exp(-1im * k * R) / (4π * R)
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
        return (singular=1.0/(4π*R), regular=complex(0.0))
    end
    
    singular_part = 1.0 / (4π * R)
    regular_part = (exp(-1im * k * R) - 1.0) / (4π * R)
    
    return (singular=singular_part, regular=regular_part)
end

end
