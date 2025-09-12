#!/usr/bin/env julia

# Test script to verify that optimized functions produce identical results
using LinearAlgebra
using StaticArrays

# Original green_3d function
function green_3d_original(r_obs::SVector{3, Float64}, r_src::SVector{3, Float64}, k::Float64)
    R = norm(r_obs - r_src)
    if R < 1e-12
        return complex(0.0)
    end
    return exp(-1im * k * R) / (4π * R)
end

# Optimized green_3d function
const FOUR_PI = 4π
function green_3d_optimized(r_obs::SVector{3, Float64}, r_src::SVector{3, Float64}, k::Float64)
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

# Test cases
test_points = [
    (SVector(1.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)),
    (SVector(0.5, 0.5, 0.5), SVector(0.1, 0.1, 0.1)),
    (SVector(2.0, 1.0, 3.0), SVector(1.0, 2.0, 0.5)),
    (SVector(0.0, 0.0, 0.0), SVector(1e-13, 0.0, 0.0)),  # Near-singular case
    (SVector(10.0, 5.0, 2.0), SVector(3.0, 1.0, 8.0))   # Far field case
]

k_values = [0.1, 1.0, 10.0, 100.0]

println("Testing Green's function optimizations...")
println("=" ^ 50)

max_error = 0.0
total_tests = 0

for k in k_values
    for (r_obs, r_src) in test_points
        original = green_3d_original(r_obs, r_src, k)
        optimized = green_3d_optimized(r_obs, r_src, k)
        
        error = abs(original - optimized)
        rel_error = abs(original) > 1e-12 ? error / abs(original) : error
        
        global total_tests += 1
        global max_error = max(max_error, rel_error)
        
        if rel_error > 1e-12
            println("MISMATCH: k=$k, r_obs=$r_obs, r_src=$r_src")
            println("  Original: $original")
            println("  Optimized: $optimized")
            println("  Relative error: $rel_error")
        end
    end
end

println("Completed $total_tests tests")
println("Maximum relative error: $max_error")

if max_error < 1e-12
    println("✓ All tests passed! Optimized functions are mathematically equivalent.")
else
    println("✗ Some tests failed! Check the implementation.")
end

# Test manual dot product optimization
println("\nTesting manual dot product optimization...")
test_vectors = [
    (SVector(1.0, 2.0, 3.0), SVector(4.0, 5.0, 6.0)),
    (SVector(0.1, 0.2, 0.3), SVector(0.4, 0.5, 0.6)),
    (SVector(-1.0, 2.0, -3.0), SVector(4.0, -5.0, 6.0))
]

dot_max_error = 0.0
for (v1, v2) in test_vectors
    original_dot = dot(v1, v2)
    manual_dot = v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3]
    
    error = abs(original_dot - manual_dot)
    global dot_max_error = max(dot_max_error, error)
    
    if error > 1e-15
        println("DOT MISMATCH: v1=$v1, v2=$v2")
        println("  Original: $original_dot")
        println("  Manual: $manual_dot")
        println("  Error: $error")
    end
end

println("Maximum dot product error: $dot_max_error")
if dot_max_error < 1e-15
    println("✓ Manual dot product optimization is correct.")
else
    println("✗ Manual dot product has errors.")
end
