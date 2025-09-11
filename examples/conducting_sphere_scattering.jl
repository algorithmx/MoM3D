"""
Conducting Sphere Scattering Example

This example demonstrates the Method of Moments solution for electromagnetic scattering
from a perfectly conducting sphere, based on the theoretical framework described in
Gibson's "The Method of Moments in Electromagnetics" Sections 6.6.2 and 7.7.2.

Theoretical Background:
- The conducting sphere is a canonical scattering problem with known analytical solution (Mie theory)
- Gibson Section 6.6.2 describes the EFIE formulation for conducting objects
- Gibson Section 7.7.2 provides validation examples with 2-meter diameter spheres
- Expected accuracy: < 0.1 dB error compared to Mie theory for well-conditioned meshes

Key Requirements (Gibson Section 7.6.1):
- Triangle aspect ratios should be reasonable (< 5.0 preferred)
- Minimum angles should be > 10° to avoid poorly conditioned matrices
- Mesh should be watertight for MFIE applications (Gibson Section 7.6.2)

Applications:
- Radar cross section prediction
- Antenna scattering analysis
- Electromagnetic compatibility studies
- Validation of numerical methods
"""

using MoM3D
using LinearAlgebra
using StaticArrays
using Printf

function create_sphere_mesh(radius::Float64, n_theta::Int, n_phi::Int)
    """
    Create a triangular mesh for a sphere using spherical coordinates.
    
    Parameters:
    - radius: Sphere radius in meters
    - n_theta: Number of divisions in polar direction (θ)
    - n_phi: Number of divisions in azimuthal direction (φ)
    
    Returns:
    - vertices: Array of 3D vertex coordinates
    - triangle_indices: Array of triangle connectivity
    
    Note: Mesh quality depends on n_theta and n_phi values.
    Gibson Section 7.6.1 recommends avoiding thin triangles near poles.
    """
    vertices = SVector{3, Float64}[]
    triangle_indices = SVector{3, Int}[]
    
    # Generate vertices using spherical coordinates
    for i in 0:n_theta
        theta = π * i / n_theta  # Polar angle: 0 to π
        for j in 0:n_phi-1
            phi = 2π * j / n_phi  # Azimuthal angle: 0 to 2π
            
            # Convert to Cartesian coordinates
            x = radius * sin(theta) * cos(phi)
            y = radius * sin(theta) * sin(phi)
            z = radius * cos(theta)
            
            push!(vertices, SVector(x, y, z))
        end
    end
    
    # Generate triangular connectivity
    for i in 0:n_theta-1
        for j in 0:n_phi-1
            v1 = i * n_phi + j + 1
            v2 = i * n_phi + ((j + 1) % n_phi) + 1
            v3 = (i + 1) * n_phi + j + 1
            v4 = (i + 1) * n_phi + ((j + 1) % n_phi) + 1
            
            # Handle pole regions to avoid degenerate triangles
            if i == 0  # North pole
                push!(triangle_indices, SVector(v1, v3, v4))
            elseif i == n_theta - 1  # South pole
                push!(triangle_indices, SVector(v1, v2, v3))
            else  # Regular region
                push!(triangle_indices, SVector(v1, v2, v4))
                push!(triangle_indices, SVector(v1, v4, v3))
            end
        end
    end
    
    return vertices, triangle_indices
end

function analyze_sphere_mesh_quality(mesh::Mesh3D, radius::Float64)
    """
    Analyze mesh quality according to Gibson Section 7.6.1 requirements.
    
    Quality Metrics:
    - Aspect ratios: Should be < 5.0 for good conditioning
    - Minimum angles: Should be > 10° to avoid numerical issues
    - Watertight: Required for MFIE formulation (Gibson Section 7.6.2)
    - T-junctions: Should be absent for proper current flow
    """
    println("=== Sphere Mesh Quality Analysis ===")
    println("Sphere radius: $(radius) meters")
    println("Number of vertices: $(length(mesh.vertices))")
    println("Number of triangles: $(length(mesh.triangles))")
    println("Number of edges: $(mesh.num_edges)")
    
    # Mesh topology
    println("\nMesh Topology:")
    println("  Watertight: $(mesh.is_watertight)")
    println("  Boundary edges: $(length(mesh.boundary_edges))")
    println("  Interior edges: $(length(mesh.interior_edges))")
    
    # Quality metrics (Gibson Section 7.6.1)
    quality = mesh.mesh_quality
    println("\nQuality Metrics (Gibson Section 7.6.1):")
    @printf "  Min aspect ratio: %.3f\n" quality.min_aspect_ratio
    @printf "  Max aspect ratio: %.3f\n" quality.max_aspect_ratio
    @printf "  Mean aspect ratio: %.3f\n" quality.mean_aspect_ratio
    @printf "  Min angle: %.1f°\n" quality.min_angle_deg
    @printf "  Max angle: %.1f°\n" quality.max_angle_deg
    @printf "  Degenerate triangles: %d\n" quality.degenerate_triangles
    @printf "  Thin triangles (AR > 10): %d\n" quality.thin_triangles
    
    # Quality assessment
    println("\nQuality Assessment:")
    if quality.max_aspect_ratio < 5.0
        println("  ✓ Aspect ratios are good (< 5.0)")
    else
        println("  ⚠ Some triangles have poor aspect ratios (> 5.0)")
    end
    
    if quality.min_angle_deg > 10.0
        println("  ✓ Minimum angles are acceptable (> 10°)")
    else
        println("  ⚠ Some triangles have very small angles (< 10°)")
    end
    
    if mesh.is_watertight
        println("  ✓ Mesh is watertight (suitable for MFIE)")
    else
        println("  ⚠ Mesh is not watertight (MFIE may be inaccurate)")
    end
    
    # T-junction detection (Gibson Section 7.6.2)
    t_junctions = detect_t_junctions(mesh)
    if isempty(t_junctions)
        println("  ✓ No T-junctions detected")
    else
        println("  ⚠ $(length(t_junctions)) T-junctions detected")
    end
    
    return quality
end

function estimate_electrical_size(radius::Float64, frequency::Float64)
    """
    Estimate electrical size of sphere for scattering analysis.
    
    Parameters:
    - radius: Sphere radius in meters
    - frequency: Frequency in Hz
    
    Returns:
    - ka: Electrical size parameter (k * radius)
    - wavelength: Free-space wavelength in meters
    
    Note: Gibson examples typically use ka values from 0.1 to 10.0
    """
    c = 2.998e8  # Speed of light in m/s
    wavelength = c / frequency
    k = 2π / wavelength  # Wavenumber
    ka = k * radius
    
    return ka, wavelength
end

function demonstrate_sphere_scattering()
    """
    Main demonstration of conducting sphere scattering analysis.
    
    This example follows the validation approach described in Gibson Section 7.7.2:
    1. Create sphere mesh with good quality triangles
    2. Analyze mesh quality according to Gibson requirements
    3. Estimate electrical size for different frequencies
    4. Discuss expected MoM solution accuracy
    """
    println("Conducting Sphere Scattering Analysis")
    println("Based on Gibson 'Method of Moments in Electromagnetics'")
    println("Sections 6.6.2 (Theory) and 7.7.2 (Validation)")
    println(repeat("=", 60))
    
    # Sphere parameters (Gibson Section 7.7.2 uses 2-meter diameter)
    radius = 1.0  # meters (2-meter diameter)
    
    # Create mesh with moderate resolution
    println("Creating sphere mesh...")
    n_theta = 16  # Polar divisions
    n_phi = 32    # Azimuthal divisions
    
    vertices, triangle_indices = create_sphere_mesh(radius, n_theta, n_phi)
    mesh = Mesh3D(vertices, triangle_indices)
    
    # Analyze mesh quality
    quality = analyze_sphere_mesh_quality(mesh, radius)
    
    # Electrical size analysis for different frequencies
    println("\n=== Electrical Size Analysis ===")
    frequencies = [10e6, 50e6, 100e6, 300e6, 500e6]  # 10 MHz to 500 MHz
    
    for freq in frequencies
        ka, wavelength = estimate_electrical_size(radius, freq)
        @printf "Frequency: %6.0f MHz, λ = %5.2f m, ka = %5.2f\n" freq/1e6 wavelength ka
    end
    
    # Expected accuracy discussion
    println("\n=== Expected MoM Solution Accuracy ===")
    println("Gibson Section 7.7.2 validation results:")
    println("  • RCS accuracy: < 0.1 dB error vs. Mie theory")
    println("  • Frequency range: 10-500 MHz (ka ≈ 0.02 to 10.5)")
    println("  • Mesh requirements: Good aspect ratios, watertight surface")
    println("  • EFIE formulation recommended over MFIE (Section 7.4.4)")
    
    # Mesh refinement recommendations
    println("\n=== Mesh Refinement Recommendations ===")
    if quality.max_aspect_ratio > 3.0
        println("  • Consider increasing n_phi to improve aspect ratios")
    end
    if quality.min_angle_deg < 15.0
        println("  • Consider more uniform mesh distribution")
    end
    if !mesh.is_watertight
        println("  • Repair mesh connectivity for MFIE applications")
    end
    
    println("\n=== Next Steps ===")
    println("1. Implement EFIE matrix assembly (Gibson Chapter 6)")
    println("2. Apply plane wave excitation")
    println("3. Solve linear system for surface currents")
    println("4. Compute radar cross section (Gibson Section 7.7.2)")
    println("5. Compare with Mie theory for validation")
    
    return mesh, quality
end

# Run demonstration if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    mesh, quality = demonstrate_sphere_scattering()
end
