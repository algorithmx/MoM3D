"""
Radar Cross Section (RCS) Prediction Example

This example demonstrates the Method of Moments approach for predicting radar cross
section of conducting objects, based on Gibson's "The Method of Moments in
Electromagnetics" Chapter 7.

Theoretical Background:
- RCS is a measure of how detectable an object is by radar
- Gibson Section 7.7.2 provides RCS validation for conducting spheres
- Gibson Section 7.7.3 discusses EMCC benchmark targets for RCS validation
- MoM provides surface currents, from which far-field RCS is calculated

Key Concepts:
- Bistatic and monostatic RCS calculations
- Far-field approximation and radiation integrals
- Frequency dependence and resonance effects
- Polarization considerations (VV, HH, VH, HV)

Applications:
- Stealth technology development
- Radar system design and analysis
- Electromagnetic compatibility studies
- Target identification and classification
"""

using MoM3D
using LinearAlgebra
using StaticArrays
using Printf

function create_simple_target_meshes()
    """
    Create meshes for simple canonical RCS targets.
    
    These targets are commonly used for RCS validation:
    - Sphere: Known analytical solution (Mie theory)
    - Flat plate: Simple geometry with edge diffraction
    - Corner reflector: Multiple scattering example
    
    Returns dictionaries of vertices and triangle indices for each target.
    """
    targets = Dict()
    
    # 1. Conducting Sphere (Gibson Section 7.7.2)
    sphere_radius = 1.0  # meters
    n_theta, n_phi = 16, 32
    sphere_vertices, sphere_triangles = create_sphere_mesh(sphere_radius, n_theta, n_phi)
    targets["sphere"] = (vertices=sphere_vertices, triangles=sphere_triangles, 
                        description="Conducting sphere, radius = $sphere_radius m")
    
    # 2. Flat Plate
    plate_length, plate_width = 2.0, 1.0  # meters
    n_x, n_y = 20, 10
    plate_vertices, plate_triangles = create_plate_mesh(plate_length, plate_width, n_x, n_y)
    targets["plate"] = (vertices=plate_vertices, triangles=plate_triangles,
                       description="Flat plate, $(plate_length) × $(plate_width) m")
    
    # 3. Corner Reflector (simplified)
    corner_size = 1.0  # meters
    corner_vertices, corner_triangles = create_corner_reflector_mesh(corner_size)
    targets["corner"] = (vertices=corner_vertices, triangles=corner_triangles,
                        description="Corner reflector, size = $corner_size m")
    
    return targets
end

function create_sphere_mesh(radius::Float64, n_theta::Int, n_phi::Int)
    """Create spherical mesh (same as in conducting_sphere_scattering.jl)"""
    vertices = SVector{3, Float64}[]
    triangle_indices = SVector{3, Int}[]
    
    for i in 0:n_theta
        theta = π * i / n_theta
        for j in 0:n_phi-1
            phi = 2π * j / n_phi
            
            x = radius * sin(theta) * cos(phi)
            y = radius * sin(theta) * sin(phi)
            z = radius * cos(theta)
            
            push!(vertices, SVector(x, y, z))
        end
    end
    
    for i in 0:n_theta-1
        for j in 0:n_phi-1
            v1 = i * n_phi + j + 1
            v2 = i * n_phi + ((j + 1) % n_phi) + 1
            v3 = (i + 1) * n_phi + j + 1
            v4 = (i + 1) * n_phi + ((j + 1) % n_phi) + 1
            
            if i == 0
                push!(triangle_indices, SVector(v1, v3, v4))
            elseif i == n_theta - 1
                push!(triangle_indices, SVector(v1, v2, v3))
            else
                push!(triangle_indices, SVector(v1, v2, v4))
                push!(triangle_indices, SVector(v1, v4, v3))
            end
        end
    end
    
    return vertices, triangle_indices
end

function create_plate_mesh(length::Float64, width::Float64, n_x::Int, n_y::Int)
    """Create flat plate mesh"""
    vertices = SVector{3, Float64}[]
    triangle_indices = SVector{3, Int}[]
    
    for i in 0:n_y
        for j in 0:n_x
            x = (j / n_x - 0.5) * length
            y = (i / n_y - 0.5) * width
            z = 0.0
            push!(vertices, SVector(x, y, z))
        end
    end
    
    for i in 0:n_y-1
        for j in 0:n_x-1
            v1 = i * (n_x + 1) + j + 1
            v2 = i * (n_x + 1) + j + 2
            v3 = (i + 1) * (n_x + 1) + j + 1
            v4 = (i + 1) * (n_x + 1) + j + 2
            
            push!(triangle_indices, SVector(v1, v2, v3))
            push!(triangle_indices, SVector(v2, v4, v3))
        end
    end
    
    return vertices, triangle_indices
end

function create_corner_reflector_mesh(size::Float64)
    """
    Create a simple corner reflector (two perpendicular plates).
    
    This is a simplified 90° corner reflector consisting of two
    perpendicular square plates meeting at a right angle.
    """
    vertices = SVector{3, Float64}[]
    triangle_indices = SVector{3, Int}[]
    
    n_div = 10  # Divisions per plate edge
    
    # Plate 1: YZ plane (x = 0)
    for i in 0:n_div
        for j in 0:n_div
            x = 0.0
            y = i * size / n_div
            z = j * size / n_div
            push!(vertices, SVector(x, y, z))
        end
    end
    
    # Plate 2: XZ plane (y = 0)
    for i in 0:n_div
        for j in 0:n_div
            x = i * size / n_div
            y = 0.0
            z = j * size / n_div
            push!(vertices, SVector(x, y, z))
        end
    end
    
    # Triangulate Plate 1
    for i in 0:n_div-1
        for j in 0:n_div-1
            v1 = i * (n_div + 1) + j + 1
            v2 = i * (n_div + 1) + j + 2
            v3 = (i + 1) * (n_div + 1) + j + 1
            v4 = (i + 1) * (n_div + 1) + j + 2
            
            push!(triangle_indices, SVector(v1, v2, v3))
            push!(triangle_indices, SVector(v2, v4, v3))
        end
    end
    
    # Triangulate Plate 2 (offset vertex indices)
    offset = (n_div + 1)^2
    for i in 0:n_div-1
        for j in 0:n_div-1
            v1 = offset + i * (n_div + 1) + j + 1
            v2 = offset + i * (n_div + 1) + j + 2
            v3 = offset + (i + 1) * (n_div + 1) + j + 1
            v4 = offset + (i + 1) * (n_div + 1) + j + 2
            
            push!(triangle_indices, SVector(v1, v3, v2))  # Reverse orientation
            push!(triangle_indices, SVector(v2, v3, v4))
        end
    end
    
    return vertices, triangle_indices
end

function analyze_rcs_target(mesh::Mesh3D, target_name::String, frequency::Float64)
    """
    Analyze RCS target mesh and estimate electrical properties.
    
    Parameters:
    - mesh: Target mesh
    - target_name: Descriptive name
    - frequency: Analysis frequency in Hz
    
    This function provides mesh quality analysis and electrical size
    estimation for RCS prediction accuracy assessment.
    """
    println("=== RCS Target Analysis: $target_name ===")
    
    # Basic mesh properties
    println("Mesh Properties:")
    println("  Vertices: $(length(mesh.vertices))")
    println("  Triangles: $(length(mesh.triangles))")
    println("  Edges: $(mesh.num_edges)")
    println("  Watertight: $(mesh.is_watertight)")
    
    # Electrical size analysis
    c = 2.998e8  # Speed of light
    wavelength = c / frequency
    k = 2π / wavelength
    
    # Estimate target size
    x_coords = [v[1] for v in mesh.vertices]
    y_coords = [v[2] for v in mesh.vertices]
    z_coords = [v[3] for v in mesh.vertices]
    
    x_size = maximum(x_coords) - minimum(x_coords)
    y_size = maximum(y_coords) - minimum(y_coords)
    z_size = maximum(z_coords) - minimum(z_coords)
    max_size = max(x_size, y_size, z_size)
    
    electrical_size = k * max_size
    
    println("\nElectrical Properties:")
    @printf "  Frequency: %.1f MHz\n" frequency/1e6
    @printf "  Wavelength: %.3f m\n" wavelength
    @printf "  Target size: %.3f × %.3f × %.3f m\n" x_size y_size z_size
    @printf "  Electrical size (k·L): %.2f\n" electrical_size
    
    # Mesh density analysis
    avg_edge_length = sum(edge.length for edge in mesh.edges) / length(mesh.edges)
    mesh_density = wavelength / avg_edge_length
    
    println("\nMesh Density:")
    @printf "  Average edge length: %.4f m\n" avg_edge_length
    @printf "  Edges per wavelength: %.1f\n" mesh_density
    
    if mesh_density > 10.0
        println("  ✓ Mesh density adequate for RCS analysis")
    else
        println("  ⚠ Consider mesh refinement for better accuracy")
    end
    
    # Quality assessment
    quality = mesh.mesh_quality
    println("\nMesh Quality:")
    @printf "  Max aspect ratio: %.3f\n" quality.max_aspect_ratio
    @printf "  Min angle: %.1f°\n" quality.min_angle_deg
    
    if quality.max_aspect_ratio < 5.0 && quality.min_angle_deg > 10.0
        println("  ✓ Mesh quality suitable for RCS analysis")
    else
        println("  ⚠ Mesh quality may affect RCS accuracy")
    end
    
    return electrical_size, mesh_density
end

function estimate_rcs_characteristics(target_name::String, electrical_size::Float64)
    """
    Estimate expected RCS characteristics for different target types.
    
    This provides theoretical expectations for validation purposes.
    """
    println("\n=== Expected RCS Characteristics ===")
    
    if target_name == "sphere"
        println("Conducting Sphere (Mie Theory):")
        println("  • Low frequency (ka << 1): σ ∝ k⁴ (Rayleigh scattering)")
        println("  • Resonance region (ka ≈ 1): Complex oscillatory behavior")
        println("  • High frequency (ka >> 1): σ ≈ πa² (optical limit)")
        println("  • Gibson Section 7.7.2: Accuracy < 0.1 dB vs. theory")
        
        if electrical_size < 0.5
            println("  → Current case: Rayleigh scattering regime")
        elseif electrical_size < 5.0
            println("  → Current case: Resonance/Mie scattering regime")
        else
            println("  → Current case: Optical scattering regime")
        end
        
    elseif target_name == "plate"
        println("Flat Plate:")
        println("  • Normal incidence: σ = 4πA²/λ² (specular reflection)")
        println("  • Edge diffraction: Additional contributions at oblique angles")
        println("  • Polarization dependence: Different for E∥ and E⊥")
        println("  • Frequency dependence: σ ∝ f² for fixed physical size")
        
    elseif target_name == "corner"
        println("Corner Reflector:")
        println("  • High RCS due to multiple reflections")
        println("  • Retroreflective properties")
        println("  • Angular dependence: Peak at normal incidence")
        println("  • Frequency dependence: Complex due to multiple paths")
        
    end
end

function demonstrate_rcs_prediction_workflow()
    """
    Main demonstration of RCS prediction workflow using Method of Moments.
    
    This example follows Gibson's approach for RCS calculation and validation.
    """
    println("Radar Cross Section (RCS) Prediction using Method of Moments")
    println("Based on Gibson 'Method of Moments in Electromagnetics' Chapter 7")
    println("=" ^ 70)
    
    # Create target meshes
    println("Creating canonical RCS targets...")
    targets = create_simple_target_meshes()
    
    # Analysis frequency
    frequency = 300e6  # 300 MHz
    
    # Analyze each target
    for (target_name, target_data) in targets
        println("\n" * "="^60)
        
        # Create mesh
        mesh = Mesh3D(target_data.vertices, target_data.triangles)
        
        # Analyze target
        electrical_size, mesh_density = analyze_rcs_target(mesh, target_name, frequency)
        
        # Estimate RCS characteristics
        estimate_rcs_characteristics(target_name, electrical_size)
    end
    
    # General MoM RCS workflow
    println("\n" * "="^60)
    println("=== MoM RCS Calculation Workflow ===")
    
    println("\n1. Problem Setup:")
    println("   • Define target geometry and mesh")
    println("   • Specify incident plane wave (direction, polarization)")
    println("   • Choose analysis frequency range")
    
    println("\n2. MoM Solution:")
    println("   • Assemble EFIE impedance matrix")
    println("   • Compute incident field excitation vector")
    println("   • Solve for surface current distribution")
    
    println("\n3. Far-Field Calculation:")
    println("   • Apply far-field approximation")
    println("   • Integrate surface currents for scattered field")
    println("   • Calculate RCS: σ = lim(4πr²|E_s|²/|E_i|²)")
    
    println("\n4. Post-Processing:")
    println("   • Plot RCS vs. frequency")
    println("   • Generate RCS patterns (θ, φ dependence)")
    println("   • Analyze polarization effects")
    println("   • Compare with analytical solutions")
    
    # Validation approaches
    println("\n=== Validation Approaches ===")
    println("Gibson Section 7.7.2 - Conducting Sphere:")
    println("  • Compare with Mie theory")
    println("  • Expected accuracy: < 0.1 dB")
    println("  • Frequency range: 10-500 MHz")
    
    println("\nGibson Section 7.7.3 - EMCC Benchmark Targets:")
    println("  • Wedge cylinder")
    println("  • Plate cylinder") 
    println("  • Business card geometry")
    println("  • Measured vs. computed comparisons")
    
    println("\nGeneral Validation:")
    println("  • Reciprocity checks")
    println("  • Convergence with mesh refinement")
    println("  • Energy conservation")
    println("  • Cross-validation with other codes")
    
    # Expected challenges
    println("\n=== Common Challenges ===")
    println("• Edge diffraction modeling")
    println("• Multiple scattering effects")
    println("• Resonance phenomena")
    println("• Computational complexity for large targets")
    println("• Mesh quality requirements")
    
    return targets
end

# Run demonstration if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    targets = demonstrate_rcs_prediction_workflow()
end
