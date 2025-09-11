"""
Antenna Impedance Calculation Example

This example demonstrates the Method of Moments approach for calculating antenna
input impedance, based on the theoretical framework in Gibson's "The Method of
Moments in Electromagnetics."

Theoretical Background:
- Antenna impedance calculation is a fundamental MoM application
- Gibson discusses antenna problems throughout the text, particularly wire antennas
- Input impedance Z = V/I where V is applied voltage and I is input current
- MoM provides surface current distribution, from which impedance is calculated

Key Concepts:
- Delta-gap excitation for voltage sources
- Current integration for input impedance
- Resonance behavior and bandwidth analysis
- Validation against analytical solutions (thin wire theory)

Applications:
- Antenna design and optimization
- Impedance matching network design
- Bandwidth and efficiency analysis
- Electromagnetic compatibility studies
"""

using MoM3D
using LinearAlgebra
using StaticArrays
using Printf

function create_dipole_mesh(length::Float64, radius::Float64, n_segments::Int)
    """
    Create a simple wire dipole mesh using cylindrical segments.
    
    Parameters:
    - length: Total dipole length in meters
    - radius: Wire radius in meters  
    - n_segments: Number of segments along the wire
    
    Returns:
    - vertices: Array of 3D vertex coordinates
    - triangle_indices: Array of triangle connectivity
    
    Note: This is a simplified representation. Real antenna modeling
    requires careful treatment of wire junctions and feed points.
    """
    vertices = SVector{3, Float64}[]
    triangle_indices = SVector{3, Int}[]
    
    # Create vertices along the wire
    half_length = length / 2.0
    segment_length = length / n_segments
    
    # Generate vertices for cylindrical wire representation
    n_circumferential = 8  # Points around wire circumference
    
    for i in 0:n_segments
        z = -half_length + i * segment_length
        
        for j in 0:n_circumferential-1
            theta = 2π * j / n_circumferential
            x = radius * cos(theta)
            y = radius * sin(theta)
            
            push!(vertices, SVector(x, y, z))
        end
    end
    
    # Generate triangular connectivity for wire surface
    for i in 0:n_segments-1
        for j in 0:n_circumferential-1
            # Current ring vertices
            v1 = i * n_circumferential + j + 1
            v2 = i * n_circumferential + ((j + 1) % n_circumferential) + 1
            
            # Next ring vertices
            v3 = (i + 1) * n_circumferential + j + 1
            v4 = (i + 1) * n_circumferential + ((j + 1) % n_circumferential) + 1
            
            # Two triangles per quadrilateral
            push!(triangle_indices, SVector(v1, v2, v3))
            push!(triangle_indices, SVector(v2, v4, v3))
        end
    end
    
    return vertices, triangle_indices
end

function create_patch_antenna_mesh(length::Float64, width::Float64, n_x::Int, n_y::Int)
    """
    Create a rectangular patch antenna mesh.
    
    Parameters:
    - length: Patch length in meters (resonant dimension)
    - width: Patch width in meters
    - n_x: Number of divisions along length
    - n_y: Number of divisions along width
    
    Returns:
    - vertices: Array of 3D vertex coordinates
    - triangle_indices: Array of triangle connectivity
    
    Note: This creates only the top surface. A complete model would
    include ground plane and substrate effects.
    """
    vertices = SVector{3, Float64}[]
    triangle_indices = SVector{3, Int}[]
    
    # Generate vertices on rectangular patch
    for i in 0:n_y
        for j in 0:n_x
            x = (j / n_x - 0.5) * length
            y = (i / n_y - 0.5) * width
            z = 0.0  # Patch at z = 0
            
            push!(vertices, SVector(x, y, z))
        end
    end
    
    # Generate triangular connectivity
    for i in 0:n_y-1
        for j in 0:n_x-1
            v1 = i * (n_x + 1) + j + 1
            v2 = i * (n_x + 1) + j + 2
            v3 = (i + 1) * (n_x + 1) + j + 1
            v4 = (i + 1) * (n_x + 1) + j + 2
            
            # Two triangles per rectangle
            push!(triangle_indices, SVector(v1, v2, v3))
            push!(triangle_indices, SVector(v2, v4, v3))
        end
    end
    
    return vertices, triangle_indices
end

function analyze_antenna_mesh(mesh::Mesh3D, antenna_type::String)
    """
    Analyze antenna mesh quality for MoM impedance calculation.
    
    Quality requirements for antenna analysis:
    - Good aspect ratios for matrix conditioning
    - Adequate mesh density (typically λ/10 to λ/20)
    - Proper treatment of feed regions
    - Watertight mesh for current continuity
    """
    println("=== Antenna Mesh Analysis: $antenna_type ===")
    println("Number of vertices: $(length(mesh.vertices))")
    println("Number of triangles: $(length(mesh.triangles))")
    println("Number of edges: $(mesh.num_edges)")
    
    # Quality metrics
    quality = mesh.mesh_quality
    println("\nMesh Quality:")
    @printf "  Max aspect ratio: %.3f\n" quality.max_aspect_ratio
    @printf "  Mean aspect ratio: %.3f\n" quality.mean_aspect_ratio
    @printf "  Min angle: %.1f°\n" quality.min_angle_deg
    @printf "  Max angle: %.1f°\n" quality.max_angle_deg
    
    # Quality assessment for antenna applications
    println("\nQuality Assessment for Antenna Analysis:")
    if quality.max_aspect_ratio < 5.0
        println("  ✓ Aspect ratios suitable for MoM analysis")
    else
        println("  ⚠ Poor aspect ratios may affect convergence")
    end
    
    if quality.min_angle_deg > 15.0
        println("  ✓ Triangle angles are well-conditioned")
    else
        println("  ⚠ Small angles may cause numerical issues")
    end
    
    return quality
end

function estimate_resonant_frequency(length::Float64, antenna_type::String)
    """
    Estimate resonant frequency for common antenna types.
    
    Parameters:
    - length: Characteristic antenna dimension in meters
    - antenna_type: "dipole" or "patch"
    
    Returns:
    - f_resonant: Estimated resonant frequency in Hz
    - wavelength: Corresponding wavelength in meters
    
    Note: These are approximate values. Actual resonance depends on
    wire radius, substrate properties, and other factors.
    """
    c = 2.998e8  # Speed of light in m/s
    
    if antenna_type == "dipole"
        # Half-wave dipole: length ≈ λ/2
        wavelength = 2.0 * length
    elseif antenna_type == "patch"
        # Rectangular patch: length ≈ λ/2 (with fringing effects)
        # Effective length is slightly longer due to fringing
        effective_length = length * 1.05  # Approximate fringing correction
        wavelength = 2.0 * effective_length
    else
        error("Unknown antenna type: $antenna_type")
    end
    
    f_resonant = c / wavelength
    
    return f_resonant, wavelength
end

function demonstrate_dipole_impedance()
    """
    Demonstrate dipole antenna impedance calculation setup.
    
    This example shows the mesh preparation and analysis steps
    for a wire dipole antenna impedance calculation.
    """
    println("Wire Dipole Antenna Impedance Analysis")
    println(repeat("=", 45))
    
    # Dipole parameters
    dipole_length = 0.5  # meters (approximately λ/2 at 300 MHz)
    radius = 0.001  # meters (thin wire)
    n_segments = 20
    
    println("Dipole parameters:")
    @printf "  Length: %.3f m\n" dipole_length
    @printf "  Radius: %.3f m\n" radius
    @printf "  Segments: %d\n" n_segments
    
    # Estimate resonant frequency
    f_res, wavelength = estimate_resonant_frequency(dipole_length, "dipole")
    @printf "  Estimated resonant frequency: %.1f MHz\n" f_res/1e6
    @printf "  Corresponding wavelength: %.3f m\n" wavelength
    
    # Create mesh
    println("\nCreating dipole mesh...")
    vertices, triangle_indices = create_dipole_mesh(dipole_length, radius, n_segments)
    mesh = Mesh3D(vertices, triangle_indices)
    
    # Analyze mesh quality
    quality = analyze_antenna_mesh(mesh, "Wire Dipole")
    
    # Mesh density analysis
    avg_edge_length = sum(edge.length for edge in mesh.edges) / length(mesh.edges)
    mesh_density = wavelength / avg_edge_length
    
    println("\nMesh Density Analysis:")
    @printf "  Average edge length: %.4f m\n" avg_edge_length
    @printf "  Mesh density: %.1f edges per wavelength\n" mesh_density
    
    if mesh_density > 10.0
        println("  ✓ Mesh density adequate for MoM analysis")
    else
        println("  ⚠ Consider mesh refinement for better accuracy")
    end
    
    return mesh, f_res
end

function demonstrate_patch_impedance()
    """
    Demonstrate patch antenna impedance calculation setup.
    
    This example shows the mesh preparation for a rectangular
    microstrip patch antenna impedance calculation.
    """
    println("\nRectangular Patch Antenna Impedance Analysis")
    println(repeat("=", 48))
    
    # Patch parameters (typical values for 2.4 GHz)
    patch_length = 0.031  # meters (resonant dimension)
    width = 0.024   # meters
    n_x = 15
    n_y = 12
    
    println("Patch parameters:")
    @printf "  Length: %.3f m\n" patch_length
    @printf "  Width: %.3f m\n" width
    @printf "  Mesh divisions: %d × %d\n" n_x n_y
    
    # Estimate resonant frequency
    f_res, wavelength = estimate_resonant_frequency(patch_length, "patch")
    @printf "  Estimated resonant frequency: %.1f GHz\n" f_res/1e9
    @printf "  Corresponding wavelength: %.3f m\n" wavelength
    
    # Create mesh
    println("\nCreating patch mesh...")
    vertices, triangle_indices = create_patch_antenna_mesh(patch_length, width, n_x, n_y)
    mesh = Mesh3D(vertices, triangle_indices)
    
    # Analyze mesh quality
    quality = analyze_antenna_mesh(mesh, "Rectangular Patch")
    
    return mesh, f_res
end

function demonstrate_impedance_calculation_workflow()
    """
    Main demonstration of antenna impedance calculation workflow.
    
    This example outlines the complete process for MoM-based
    antenna impedance calculation, following Gibson's methodology.
    """
    println("Method of Moments Antenna Impedance Calculation")
    println("Based on Gibson 'Method of Moments in Electromagnetics'")
    println(repeat("=", 60))
    
    # Demonstrate different antenna types
    dipole_mesh, dipole_freq = demonstrate_dipole_impedance()
    patch_mesh, patch_freq = demonstrate_patch_impedance()
    
    # General MoM workflow for impedance calculation
    println("\n=== MoM Impedance Calculation Workflow ===")
    println("1. Mesh Generation:")
    println("   • Create surface triangulation of antenna geometry")
    println("   • Ensure adequate mesh density (λ/10 to λ/20)")
    println("   • Verify mesh quality (aspect ratios, angles)")
    
    println("\n2. Basis Function Setup:")
    println("   • Define RWG basis functions on mesh edges")
    println("   • Handle feed point excitation (delta-gap source)")
    println("   • Ensure current continuity across triangles")
    
    println("\n3. Matrix Assembly:")
    println("   • Compute EFIE impedance matrix elements")
    println("   • Apply numerical integration (singular/regular)")
    println("   • Include feed point boundary conditions")
    
    println("\n4. System Solution:")
    println("   • Solve Z·I = V for surface current coefficients")
    println("   • Extract input current at feed point")
    println("   • Calculate input impedance Z_in = V_feed / I_feed")
    
    println("\n5. Post-Processing:")
    println("   • Plot impedance vs. frequency")
    println("   • Identify resonant frequencies")
    println("   • Calculate bandwidth and efficiency")
    println("   • Validate against analytical solutions")
    
    # Expected results
    println("\n=== Expected Results ===")
    println("Wire Dipole (λ/2):")
    println("  • Input resistance: ~73 Ω at resonance")
    println("  • Input reactance: ~0 Ω at resonance")
    println("  • Bandwidth: ~10% for VSWR < 2:1")
    
    println("\nRectangular Patch:")
    println("  • Input resistance: ~100-300 Ω (depends on feed location)")
    println("  • Input reactance: ~0 Ω at resonance")
    println("  • Bandwidth: ~2-5% for VSWR < 2:1")
    
    println("\n=== Validation Approaches ===")
    println("• Compare with analytical solutions (thin wire theory)")
    println("• Cross-validate with commercial EM simulators")
    println("• Verify convergence with mesh refinement")
    println("• Check reciprocity and energy conservation")
    
    return dipole_mesh, patch_mesh
end

# Run demonstration if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    dipole_mesh, patch_mesh = demonstrate_impedance_calculation_workflow()
end
