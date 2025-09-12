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

function compute_input_impedance(mesh::Mesh3D, frequency::Float64; V_feed::ComplexF64 = 1.0+0im)
    # Assemble EFIE matrix
    println("Assembling EFIE matrix at $(frequency/1e6) MHz...")
    Z = assemble_efie_matrix(mesh, frequency; progress=true)

    # Build a delta-gap excitation vector: apply unit voltage across the central edge
    # Find a candidate feed edge: interior edge closest to origin (0,0,0)
    feed_edge_idx = 0
    min_dist = Inf
    for (i, e) in enumerate(mesh.edges)
        d = norm(e.center)
        if d < min_dist
            min_dist = d
            feed_edge_idx = i
        end
    end

    if feed_edge_idx == 0
        error("No feed edge found in mesh")
    end

    # Incident field function for delta-gap: create right-hand side vector b such that
    # the voltage across the feed edge equals V_feed. A simple approach is to set b
    # equal to zeros except a driving entry on the feed edge. In a rigorous delta-gap
    # formulation b is constructed from the excitation integral; here we approximate
    # by placing the voltage on the RWG associated with the feed edge.
    n = mesh.num_edges
    b = zeros(ComplexF64, n)
    b[feed_edge_idx] = V_feed

    # Solve the MoM system Z * I = b
    println("Solving linear system (n = $(n))...")
    I, report = solve_mom_system(Z, b; method=:direct, return_report=true)
    println("Solver report: method=$(report.method) converged=$(report.converged) residual=$(report.residual_norm)")

    # Estimate input current: for a delta-gap feed the input current is the coefficient
    # on the driven RWG; use that coefficient as I_feed.
    I_feed = I[feed_edge_idx]

    Z_in = V_feed / I_feed

    println("Computed input impedance at $(frequency/1e6) MHz: $(Z_in)")

    # Return solver report and frequency for later summarization
    return Z_in, I, feed_edge_idx, report, frequency
end


# ======================================================

using Profile
using ProfileView
Profile.clear()

# Demonstrate different antenna types and compute impedances
dipole_mesh, dipole_freq = demonstrate_dipole_impedance() ;
patch_mesh, patch_freq = demonstrate_patch_impedance() ;

# Compute input impedances (unit voltage feed)


Profile.clear()

@profile Z_dipole, I_dipole, feed_idx_d, report_d, freq_d = compute_input_impedance(dipole_mesh, dipole_freq);

ProfileView.view()


Z_patch, I_patch, feed_idx_p, report_p, freq_p = compute_input_impedance(patch_mesh, patch_freq)
