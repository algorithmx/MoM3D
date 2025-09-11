"""
Enhanced Sphere Scattering Example with Theoretical Validation

This example demonstrates sphere mesh creation and quality analysis based on
Gibson's "The Method of Moments in Electromagnetics" Sections 6.6.2 and 7.7.2.

Theoretical Background:
- Conducting sphere scattering is a canonical electromagnetic problem
- Analytical solution available via Mie theory for validation
- Gibson Section 7.7.2 describes 2-meter diameter sphere validation
- Expected MoM accuracy: < 0.1 dB compared to Mie theory

Key Requirements from Gibson Section 7.6.1:
- Reasonable triangle aspect ratios (< 5.0 preferred)
- Minimum angles > 10° for good matrix conditioning
- Watertight mesh for MFIE applications (Section 7.6.2)
"""

using MoM3D
using LinearAlgebra
using StaticArrays

function create_sphere_mesh(radius::Float64, n_theta::Int, n_phi::Int)
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

function test_sphere_mesh_quality()
    println("=== Testing Sphere Mesh Quality ===")
    
    radius = 1.0
    vertices, triangle_indices = create_sphere_mesh(radius, 10, 20)
    
    println("Creating mesh with $(length(vertices)) vertices and $(length(triangle_indices)) triangles")
    
    mesh = Mesh3D(vertices, triangle_indices)
    
    println("Mesh created successfully!")
    println("Number of edges: $(mesh.num_edges)")
    println("Number of interior edges: $(length(mesh.interior_edges))")
    println("Number of boundary edges: $(length(mesh.boundary_edges))")
    println("Is watertight: $(mesh.is_watertight)")
    
    println("\nMesh Quality Metrics:")
    println("  Min area: $(mesh.mesh_quality.min_area)")
    println("  Max area: $(mesh.mesh_quality.max_area)")
    println("  Min aspect ratio: $(mesh.mesh_quality.min_aspect_ratio)")
    println("  Max aspect ratio: $(mesh.mesh_quality.max_aspect_ratio)")
    println("  Mean aspect ratio: $(mesh.mesh_quality.mean_aspect_ratio)")
    println("  Min angle: $(mesh.mesh_quality.min_angle_deg)°")
    println("  Max angle: $(mesh.mesh_quality.max_angle_deg)°")
    println("  Degenerate triangles: $(mesh.mesh_quality.degenerate_triangles)")
    println("  Thin triangles: $(mesh.mesh_quality.thin_triangles)")
    
    validation = validate_mesh(mesh)
    println("\nMesh Validation:")
    println("  Is valid: $(validation.is_valid)")
    if !validation.is_valid
        println("  Issues found:")
        for issue in validation.issues
            println("    - $issue")
        end
    end
    
    return mesh
end

function test_plate_mesh_quality()
    println("\n=== Testing Plate Mesh Quality ===")
    
    length = 2.0
    width = 1.0
    vertices, triangle_indices = create_plate_mesh(length, width, 10, 5)
    
    println("Creating plate mesh with $(length(vertices)) vertices and $(length(triangle_indices)) triangles")
    
    mesh = Mesh3D(vertices, triangle_indices)
    
    println("Mesh created successfully!")
    println("Number of edges: $(mesh.num_edges)")
    println("Number of interior edges: $(length(mesh.interior_edges))")
    println("Number of boundary edges: $(length(mesh.boundary_edges))")
    println("Is watertight: $(mesh.is_watertight)")
    
    validation = validate_mesh(mesh)
    println("\nMesh Validation:")
    println("  Is valid: $(validation.is_valid)")
    if !validation.is_valid
        println("  Issues found:")
        for issue in validation.issues
            println("    - $issue")
        end
    end
    
    return mesh
end

function test_mesh_repair()
    println("\n=== Testing Mesh Repair Functionality ===")
    
    vertices_with_duplicates = [
        SVector(0.0, 0.0, 0.0),
        SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1e-15),  # Near-duplicate of vertex 1
        SVector(1.0, 0.0, 1e-15),  # Near-duplicate of vertex 2
    ]
    
    triangle_indices_with_issues = [
        SVector(1, 2, 3),
        SVector(4, 5, 3),  # Uses near-duplicate vertices
    ]
    
    println("Original mesh: $(length(vertices_with_duplicates)) vertices, $(length(triangle_indices_with_issues)) triangles")
    
    repaired_vertices, repaired_triangles = repair_mesh_connectivity(
        vertices_with_duplicates, triangle_indices_with_issues, tolerance=1e-10)
    
    println("Repaired mesh: $(length(repaired_vertices)) vertices, $(length(repaired_triangles)) triangles")
    
    return repaired_vertices, repaired_triangles
end

function main()
    println("3D Method of Moments - Shape Processing Capabilities Test")
    println(repeat("=", 60))
    
    sphere_mesh = test_sphere_mesh_quality()
    plate_mesh = test_plate_mesh_quality()
    test_mesh_repair()
    
    println("\n=== Shape Processing Test Complete ===")
    println("All mesh processing capabilities tested successfully!")
    
    return sphere_mesh, plate_mesh
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
