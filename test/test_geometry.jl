using Test
using MoM3D
using LinearAlgebra
using StaticArrays

@testset "Geometry Tests" begin
    
    @testset "Triangle Creation" begin
        v1 = SVector(0.0, 0.0, 0.0)
        v2 = SVector(1.0, 0.0, 0.0)
        v3 = SVector(0.0, 1.0, 0.0)
        
        tri = Triangle(v1, v2, v3)
        
        @test tri.area ≈ 0.5
        @test norm(tri.normal - SVector(0.0, 0.0, 1.0)) < 1e-12
    end
    
    @testset "Simple Mesh Creation" begin
        vertices = [
            SVector(0.0, 0.0, 0.0),
            SVector(1.0, 0.0, 0.0),
            SVector(0.0, 1.0, 0.0),
            SVector(1.0, 1.0, 0.0)
        ]
        
        triangle_indices = [
            SVector(1, 2, 3),
            SVector(2, 4, 3)
        ]
        
        mesh = Mesh3D(vertices, triangle_indices)
        
        @test length(mesh.vertices) == 4
        @test length(mesh.triangles) == 2
        @test mesh.num_edges > 0
        @test !mesh.is_watertight  # This is an open surface
    end
    
    @testset "Aspect Ratio Calculation" begin
        # Equilateral triangle
        v1 = SVector(0.0, 0.0, 0.0)
        v2 = SVector(1.0, 0.0, 0.0)
        v3 = SVector(0.5, sqrt(3)/2, 0.0)
        
        tri_equilateral = Triangle(v1, v2, v3)
        aspect_ratio = compute_aspect_ratio(tri_equilateral)
        
        @test aspect_ratio ≈ 1.0 atol=1e-10
        
        # Thin triangle
        v3_thin = SVector(0.5, 0.01, 0.0)
        tri_thin = Triangle(v1, v2, v3_thin)
        aspect_ratio_thin = compute_aspect_ratio(tri_thin)
        
        @test aspect_ratio_thin > 10.0
    end
    
    @testset "Angle Calculations" begin
        # Right triangle
        v1 = SVector(0.0, 0.0, 0.0)
        v2 = SVector(1.0, 0.0, 0.0)
        v3 = SVector(0.0, 1.0, 0.0)
        
        tri = Triangle(v1, v2, v3)
        
        min_angle = compute_min_angle(tri)
        max_angle = compute_max_angle(tri)
        
        @test min_angle ≈ π/4 atol=1e-10  # 45 degrees
        @test max_angle ≈ π/2 atol=1e-10  # 90 degrees
    end
    
    @testset "Mesh Quality Analysis" begin
        vertices = [
            SVector(0.0, 0.0, 0.0),
            SVector(1.0, 0.0, 0.0),
            SVector(0.0, 1.0, 0.0)
        ]
        
        triangle_indices = [SVector(1, 2, 3)]
        
        mesh = Mesh3D(vertices, triangle_indices)
        
        @test mesh.mesh_quality.min_area > 0
        @test mesh.mesh_quality.max_area > 0
        @test mesh.mesh_quality.degenerate_triangles == 0
        @test mesh.mesh_quality.min_angle_deg > 0
        @test mesh.mesh_quality.max_angle_deg < 180
    end
    
    @testset "Mesh Validation" begin
        # Create a simple valid mesh
        vertices = [
            SVector(0.0, 0.0, 0.0),
            SVector(1.0, 0.0, 0.0),
            SVector(0.0, 1.0, 0.0)
        ]
        
        triangle_indices = [SVector(1, 2, 3)]
        mesh = Mesh3D(vertices, triangle_indices)
        
        validation = validate_mesh(mesh)
        
        @test !validation.watertight  # Single triangle is not watertight
        @test validation.quality_metrics.degenerate_triangles == 0
    end
    
    @testset "Mesh Repair" begin
        # Create vertices with near-duplicates
        vertices_with_duplicates = [
            SVector(0.0, 0.0, 0.0),
            SVector(1.0, 0.0, 0.0),
            SVector(0.0, 1.0, 0.0),
            SVector(1e-15, 1e-15, 1e-15),  # Near-duplicate of first vertex
        ]
        
        triangle_indices = [
            SVector(1, 2, 3),
            SVector(4, 2, 3)
        ]
        
        repaired_vertices, repaired_triangles = repair_mesh_connectivity(
            vertices_with_duplicates, triangle_indices, tolerance=1e-10)
        
        @test length(repaired_vertices) < length(vertices_with_duplicates)
        @test length(repaired_triangles) <= length(triangle_indices)
    end
end
