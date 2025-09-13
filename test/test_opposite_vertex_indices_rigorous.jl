using Test
using MoM3D
using StaticArrays

@testset "Opposite vertex indices - rigorous" begin

    @testset "Interior edge opposite indices (unordered)" begin
        verts = [SVector(0.0,0.0,0.0), SVector(1.0,0.0,0.0), SVector(0.0,1.0,0.0), SVector(1.0,1.0,0.0)]
        tris = [SVector(1,2,3), SVector(2,1,4)]
        mesh = Mesh3D(verts, tris)

        # find edge shared by 1 and 2
    eidx = findfirst(e -> Set([e.vertex1, e.vertex2]) == Set([1,2]), mesh.edges)
        @test eidx !== nothing
        edge = mesh.edges[eidx]

        # opposites should be 3 and 4 in some order
        opps = Set([edge.opp_plus, edge.opp_minus])
        @test opps == Set([3,4])
    end

    @testset "Boundary edges not created as Edge objects" begin
        verts = [SVector(0.0,0.0,0.0), SVector(1.0,0.0,0.0), SVector(0.0,1.0,0.0)]
        tris = [SVector(1,2,3)]
        mesh = Mesh3D(verts, tris)

        # find_edges exposes boundary edge map via second return value
        edges, edge_map, node_conn, facet_conn = find_edges(verts, tris)
        @test haskey(edge_map, (1,2))
        @test length(edge_map[(1,2)]) == 1

        # edges list should be empty as there are no interior edges
        @test isempty(edges)
    end

    @testset "Degenerate triangles (colinear) compute indices correctly" begin
        # colinear points
        verts = [SVector(0.0,0.0,0.0), SVector(1.0,0.0,0.0), SVector(2.0,0.0,0.0), SVector(0.0,1.0,0.0)]
        tris = [SVector(1,2,3), SVector(2,1,4)]
        mesh = Mesh3D(verts, tris)

        # ensure opposites are taken from triangle indices, even if area is zero
    eidx = findfirst(e -> Set([e.vertex1, e.vertex2]) == Set([1,2]), mesh.edges)
        @test eidx !== nothing
        edge = mesh.edges[eidx]
        opps = Set([edge.opp_plus, edge.opp_minus])
        @test opps == Set([3,4])
    end

    @testset "Non-manifold edge (more than two adjacent triangles)" begin
        verts = [SVector(0.0,0.0,0.0), SVector(1.0,0.0,0.0), SVector(0.0,1.0,0.0), SVector(-1.0,1.0,0.0), SVector(2.0,1.0,0.0)]
        # triangles 1,2,3 ; 1,2,4 ; 1,2,5 share edge (1,2)
        tris = [SVector(1,2,3), SVector(1,2,4), SVector(1,2,5)]
        edges, edge_map, node_conn, facet_conn = find_edges(verts, tris)

        @test haskey(edge_map, (1,2))
        @test length(edge_map[(1,2)]) == 3
        # find_edges should not create an Edge object for non-manifold shared-by-3
        eidx = findfirst(e -> Set([e.vertex1, e.vertex2]) == Set([1,2]), edges)
        @test eidx === nothing
    end

    @testset "update_mesh! recomputes opp indices correctly" begin
        verts = [SVector(0.0,0.0,0.0), SVector(1.0,0.0,0.0), SVector(0.0,1.0,0.0), SVector(1.0,1.0,0.0)]
        tris = [SVector(1,2,3), SVector(2,1,4)]
        mesh = Mesh3D(verts, tris)
        eidx = findfirst(e -> Set([e.vertex1, e.vertex2]) == Set([1,2]), mesh.edges)
        @test eidx !== nothing
        edge = mesh.edges[eidx]
        @test Set([edge.opp_plus, edge.opp_minus]) == Set([3,4])

        # modify triangles: replace second triangle with a new one sharing edge 1-2 but opposite vertex 5
        push!(mesh.vertices, SVector(2.0,2.0,0.0))
        tri_new = SVector(2,1,5)
        update_mesh!(mesh, mesh.vertices, [tris[1], tri_new])

        # find updated edge
        eidx2 = findfirst(e -> Set([e.vertex1, e.vertex2]) == Set([1,2]), mesh.edges)
        @test eidx2 !== nothing
        edge2 = mesh.edges[eidx2]
        @test Set([edge2.opp_plus, edge2.opp_minus]) == Set([3,5])
    end

end
