using Test
using MoM3D
using StaticArrays

@testset "recompute_edge_opposites! and RWG cache invalidation" begin
    verts = [SVector(0.0,0.0,0.0), SVector(1.0,0.0,0.0), SVector(0.0,1.0,0.0), SVector(1.0,1.0,0.0)]
    tris = [SVector(1,2,3), SVector(2,1,4)]
    mesh = Mesh3D(verts, tris)

    eidx = findfirst(e -> Set([e.vertex1, e.vertex2]) == Set([1,2]), mesh.edges)
    @test eidx !== nothing
    edge = mesh.edges[eidx]
    @test Set([edge.opp_plus, edge.opp_minus]) == Set([3,4])

    # Modify second triangle in-place to change opposite vertex from 4 -> 5
    push!(mesh.vertices, SVector(2.0,2.0,0.0))
    mesh.triangles[2] = Triangle(mesh.vertices[2], mesh.vertices[1], mesh.vertices[5])

    # At this point the edges still have old opp indices
    edge_old = mesh.edges[eidx]
    @test Set([edge_old.opp_plus, edge_old.opp_minus]) == Set([3,4])

    # Recompute opposites
    recompute_edge_opposites!(mesh)

    edge_new = mesh.edges[eidx]
    @test Set([edge_new.opp_plus, edge_new.opp_minus]) == Set([3,5])

    # Ensure cached RWGFunctions will pick up new opps (build cache first)
    rwgs = get_rwgs(mesh)
    rwg = rwgs[eidx]

    # centroid of triangle plus
    tri_plus = mesh.triangles[rwg.triangle_plus]
    v1,v2,v3 = tri_plus.vertices
    p = (v1 + v2 + v3) / 3.0

    # should not error now (opp indices present)
    f = evaluate_rwg(rwg, p, mesh)
    @test f != SVector(0.0,0.0,0.0)
end
