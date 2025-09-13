using Test
using MoM3D
using StaticArrays

@testset "Edge opposite indices and RWG fast path" begin
    # Simple mesh: two triangles sharing edge (1,2)
    verts = [SVector(0.0,0.0,0.0), SVector(1.0,0.0,0.0), SVector(0.0,1.0,0.0)]
    # triangles: (1,2,3) and (2,1,3) share edge (1,2)
    tris = [SVector(1,2,3), SVector(2,1,3)]

    mesh = Mesh3D(verts, tris)

    @test length(mesh.edges) >= 1

    # find edge with vertices 1 and 2
    eidx = findfirst(e -> (e.vertex1 == 1 && e.vertex2 == 2) || (e.vertex1 == 2 && e.vertex2 == 1), mesh.edges)
    @test eidx !== nothing
    edge = mesh.edges[eidx]

    @test edge.opp_plus != 0
    @test edge.opp_minus != 0

    # Build RWG and evaluate at a point inside tri_plus
    rwgs = MoM3D.get_rwgs(mesh)
    rwg = rwgs[eidx]

    tri_plus = mesh.triangles[rwg.triangle_plus]
    # pick centroid of tri_plus
    v1,v2,v3 = tri_plus.vertices
    p = (v1 + v2 + v3) / 3.0

    # Evaluate using current fast path
    f_fast = evaluate_rwg(rwg, p, mesh)

    # Force error by zeroing opp_plus temporarily and evaluating: no fallback path exists
    saved_opp = edge.opp_plus
    mesh.edges[eidx] = Edge(edge.triangle_plus, edge.triangle_minus, edge.vertex1, edge.vertex2, mesh.vertices[edge.vertex1], mesh.vertices[edge.vertex2], 0, edge.opp_minus)
    rwgs2 = MoM3D.get_rwgs(mesh)
    rwg2 = rwgs2[eidx]

    @test_throws ErrorException evaluate_rwg(rwg2, p, mesh)

    # restore
    mesh.edges[eidx] = Edge(edge.triangle_plus, edge.triangle_minus, edge.vertex1, edge.vertex2, mesh.vertices[edge.vertex1], mesh.vertices[edge.vertex2], saved_opp, edge.opp_minus)
end
