using Test
using StaticArrays
using LinearAlgebra
using MoM3D

const G = MoM3D.Geometry
const IE = MoM3D.IntegralEquations
const BF = MoM3D.BasisFunctions

# Helper: build a small mesh (two triangles sharing an edge)
function make_simple_mesh()
    verts = [
        SVector(0.0, 0.0, 0.0),
        SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(1.0, 1.0, 0.0)
    ]
    tris = [SVector(1,2,3), SVector(2,4,3)]
    return G.Mesh3D(verts, tris)
end

@testset "IntegrationCache correctness (enriched)" begin
    mesh = make_simple_mesh()
    ntri = length(mesh.triangles)

    # Build caches for a few orders
    orders = (3, 5, 7)
    caches = IE.ensure_integration_caches(mesh; orders=orders)
    @test isa(caches, G.IntegrationCache)

    # Validate order mapping/shape
    @test length(caches.orders) == length(orders)
    for ord in caches.orders
        @test haskey(caches.order_index, ord)
        oi = caches.order_index[ord]
        @test length(caches.tri_cart_pts[oi]) == ntri
        bpts, w = caches.gauss_rules[oi]
        @test length(caches.tri_cart_pts[oi][1]) == length(bpts)
    end

    # For each triangle and each order: cached RWG vectors must equal direct evaluation
    rwgs = BF.get_rwgs(mesh)
    @test !isempty(rwgs)

    for (oi, ord) in enumerate(caches.orders)
        tri_cart = caches.tri_cart_pts[oi]
        tri_rwg = caches.tri_rwg_vals[oi]
        for tid in 1:ntri
            pts = tri_cart[tid]
            rwg_dict = tri_rwg[tid]

            # Every RWG adjacent to tid (via mesh edges) should have an entry
            # find adjacency via edges stored in mesh
            adjacent_rwgs = Int[]
            for (i, r) in enumerate(rwgs)
                if r.triangle_plus == tid || r.triangle_minus == tid
                    push!(adjacent_rwgs, i)
                end
            end

            # If there are adjacent rwgs, ensure cache contains them and values match eval
            for rglobal in adjacent_rwgs
                r = rwgs[rglobal]
                @test haskey(rwg_dict, r.edge_index)
                vals = rwg_dict[r.edge_index]
                @test length(vals) == length(pts)
                for (i, p) in enumerate(pts)
                    direct = BF.evaluate_rwg(r, p, mesh)
                    cached = vals[i]
                    # Use a slightly relaxed tolerance for geometry rounding
                    @test isapprox(direct, cached; atol=1e-12, rtol=0)
                end
            end
        end
    end

    # Numeric parity: pick RWG pair (use same RWG twice if mesh is very small)
    rwg_m = rwgs[1]
    rwg_n = length(rwgs) >= 2 ? rwgs[2] : rwgs[1]
    tri1 = mesh.triangles[1]
    tri2 = mesh.triangles[2]
    k = 2π

    el_cached = IE.compute_regular_efie_element(1, 2, tri1, tri2, rwg_m, rwg_n, mesh, k; caches=caches)
    el_direct = IE.compute_regular_efie_element(1, 2, tri1, tri2, rwg_m, rwg_n, mesh, k; caches=nothing)

    # Relax tolerance for potential numerical extraction differences
    @test isapprox(el_cached, el_direct; atol=1e-8, rtol=1e-8)
end
