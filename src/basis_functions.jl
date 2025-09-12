module BasisFunctions

using LinearAlgebra
using StaticArrays
using ..Geometry

export RWGFunction, evaluate_rwg, evaluate_rwg_divergence, get_rwgs, invalidate_rwg_cache!

struct RWGFunction
    edge_index::Int
    triangle_plus::Int
    triangle_minus::Int
    edge_length::Float64
    area_plus::Float64
    area_minus::Float64

    function RWGFunction(edge::Edge, mesh::Mesh3D)
        # Determine the index of this edge within the mesh.edges array
        edge_idx = findfirst(e -> (e.vertex1 == edge.vertex1 && e.vertex2 == edge.vertex2 && e.triangle_plus == edge.triangle_plus && e.triangle_minus == edge.triangle_minus), mesh.edges)
        if edge_idx === nothing
            error("Edge not found in mesh when constructing RWGFunction")
        end
        tri_plus = mesh.triangles[edge.triangle_plus]
        tri_minus = mesh.triangles[edge.triangle_minus]
        new(edge_idx, edge.triangle_plus, edge.triangle_minus, edge.length, tri_plus.area, tri_minus.area)
    end
end

"""
get_rwgs(mesh::Mesh3D) -> Vector{RWGFunction}

Return a cached vector of `RWGFunction` instances for `mesh`.
If the cache does not exist it is built and stored on the mesh, and
the same instances are returned on subsequent calls. This avoids
repeated construction and allocations when multiple routines (assembly,
postprocessing, etc.) need the same RWG basis.

Note: If you modify the mesh geometry (vertices, triangles, edges),
you must invalidate the cache via `invalidate_rwg_cache!(mesh)` or
call `update_mesh!` which clears the cache automatically.
"""
function get_rwgs(mesh::Mesh3D)
    # Cache RWGFunction instances on the mesh to avoid repeated construction
    if mesh.rwg_cache !== nothing
        return mesh.rwg_cache
    end

    n_edges = mesh.num_edges
    rwgs = Vector{RWGFunction}(undef, n_edges)
    for i in 1:n_edges
        rwgs[i] = RWGFunction(mesh.edges[i], mesh)
    end

    mesh.rwg_cache = rwgs
    return rwgs
end

"""
invalidate_rwg_cache!(mesh::Mesh3D)

Clear any cached RWG basis stored on `mesh`. Call this after performing
in-place changes to mesh geometry so subsequent calls to `get_rwgs`
reconstruct the basis from the updated mesh data.
"""
function invalidate_rwg_cache!(mesh::Mesh3D)
    mesh.rwg_cache = nothing
    return nothing
end

function evaluate_rwg(rwg::RWGFunction, point::SVector{3,Float64}, mesh::Mesh3D)
    edge = mesh.edges[rwg.edge_index]

    if point_in_triangle(point, mesh.triangles[rwg.triangle_plus])
        vertex_opposite = find_opposite_vertex(edge, mesh.triangles[rwg.triangle_plus], mesh)
        return (rwg.edge_length / (2.0 * rwg.area_plus)) * (point - vertex_opposite)
    elseif point_in_triangle(point, mesh.triangles[rwg.triangle_minus])
        vertex_opposite = find_opposite_vertex(edge, mesh.triangles[rwg.triangle_minus], mesh)
        return -(rwg.edge_length / (2.0 * rwg.area_minus)) * (point - vertex_opposite)
    else
        return SVector(0.0, 0.0, 0.0)
    end
end

function evaluate_rwg_divergence(rwg::RWGFunction)
    return rwg.edge_length / rwg.area_plus + rwg.edge_length / rwg.area_minus
end

function point_in_triangle(point::SVector{3,Float64}, triangle::Triangle)
    v1, v2, v3 = triangle.vertices

    v0 = v3 - v1
    v1_new = v2 - v1
    v2_new = point - v1

    dot00 = dot(v0, v0)
    dot01 = dot(v0, v1_new)
    dot02 = dot(v0, v2_new)
    dot11 = dot(v1_new, v1_new)
    dot12 = dot(v1_new, v2_new)

    inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * inv_denom
    v = (dot00 * dot12 - dot01 * dot02) * inv_denom

    return (u >= 0) && (v >= 0) && (u + v <= 1)
end

function find_opposite_vertex(edge::Edge, triangle::Triangle, mesh::Mesh3D)
    vertices = triangle.vertices
    for vertex in vertices
        v1 = mesh.vertices[edge.vertex1]
        v2 = mesh.vertices[edge.vertex2]
        if norm(vertex - v1) > 1e-12 && norm(vertex - v2) > 1e-12
            return vertex
        end
    end
    error("Could not find opposite vertex")
end

end
