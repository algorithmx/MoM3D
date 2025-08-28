module BasisFunctions

using LinearAlgebra
using StaticArrays
using ..Geometry

export RWGFunction, evaluate_rwg, evaluate_rwg_divergence

struct RWGFunction
    edge_index::Int
    triangle_plus::Int
    triangle_minus::Int
    edge_length::Float64
    area_plus::Float64
    area_minus::Float64
    
    function RWGFunction(edge::Edge, mesh::Mesh3D)
        tri_plus = mesh.triangles[edge.triangle_plus]
        tri_minus = mesh.triangles[edge.triangle_minus]
        new(edge.triangle_plus, edge.triangle_plus, edge.triangle_minus, 
            edge.length, tri_plus.area, tri_minus.area)
    end
end

function evaluate_rwg(rwg::RWGFunction, point::SVector{3, Float64}, mesh::Mesh3D)
    edge = mesh.edges[rwg.edge_index]
    
    if point_in_triangle(point, mesh.triangles[rwg.triangle_plus])
        vertex_opposite = find_opposite_vertex(edge, mesh.triangles[rwg.triangle_plus])
        return (rwg.edge_length / (2.0 * rwg.area_plus)) * (point - vertex_opposite)
    elseif point_in_triangle(point, mesh.triangles[rwg.triangle_minus])
        vertex_opposite = find_opposite_vertex(edge, mesh.triangles[rwg.triangle_minus])
        return -(rwg.edge_length / (2.0 * rwg.area_minus)) * (point - vertex_opposite)
    else
        return SVector(0.0, 0.0, 0.0)
    end
end

function evaluate_rwg_divergence(rwg::RWGFunction)
    return rwg.edge_length / rwg.area_plus + rwg.edge_length / rwg.area_minus
end

function point_in_triangle(point::SVector{3, Float64}, triangle::Triangle)
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

function find_opposite_vertex(edge::Edge, triangle::Triangle)
    vertices = triangle.vertices
    for vertex in vertices
        if vertex != mesh.vertices[edge.vertex1] && vertex != mesh.vertices[edge.vertex2]
            return vertex
        end
    end
    error("Could not find opposite vertex")
end

end
