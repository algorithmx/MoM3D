module Geometry

using LinearAlgebra
using StaticArrays
using GeometryBasics

export Triangle, Edge, Mesh3D, find_edges, validate_mesh, analyze_mesh_quality, detect_t_junctions, repair_mesh_connectivity, compute_aspect_ratio, compute_min_angle, compute_max_angle

struct Triangle
    vertices::SVector{3,SVector{3,Float64}}
    area::Float64
    normal::SVector{3,Float64}

    function Triangle(v1::SVector{3,Float64}, v2::SVector{3,Float64}, v3::SVector{3,Float64})
        vertices = SVector(v1, v2, v3)
        edge1 = v2 - v1
        edge2 = v3 - v1
        normal = normalize(cross(edge1, edge2))
        area = 0.5 * norm(cross(edge1, edge2))
        new(vertices, area, normal)
    end
end

struct Edge
    triangle_plus::Int
    triangle_minus::Int
    vertex1::Int
    vertex2::Int
    length::Float64
    center::SVector{3,Float64}

    function Edge(tri_plus::Int, tri_minus::Int, v1_idx::Int, v2_idx::Int,
        v1::SVector{3,Float64}, v2::SVector{3,Float64})
        length = norm(v2 - v1)
        center = 0.5 * (v1 + v2)
        new(tri_plus, tri_minus, v1_idx, v2_idx, length, center)
    end
end

"""
Mesh3D

Fields:
- `vertices`, `triangles`, `edges`, ... : core mesh data
- `rwg_cache::Union{Nothing, Vector{Any}}` : lazy cache for RWG basis functions.

The `rwg_cache` field stores the vector returned by `get_rwgs(mesh)` and
is populated lazily. If you mutate the mesh in-place, call
`invalidate_rwg_cache!(mesh)` or use `update_mesh!` which clears the cache
automatically.
"""
mutable struct Mesh3D
    vertices::Vector{SVector{3,Float64}}
    triangles::Vector{Triangle}
    edges::Vector{Edge}
    boundary_edges::Vector{Tuple{Int,Int}}
    interior_edges::Vector{Int}
    num_edges::Int
    is_watertight::Bool
    mesh_quality::NamedTuple
    rwg_cache::Union{Nothing,Vector{Any}}
    integration_cache::Union{Nothing,Dict{Symbol,Any}}

    function Mesh3D(vertices::Vector{SVector{3,Float64}}, triangle_indices::Vector{SVector{3,Int}})
        triangles = [Triangle(vertices[idx[1]], vertices[idx[2]], vertices[idx[3]])
                     for idx in triangle_indices]

        edges, edge_map, node_conn, facet_conn = find_edges(vertices, triangle_indices)

        boundary_edges = Tuple{Int,Int}[]
        interior_edges = Int[]

        for ((v1, v2), facets) in edge_map
            if length(facets) == 1
                push!(boundary_edges, (v1, v2))
            elseif length(facets) == 2
                edge_idx = findfirst(e -> (e.vertex1 == v1 && e.vertex2 == v2), edges)
                if edge_idx !== nothing
                    push!(interior_edges, edge_idx)
                end
            end
        end

        is_watertight = length(boundary_edges) == 0
        mesh_quality = analyze_mesh_quality(triangles, vertices)

        new(vertices, triangles, edges, boundary_edges, interior_edges,
            length(edges), is_watertight, mesh_quality, nothing, nothing)
    end
end

"""
update_mesh!(mesh::Mesh3D, vertices, triangle_indices)

Update `mesh` in-place with new `vertices` and `triangle_indices`. This
recomputes triangles, edges and connectivity and automatically clears the
`rwg_cache` so subsequent calls to `get_rwgs(mesh)` will rebuild the
basis for the updated geometry.
"""
function update_mesh!(mesh::Mesh3D, vertices::Vector{SVector{3,Float64}}, triangle_indices::Vector{SVector{3,Int}})
    # Recompute triangles, edges and connectivity in-place and invalidate RWG cache
    triangles = [Triangle(vertices[idx[1]], vertices[idx[2]], vertices[idx[3]])
                 for idx in triangle_indices]
    edges, edge_map, node_conn, facet_conn = find_edges(vertices, triangle_indices)

    boundary_edges = Tuple{Int,Int}[]
    interior_edges = Int[]

    for ((v1, v2), facets) in edge_map
        if length(facets) == 1
            push!(boundary_edges, (v1, v2))
        elseif length(facets) == 2
            edge_idx = findfirst(e -> (e.vertex1 == v1 && e.vertex2 == v2), edges)
            if edge_idx !== nothing
                push!(interior_edges, edge_idx)
            end
        end
    end

    is_watertight = length(boundary_edges) == 0
    mesh_quality = analyze_mesh_quality(triangles, vertices)

    mesh.vertices = vertices
    mesh.triangles = triangles
    mesh.edges = edges
    mesh.boundary_edges = boundary_edges
    mesh.interior_edges = interior_edges
    mesh.num_edges = length(edges)
    mesh.is_watertight = is_watertight
    mesh.mesh_quality = mesh_quality
    mesh.rwg_cache = nothing
    mesh.integration_cache = nothing

    return mesh
end

function find_edges(vertices::Vector{SVector{3,Float64}}, triangle_indices::Vector{SVector{3,Int}})
    n_vertices = length(vertices)
    n_triangles = length(triangle_indices)

    node_connectivity = [Int[] for _ in 1:n_vertices]
    facet_connectivity = [Int[] for _ in 1:n_vertices]

    for (tri_idx, tri) in enumerate(triangle_indices)
        sorted_nodes = sort([tri[1], tri[2], tri[3]])

        for node in sorted_nodes
            push!(facet_connectivity[node], tri_idx)
        end

        if !(sorted_nodes[2] in node_connectivity[sorted_nodes[1]])
            push!(node_connectivity[sorted_nodes[1]], sorted_nodes[2])
        end
        if !(sorted_nodes[3] in node_connectivity[sorted_nodes[1]])
            push!(node_connectivity[sorted_nodes[1]], sorted_nodes[3])
        end
        if !(sorted_nodes[3] in node_connectivity[sorted_nodes[2]])
            push!(node_connectivity[sorted_nodes[2]], sorted_nodes[3])
        end
    end

    edges = Edge[]
    edge_map = Dict{Tuple{Int,Int},Vector{Int}}()

    for node1 in 1:n_vertices
        for node2 in node_connectivity[node1]
            common_facets = intersect(facet_connectivity[node1], facet_connectivity[node2])
            edge_key = (node1, node2)
            edge_map[edge_key] = common_facets

            if length(common_facets) == 2
                edge = Edge(common_facets[1], common_facets[2], node1, node2,
                    vertices[node1], vertices[node2])
                push!(edges, edge)
            elseif length(common_facets) > 2
                @warn "Edge ($node1, $node2) touches $(length(common_facets)) triangles - non-manifold geometry detected"
            end
        end
    end

    return edges, edge_map, node_connectivity, facet_connectivity
end

function analyze_mesh_quality(triangles::Vector{Triangle}, vertices::Vector{SVector{3,Float64}})
    if isempty(triangles)
        return (min_area=0.0, max_area=0.0, min_aspect_ratio=0.0, max_aspect_ratio=0.0,
            mean_aspect_ratio=0.0, min_angle_deg=0.0, max_angle_deg=0.0,
            degenerate_triangles=0, thin_triangles=0)
    end

    areas = [tri.area for tri in triangles]
    aspect_ratios = [compute_aspect_ratio(tri) for tri in triangles]
    min_angles = [compute_min_angle(tri) for tri in triangles]
    max_angles = [compute_max_angle(tri) for tri in triangles]

    degenerate_count = count(area -> area < 1e-12, areas)
    thin_count = count(ratio -> ratio > 10.0, aspect_ratios)

    return (
        min_area=minimum(areas),
        max_area=maximum(areas),
        min_aspect_ratio=minimum(aspect_ratios),
        max_aspect_ratio=maximum(aspect_ratios),
        mean_aspect_ratio=sum(aspect_ratios) / length(aspect_ratios),
        min_angle_deg=minimum(min_angles) * 180.0 / π,
        max_angle_deg=maximum(max_angles) * 180.0 / π,
        degenerate_triangles=degenerate_count,
        thin_triangles=thin_count
    )
end

function validate_mesh(mesh::Mesh3D; aspect_ratio_threshold::Float64=10.0,
    min_angle_threshold::Float64=5.0)
    issues = String[]

    if !mesh.is_watertight
        push!(issues, "Mesh is not watertight: $(length(mesh.boundary_edges)) boundary edges found")
    end

    if mesh.mesh_quality.degenerate_triangles > 0
        push!(issues, "$(mesh.mesh_quality.degenerate_triangles) degenerate triangles with area < 1e-12")
    end

    if mesh.mesh_quality.max_aspect_ratio > aspect_ratio_threshold
        push!(issues, "Poor aspect ratios detected: max = $(mesh.mesh_quality.max_aspect_ratio)")
    end

    if mesh.mesh_quality.min_angle_deg < min_angle_threshold
        push!(issues, "Small angles detected: min = $(mesh.mesh_quality.min_angle_deg)°")
    end

    if mesh.mesh_quality.thin_triangles > 0
        push!(issues, "$(mesh.mesh_quality.thin_triangles) thin triangles with aspect ratio > 10")
    end

    t_junctions = detect_t_junctions(mesh)
    if !isempty(t_junctions)
        push!(issues, "$(length(t_junctions)) T-junctions detected")
    end

    return (
        is_valid=isempty(issues),
        issues=issues,
        watertight=mesh.is_watertight,
        quality_metrics=mesh.mesh_quality,
        t_junctions=t_junctions
    )
end

function compute_aspect_ratio(triangle::Triangle)
    v1, v2, v3 = triangle.vertices
    edge_lengths = [norm(v2 - v1), norm(v3 - v2), norm(v1 - v3)]
    longest_edge = maximum(edge_lengths)
    shortest_edge = minimum(edge_lengths)

    if shortest_edge < 1e-12
        return Inf
    end

    return longest_edge / shortest_edge
end

function compute_min_angle(triangle::Triangle)
    v1, v2, v3 = triangle.vertices

    a = norm(v3 - v2)  # opposite to v1
    b = norm(v3 - v1)  # opposite to v2  
    c = norm(v2 - v1)  # opposite to v3

    if a < 1e-12 || b < 1e-12 || c < 1e-12
        return 0.0
    end

    angle1 = acos(clamp((b^2 + c^2 - a^2) / (2 * b * c), -1.0, 1.0))
    angle2 = acos(clamp((a^2 + c^2 - b^2) / (2 * a * c), -1.0, 1.0))
    angle3 = acos(clamp((a^2 + b^2 - c^2) / (2 * a * b), -1.0, 1.0))

    return minimum([angle1, angle2, angle3])
end

function compute_max_angle(triangle::Triangle)
    v1, v2, v3 = triangle.vertices

    a = norm(v3 - v2)  # opposite to v1
    b = norm(v3 - v1)  # opposite to v2  
    c = norm(v2 - v1)  # opposite to v3

    if a < 1e-12 || b < 1e-12 || c < 1e-12
        return π
    end

    angle1 = acos(clamp((b^2 + c^2 - a^2) / (2 * b * c), -1.0, 1.0))
    angle2 = acos(clamp((a^2 + c^2 - b^2) / (2 * a * c), -1.0, 1.0))
    angle3 = acos(clamp((a^2 + b^2 - c^2) / (2 * a * b), -1.0, 1.0))

    return maximum([angle1, angle2, angle3])
end

function detect_t_junctions(mesh::Mesh3D; tolerance::Float64=1e-10)
    t_junctions = Tuple{Int,Int,Int}[]  # (triangle_idx, edge_vertex1, edge_vertex2)

    for (tri_idx, triangle) in enumerate(mesh.triangles)
        vertices = triangle.vertices

        for i in 1:3
            v1 = vertices[i]
            v2 = vertices[i%3+1]
            edge_vec = v2 - v1
            edge_length = norm(edge_vec)

            if edge_length < tolerance
                continue
            end

            for (other_idx, other_vertex) in enumerate(mesh.vertices)
                if other_idx == tri_idx
                    continue
                end

                if norm(other_vertex - v1) < tolerance || norm(other_vertex - v2) < tolerance
                    continue
                end

                t = dot(other_vertex - v1, edge_vec) / (edge_length^2)

                if 0.0 < t < 1.0
                    closest_point = v1 + t * edge_vec
                    distance = norm(other_vertex - closest_point)

                    if distance < tolerance
                        push!(t_junctions, (tri_idx, i, i % 3 + 1))
                    end
                end
            end
        end
    end

    return t_junctions
end

function repair_mesh_connectivity(vertices::Vector{SVector{3,Float64}},
    triangle_indices::Vector{SVector{3,Int}};
    tolerance::Float64=1e-10)
    vertex_map = Dict{Int,Int}()
    unique_vertices = SVector{3,Float64}[]

    for (i, vertex) in enumerate(vertices)
        merged = false
        for (j, unique_vertex) in enumerate(unique_vertices)
            if norm(vertex - unique_vertex) < tolerance
                vertex_map[i] = j
                merged = true
                break
            end
        end

        if !merged
            push!(unique_vertices, vertex)
            vertex_map[i] = length(unique_vertices)
        end
    end

    repaired_triangles = SVector{3,Int}[]
    for tri in triangle_indices
        new_tri = SVector(vertex_map[tri[1]], vertex_map[tri[2]], vertex_map[tri[3]])

        if new_tri[1] != new_tri[2] && new_tri[2] != new_tri[3] && new_tri[1] != new_tri[3]
            push!(repaired_triangles, new_tri)
        end
    end

    return unique_vertices, repaired_triangles
end

end
