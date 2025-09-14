module Geometry

using LinearAlgebra
using StaticArrays
using GeometryBasics

export Triangle, Edge, Mesh3D, find_edges, update_mesh!, validate_mesh, analyze_mesh_quality, detect_t_junctions, repair_mesh_connectivity, compute_aspect_ratio, compute_min_angle, compute_max_angle, point_in_triangle, find_opposite_vertex, find_opposite_vertex_index, recompute_edge_opposites!

# Typed integration cache used by IntegralEquations to store precomputed quadrature points and RWG values
struct IntegrationCache
    orders::Vector{Int}
    order_index::Dict{Int,Int}
    gauss_rules::Vector{Tuple{Vector{SVector{2,Float64}},Vector{Float64}}}
    tri_cart_pts::Vector{Vector{Vector{SVector{3,Float64}}}}
    tri_rwg_vals::Vector{Vector{Dict{Int,Vector{SVector{3,Float64}}}}}
end

struct Triangle
    vertices::SVector{3,SVector{3,Float64}}
    area::Float64
    normal::SVector{3,Float64}
    # Precomputed for fast point-in-triangle tests
    e0::SVector{3,Float64}    # v3 - v1
    e1::SVector{3,Float64}    # v2 - v1
    dot00::Float64
    dot01::Float64
    dot11::Float64
    inv_denom::Float64

    function Triangle(v1::SVector{3,Float64}, v2::SVector{3,Float64}, v3::SVector{3,Float64})
        vertices = SVector(v1, v2, v3)
        edge1 = v2 - v1
        edge2 = v3 - v1
        normal = normalize(cross(edge1, edge2))
        area = 0.5 * norm(cross(edge1, edge2))

        e0 = v3 - v1
        e1 = v2 - v1
        dot00 = dot(e0, e0)
        dot01 = dot(e0, e1)
        dot11 = dot(e1, e1)
        denom = dot00 * dot11 - dot01 * dot01
        inv_denom = abs(denom) < 1e-16 ? 0.0 : 1.0 / denom

        new(vertices, area, normal, e0, e1, dot00, dot01, dot11, inv_denom)
    end
end

struct Edge
    triangle_plus::Int
    triangle_minus::Int
    vertex1::Int
    vertex2::Int
    length::Float64
    center::SVector{3,Float64}
    opp_plus::Int
    opp_minus::Int

    function Edge(tri_plus::Int, tri_minus::Int, v1_idx::Int, v2_idx::Int,
        v1::SVector{3,Float64}, v2::SVector{3,Float64}, opp_plus::Int, opp_minus::Int)
        length = norm(v2 - v1)
        center = 0.5 * (v1 + v2)
        new(tri_plus, tri_minus, v1_idx, v2_idx, length, center, opp_plus, opp_minus)
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
    integration_cache::Union{Nothing,IntegrationCache}

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
                t1, t2 = common_facets[1], common_facets[2]
                tri1 = triangle_indices[t1]
                tri2 = triangle_indices[t2]

                # find the opposite vertex index in each triangle (the vertex that is not node1 or node2)
                opp1 = (tri1[1] != node1 && tri1[1] != node2) ? tri1[1] : ((tri1[2] != node1 && tri1[2] != node2) ? tri1[2] : tri1[3])
                opp2 = (tri2[1] != node1 && tri2[1] != node2) ? tri2[1] : ((tri2[2] != node1 && tri2[2] != node2) ? tri2[2] : tri2[3])

                edge = Edge(t1, t2, node1, node2, vertices[node1], vertices[node2], opp1, opp2)
                push!(edges, edge)
            elseif length(common_facets) > 2
                @warn "Edge ($node1, $node2) touches $(length(common_facets)) triangles - non-manifold geometry detected"
            end
        end
    end

    return edges, edge_map, node_connectivity, facet_connectivity
end

"""
    recompute_edge_opposites!(mesh::Mesh3D)

Recompute the `opp_plus` and `opp_minus` fields for every `Edge` in `mesh`
based on the current `mesh.triangles` and `mesh.vertices`. This is useful
when connectivity/triangles have changed in-place but the `edges` vector
was not rebuilt via `update_mesh!`.

This function will also call `invalidate_rwg_cache!(mesh)` to ensure any
cached RWGFunctions observe the new opposite indices.
"""
function recompute_edge_opposites!(mesh::Mesh3D)
    # Build a fast map from vertex coordinate -> index for O(1) lookup
    vmap = Dict{SVector{3,Float64},Int}()
    for (i, v) in enumerate(mesh.vertices)
        vmap[v] = i
    end

    # Build edge -> facet list using the vertex index map
    edge_to_facets = Dict{Tuple{Int,Int},Vector{Int}}()
    for (tidx, tri) in enumerate(mesh.triangles)
        tri_nodes = Int[]
        for v in tri.vertices
            idx = get(vmap, v, nothing)
            if idx === nothing
                # Fallback to tolerant search if direct equality failed
                idx = findfirst(x -> isapprox(x, v; atol=1e-12), mesh.vertices)
                if idx === nothing
                    error("Vertex in triangle not found in mesh.vertices while recomputing opposites")
                end
            end
            push!(tri_nodes, idx)
        end

        pairs = ((min(tri_nodes[1], tri_nodes[2]), max(tri_nodes[1], tri_nodes[2])),
            (min(tri_nodes[1], tri_nodes[3]), max(tri_nodes[1], tri_nodes[3])),
            (min(tri_nodes[2], tri_nodes[3]), max(tri_nodes[2], tri_nodes[3])))

        for p in pairs
            push!(get!(edge_to_facets, p, Int[]), tidx)
        end
    end

    # Update each Edge's opposite indices using the prebuilt facet lists
    for (i, edge) in enumerate(mesh.edges)
        v1, v2 = edge.vertex1, edge.vertex2
        key = (min(v1, v2), max(v1, v2))
        facets = get(edge_to_facets, key, Int[])
        if length(facets) == 2
            t1, t2 = facets[1], facets[2]

            # Get triangle node indices efficiently using the vmap
            tri1 = mesh.triangles[t1]
            tri2 = mesh.triangles[t2]
            tri1_nodes = [get(vmap, v, findfirst(x -> isapprox(x, v; atol=1e-12), mesh.vertices)) for v in tri1.vertices]
            tri2_nodes = [get(vmap, v, findfirst(x -> isapprox(x, v; atol=1e-12), mesh.vertices)) for v in tri2.vertices]

            opp1 = (tri1_nodes[1] != v1 && tri1_nodes[1] != v2) ? tri1_nodes[1] : ((tri1_nodes[2] != v1 && tri1_nodes[2] != v2) ? tri1_nodes[2] : tri1_nodes[3])
            opp2 = (tri2_nodes[1] != v1 && tri2_nodes[1] != v2) ? tri2_nodes[1] : ((tri2_nodes[2] != v1 && tri2_nodes[2] != v2) ? tri2_nodes[2] : tri2_nodes[3])

            new_opp_plus = edge.triangle_plus == t1 ? opp1 : (edge.triangle_plus == t2 ? opp2 : 0)
            new_opp_minus = edge.triangle_minus == t1 ? opp1 : (edge.triangle_minus == t2 ? opp2 : 0)

            mesh.edges[i] = Edge(edge.triangle_plus, edge.triangle_minus, v1, v2, mesh.vertices[v1], mesh.vertices[v2], new_opp_plus, new_opp_minus)
        end
    end

    # Invalidate RWG cache so cached RWGFunctions pick up new opp indices
    invalidate_rwg_cache!(mesh)

    return nothing
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

"""
        point_in_triangle(point::SVector{3,Float64}, triangle::Triangle) -> Bool

Check whether a 3D point lies inside a triangle using barycentric
coordinates computed in the triangle's plane.

Arguments
- `point` : An `SVector{3,Float64}` representing the point to test.
- `triangle` : A `Triangle` containing three vertices (SVector{3,Float64}).

Returns
- `Bool` : `true` when the projected barycentric coordinates are inside
    the triangle (including edges), `false` otherwise.

Notes
- This routine computes barycentric coordinates in the plane defined by
    the triangle; the caller should ensure the point is intended to be
    tested against that plane (or accept the implicit projection). Small
    numerical tolerances are not handled here — callers may want to apply
    an epsilon when testing points near edges.
"""
function point_in_triangle(point::SVector{3,Float64}, triangle::Triangle; tol::Float64=1e-8)
    # Unpack triangle vertices and precomputed normal
    v1, v2, v3 = triangle.vertices
    n = triangle.normal

    # Quick reject: check distance from point to triangle plane
    # signed distance = dot(n, point - v1)
    dist = dot(n, point - v1)
    if abs(dist) > tol
        return false
    end

    # Project point onto plane (robust for near-coplanar points)
    p_proj = point - dist * n

    # Use precomputed edge vectors and dot-products stored in Triangle
    v2v = p_proj - v1
    dot02 = dot(triangle.e0, v2v)
    dot12 = dot(triangle.e1, v2v)

    if triangle.inv_denom == 0.0
        # Degenerate triangle
        return false
    end

    u = (triangle.dot11 * dot02 - triangle.dot01 * dot12) * triangle.inv_denom
    v = (triangle.dot00 * dot12 - triangle.dot01 * dot02) * triangle.inv_denom

    eps = -tol
    return (u >= eps) && (v >= eps) && (u + v <= 1.0 - eps)
end


"""
        find_opposite_vertex(edge::Edge, triangle::Triangle, mesh::Mesh3D) -> SVector{3,Float64}

Return the vertex coordinate of `triangle` that is opposite the given `edge`.

Arguments
- `edge::Edge` : an `Edge` whose `vertex1` and `vertex2` are indices into `mesh.vertices`.
- `triangle::Triangle` : a triangle storing three vertex coordinates (`SVector{3,Float64}`).
- `mesh::Mesh3D` : the containing mesh used to look up the coordinates of the edge endpoints.

Returns
- `SVector{3,Float64}` : the coordinate of the triangle vertex that is not equal (within a small tolerance) to
    either `mesh.vertices[edge.vertex1]` or `mesh.vertices[edge.vertex2]`.

Notes
- Comparison is done by coordinate distance with a tiny tolerance (1e-12) to avoid strict floating-point equality.
- The function returns the coordinate vector, not the vertex index; if an index is required prefer comparing
    indices or add a variant that returns the index as well.
- If no distinct opposite vertex is found (e.g. degenerate triangle or inconsistent mesh), the function throws an error.
"""
function find_opposite_vertex(edge::Edge, triangle::Triangle, mesh::Mesh3D)
    # Rely on precomputed opposite indices stored on the edge. The triangle
    # argument must be one of the two triangles adjacent to the edge.
    if edge.triangle_plus > 0 && triangle === mesh.triangles[edge.triangle_plus]
        opp = edge.opp_plus
    elseif edge.triangle_minus > 0 && triangle === mesh.triangles[edge.triangle_minus]
        opp = edge.opp_minus
    else
        error("Provided triangle is not adjacent to the edge")
    end

    if opp == 0
        error("Opposite vertex index not set for this edge side")
    end

    return mesh.vertices[opp]
end

@inline dist2(a::SVector{3,Float64}, b::SVector{3,Float64}) = dot(a - b, a - b)


"""
    find_opposite_vertex_index(edge::Edge; side::Symbol = :plus) -> Int

Return the precomputed opposite vertex index for `edge`.

Arguments
- `edge::Edge` : the edge with `opp_plus` and `opp_minus` fields set during `find_edges`.
- `side::Symbol` : either `:plus` or `:minus` to select the corresponding opposite vertex.

Returns
- `Int` : the vertex index of the opposite vertex for the requested side (0 if missing).
"""
function find_opposite_vertex_index(edge::Edge; side::Symbol=:plus)
    if side === :plus
        return edge.opp_plus
    elseif side === :minus
        return edge.opp_minus
    else
        error("side must be :plus or :minus")
    end
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
