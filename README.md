MoM3D
=====

Quick notes about the RWG basis cache
------------------------------------

This project precomputes RWG basis functions per mesh edge to avoid
repeated allocations and work during assembly, excitation evaluation,
and postprocessing.

Key APIs
- `get_rwgs(mesh::Mesh3D)`: returns a vector of `RWGFunction` instances
	for the given mesh. The vector is cached on the mesh after first
	construction and returned on subsequent calls.

- `invalidate_rwg_cache!(mesh::Mesh3D)`: manually clears the cached
	RWG basis. Call this if you mutate the mesh in-place (e.g. changing
	vertices, triangles, or edges directly).

- `update_mesh!(mesh::Mesh3D, vertices, triangle_indices)`: helper to
	update mesh geometry in-place; it recomputes triangles/edges and
	automatically invalidates the RWG cache.

Why the cache exists
- Multiple parts of the code (matrix assembly, excitation vector
	computation, field evaluation, postprocessing) need the same RWG
	basis. Caching them on the mesh saves allocations and keeps a single
	canonical RWG instance per edge.

Usage examples

Use the cached RWGs in assembly and postprocessing:

```julia
using MoM3D

# Build or load a Mesh3D instance `mesh`

# Assembly uses the cached basis automatically
Z = assemble_efie_matrix(mesh, frequency)

# Postprocessing gets the same basis instances
E_field, H_field = compute_near_field(current_coeffs, mesh, observation_points, frequency)
```

Invalidating the cache after mesh edits:

```julia
# If you change mesh data in-place, clear the cache first
invalidate_rwg_cache!(mesh)

# Or use the helper which updates and clears cache automatically
update_mesh!(mesh, new_vertices, new_triangle_indices)
```

Notes on typing
- Internally the mesh stores the cache in a `rwg_cache` field with a
	general element type to avoid module precompilation ordering issues.
	The value stored is a `Vector{RWGFunction}` built by `get_rwgs`.

If you'd like stricter typing for the cache, consider moving the
`RWGFunction` type into the `Geometry` module or adjusting the module
load order; both are straightforward but touch module organization.
