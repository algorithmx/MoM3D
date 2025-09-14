Integration cache design — MoM3D
================================

Goal
----
Provide a typed, compact, and thread-safe cache for quadrature rules, precomputed cartesian quadrature points per triangle, and pre-evaluated RWG basis values per triangle. The cache should minimize runtime hash lookups and allocations during matrix assembly and make common code paths allocation-light and fast.

Location
--------
- Implementation: `src/integral_equations.jl` (`ensure_integration_caches`, consumers) and `src/geometry.jl` (typed cache and `TriRWGVals` definitions).
- Tests: `test/test_integration_cache.jl` and `test/test_get_rwg_vals.jl`

Core types
----------
- Geometry.IntegrationCache
  - orders::Vector{Int}
  - order_index::Dict{Int,Int}  # maps quadrature order -> small integer index (1..n_orders)
  - gauss_rules::Vector{Tuple{Vector{SVector{2,Float64}},Vector{Float64}}}  # per-order barycentric pts and weights
  - tri_cart_pts::Vector{Vector{Vector{SVector{3,Float64}}}}  # tri_cart_pts[oi][tri_idx] -> Vector of cartesian points
  - tri_rwg_vals::Vector{Vector{TriRWGVals}}  # tri_rwg_vals[oi][tri_idx] -> TriRWGVals for that triangle and order

- Geometry.TriRWGVals
  - rwg_indices::Vector{Int}  # global RWG edge indices adjacent to this triangle
  - rwg_values::Vector{Vector{SVector{3,Float64}}}  # aligned pre-evaluated RWG vectors for each index

Design rationale
----------------
- Use an `orders` vector + `order_index` map to translate sparse quadrature order integers (3,5,7,9...) into dense small integer indices. This avoids repeated dictionary lookups keyed by order across inner loops.
- Store per-order arrays (`gauss_rules`, `tri_cart_pts`, `tri_rwg_vals`) indexed by the small integer `oi` for fast array indexing.
- Compact per-triangle RWG storage (`TriRWGVals`) stores only the RWG indices adjacent to the triangle and aligned vectors. This removes the need for a global Dict lookup and reduces memory overhead and indirection compared with per-triangle Dicts.
- Accessor `get_rwg_vals(tri_vals::TriRWGVals, edge_index::Int)` provides a tiny fast-path for the very common adjacency sizes (0..3) and falls back to a small linear search for larger adjacency counts.

Concurrency
-----------
- Building and storing the cache on a `Mesh3D` is guarded by a module-level `ReentrantLock` (`INTEGRATION_CACHE_LOCK`) and follows double-checked locking:
  1. If `mesh.integration_cache !== nothing`, return it.
  2. Lock, check again, build the cache and store it on `mesh.integration_cache`.
  3. Unlock and return.
- Consumers read the cache without locks (the cache is immutable after construction) so reads are lock-free and fast.

Cache invalidation
------------------
- Any mesh mutation that changes triangle connectivity, vertex positions, or edge opposites must set `mesh.integration_cache = nothing` so the next assembly rebuilds the cache. Current places that invalidate caches include `update_mesh!`, `recompute_edge_opposites!` and other mesh-modifying routines.

APIs
----
- `ensure_integration_caches(mesh::Mesh3D; orders=(3,5,7,9)) -> IntegrationCache`
  - Build if missing, otherwise return the cached object.
  - Precomputes *per-order* Gauss rules, per-triangle cartesian quadrature points (barycentric->cartesian hoisted), and per-triangle RWG evaluations for each adjacent RWG.

- `get_rwg_vals(tri_vals::TriRWGVals, edge_index::Int) -> Union{Vector{SVector{3,Float64}}, Nothing}`
  - Fast-path for rwg index counts 0..3 (direct comparisons), fallback small linear search for larger counts.

Testing
-------
- `test/test_integration_cache.jl` verifies cache shapes, per-triangle cached RWG values equal direct `evaluate_rwg` evaluations, and numeric parity between cached and direct integrals for a pair of triangles.
- `test/test_get_rwg_vals.jl` verifies correctness across adjacency sizes 0..5 and performs an allocation-sensitive check for small sizes. The allocation check is warmed-up and allows a small budget to avoid brittle failures across platforms.

Performance notes & next steps
------------------------------
- The current `get_rwg_vals` micro-optimization reduces the common-case looping overhead for triangles with 0–3 adjacent RWGs. This covers the large majority of triangular meshes where edges are shared by 2 triangles (so adjacency per triangle is typically small).
- If more performance is needed:
  - Replace `rwg_indices::Vector{Int}` and `rwg_values::Vector{Vector{SVector{3,Float64}}}` with small fixed-size tuple storage for common sizes (0/1/2/3) to eliminate vector allocation and indirection. This increases code complexity but can further reduce allocations.
  - Add micro-benchmarks (BenchmarkTools.jl) to measure per-element assembly time and allocations on representative meshes.
  - Implement a specialized `TriRWGVals` layout where very small adjacency cases inline values in the struct (Union{Empty, FixedSize2, FixedSize3}) to avoid heap allocations for the `rwg_values` outer vector.
  - Profile assembly with `@profile` and `ProfileView`/`TraceView` to find remaining hotspots.

Notes
-----
- The cache is intentionally not backward-compatible with older dict-based caches; all call sites were migrated.
- The tests include tolerances for numerical extraction differences; small differences may still appear due to quadrature choices.

Contact
-------
For questions or to request the micro-benchmarking patch, open an issue or ask for an implementation in the `augment-improve` branch.
