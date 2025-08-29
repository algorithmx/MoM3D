# MoM3D – Method of Moments (MoM) Code Review Progress Report

Date: 2025-08-29
Scope: Review of src/ for current MoM implementation (no code changes).

## Executive Summary
- The project scaffolding is solid: clear modular separation (Geometry, BasisFunctions, GreenFunctions, Integration, IntegralEquations, Solvers, PostProcess) with a top-level MoM3D module wiring them together.
- Core numerical pieces are drafted, but several critical correctness and scoping issues will prevent the code from running as-is and/or producing valid results.
- Highest-priority blockers include: inconsistent/incorrect use of `mesh` in several functions, incorrect RWG basis constructor/usage, misuse of near-field/far-field integration (double vs single integration), undefined constants (μ0, ε0) due to export/scope, and an invalid near-singular integration strategy.
- MFIE/CFIE are declared in exports but not implemented.

Immediate recommendation: address the “Blockers” below before further feature work or performance tuning.

---

## Module-by-Module Review

### 1) src/geometry.jl (module Geometry)
Strengths:
- Triangle/Edge/Mesh3D types with derived geometry (area, normal, edge length/center).
- Mesh validation and quality metrics (aspect ratio, angles, degeneracy) and helper utilities (T-junction detection, basic connectivity repair).

Concerns/Issues:
- Edge orientation for RWG: `find_edges` sets `(triangle_plus, triangle_minus)` as `common_facets[1]`/`[2]` without ensuring orientation conforms to RWG conventions (plus minus associated to edge direction relative to local triangle orientation). This can flip signs downstream.
- `detect_t_junctions`: logic compares a vertex loop index (`other_idx`) to a triangle index (`tri_idx`), which are different kinds of indices; this likely yields false positives/negatives.
- Performance: O(N^2) scans (e.g., in `repair_mesh_connectivity`); acceptable for small meshes, but may need optimization for larger ones.

Status: Foundational but requires orientation correctness review and bug fixes in T-junction detection.

### 2) src/basis_functions.jl (module BasisFunctions)
Strengths:
- RWG basis structure and evaluation stubs (value and divergence) are present.

Blockers/Issues:
- RWGFunction constructor field mapping appears wrong:
  - Struct fields: `(edge_index, triangle_plus, triangle_minus, edge_length, area_plus, area_minus)`.
  - Constructor uses: `new(edge.triangle_plus, edge.triangle_plus, edge.triangle_minus, edge.length, ...)`.
  - This sets `edge_index` to `edge.triangle_plus` (a triangle index!), and sets `triangle_plus` to the same value. This will corrupt indexing and break `evaluate_rwg`, which does `edge = mesh.edges[rwg.edge_index]`.
- `find_opposite_vertex(edge, triangle)` references `mesh.vertices[...]` but `mesh` is not in scope (function signature lacks `mesh`). This is a hard runtime error.
- `evaluate_rwg_divergence(rwg)`: returns `edge_length/area_plus + edge_length/area_minus` (always positive), whereas RWG divergence is piecewise constant with opposite signs in adjacent triangles; treating it as a single global constant is not correct for integrands that depend on observation/source locations.

Status: Needs constructor and API fixes; divergence handling should be location-aware (or two separate per-triangle divergences).

### 3) src/green_functions.jl (module GreenFunctions)
Strengths:
- Free-space scalar Green’s function `green_3d` and gradient; wavenumber helper; singular-extraction helper (`green_3d_singular_extraction`).

Concerns/Issues:
- Constants `μ0`, `ε0` are defined but not exported. Other modules use them unqualified (e.g., `sqrt(μ0/ε0)`), which will fail unless fully-qualified (e.g., `GreenFunctions.μ0`).

Status: Functionality OK; ensure proper export or use qualified references in other modules.

### 4) src/integration.jl (module Integration)
Strengths:
- Triangle quadrature rules with weights normalized to 1.0; conversion `barycentric_to_cartesian`.
- Regular and singular integral routines for triangle integrals.

Blockers/Issues:
- `integrate_near_singular` uses `QuadGK.quadgk` with vector bounds to attempt multi-dimensional integration. `QuadGK` is 1D; passing vector bounds is invalid. A proper near-singular treatment should use specialized triangle subdivision or a true multi-D cubature (e.g., HCubature/Cubature) and geometry-aware mapping.
- `integrate_regular` performs a double integral (obs × src) and multiplies by `tri_src.area * tri_obs.area`. This is fine for true double integrals, but is later misused in PostProcess for single-surface integrals (see below).

Status: Regular/singular routines acceptable as drafts; near-singular routine must be reworked.

### 5) src/integral_equations.jl (module IntegralEquations)
Strengths:
- EFIE assembly skeleton with path-dependent treatment (overlap, adjacent, far) and singular-extraction option.

Blockers/Issues:
- Unscoped `mesh` usage: helper functions `compute_regular_efie_element`, `compute_singular_efie_element`, and `compute_near_singular_efie_element` use `mesh` inside integrands but do not accept it as a parameter; they are called without `mesh`. This is a hard runtime error.
- Reliance on `evaluate_rwg_divergence(rwg)` as a global constant (see BasisFunctions notes). Scalar term in EFIE should reflect piecewise divergence and correct testing/expansion.
- Adjacency/overlap detection: very coarse. `triangles_overlap` checks identity only; adjacency is based on coordinate equality within 1e-12; acceptable as a first pass but may be brittle for noisy meshes.
- Exports declare `assemble_mfie_matrix, assemble_cfie_matrix` but no implementations exist.
- Material constants: Uses `η = sqrt(μ0/ε0)` where `μ0, ε0` are not imported/exported correctly (see GreenFunctions). Will be undefined here.

Status: Conceptually aligned but currently not executable due to missing parameters/scope and missing MFIE/CFIE.

### 6) src/solvers.jl (module Solvers)
Strengths:
- Direct and iterative solvers (GMRES, BiCGSTAB(l)) via IterativeSolvers.jl; basic regularization fallback for direct solve; condition number estimate helper; diagonal preconditioner helper.

Concerns/Issues:
- Preconditioner is defined but not wired into iterative solves.
- No stopping/reporting abstraction (everything verbose=true) and no residual/history exposure.

Status: Usable baseline once upstream matrix assembly is fixed; could be enhanced with preconditioning and solver options.

### 7) src/postprocess.jl (module PostProcess)
Strengths:
- API stubs for near-field, far-field, RCS, and current density evaluation.

Blockers/Issues:
- Unscoped constants and functions: uses `wavenumber`, `μ0`, `ε0` unqualified; only `wavenumber` is exported; `μ0`, `ε0` are not—will be undefined.
- Near-field/magnetic-field integrals misuse double integral routine:
  - Calls `integrate_regular(integrand, tri_src, Triangle(r_obs, r_obs, r_obs), 3)`; the observation “triangle” is degenerate (area=0), so Jacobian product `tri_src.area * tri_obs.area` is zero → returns zero fields.
  - For these tasks, one needs a single-surface integral over source triangle with fixed observation point (or specialized kernels), not a double integral over obs × src.
- Far-field integrals similarly use `integrate_regular` with both triangles, effectively performing a double integral when only a single surface integral over `r_src` is intended.

Status: Non-functional as written; needs single-integral routines and kernel-specific formulations.

---

## Cross-Cutting Issues (Blockers)
1) Scope/Namespace bugs
   - Unqualified use of `μ0`, `ε0` in IntegralEquations and PostProcess (not exported). Use `GreenFunctions.μ0`, `GreenFunctions.ε0` or export them.
   - Multiple helpers reference `mesh` without it being passed (BasisFunctions.find_opposite_vertex; IntegralEquations compute_* helpers). These must accept `mesh` as a parameter.

2) RWG basis implementation
   - Constructor fills fields incorrectly; `edge_index` should refer to the index into `mesh.edges`, not a triangle index. Orientation and opposite vertex retrieval must be index-based and consistent with RWG definitions.
   - Divergence should be piecewise ±l/A per adjacent triangle, not a single global sum.

3) Near-/far-field integration design
   - Current code uses a double-integral routine with degenerate observation geometry, producing zero results. Need dedicated single-integral quadrature over the source triangle(s) with kernels evaluated at fixed observation directions/points.

4) Near-singular integration
   - `QuadGK` cannot be used for multi-D; near-singular handling needs proper subdivision or multi-D cubature with adapted transforms.

5) Missing features
   - MFIE and CFIE assembly functions are exported but not implemented.
   - No tests or example scripts exercising the pipeline end-to-end from mesh → Z, b → solve → fields/RCS.

---

## What Appears Complete/Usable After Fixes
- Geometry/Mesh: suitable for small problems once T-junction detection check is corrected and RWG orientation strategy is verified.
- Green’s functions: scalar/gradient routines usable.
- Regular/singular integrals: usable as building blocks; near-singular requires redesign.
- Solvers: acceptable baseline; can solve once Z, b are valid.

---

## Prioritized Next Steps (no code changes done yet)
1) Correct scoping/namespace
   - Qualify or export `μ0`, `ε0` and pass `mesh` explicitly where needed.
2) Fix RWG basis
   - Correct RWGFunction constructor field mapping; implement `find_opposite_vertex` using vertex indices (and provide it `mesh`).
   - Make divergence evaluation piecewise-aware (or provide two functions for plus/minus triangles).
3) Make EFIE assembly executable
   - Add `mesh` parameter to compute_* EFIE element helpers; review the scalar term formulation using correct divergence handling.
   - Validate overlap/adjacency logic; add robust shared-edge/vertex checks.
4) Redesign near-/far-field evaluation
   - Implement single-surface quadrature routines for fixed observation points/directions; remove degenerate-triangle pattern.
5) Replace near-singular integration
   - Use triangle subdivision strategies or a proper multi-D cubature package; consider Duffy transform for weakly singular integrals.
6) Implement MFIE/CFIE (optional after EFIE is stable)
   - Provide assembly for MFIE and combined-field (with coupling parameter α) and validate on canonical scatterers.
7) Add tests and example validations
   - Unit tests for geometry and RWG evaluation; integration tests assembling small Z for simple meshes; comparison against known results (e.g., sphere/Mie for RCS).

---

## Risks and Validation Plan
- Risk: Orientation and divergence mistakes lead to non-physical impedance matrices (loss of symmetry properties, poor conditioning).
- Risk: Inadequate near-singular treatment will cause large discretization errors for touching/adjacent triangles.
- Validation:
  - Mesh quality checks on input meshes; refine and repair before solve.
  - Assemble EFIE for a tiny mesh (2 triangles) and compare against analytical/closed-form expectations where possible.
  - End-to-end sphere scattering vs. Mie series once EFIE works; perform mesh refinement study.

---

## Closing
The codebase is a good foundation but currently non-executable for key workflows due to several correctness and API issues. Addressing the listed blockers will unlock assembly and solution phases, after which postprocessing and performance work can proceed.

