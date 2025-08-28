module MoM3D

using LinearAlgebra
using StaticArrays

include("geometry.jl")
include("basis_functions.jl")
include("green_functions.jl")
include("integration.jl")
include("integral_equations.jl")
include("solvers.jl")
include("postprocess.jl")

using .Geometry
using .BasisFunctions
using .GreenFunctions
using .Integration
using .IntegralEquations
using .Solvers
using .PostProcess

export Mesh3D, Triangle, Edge, RWGFunction
export green_3d, wavenumber
export gauss_triangle, integrate_regular, integrate_singular
export assemble_efie_matrix, compute_excitation_vector
export solve_mom_system
export compute_near_field, compute_far_field, compute_rcs
export validate_mesh, analyze_mesh_quality, detect_t_junctions, repair_mesh_connectivity
export compute_aspect_ratio, compute_min_angle, compute_max_angle

end # module MoM3D
