using Test
using MoM3D

println("Running MoM3D Test Suite")
println("=" ^ 40)

@testset "MoM3D.jl Tests" begin
    include("test_geometry.jl")
    include("test_edge_opposite_indices.jl")
    include("test_opposite_vertex_indices_rigorous.jl")
    include("test_theoretical_validation.jl")
    include("test_solvers.jl")
    include("test_green_functions.jl")
    include("test_integration.jl")
    include("test_integration_cache.jl")
end

println("\nAll tests completed!")
