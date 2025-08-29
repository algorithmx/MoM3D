using Test
using MoM3D

println("Running MoM3D Test Suite")
println("=" ^ 40)

@testset "MoM3D.jl Tests" begin
    include("test_geometry.jl")
    include("test_theoretical_validation.jl")
    include("test_solvers.jl")
end

println("\nAll tests completed!")
