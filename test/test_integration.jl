using Test
using MoM3D
using StaticArrays

@testset "Integration Tests" begin
    using MoM3D.Integration: _fallback_gauss_triangle, gauss_triangle, barycentric_to_cartesian, integrate_regular, integrate_singular

    # Helper: integrate monomial x^i y^j over unit triangle (barycentric coords mapped)
    function integrate_monomial_over_triangle(i, j)
        # Integral over unit triangle (vertices (0,0),(1,0),(0,1)) of x^i y^j = i! j! / ( (i+j+2)! )
        return factorial(i) * factorial(j) / factorial(i + j + 2)
    end

    @testset "Fallback rule correctness (low-order)" begin
        # Test that fallback 1-point integrates constant exactly
        pts, w = _fallback_gauss_triangle(1)
        @test length(pts) == 1
        @test isapprox(sum(w), 1.0; atol=1e-12)

        # 3-point rule should integrate linear polynomials exactly (degree 1)
    pts3, w3 = _fallback_gauss_triangle(3)
    # Convert barycentric (ξ1,ξ2) to cartesian on unit triangle with vertices
    # (0,0),(1,0),(0,1): x = ξ2, y = 1 - ξ1 - ξ2
    fx = sum(w3[i] * (pts3[i][2]) for i in eachindex(w3)) * 0.5
    fy = sum(w3[i] * (1.0 - pts3[i][1] - pts3[i][2]) for i in eachindex(w3)) * 0.5
    @test isapprox(fx, integrate_monomial_over_triangle(1,0); atol=1e-12)
    @test isapprox(fy, integrate_monomial_over_triangle(0,1); atol=1e-12)

        # 7-point Dunavant should integrate up to degree 5 polynomials; test degree 2
    pts7, w7 = _fallback_gauss_triangle(7)
        # test ∫ x^2, ∫ x*y, ∫ y^2
    # Convert barycentric to cartesian and include area factor (0.5)
    fx2 = sum(w7[i] * (pts7[i][2]^2) for i in eachindex(w7)) * 0.5
    fxy = sum(w7[i] * (pts7[i][2] * (1.0 - pts7[i][1] - pts7[i][2])) for i in eachindex(w7)) * 0.5
    fy2 = sum(w7[i] * ((1.0 - pts7[i][1] - pts7[i][2])^2) for i in eachindex(w7)) * 0.5
        @test isapprox(fx2, integrate_monomial_over_triangle(2,0); atol=1e-10)
        @test isapprox(fxy, integrate_monomial_over_triangle(1,1); atol=1e-10)
        @test isapprox(fy2, integrate_monomial_over_triangle(0,2); atol=1e-10)
    end

    @testset "gauss_triangle (SimplexQuad path) matches fallback for low orders" begin
        for n in (1,3,7)
            pts_lib, w_lib = gauss_triangle(n)
            pts_fb, w_fb = _fallback_gauss_triangle(n)
            # Convert both to cartesian averages (include area factor) and compare
            mx_lib = sum(w_lib[i] * pts_lib[i][1] for i in eachindex(w_lib)) * 0.5
            my_lib = sum(w_lib[i] * pts_lib[i][2] for i in eachindex(w_lib)) * 0.5
            mx_fb = sum(w_fb[i] * pts_fb[i][2] for i in eachindex(w_fb)) * 0.5
            my_fb = sum(w_fb[i] * (1.0 - pts_fb[i][1] - pts_fb[i][2]) for i in eachindex(w_fb)) * 0.5
            @test isapprox(mx_lib, mx_fb; atol=1e-8)
            @test isapprox(my_lib, my_fb; atol=1e-8)
        end
    end

    @testset "Quadrature weight normalization" begin
        for n in (1,3,4,5,7,9)
            pts, w = gauss_triangle(n)
            @test isapprox(sum(w), 1.0; atol=1e-12)
            ptsf, wf = _fallback_gauss_triangle(n)
            @test isapprox(sum(wf), 1.0; atol=1e-12)
        end
    end

    @testset "Harder integration tests (regular and singular)" begin
        # Construct a unit right triangle in 3D for Geometry.Triangle
        v1 = SVector(0.0, 0.0, 0.0)
        v2 = SVector(1.0, 0.0, 0.0)
        v3 = SVector(0.0, 1.0, 0.0)
    tri = MoM3D.Geometry.Triangle(v1, v2, v3)

        # Test integrate_singular with a smooth integrand (should equal area * value at centroid)
        integrand_constant = r -> 2.5 + 0.1* r[1] - 0.2 * r[2]
    res = integrate_singular(integrand_constant, tri, 7)
        # Analytical integral over triangle of integrand = area * average. For linear function average equals value at centroid (1/3,1/3)
        centroid = SVector(1/3,1/3,0.0)
        expected = tri.area * integrand_constant(centroid)
        @test isapprox(res, expected; atol=1e-6)

        # Test integrate_regular for two identical triangles and a separable kernel f(r_obs,r_src)=1
        tri_src = tri
        tri_obs = tri
        res2 = integrate_regular((ro, rs) -> 1.0, tri_src, tri_obs, quad_order=3)
        # Integral over source and obs of 1 = area_src * area_obs
        @test isapprox(res2, tri.area * tri.area; atol=1e-10)
    end
end
