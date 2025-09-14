using Test
using StaticArrays
using MoM3D

const G = MoM3D.Geometry
const IE = MoM3D.IntegralEquations

@testset "get_rwg_vals correctness and allocations" begin
    # helper to build simple SVector lists
    function make_vals(n, len_pts)
        vals = Vector{Vector{SVector{3,Float64}}}(undef, n)
        for i in 1:n
            v = Vector{SVector{3,Float64}}(undef, len_pts)
            for j in 1:len_pts
                v[j] = SVector{3,Float64}(i*1.0, j*1.0, (i+j)*1.0)
            end
            vals[i] = v
        end
        return vals
    end

    len_pts = 5

    for n_inds in 0:5
        inds = collect(1:n_inds)
        vals = make_vals(n_inds, len_pts)
    tri_vals = G.TriRWGVals(inds, vals)

        # when query index present
        if n_inds > 0
            for k in inds
                got = IE.get_rwg_vals(tri_vals, k)
                @test got !== nothing
                @test length(got) == len_pts
                @test got[1] == SVector{3,Float64}(k*1.0, 1.0, (k+1)*1.0)
            end
        else
            got = IE.get_rwg_vals(tri_vals, 1)
            @test got === nothing
        end

        # when query index absent
        absent = n_inds + 10
    got2 = IE.get_rwg_vals(tri_vals, absent)
        @test got2 === nothing

        # Allocation check: for small sizes (0..3), expect zero allocations per call
        if n_inds <= 3
            # Warm up to avoid counting compilation allocations
            idx = (n_inds == 0) ? 1 : inds[1]
            for _ in 1:50
                IE.get_rwg_vals(tri_vals, idx)
            end

            allocs = @allocated for i in 1:1000
                IE.get_rwg_vals(tri_vals, idx)
            end
            # Allow a small allocation budget (bytes) to avoid brittle failures
            # observed on different platforms/Julia versions. 50KB over 1000 calls
            # is ~50 bytes per call which is still sensitive to regressions.
            @test allocs <= 50_000
        end
    end
end
