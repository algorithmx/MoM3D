using Pkg
Pkg.activate(".")
using MoM3D
using MoM3D.Integration: gauss_triangle, _fallback_gauss_triangle

for n in (1,3,4,5,7,9)
    pts, w = gauss_triangle(n)
    println("gauss_triangle($n): count=$(length(w)), sum_weights=$(sum(w))")
    ptsf, wf = _fallback_gauss_triangle(n)
    println("fallback($n): count=$(length(wf)), sum_weights=$(sum(wf))\n")
end
