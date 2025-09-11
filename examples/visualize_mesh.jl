"""
Visualize antenna meshes (dipole + patch) from
`antenna_impedance_calculation.jl` using GLMakie.

Usage:
  julia --project=. examples/visualize_mesh.jl

Install dependencies first (see README below):
  julia --project=. -e 'using Pkg; Pkg.add(["GLMakie"])'

If you run this on a headless server, the script will still save PNGs and an OBJ file.
"""

using GLMakie

# Reuse functions defined in the example file
include(joinpath(@__DIR__, "antenna_impedance_calculation.jl"))

# Helper: convert your SVector vertices and triangle indices into Makie-friendly structures
function makie_points(vertices)
    # returns a Vector{Point3f0}
    return Point3f0.([(v[1], v[2], v[3]) for v in vertices])
end

function makie_faces(triangle_indices)
    # returns Vector{NTuple{3,Int}}
    return Tuple.(triangle_indices)
end

function save_obj(path::AbstractString, vertices, triangles)
    open(path, "w") do io
        for v in vertices
            println(io, "v ", v[1], " ", v[2], " ", v[3])
        end
        for t in triangles
            # OBJ is 1-based already
            println(io, "f ", t[1], " ", t[2], " ", t[3])
        end
    end
end

function plot_and_save(vertices, triangles; title="mesh", save_png=true, save_obj_file=true)
    pts = makie_points(vertices)
    faces = makie_faces(triangles)

    scene = mesh(pts, faces; color = :lightsteelblue, shading = true)
    scene.center = Observable(Point3f0(0,0,0))
    axis3d = scene

    if save_png
        pngfile = "$(title).png"
        save(pngfile, scene)
        println("Saved PNG: $pngfile")
    end

    if save_obj_file
        objfile = "$(title).obj"
        save_obj(objfile, vertices, triangles)
        println("Saved OBJ: $objfile")
    end

    return scene
end

function main()
    # Dipole example (reuses create_dipole_mesh from the included file)
    length = 0.5
    radius = 0.001
    n_segments = 20
    println("Creating dipole mesh...")
    vertices_dipole, triangles_dipole = create_dipole_mesh(length, radius, n_segments)
    println("Dipole vertices: ", length(vertices_dipole), ", triangles: ", length(triangles_dipole))

    println("Plotting dipole mesh and saving files...")
    scene1 = plot_and_save(vertices_dipole, triangles_dipole; title="dipole_mesh")

    # Patch example
    length_p = 0.031
    width_p = 0.024
    n_x = 15
    n_y = 12
    println("Creating patch mesh...")
    vertices_patch, triangles_patch = create_patch_antenna_mesh(length_p, width_p, n_x, n_y)
    println("Patch vertices: ", length(vertices_patch), ", triangles: ", length(triangles_patch))

    println("Plotting patch mesh and saving files...")
    scene2 = plot_and_save(vertices_patch, triangles_patch; title="patch_mesh")

    println("Done. If a display is available, two windows should appear; PNGs/OBJs saved to the current directory.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
