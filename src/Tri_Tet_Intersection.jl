__precompile__(true)

module Tri_Tet_Intersections

using StaticArrays
using LinearAlgebra

# import

include("types.jl")
include("utility.jl")

export
    # types.jl
    Triangle,
    Tetrahedron,

    # utility.jl
    volume,
    asMatOnePad,
    triangleCross,
    area,
    centroid,
    verifyVolume

end
