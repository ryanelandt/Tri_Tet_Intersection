__precompile__(true)

module Tri_Tet_Intersection

using StaticArrays
# using 

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
