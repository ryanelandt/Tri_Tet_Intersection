__precompile__(true)

module Tri_Tet_Intersections

using StaticArrays
using LinearAlgebra


include("types.jl")
include("utility.jl")
include("clipped_polygon.jl")

export
    # types.jl
    Triangle,
    Tetrahedron,

    # utility.jl
    volume,
    asMatOnePad,
    asMat,
    getTop,
    triangleCross,
    area,
    centroid,
    triangleNormal,

    # clipped_polygon.jl
    ClippedPolygon,
    weightPoly,
    add!,
    setTriangleTetMutable,
    modularTriTetClip,
    clipTriangleWithOneTetPlane,
    tet_clip_poly_to_cartesian!

end
