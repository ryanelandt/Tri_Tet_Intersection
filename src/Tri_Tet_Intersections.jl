__precompile__(true)

module Tri_Tet_Intersections

using StaticArrays
using LinearAlgebra
using RigidBodyDynamics.Spatial


include("types.jl")
include("utility.jl")
include("clipped_polygon.jl")
include("quadrature.jl")
include("static_clip.jl")
include("plane_tet_intersection.jl")

export
    # types.jl
    Triangle,
    Tetrahedron,
    Triangle3D,
    Tetrahedron3D,

    # utility.jl
    unPad,
    onePad,
    zeroPad,
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
    tet_clip_poly_to_cartesian!,

    # quadrature.jl
    TriTetQuadRule,
    getTriQuadRule,
    getTetQuadRule,

    # static_clip.jl
    clip_in_tet_coordinates,
    # clip,

    # plane_tet_intersection.jl
    clip_plane_tet

end
