__precompile__(true)

module Tri_Tet_Intersections

using StaticArrays
using LinearAlgebra
using RigidBodyDynamics.Spatial


include("types.jl")
include("utility.jl")
include("clipped_polygon.jl")
include("quadrature.jl")

export
    # types.jl
    Triangle,
    Tetrahedron,
    Triangle3D,
    Tetrahedron3D,

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
    tet_clip_poly_to_cartesian!,

    # quadrature.jl
    TriTetQuadRule,
    getTriQuadRule,
    getTetQuadRule


end
