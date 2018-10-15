__precompile__(true)

module Tri_Tet_Intersections

using StaticArrays
using LinearAlgebra
using RigidBodyDynamics.Spatial

include("poly_eight.jl")
include("types.jl")
include("utility.jl")
include("clipped_polygon.jl")
include("quadrature.jl")
include("static_clip.jl")
include("plane_tet_intersection.jl")
include("test_utility.jl")

export
    # poly_eight.jl
    poly_eight,
    mul_then_un_pad,
    one_pad_then_mul,

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

    # # clipped_polygon.jl
    # ClippedPolygon,
    # weightPoly,
    # add!,
    # setTriangleTetMutable,
    # modularTriTetClip,
    # clipTriangleWithOneTetPlane,
    # tet_clip_poly_to_cartesian!,

    # quadrature.jl
    TriTetQuadRule,
    getTriQuadRule,
    getTetQuadRule,

    # static_clip.jl
    poly_eight,
    clip_in_tet_coordinates,
    # clip,

    # plane_tet_intersection.jl
    clip_plane_tet,

    # test_utility.jl
    roll_non_degenerate_tet,
    make_4_sided,
    normal,
    dist_from_plane,
    project_into_plane,
    signed_area

end
