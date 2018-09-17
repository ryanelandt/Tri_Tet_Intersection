verifyVolume(t::Tetrahedron{T}) where {T} = LinearAlgebra.det(asMatOnePad(t)) * Float64(-1/6)  # TODO: move to test/test_util.jl

v1 = SVector{3,Float64}(0.0, 0.0, 0.0)
v2 = SVector{3,Float64}(1.0, 0.0, 0.0)
v3 = SVector{3,Float64}(0.0, 1.0, 0.0)

@testset "triangle" begin
    t = Triangle(v1, v2, v3)
    @test area(t) ≈ 0.5
    @test centroid(t) ≈ SVector{3,Float64}(1/3,1/3,0)
    @test triangleCross(t) ≈ SVector{3,Float64}(0,0,1)

    sv_t = SVector{3,SVector{3,Float64}}(v1, v2, v3)
    @test area(sv_t) ≈ 0.5
    @test centroid(sv_t) ≈ SVector{3,Float64}(1/3,1/3,0)
    @test triangleCross(sv_t) ≈ SVector{3,Float64}(0,0,1)
end
