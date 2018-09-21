verifyVolume(t::Tetrahedron{T}) where {T} = LinearAlgebra.det(asMatOnePad(t)) * Float64(-1/6)

v1 = SVector{3,Float64}(0.0, 0.0, 0.0)
v2 = SVector{3,Float64}(1.0, 0.0, 0.0)
v3 = SVector{3,Float64}(0.0, 1.0, 0.0)
v4 = SVector{3,Float64}(0.0, 0.0, 1.0)

@testset "tetrahedron" begin
    tet_Tup = (v1, v2, v3, v4)
    tet_Tet = (Tetrahedron(tet_Tup...),)
    tet_SV  = (SVector{4,SVector{3,Float64}}(tet_Tup...),)
    for tet_rep = (tet_Tet, tet_SV, tet_Tup)
        @test volume(tet_rep...) ≈ 1/6
        @test centroid(tet_rep...) ≈ SVector{3,Float64}(1,1,1)/4
    end
end
