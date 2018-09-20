v1 = SVector{3,Float64}(0.0, 0.0, 0.0)
v2 = SVector{3,Float64}(1.0, 0.0, 0.0)
v3 = SVector{3,Float64}(0.0, 1.0, 0.0)

@testset "triangle" begin
    tri_Tup = (v1, v2, v3)
    tri_Tet = (Triangle(tri_Tup...),)
    tri_SV  = (SVector{3,SVector{3,Float64}}(tri_Tup...),)
    for tri_rep = (tri_Tet, tri_SV, tri_Tup)
        @test area(tri_rep...) ≈ 0.5
        @test centroid(tri_rep...) ≈ SVector{3,Float64}(1/3,1/3,0)
        @test triangleCross(tri_rep...) ≈ SVector{3,Float64}(0,0,1)
    end
end
