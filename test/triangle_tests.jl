v1 = SVector{3,Float64}(0.0, 0.0, 0.0)
v2 = SVector{3,Float64}(1.0, 0.0, 0.0)
v3 = SVector{3,Float64}(0.0, 1.0, 0.0)

@testset "triangle" begin
    tri_Sep = (v1, v2, v3)
    tri_SV  = (SVector{3,SVector{3,Float64}}(tri_Sep...),)
    tri_Tup  = ((tri_Sep...),)
    for tri_rep = (tri_Sep, tri_SV, tri_Tup)
        @test area(tri_rep...) ≈ 0.5
        @test centroid(tri_rep...) ≈ SVector{3,Float64}(1/3,1/3,0)
        @test triangleCross(tri_rep...) ≈ SVector{3,Float64}(0,0,1)
    end
end
