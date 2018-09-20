@testset "exports" begin
    # Ensure that every exported name is actually defined
    for name in names(Tri_Tet_Intersections)
        @test isdefined(Tri_Tet_Intersections, name)
    end
end
