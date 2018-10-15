
function minimum_area(p::poly_eight{3,Float64}, n̂::SVector{3,Float64}, r::SVector{3,Float64})
    isempty(p) && (return -Inf)
    min_area = Inf
    length_clip = length(p)
    for k = 1:length_clip
        area_k = signed_area(n̂, p[k], p[mod1(k + 1, length_clip)], r)
        min_area = min(area_k, min_area)
    end
    return min_area
end

@testset "static clip" begin
    tol = 1.0e-13
    n_empty = 0
    n_hits = zeros(MVector{8,Int64})
    for k_tet = 1:500000
        if n_hits[8] <= 2
            r_orig = make_4_sided()

            n̂ = triangleNormal(r_orig[1], r_orig[2], r_orig[3])
            n̂_check = triangleNormal(r_orig[2], r_orig[3], r_orig[4])
            (n̂ ≈ n̂_check) || error("quad points do not lie on a plane")
            plane = SMatrix{1,4,Float64,4}(n̂[1], n̂[2], n̂[3], -dot(n̂, r_orig[1]))

            tet = roll_non_degenerate_tet()
            A = asMatOnePad(tet[1], tet[2], tet[3], tet[4])
            inv_A = inv(A)
            ζ_orig = one_pad_then_mul(inv_A, r_orig)
            ζ_clip = clip_in_tet_coordinates(ζ_orig)
            r_clip = mul_then_un_pad(A, ζ_clip)

            for k = 1:length(r_clip)  # confirm that all vertices on plane
                # this loose tolerance is necessary amazingly
                # TODO: figure out why this is as it seems excessive
                @test abs(dist_from_plane(plane, r_clip[k])) < (100 * tol)
            end

            for k = 1:30  # test 30 random points
                r_in_plane = project_into_plane(plane, randn(SVector{3,Float64}))
                min_ζ = minimum(inv_A * onePad(r_in_plane))
                min_area_orig = minimum_area(r_orig, n̂, r_in_plane)  # is_inside origional
                min_area_clip = minimum_area(r_clip, n̂, r_in_plane)  # is_inside clip
                if tol < min_area_clip  # inside clipping
                    point_inside_ζ = -tol < min_ζ
                    point_inside_orig = -tol < min_area_orig
                    @test point_inside_ζ && point_inside_orig
                else
                    point_outside_ζ = min_ζ < tol
                    point_outside_orig = min_area_orig < tol
                    @test point_outside_ζ || point_outside_orig
                end
            end

            if isempty(r_clip)
                n_empty += 1
            else
                n_hits[length(r_clip)] += 1
            end
        end
    end
    @test 2 <= n_hits[8]  # sanity check
    @test 1000 < n_empty  # sanity check
end



#
# # TODO: write good test for both triangle and convex quadlateral cases
# # 1. confirm that all clipped vertices lie in a plane
# # 2. for each of 1000 random points...
# #   if lies in clipped region
# #        lies inside both origional geometry and clipping tet
# #   else
# #        lies outside either origional geometry or clipping tet
# #
#
# function calcSignedArea(A_top, r̃, n̂, i_1::Int64, i_2, the_out::NTuple{N,T}) where {N,T}
#     v_1 = A_top * the_out[i_1]
#     v_2 = A_top * the_out[i_2]
#     vector_area = triangleCross(v_1, v_2, r̃)
#     return dot(vector_area, n̂)
# end
#
# function findMinSignedArea(A_top, r̃, n̂, the_out::NTuple{N,T}) where {N,T}
#     min_signed_area = Inf
#     for k = 1:N
#         signed_area_k = calcSignedArea(A_top, r̃, n̂, k, mod1(k+1, N), the_out)
#         min_signed_area = min(signed_area_k, min_signed_area)
#     end
#     return min_signed_area
# end
#
# function test_phi(ϕ, inv_A, A3, A_top, n̂, tup_p123, the_out::NTuple{N,T}) where {N,T}
#     (N == 0) && (return true)  # TODO: this case is not handled correctly
#
#     ϕ /= sum(ϕ)  # phi is in triangular coordinates and needs to sum to 1
#     r̃ = A3 * ϕ
#     p1, p2, p3 = tup_p123
#     (1.0e-14 < abs(volume(p1, p2, p3, r̃))) && error("cart_pick not in triangle plane")
#     ζ = inv_A * onePad(r̃)
#     (sum(ζ) ≈ 1) || error("tet point doesn't sum to 1")
#     s_min = findMinSignedArea(A_top, r̃, n̂, the_out)
#     all(0.0 .< ζ) && (return -1.0e-14 < s_min)
#     return s_min < 1.0e-14
# end
#
# function test_coplanar(A_top, tup_p123, the_out::NTuple{N, T}) where {N,T}
#     p1, p2, p3 = tup_p123
#     for k = 1:N
#         p_clip = A_top * the_out[k]
#         (volume(p1, p2, p3, p_clip) < -1.0e-14) && (return false)
#     end
#     return true
# end
#
#
# function test_clip_inside_tri(A_top, n̂, tup_p123, the_out::NTuple{N,T}) where {N,T}
#     p1, p2, p3 = tup_p123
#     for k = 1:N
#         p_clip = A_top * the_out[k]
#         s_1 = signed_area(n̂, p1, p2, p_clip)
#         s_2 = signed_area(n̂, p2, p3, p_clip)
#         s_3 = signed_area(n̂, p3, p1, p_clip)
#         (min(s_1, s_2, s_3) < -1.0e-14) && (return false)
#     end
#     return true
# end
#
# function test_clip_sign(the_out::NTuple{N,SVector{4,T}}) where {N,T}
#     for k = 1:N
#         ζ = the_out[k]
#         any(ζ .< -1.0e-14) && (return false)
#         (sum(ζ) ≈ 1.0) || (return false)
#     end
#     return true
# end
#
# function signed_area(n̂::SVector{3,T}, p1::SVector{3,T}, p2::SVector{3,T}, p3::SVector{3,T}) where {T}
#     vec_area = triangleCross(p1, p2, p3)
#     return dot(vec_area, n̂)
# end
#
# using LinearAlgebra
# using Test
#
#
# @testset "Tri Tet Clipping" begin
#
#     v1 = SVector{3,Float64}(1.0,0.0,0.0)
#     v2 = SVector{3,Float64}(0.0,1.0,0.0)
#     v3 = SVector{3,Float64}(0.0,0.0,1.0)
#     v4 = SVector{3,Float64}(1.0,1.0,1.0)
#     A = asMatOnePad(v1,v2,v3,v4)
#     inv_A = inv(A)
#
#     n_seven_need = 2
#     n_hit_by_side = zeros(Int64, 7)
#     n_fail = zeros(Int64, 1)
#
#     for k = 1:100000
#         if n_hit_by_side[7] < n_seven_need
#             p1 = randn(SVector{3,Float64})
#             p2 = randn(SVector{3,Float64})
#             p3 = randn(SVector{3,Float64})
#             tup_p123 = (p1, p2, p3)
#
#             n̂ = triangleNormal(p1, p2, p3)
#             z1 = inv_A * onePad(p1)
#             z2 = inv_A * onePad(p2)
#             z3 = inv_A * onePad(p3)
#
#             A_top = getTop(A)
#             A3 = asMat(p1, p2, p3)
#             # is_exist, the_out = clip(z1, z2, z3, 1)
#             is_exist, the_out = clip_in_tet_coordinates(z1, z2, z3)
#
#             n_clip_vert = length(the_out)
#             if n_clip_vert == 0
#                 n_fail[1] += 1
#                 @test !is_exist
#             else
#                 n_hit_by_side[n_clip_vert] += 1
#                 @test is_exist
#             end
#
#             @test test_clip_sign(the_out)
#             @test test_coplanar(A_top, tup_p123, the_out)
#             @test test_clip_inside_tri(A_top, n̂, tup_p123, the_out)
#
#             for k = 1:10
#                 ϕ = rand(SVector{3,Float64})
#                 @test test_phi(ϕ, inv_A, A3, A_top, n̂, tup_p123, the_out)
#             end
#             for k = 1:3
#                 ϕ = MVector{3,Float64}(0,0,0)
#                 ϕ[k] += 1
#                 @test test_phi(SVector{3,Float64}(ϕ), inv_A, A3, A_top, n̂, tup_p123, the_out)
#             end
#         end
#     end
#     @test (n_seven_need <= n_hit_by_side[7])
# end
