
# clipping functions assume that the triangle does not lie on a side of the tet

clip_in_tet_coordinates(z::NTuple{N,SVector{4,T}}) where {N,T} = clip_by_tet_plane(z, 1)

function clip_by_tet_plane(z1::SVector{4,T}, z2::NTuple{N,SVector{4,T}}, z3::SVector{4,T}, i::Int64) where {N,T}
    z = (z1, z2..., z3)
    return clip_by_tet_plane(z, i)
end

function tuple_forward(z::NTuple{N,SVector{4,T}}, i::Int64) where {N,T}
    return Tuple(z[mod1(k + i - 1, N)] for k = 1:N)
end

function clip_by_tet_plane(z::NTuple{N,SVector{4,T}}, i::Int64) where {N,T}
    (i == 5) && (return true, z)
    s = Tuple(z[k][i] for k = 1:N)
    is_non_pos = (s .<= 0.0)
    all(is_non_pos) && (return false, NTuple{0,SVector{4,T}}())
    if all(0.0 .<= s)  # non_negative
        return clip_by_tet_plane(z, i + 1)
    else
        for k = 1:N
            if is_non_pos[k] && !is_non_pos[mod1(k + 1, N)]
                z = tuple_forward(z, k)
                return cut_in_tet_coordinates(z, i)
            end
        end
    end
    error("this code should be unreachable")
end

function cut_in_tet_coordinates(z::NTuple{N,SVector{4,T}}, i::Int64) where {N,T}
    (z[N-1][i] <= 0.0) && (return clip_by_tet_plane(z[1:(N - 1)], i))
    z_start = clip_node(z[1], z[2], i)
    if 0.0 < z[N][i]
        z_end = clip_node(z[1], z[N], i)
        return clip_by_tet_plane(z_start, z[2:N], z_end, i + 1)
    else
        z_end = clip_node(z[N], z[N-1], i)
        return clip_by_tet_plane(z_start, z[2:(N-1)], z_end, i + 1)
    end
end

function clip_node(z_non::SVector{4,T}, z_pos::SVector{4,T}, i::Int64) where {T}
    v_non = z_non[i]
    v_pos = z_pos[i]
    return weightPoly(z_non, z_pos, v_non, v_pos)
end



# function clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, i::Int64) where {T}
#     (i == 5) && (return true, (z1, z2, z3))
#
#     s = SVector{3,T}(z1[i], z2[i], z3[i])
#     is_non_pos = (s .<= 0.0)
#     all(is_non_pos) && (return false, NTuple{0,SVector{4,T}}())
#     if all(0.0 .<= s)  # non_negative
#         return clip(z1, z2, z3, i+1)
#     else
#         is_non_pos[1] && !is_non_pos[2] && (return cut_clip(z1, z2, z3, i))
#         is_non_pos[2] && !is_non_pos[3] && (return cut_clip(z2, z3, z1, i))
#         is_non_pos[3] && !is_non_pos[1] && (return cut_clip(z3, z1, z2, i))
#     end
#     error("this code should be unreachable")
# end
#
# function clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}, i::Int64) where {T}
#     (i == 5) && (return true, (z1, z2, z3, z4))
#
#     s = SVector{4,T}(z1[i], z2[i], z3[i], z4[i])
#     is_non_pos = (s .<= 0.0)
#     all(is_non_pos) && (return false, NTuple{0,SVector{4,T}}())
#     if all(0.0 .<= s)  # non_negative
#         return clip(z1, z2, z3, z4, i + 1)
#     else
#         is_non_pos[1] && !is_non_pos[2] && (return cut_clip(z1, z2, z3, z4, i))
#         is_non_pos[2] && !is_non_pos[3] && (return cut_clip(z2, z3, z4, z1, i))
#         is_non_pos[3] && !is_non_pos[4] && (return cut_clip(z3, z4, z1, z2, i))
#         is_non_pos[4] && !is_non_pos[1] && (return cut_clip(z4, z1, z2, z3, i))
#     end
#     error("this code should be unreachable")
# end
#
# function clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}, z5::SVector{4,T}, i::Int64) where {T}
#     (i == 5) && (return true, (z1, z2, z3, z4, z5))
#
#     s = SVector{5,T}(z1[i], z2[i], z3[i], z4[i], z5[i])
#     is_non_pos = (s .<= 0.0)
#     all(is_non_pos) && (return false, NTuple{0,SVector{4,T}}())
#     if all(0.0 .<= s)  # non_negative
#         return clip(z1, z2, z3, z4, z5, i + 1)
#     else
#         is_non_pos[1] && !is_non_pos[2] && (return cut_clip(z1, z2, z3, z4, z5, i))
#         is_non_pos[2] && !is_non_pos[3] && (return cut_clip(z2, z3, z4, z5, z1, i))
#         is_non_pos[3] && !is_non_pos[4] && (return cut_clip(z3, z4, z5, z1, z2, i))
#         is_non_pos[4] && !is_non_pos[5] && (return cut_clip(z4, z5, z1, z2, z3, i))
#         is_non_pos[5] && !is_non_pos[1] && (return cut_clip(z5, z1, z2, z3, z4, i))
#     end
#     error("this code should be unreachable")
# end
#
# function clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}, z5::SVector{4,T}, z6::SVector{4,T}, i::Int64) where {T}
#     (i == 5) && (return true, (z1, z2, z3, z4, z5, z6))
#
#     s = SVector{6,T}(z1[i], z2[i], z3[i], z4[i], z5[i], z6[i])
#     is_non_pos = (s .<= 0.0)
#     all(is_non_pos) && (return false, NTuple{0,SVector{4,T}}())
#     if all(0.0 .<= s)  # non_negative
#         return clip(z1, z2, z3, z4, z5, z6, i + 1)
#     else
#         is_non_pos[1] && !is_non_pos[2] && (return cut_clip(z1, z2, z3, z4, z5, z6, i))
#         is_non_pos[2] && !is_non_pos[3] && (return cut_clip(z2, z3, z4, z5, z6, z1, i))
#         is_non_pos[3] && !is_non_pos[4] && (return cut_clip(z3, z4, z5, z6, z1, z2, i))
#         is_non_pos[4] && !is_non_pos[5] && (return cut_clip(z4, z5, z6, z1, z2, z3, i))
#         is_non_pos[5] && !is_non_pos[6] && (return cut_clip(z5, z6, z1, z2, z3, z4, i))
#         is_non_pos[6] && !is_non_pos[1] && (return cut_clip(z6, z1, z2, z3, z4, z5, i))
#     end
#     error("this code should be unreachable")
# end
#
# # z1 is gauranteed to be non-positive
# # z2 is gauranteed to be positive
# function cut_clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, i::Int64) where {T}
#     z_start = clip_node(z1, z2, i)
#     if 0.0 < z3[i]
#         z_end = clip_node(z1, z3, i)
#         return clip(z_start, z2, z3, z_end, i + 1)
#     else
#         z_end = clip_node(z3, z2, i)
#         return clip(z_start, z2, z_end, i + 1)
#     end
# end
#
# function cut_clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}, i::Int64) where {T}
#     (z3[i] <= 0.0) && (return cut_clip(z1, z2, z3, i))
#     z_start = clip_node(z1, z2, i)
#     if 0.0 < z4[i]
#         z_end = clip_node(z1, z4, i)
#         return clip(z_start, z2, z3, z4, z_end, i + 1)
#     else
#         z_end = clip_node(z4, z3, i)
#         return clip(z_start, z2, z3, z_end, i + 1)
#     end
# end
#
# function cut_clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}, z5::SVector{4,T}, i::Int64) where {T}
#     (z4[i] <= 0.0) && (return cut_clip(z1, z2, z3, z4, i))
#     z_start = clip_node(z1, z2, i)
#     if 0.0 < z5[i]
#         z_end = clip_node(z1, z5, i)
#         return clip(z_start, z2, z3, z4, z5, z_end, i + 1)
#     else
#         z_end = clip_node(z5, z4, i)
#         return clip(z_start, z2, z3, z4, z_end, i + 1)
#     end
# end
#
# function cut_clip(z1::SVector{4,T}, z2::SVector{4,T}, z3::SVector{4,T}, z4::SVector{4,T}, z5::SVector{4,T}, z6::SVector{4,T}, i::Int64) where {T}
#     (z5[i] <= 0.0) && (return cut_clip(z1, z2, z3, z4, z5, i))
#     z_start = clip_node(z1, z2, i)
#     if 0.0 <= z6[i]
#         z_end = clip_node(z1, z6, i)
#         return true, (z_start, z2, z3, z4, z5, z6, z_end)
#     else
#         z_end = clip_node(z6, z5, i)
#         return true, (z_start, z2, z3, z4, z5, z_end)
#     end
# end
#
# function clip_node(z_non::SVector{4,T}, z_pos::SVector{4,T}, i::Int64) where {T}
#     v_non = z_non[i]
#     v_pos = z_pos[i]
#     return weightPoly(z_non, z_pos, v_non, v_pos)
# end
