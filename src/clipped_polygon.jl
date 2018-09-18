struct ClippedPolygon{N,T}
    v::MVector{7,SVector{N,T}}
    i::Int64
    function ClippedPolygon{N,T}() where {N,T}
        v = MVector{7,SVector{N,T}}(undef)
        i = 0
        return new{N,T}(v, i)
    end
end

Base.length(c::ClippedPolygon) = c.i
function Base.empty!(c::ClippedPolygon)
    c.i = 0
    return nothing
end
Base.isempty(c::ClippedPolygon) = (c.i == 0)
Base.getindex(c::ClippedPolygon, k::Int64) = c.v[k]
function weightPoly(p1::SVector{N,T}, p2::SVector{N,T}, w1::T, w2::T) where {N,T}
  sum_weight = w1 - w2
  c1 = w1 / sum_weight
  c2 = w2 / sum_weight
  return c1 * p2 - c2 * p1
end
@inline function add!(c::ClippedPolygon{N,T}, a::SVector{N,T}) where {N,T}
  c.i += 1
  c.v[c.i] = a
  return nothing
end
function centroid(p_new::ClippedPolygon{3,T}) where {T}
  cart_a = p_new.v[1]
  cart_c = p_new.v[2]  # because c becomes b
  cum_sum  = zero(T)
  cum_prod = zeros(SVector{3,T})
  for k = 3:p_new.i
    cart_b = cart_c
    cart_c = p_new.v[k]
    area_ = area(cart_a, cart_b, cart_c)
    if 0.0 < area_  # when area is 0.0 there in a nan in the partials
      cum_prod += area_ * centroid(cart_a, cart_b, cart_c)
      cum_sum  += area_
    end
  end
  if cum_sum == 0.0  # just return an arbitrary vertex if the area is zero
    return cum_sum, cart_a
  else
    return cum_sum, cum_prod / cum_sum
  end
end

function setTriangleTetMutable(p_tet_1::ClippedPolygon{N,T}, p1::SVector{4,T}, p2::SVector{4,T}, p3::SVector{4,T}) where {N,T}
  @inbounds p_tet_1.v[1] = p1
  @inbounds p_tet_1.v[2] = p2
  @inbounds p_tet_1.v[3] = p3
  p_tet_1.i = 3
  return nothing
end
function modularTriTetClip(p_tet_1::ClippedPolygon{N,T}, p_tet_2::ClippedPolygon{N,T}) where {N,T}
  clipTriangleWithOneTetPlane(p_tet_1, p_tet_2, 1)  # info ends up in 2
  if p_tet_2.i <= 2 # exit early if the area is gauranteed to be zero
    p_tet_1.i = 0
    return nothing
  end

  clipTriangleWithOneTetPlane(p_tet_2, p_tet_1, 2)  # info ends up in 1
  (p_tet_1.i <= 2) && (return nothing)

  clipTriangleWithOneTetPlane(p_tet_1, p_tet_2, 3)  # info ends up in 2
  if p_tet_2.i <= 2
    p_tet_1.i = 0
    return nothing
  end

  clipTriangleWithOneTetPlane(p_tet_2, p_tet_1, 4)  # info ends up in 1
  return nothing
end
function clipTriangleWithOneTetPlane(p_tet_1::ClippedPolygon{N,T}, p_tet_2::ClippedPolygon{N,T}, k_cut_plane::Int64) where {N,T}
  # this takes prefilled p_tet_1 and clips the selected clipped plane into p_tet_2
  p_tet_2.i = 0
  zeta_prev = p_tet_1.v[p_tet_1.i]
  val_prev = zeta_prev[k_cut_plane]
  prev_val_zero = (0.0 == val_prev)
  prev_val_pos = (0.0 < val_prev)
  for k = 1:p_tet_1.i
    @inbounds zeta_next = p_tet_1.v[k]
    @inbounds val_next = zeta_next[k_cut_plane]
    next_val_zero = (0.0 == val_next)
    next_val_pos = (0.0 < val_next)
    if prev_val_pos || prev_val_zero
      add!(p_tet_2, zeta_prev)
    end
    if xor(prev_val_pos, next_val_pos) && (!prev_val_zero) && (!next_val_zero)
      zeta_ = weightPoly(zeta_prev, zeta_next, val_prev, val_next)
      add!(p_tet_2, zeta_)
    end
    zeta_prev = zeta_next
    val_prev = val_next
    prev_val_zero = next_val_zero
    prev_val_pos = next_val_pos
  end
  return nothing
end
function tet_clip_poly_to_cartesian!(p_new::ClippedPolygon{3,T}, p_tet_1::ClippedPolygon{4,T}, A_w_zeta_top::SMatrix{3,4,T,12}) where {T}
  for k = 1:p_tet_1.i
    p_new.v[k] = A_w_zeta_top * p_tet_1.v[k]
  end
  p_new.i = p_tet_1.i
  return nothing
end
