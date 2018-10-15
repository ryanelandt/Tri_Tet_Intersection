
struct poly_eight{N,T}
    n::Int64
    v::NTuple{8,SVector{N,T}}
    function poly_eight{N,T}() where {N,T}
        return new{N,T}(0)
    end
    function poly_eight(n::Int64, v::NTuple{8,SVector{N,T}}) where {N,T}
        return new{N,T}(n, v)
    end
end

@inline Base.isempty(p::poly_eight) = (p.n == 0)
@inline Base.length(p::poly_eight) = p.n
@inline Base.getindex(p::poly_eight, k::Int64) = p.v[k]

function centroid(p_new::poly_eight{3,T}) where {T}
  cart_a = p_new.v[1]
  cart_c = p_new.v[2]  # because c becomes b
  cum_sum  = zero(T)
  cum_prod = zeros(SVector{3,T})
  for k = 3:p_new.n
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
