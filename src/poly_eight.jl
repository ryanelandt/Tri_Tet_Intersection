
struct poly_eight{N,T}
    n::Int64
    v::NTuple{8,SVector{N,T}}
    function poly_eight{N,T}() where {N,T}
        return new{N,T}(0)
    end
    function poly_eight(n::Int64, v::NTuple{8,SVector{N,T}}) where {N,T}
        return new{N,T}(n, v)
    end
    function poly_eight(v::NTuple{3,SVector{N,T}}) where {N,T}
        return new{N,T}(3, (v[1], v[2], v[3], v[1], v[1], v[1], v[1], v[1]))
    end
    function poly_eight(v::NTuple{4,SVector{N,T}}) where {N,T}
        return new{N,T}(4, (v[1], v[2], v[3], v[4], v[1], v[1], v[1], v[1]))
    end
    function poly_eight(v::SVector{N2,SVector{N,T}}) where {N2,N,T}
        return poly_eight(v.data)
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

function one_pad_then_mul(m::SMatrix{4,4,T1,16}, p::poly_eight{3,T2}) where {T1,T2}
    @inline one_pad_then_mul_(m::SMatrix{4,4,T1,16}, v::SVector{3,T2}) where {T1,T2} = m * onePad(v)
    t1 = one_pad_then_mul_(m, p[1])
    t2 = one_pad_then_mul_(m, p[2])
    t3 = one_pad_then_mul_(m, p[3])
    t4 = one_pad_then_mul_(m, p[4])
    length_p = length(p)
    if length_p <= 4
        return poly_eight(length_p, (t1, t2, t3, t4, t1, t1, t1, t1))
    else
        t5 = one_pad_then_mul_(m, p[5])
        t6 = one_pad_then_mul_(m, p[6])
        t7 = one_pad_then_mul_(m, p[7])
        t8 = one_pad_then_mul_(m, p[8])
        return poly_eight(length_p, (t1, t2, t3, t4, t5, t6, t7, t8))
    end
end

function mul_then_un_pad(m::SMatrix{4,4,T1,16}, p::poly_eight{4,T2}) where {T1,T2}
    @inline mul_then_un_pad_(m::SMatrix{4,4,T1,16}, v::SVector{4,T2}) where {T1,T2} = unPad(m * v)
    t1 = mul_then_un_pad_(m, p[1])
    t2 = mul_then_un_pad_(m, p[2])
    t3 = mul_then_un_pad_(m, p[3])
    t4 = mul_then_un_pad_(m, p[4])
    length_p = length(p)
    if length_p <= 4
        return poly_eight(length_p, (t1, t2, t3, t4, t1, t1, t1, t1))
    else
        t5 = mul_then_un_pad_(m, p[5])
        t6 = mul_then_un_pad_(m, p[6])
        t7 = mul_then_un_pad_(m, p[7])
        t8 = mul_then_un_pad_(m, p[8])
        return poly_eight(length_p, (t1, t2, t3, t4, t5, t6, t7, t8))
    end
end
