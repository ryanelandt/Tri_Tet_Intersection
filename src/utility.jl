
function weightPoly(p1::SVector{N,T}, p2::SVector{N,T}, w1::T, w2::T) where {N,T}
  sum_weight = w1 - w2
  c1 = w1 / sum_weight
  c2 = w2 / sum_weight
  return c1 * p2 - c2 * p1
end

unPad(a::SVector{4,T}) where {T} = SVector{3,T}(a[1], a[2], a[3])
unPad(a::SMatrix{1,4,T,4}) where {T} = SVector{3,T}(a[1], a[2], a[3])
onePad(a::SVector{3,T}) where {T} = SVector{4,T}(a[1], a[2], a[3], one(T))
zeroPad(a::SVector{3,T}) where {T} = SVector{4,T}(a[1], a[2], a[3], zero(T))

### Triangle ###
triangleCross(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T} = cross(v2 - v1, v3 - v2)
area(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T} = LinearAlgebra.norm(triangleCross(v1, v2, v3)) * 0.5
centroid(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T} = (v1 + v2 + v3) * Float64(1/3)  # 4 times faster than dividing by 3
triangleNormal(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T} = normalize(triangleCross(v1, v2, v3))

function asMatOnePad(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T}
    A = @SMatrix [
    v1[1] v2[1] v3[1];
    v1[2] v2[2] v3[2];
    v1[3] v2[3] v3[3];
    one(T)  one(T)  one(T)
    ]
    return A
end

function asMat(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T}
    A = @SMatrix [
    v1[1] v2[1] v3[1];
    v1[2] v2[2] v3[2];
    v1[3] v2[3] v3[3]
    ]
    return A
end

for funName in (:triangleCross, :area, :centroid, :triangleNormal, :asMatOnePad, :asMat)
    @eval begin
        function $funName(sv::SVector{3,SVector{3,T}}) where {T}
            return $funName(sv[1], sv[2], sv[3])
        end
        function $funName(sv::NTuple{3,SVector{3,T}}) where {T}
            return $funName(sv[1], sv[2], sv[3])
        end
    end
end

### Tetrahedron ###
function volume(t::Tetrahedron{T}) where {T}
    # NOTE: This is an algebraic refactorization of a symbolic answer
    a1, a2, a3 = t.v1[1], t.v1[2], t.v1[3]
    b1, b2, b3 = t.v2[1], t.v2[2], t.v2[3]
    c1, c2, c3 = t.v3[1], t.v3[2], t.v3[3]
    d1, d2, d3 = t.v4[1], t.v4[2], t.v4[3]

    V =       (b1 - a1) * (c2 * d3 - c3 * d2)
    V = muladd(b2 - a2,    c3 * d1 - c1 * d3, V)
    V = muladd(b3 - a3,    c1 * d2 - c2 * d1, V)

    V = muladd(c1 - d1,    a3 * b2 - a2 * b3, V)
    V = muladd(c2 - d2,    a1 * b3 - a3 * b1, V)
    V = muladd(c3 - d3,    a2 * b1 - a1 * b2, V)
    return V * Float64(1/6)
end

centroid(t::Tetrahedron{T}) where {T} = (t.v1 + t.v2 + t.v3 + t.v4) * 0.25

function asMatOnePad(t::Tetrahedron{T}) where {T}
    A = @SMatrix [
    t.v1[1] t.v2[1] t.v3[1] t.v4[1];
    t.v1[2] t.v2[2] t.v3[2] t.v4[2];
    t.v1[3] t.v2[3] t.v3[3] t.v4[3];
    one(T)  one(T)  one(T)  one(T)
    ]
    return A
end

function asMat(t::Tetrahedron{T}) where {T}
    A = @SMatrix [
    t.v1[1] t.v2[1] t.v3[1] t.v4[1];
    t.v1[2] t.v2[2] t.v3[2] t.v4[2];
    t.v1[3] t.v2[3] t.v3[3] t.v4[3]
    ]
    return A
end

for funName in (:volume, :centroid, :asMatOnePad, :asMat)
    @eval begin
        function $funName(sv::SVector{4,SVector{3,T}}) where {T}
            t = Tetrahedron(sv[1], sv[2], sv[3], sv[4])
            return $funName(t)
        end
        function $funName(sv_1::SVector{3,T}, sv_2::SVector{3,T}, sv_3::SVector{3,T}, sv_4::SVector{3,T}) where {T}
            t = Tetrahedron(sv_1, sv_2, sv_3, sv_4)
            return $funName(t)
        end
    end
end
