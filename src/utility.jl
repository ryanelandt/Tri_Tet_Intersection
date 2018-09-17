### Triangle ###
triangleCross(t::Triangle{T}) where {T} = cross(t.v2 - t.v1, t.v3 - t.v2)
function triangleCross(sv::SVector{3,SVector{3,Float64}})
    t = Triangle(sv[1], sv[2], sv[3])
    return triangleCross(t)
end

area(t::Triangle{T}) where {T} = LinearAlgebra.norm(triangleCross(t)) * 0.5
function area(sv::SVector{3,SVector{3,Float64}})
    t = Triangle(sv[1], sv[2], sv[3])
    return area(t)
end

centroid(t::Triangle{T}) where {T} = (t.v1 + t.v2 + t.v3) * Float64(1/3)  # 4 times faster than dividing by 3
function centroid(sv::SVector{3,SVector{3,Float64}})
    t = Triangle(sv[1], sv[2], sv[3])
    return centroid(t)
end

triangleNormal(t::Triangle) = normalize(triangleCross(t))
function triangleNormal(sv::SVector{3,SVector{3,Float64}})
    t = Triangle(sv[1], sv[2], sv[3])
    return triangleNormal(t)
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
function volume(sv::SVector{4,SVector{3,Float64}})
    t = Tetrahedron(sv[1], sv[2], sv[3], sv[4])
    return volume(t)
end

centroid(t::Tetrahedron{T}) where {T} = (t.v1 + t.v2 + t.v3 + t.v4) * 0.25
function centroid(sv::SVector{4,SVector{3,Float64}})
    t = Tetrahedron(sv[1], sv[2], sv[3], sv[4])
    return centroid(t)
end

function asMatOnePad(t::Tetrahedron{T}) where {T}
    A = @SMatrix [
    t.v1[1] t.v2[1] t.v3[1] t.v4[1];
    t.v1[2] t.v2[2] t.v3[2] t.v4[2];
    t.v1[3] t.v2[3] t.v3[3] t.v4[3];
    one(T)  one(T)  one(T)  one(T)
    ]
    return A
end
