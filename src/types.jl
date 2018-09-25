struct Triangle{T}
    v1::SVector{3,T}
    v2::SVector{3,T}
    v3::SVector{3,T}
    function Triangle(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T}
        return new{T}(v1, v2, v3)
    end
end

struct Tetrahedron{T}
    v1::SVector{3,T}
    v2::SVector{3,T}
    v3::SVector{3,T}
    v4::SVector{3,T}
    function Tetrahedron(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}, v4::SVector{3,T}) where {T}
        return new{T}(v1, v2, v3, v4)
    end
end

struct Triangle3D{T}
    frame::CartesianFrame3D
    v1::SVector{3,T}
    v2::SVector{3,T}
    v3::SVector{3,T}
    function Triangle3D(frame, v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T}
        return new{T}(frame, v1, v2, v3)
    end
end

struct Tetrahedron3D{T}
    frame::CartesianFrame3D
    v1::SVector{3,T}
    v2::SVector{3,T}
    v3::SVector{3,T}
    v4::SVector{3,T}
    function Tetrahedron3D(frame, v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}, v4::SVector{3,T}) where {T}
        return new{T}(frame, v1, v2, v3, v4)
    end
end
