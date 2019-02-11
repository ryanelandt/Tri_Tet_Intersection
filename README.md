# Tri_Tet_Intersections.jl

[![Build Status](https://travis-ci.com/ryanelandt/Tri_Tet_Intersections.jl.svg?branch=master)](https://travis-ci.com/ryanelandt/Tri_Tet_Intersections.jl)
[![codecov.io](https://codecov.io/github/ryanelandt/Tri_Tet_Intersections.jl/coverage.svg?branch=master)](https://codecov.io/github/ryanelandt/Tri_Tet_Intersections.jl?branch=master)

This packages computes intersections between triangles and tetrahedrons.
<!---
Basic geometry functions for triangles and tetrahedrons are also included.
-->
The intersection algorithm is the Sutherland–Hodgman algorithm implemented in tetrahedral coordinates.
This algorithm clips the triangle with the planes of the tetrahedron leaving only the parts of the triangle that lie within the tetrahedron.
Tetrahedral coordinates...

Points have three components when represented in Cartesian coordinates.
For clipping it is expedient to represent points in tetrahedral coordinates.
Tetrahedral coordinates give points as a weighted average of the position of a tetrahedron vertices.
Specifically, suppose that a tetrahedron has vertices $a$, $b$, $c$ and $d$.

```math
\begin{bmatrix}
a & b \\
c & d
\end{bmatrix}
```




<!---
end of README
-->
