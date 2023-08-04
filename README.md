## AnnuliOrthogonalPolynomials.jl

[![Build Status](https://github.com/ioannisPApapadopoulos/AnnuliOrthogonalPolynomials.jl/workflows/CI/badge.svg)](https://github.com/JuliaApproximation/MultivariateOrthogonalPolynomials.jl/actions)
[![codecov](https://github.com/ioannisPApapadopoulos/AnnuliOrthogonalPolynomials.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaApproximation/MultivariateOrthogonalPolynomials.jl)


A Julia package for orthogonal polynomials on an annulus.

This code used to be part of [AlgebraicCurveOrthogonalPolynomials.jl](https://github.com/JuliaApproximation/AlgebraicCurveOrthogonalPolynomials.jl).


This experimental package implements `ZernikeAnnulus` and `ComplexZernikeAnnulus`,
families of orthogonal polynomials on the annulus `{ρ² ≤ x² + y² ≤ 1}`
with respect to the weight `(r² - ρ²)ᵃ * (1-r²)ᵇ`. This builds on top of
[SemiclassicalOrthogonalPolynomials.jl](https://github.com/JuliaApproximation/SemiclassicalOrthogonalPolynomials.jl) and [MultivariateOrthogonalPolynomials.jl](https://github.com/JuliaApproximation/MultivariateOrthogonalPolynomials.jl).

```julia
julia> using AnnuliOrthogonalPolynomials, StaticArrays

julia> Z = ZernikeAnnulus(0.5, 0, 0)
ZernikeAnnulus{Float64}(0.5, 0.0, 0.0)

julia> Z[SVector(0.5,0.2), Block.(1:3)] # Blocked by degree
3-blocked 6-element BlockVector{Float64}:
 1.0
 ───────────────────
 0.20000000000000004
 0.5000000000000001
 ───────────────────
 1.547298721428197
 0.20000000000000007
 0.21000000000000008

julia> x,y = first.(axes(Z,1)), last.(axes(Z,1));

julia> u = Z * (Z \ (exp.(x) .* cos.(y))) # Expand in annulus OPs
ZernikeAnnulus{Float64}(0.5, 0.0, 0.0) * [1.0, 9.352800883640776e-18, 0.9999999999999999, -5.2512913921638607e-17, -2.337754850841795e-18, 0.5, 2.153350906504162e-18, 1.1925107582660448e-17, -5.513529126039963e-18, 0.1666666666666666  …  ]

julia> u[SVector(0.5,0.2)] # Evaluate expansion
1.6158566135891377
```