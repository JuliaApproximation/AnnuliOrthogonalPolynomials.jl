module AnnuliOrthogonalPolynomials
using BlockArrays, BlockBandedMatrices, ClassicalOrthogonalPolynomials, ContinuumArrays, DomainSets, FastTransforms, FillArrays,
    HarmonicOrthogonalPolynomials, LazyArrays, LinearAlgebra,
    MultivariateOrthogonalPolynomials, QuasiArrays, SemiclassicalOrthogonalPolynomials, StaticArrays 

import Base: in, axes, getindex, broadcasted, tail, +, -, *, /, \, convert, OneTo, show, summary, ==, oneto, diff
import BlockArrays: block, blockindex, _BlockedUnitRange, blockcolsupport
import ClassicalOrthogonalPolynomials: checkpoints, ShuffledR2HC, TransformFactorization, ldiv, paddeddata, jacobimatrix, orthogonalityweight, SetindexInterlace
import ContinuumArrays: Weight, weight, grid, ℵ₁, ℵ₀, unweighted, plan_transform
import HarmonicOrthogonalPolynomials: BivariateOrthogonalPolynomial, MultivariateOrthogonalPolynomial, Plan
import MultivariateOrthogonalPolynomials: BlockOneTo, ModalInterlace, laplacian, ModalTrav
import SemiclassicalOrthogonalPolynomials: divdiff, HalfWeighted

export Block, SVector, CircleCoordinate, ZernikeAnnulus, ComplexZernikeAnnulus, Laplacian

abstract type AnnuliOrthogonalPolynomial{d,T} <: MultivariateOrthogonalPolynomial{d,T} end

include("annulus.jl")

end # module AnnuliOrthogonalPolynomials
