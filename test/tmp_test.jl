using Test, AlgebraicCurveOrthogonalPolynomials, ClassicalOrthogonalPolynomials, LinearAlgebra
import ForwardDiff: derivative, hessian, gradient
import AlgebraicCurveOrthogonalPolynomials: ZernikeAnnulusTransform, ZernikeAnnulusITransform, plotgrid, plotvalues
import BlockArrays: PseudoBlockArray, blockcolsupport

a = 1
@testset "transform" begin
    for ρ in (0.1, 0.5, 2/3)
        for (a,b) in [(0,0), (1,0), (0,1), (1,1)]
            Z = ZernikeAnnulus(ρ, a, b)
            xy = axes(Z,1)
            for j = 1:10
                f = xy -> Z[xy, j]
                @test Z[:, 1:10] \ f.(xy) ≈ (1:10 .== j)
            end
        end
    end

    N = 5
    T = ZernikeAnnulusTransform{Float64}(N,0,0,0,0.5)
    Ti = ZernikeAnnulusITransform{Float64}(N,0,0,0,0.5)

    v = PseudoBlockArray(randn(sum(1:N)),1:N)
    @test T * (Ti * v) ≈ v
    @test_throws MethodError T * randn(15)

    Z = ZernikeAnnulus(0.5, 1, 1); c = [1;2;zeros(∞)]
    u = Z * c
    Bs = last(blockcolsupport(u.args[2]))
    g = plotgrid(Z, Bs)
    @test g == plotgrid(Weighted(Z), Bs)
    @test norm(plotvalues(u, g)) ≈ 13.22875655532295
    @test norm(plotvalues(Weighted(Z)*c, g)) ≈ 0.9472636114669317
end