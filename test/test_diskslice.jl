using AnnuliOrthogonalPolynomials, StaticArrays, QuadGK, Test

@testset "DiskSlice" begin
    P = JacobiDiskSlice(0.0)

    @test P[SVector(0.1,0.2),1] ≈ 1

    @test quadgk(x -> quadgk(y -> P[SVector(x,y),1], -sqrt(1-x^2), sqrt(1-x^2), atol=1e-8)[1], 0, 1, atol=1e-8)[1] ≈ π/2
    @test abs(quadgk(x -> quadgk(y -> P[SVector(x,y),2], -sqrt(1-x^2), sqrt(1-x^2), atol=1e-8)[1], 0, 1, atol=1e-8)[1]) ≤ 1E-8
    @test abs(quadgk(x -> quadgk(y -> P[SVector(x,y),3], -sqrt(1-x^2), sqrt(1-x^2), atol=1e-8)[1], 0, 1, atol=1e-8)[1]) ≤ 1E-8
    @test abs(quadgk(x -> quadgk(y -> P[SVector(x,y),4], -sqrt(1-x^2), sqrt(1-x^2), atol=1e-8)[1], 0, 1, atol=1e-8)[1]) ≤ 1E-8
    @test abs(quadgk(x -> quadgk(y -> P[SVector(x,y),5], -sqrt(1-x^2), sqrt(1-x^2), atol=1e-8)[1], 0, 1, atol=1e-8)[1]) ≤ 1E-8
    @test abs(quadgk(x -> quadgk(y -> P[SVector(x,y),3]P[SVector(x,y),2], -sqrt(1-x^2), sqrt(1-x^2), atol=1e-8)[1], 0, 1, atol=1e-8)[1]) ≤ 1E-8
    @test abs(quadgk(x -> quadgk(y -> P[SVector(x,y),4]P[SVector(x,y),2], -sqrt(1-x^2), sqrt(1-x^2), atol=1e-8)[1], 0, 1, atol=1e-8)[1]) ≤ 1E-8
    @test abs(quadgk(x -> quadgk(y -> P[SVector(x,y),5]P[SVector(x,y),2], -sqrt(1-x^2), sqrt(1-x^2), atol=1e-8)[1], 0, 1, atol=1e-8)[1]) ≤ 1E-8
    @test abs(quadgk(x -> quadgk(y -> P[SVector(x,y),4]P[SVector(x,y),3], -sqrt(1-x^2), sqrt(1-x^2), atol=1e-8)[1], 0, 1, atol=1e-8)[1]) ≤ 1E-8
    @test abs(quadgk(x -> quadgk(y -> P[SVector(x,y),5]P[SVector(x,y),3], -sqrt(1-x^2), sqrt(1-x^2), atol=1e-8)[1], 0, 1, atol=1e-8)[1]) ≤ 1E-8
    @test abs(quadgk(x -> quadgk(y -> P[SVector(x,y),5]P[SVector(x,y),4], -sqrt(1-x^2), sqrt(1-x^2), atol=1e-8)[1], 0, 1, atol=1e-8)[1]) ≤ 1E-8
end
