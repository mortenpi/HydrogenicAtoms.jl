using HydrogenicAtoms
using HydrogenicAtoms.DiracAnalytic: E, α
using Test

@testset "HydrogenicAtoms.jl" begin
    @test E(0, 1, -1) isa Real
    @test E(0, 1, -1) ≈ α^-2
    @test E(1/α, 1, -1) isa Real
    @test E(1/α, 1, -1) ≈ 0.0
    @test_throws DomainError E(138, 1, -1)

    @test E(Complex, 0, 1, -1) isa Complex
    @test E(Complex, 0, 1, -1) ≈ α^-2
    @test E(Complex, 1/α, 1, -1) isa Complex
    @test_broken E(Complex, 1/α, 1, -1) ≈ 0.0
    @test E(Complex, 138, 1, -1) isa Complex
    @test E(Complex, 138, 1, -1) ≈ 0.0 + 2231.3523399043934im
end
