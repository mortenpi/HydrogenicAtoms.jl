using HydrogenicAtoms
using HydrogenicAtoms.DiracAnalytic: E, α
using AtomicLevels: @ro_str
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

    @test E(0, ro"1s") isa Real
    @test E(0, ro"1s") ≈ α^-2
    @test E(Complex, 0, ro"1s") isa Complex
    @test E(Complex, 0, ro"1s") ≈ α^-2

    @test E(100, ro"2s") ≈ E(Complex, 100, 2, -1)
    @test E(100, ro"3p-") ≈ E(Complex, 100, 3, 1)
    @test E(100, ro"5p") ≈ E(Complex, 100, 5, -2)
end
