using BarnesDoubleGamma,
    Test,
    Documenter

# tests included in docs
DocMeta.setdocmeta!(
    BarnesDoubleGamma,
    :DocTestSetup,
    :(using BarnesDoubleGamma),
    recursive=true
)
Documenter.doctest(BarnesDoubleGamma)

@testset "BarnesDoubleGamma.jl" begin

    @testset "digamma_reg" begin
        for n in 0:-1:-5
            @test digamma_reg(n) == digamma_reg(-n+1)
        end
    end

    @testset "BDoubleGamma" begin
        bdg = BDoubleGamma(sqrt(3))
        bdg_big = BDoubleGamma(sqrt(big(3)))
        
        @test bdg(1) ≈ 1
        @test bdg_big(big"1.") ≈ 1

        @test barnesdoublegamma(1, sqrt(3)) ≈ 1
        @test barnesdoublegamma(big(1), sqrt(big(3))) ≊ 1
    end

    @testset "double gamma function" begin
        @testset "shift equations" begin
            p = precision(BigFloat)
            setprecision(BigFloat, 128)
            β = big"0.74"
            w = big"1.42"

            dg = DoubleGamma(β)
            v1 = dg(w)
            v2 = dg(w + β)
            shift = sqrt(2 * oftype(β, π)) * β^(β * w - 1 // 2) / gamma(β * w)

            @test isapprox(v2 / v1, shift, rtol=1e-8)
            setprecision(BigFloat, p)
        end

    end
end
