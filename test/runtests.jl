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

    @testset "BDoubleGamma" begin
        bdg = BDoubleGamma(sqrt(3))
        lbdg = bdg.logBDG
        bdg_big = BDoubleGamma(sqrt(big(3)))
        lbdg_big = bdg_big.logBDG
        
        @test bdg(1) ≈ 1
        @test real(lbdg(1)) < 1e-8
        @test bdg_big(big"1.") ≈ 1

        @test barnesdoublegamma(1, sqrt(3)) ≈ 1
        @test barnesdoublegamma(big(1), sqrt(big(3))) ≈ 1

        log_exact_val(τ) = (τ-1)/2 * log(2oftype(τ, π)) - log(τ)/2
        exact_val(τ) = (2oftype(τ, π))^((τ-1)/2)/sqrt(τ)

        @test bdg(sqrt(3)) ≈ exact_val(sqrt(3))
        @test bdg_big(sqrt(big(3))) ≈ exact_val(sqrt(big(3)))
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
