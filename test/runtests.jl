using BarnesDoubleGamma,
    Test,
    Documenter

# tests included in docs
# DocMeta.setdocmeta!(
#     BarnesDoubleGamma,
#     :DocTestSetup,
#     :(using BarnesDoubleGamma),
#     recursive=true
# )
# Documenter.doctest(BarnesDoubleGamma)

@testset "BarnesDoubleGamma.jl" begin
    @testset "digamma_reg" begin
        for n in 0:-1:-5
            @test digamma_reg(n) == digamma_reg(-n+1)
        end
    end

    @testset "double gamma function" begin
        @testset "shift equations" begin
            β = Double64(0.74)
            w = Double64(1.42)

            v1 = doublegamma(w, β)
            v2 = doublegamma(w + β, β)
            shift = sqrt(2*oftype(β, π)) * β^(β * w - 1 // 2) / gamma(β * w)

            @test isapprox(v2 / v1, shift, rtol=1e-8)
        end

        @testset "exact values" begin
	    @test barnesdoublegamma(1, sqrt(3)) ≈ 1
            @test isapprox(barnesdoublegamma(big(1), sqrt(big(3))), 1, rtol=1e-40)

	    @test isapprox(barnesdoublegamma(1, sqrt(3)), 1)
            @test isapprox(barnesdoublegamma(big(1), sqrt(big(3))), 1)
        end
    end
end
