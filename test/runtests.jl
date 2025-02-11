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

    @testset "double gamma function" begin
        @testset "shift equations" begin
            β = 0.74
            w = 1.42
            tol = 1e-8

            v1 = doublegamma(w, β; tol=tol)
            v2 = doublegamma(w + β, β; tol=tol)
            shift = sqrt(2π) * β^(β * w - 1 // 2) / gamma(β * w)

            @test isapprox(v2 / v1, shift, rtol=1e-8)
        end
    end
end
