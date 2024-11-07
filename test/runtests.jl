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
    @testset "shift equations" begin
        β = 0.74
        w = 1.42
        tol = 1e-8
        
        v1 = doublegamma(w, β, tol)
        v2 = doublegamma(w+β, β, tol)
        shift = sqrt(2π) * β^(β*w - 1//2)/gamma(β*w)

        @test isapprox(v2/v1, shift, rtol=1e-8)
    end

    # function Cref(β, ri, si)
    #     return prod(
    #         inv(doublegamma(
    #             (β+inv(β))/2 + β/2*abs(dot(e, ri)) + inv(β)/2*dot(e, si),
    #             β,
    #             1e-16
    #         ))
    #         for e in product((-1, 1), (-1, 1), (-1, 1))
    #     )
    # end
    # β = sqrt(3/4)
    # β2P(r, s) = (r*β^2-1)
    # println(Cref(β, (0, 0, 0), (β2P(1, 2), β2P(1, 2), β2P(1, 2))))
end
