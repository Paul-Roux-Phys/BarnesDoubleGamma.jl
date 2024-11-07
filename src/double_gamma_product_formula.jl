#=
# See arXiv:2208.13876
=#

function integrandC(x, τ)
    x = big(x)
    return exp((1-τ)*x)/(2*sinh(x)*sinh(τ*x)) - exp(-2*x)/(τ*x)*(exp(x)/(2*sinh(x))+1-τ/2)
end

function modularC(τ)
    P = precision(BigFloat) ÷ 2
    #temporarily increase precision to avoid artificial divergence around zero
    setprecision(BigFloat, Int(floor(1.3*P)))
    cutoff = big(2^(-P/3)) # to prevent artificial divergence around zero
    tol = inv(big(2))^P

    #compute integral
    value, error = quadgk(x -> integrandC(x, τ), big(cutoff), big(Inf), rtol = tol, order=21)
    C0 = (2/τ - 3//2 + τ/6)*cutoff + (5//12 - 1/τ + τ/12)*cutoff^2 + (4/(9*τ) - 2//9 + 1//54*τ - 1//270*τ^3)*cutoff^3
    setprecision(BigFloat, 2P)

    return 1/(2*τ)*log(2*oftype(τ, π)) - value - C0
end

function integrandD(x, τ)
    x = big(x)
    return x*exp((1-τ)*x)/(sinh(x)*sinh(τ*x)) - exp(-2*x)/(τ*x)
end

function modularD(τ)
    P = precision(BigFloat, base=10)
    #temporarily increase precision to avoid artificial divergence around zero
    setprecision(BigFloat, base=10, Int(floor(1.3*P)))
    cutoff = big(10^(-P/5)) # to prevent artificial divergence around zero
    tol = inv(big(10))^P
    value, error = quadgk( x -> integrandD(x, τ), big(0), big(Inf), rtol = tol, order=21)
    setprecision(BigFloat, base=10, P)
    return value
end

@memoize function modularcoeff_a(τ)
    return convert_precision(1/2*τ*log(2*oftype(τ, π)*τ) + 1/2*log(τ) - τ*modularC(τ), precision(τ))
end

@memoize function modularcoeff_b(τ)
    return convert_precision(-τ*log(τ) - τ^2*modularD(τ), precision(τ))
end

function log_Barnes_GN(N, z, τ)
    # keep only the minimum precision
    prec = min(precision(z), precision(τ))
    @convert_precision!(z, prec)
    @convert_precision!(τ, prec)

    #compute the sum
    res = 0
    res -= log(τ) + loggamma(complex(z))
    res += convert_precision(modularcoeff_a(τ)*z/τ + modularcoeff_b(τ)*z^2/(2*τ^2), prec)
    for m in 1:N
        mτ = m*τ
        d = digamma(mτ)
        t = trigamma(mτ)
        res += loggamma(complex(mτ)) - loggamma(complex(z + mτ)) + z*d + z^2 / 2 * t
    end
    return res
end

function factorial_big(n)::BigInt
    return factorial(big(n))
end

function polynomial_Pn_coeff(n, j, τ)
    if n == 1 && j == 0
        return 1//6
    elseif j >= n
        return 0
    elseif j == n-1
        return inv(factorial_big(n+2))
    else # j < n-1 and n > 1
        return -sum(((1+τ)^(k+2) - 1 - τ^(k+2))/factorial_big(k+2)/τ * polynomial_Pn_coeff(n-k, j, τ) for k in 1:n-1-j)
    end
end

function polynomial_Pn(n, τ)
    prec = precision(τ)
    return [convert_precision(polynomial_Pn_coeff(n, j, τ), prec) for j in 0:n-1]
end

function rest_RMN(M, N, z, τ)
    # keep only the minimum precision
    prec = min(precision(z), precision(τ))
    @convert_precision!(z, prec)
    @convert_precision!(τ, prec)

    mτ = -τ

    p(k::Int) = evalpoly(z, polynomial_Pn(k, mτ))

    coeffs_sum = vcat(0, [factorial_big(k-1) * p(k) for k in 1:M])

    return -inv(τ)*evalpoly(inv(-N*τ), convert_precision.(coeffs_sum, prec))
end

"""
    log_barnesdoublegamma(z, τ; tol = 10^(-precision(z)))

Logarithm of Barne's G-function ``\\log(G(z; τ))``, up to a given tolerance.
The tolerance is rather pessimistic. For efficiency, when using a high precision `z` it is
recommended to set the tolerance to something lower than the precision of `z`.

# Examples

```jldoctest
julia> z = 1; τ = sqrt(3); log_barnesdoublegamma(z, τ)
4.7036054521762405e-12

julia> z = sqrt(big"2"); τ = sqrt(big"3"); log_barnesdoublegamma(z, τ)
0.293394920968626643183216869255154162603276161914057065216580158409738193128475
```
"""
function log_barnesdoublegamma(z, τ, tol)
    d = -log(tol)/log(10)
    M = floor(Int, 0.8*log(10)/log(20)*d)
    N = 30*M
    return log_Barnes_GN(N, z, τ) + z^3*rest_RMN(M, N, z, τ)
end

log_barnesdoublegamma(z::Union{Float64, ComplexF64}, τ::Number; tol=1e-16) =
    log_barnesdoublegamma(z, τ, tol)

log_barnesdoublegamma(
    z::Union{BigFloat, Complex{BigFloat}},
    τ::Number;
    tol=1/big"10"^(precision(z, base=10)/3)
) = log_barnesdoublegamma(z, τ, tol)

log_barnesdoublegamma(z::Real, τ::Real, tol) =
    real(log_barnesdoublegamma(complex(z), complex(τ), tol))

log_barnesdoublegamma(z::Number, τ::Number) = log_barnesdoublegamma(float(z), float(τ))

"""
    barnesdoublegamma(z, τ; tol=10^(-precision(z)))

Barne's G-function ``G(z, τ)``, up to a given tolerance.
For efficiency, when using a high precision `z` it is
recommended to set the tolerance to something lower than the precision of `z`.

# Examples

```jldoctest
julia> z = big"1"; τ = sqrt(big"3"); barnesdoublegamma(z, τ)
0.9999999999999999999999999999999999999999999345543169494792740822465202042518853

julia> z = sqrt(big"2"); τ = sqrt(big"3"); barnesdoublegamma(z, τ)
1.340972263940081256497568500074283394055091669857109294097794355305731639439863

julia> z = big(sqrt(3)); τ = sqrt(big"3"); barnesdoublegamma(z, τ)
1.488928335365086422942328604671778776079655470648475974888882080914824321973775
```

"""
barnesdoublegamma(z, τ, tol) =
    exp(log_barnesdoublegamma(z, τ, tol))

barnesdoublegamma(z, τ) = exp(log_barnesdoublegamma(z, τ))

"""
    loggamma2(w, β; tol=10^(-precision(z)))

Logarithm of the ``Γ_2(w, β)`` function.
For efficiency, when using a high precision `w` it is
recommended to set the tolerance to something lower than the precision of `w`.
"""
function loggamma2(w, β; tol=missing)
    β = real(β-1/β) < 0 ? 1/β : β # change β -> 1/β if needed
    if tol !== missing
        l = log_barnesdoublegamma(w / β, 1 / β^2, tol)
    else
        l = log_barnesdoublegamma(w / β, 1 / β^2)
    end
    return w/(2*β)*log(2*oftype(β, π)) + (w/2*(w-β-1/β)+1)*log(β) - l
end

"""
    gamma2(w, β; tol=10^(-precision(w)))

``Γ_2(w, β)`` function.
For efficiency, when using a high precision `w` it is
recommended to set the tolerance to something lower than the precision of `w`.
"""
function gamma2(w, β; tol=missing)
    return exp(loggamma2(w, β; tol=tol))
end

"""
        logdoublegamma(w, β; tol=10^(-precision(w)))

Compute the logarithm of the double gamma function ``Γ_β(w, β)`` with precision tol
"""
function logdoublegamma(w, β; tol=missing)
    return loggamma2(w, β, tol=tol) - loggamma2((β+1/β)/2, β, tol=tol)
end

"""
        doublegamma(w, β; tol=10^(-precision(w)))

Compute the double gamma function ``Γ_β(w)`` with precision tol
"""
function doublegamma(w, β; tol=missing)
    exp(logdoublegamma(w, β, tol=tol))
end