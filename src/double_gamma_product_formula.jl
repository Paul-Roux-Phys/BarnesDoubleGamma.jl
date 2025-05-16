#=
# See arXiv:2208.13876
=#

function integrandC(x, τ)
    one_minus_τ = 1 - τ
    sinh_x = sinh(x)
    sinh_τx = sinh(τ * x)
    exp1 = exp(one_minus_τ * x)
    exp2 = exp(-2 * x)
    exp3 = exp(x)

    return exp1 / (2 * sinh_x * sinh_τx) -
           exp2 / (τ * x) * (exp3 / (2 * sinh_x) + 1 - τ / 2)
end

function modularC(τ::T)::Tuple{T, T} where {T}
    prec = precision(real(τ))
    U = BigFloat
    P = prec ÷ 2
    cutoff = inv(U(2)^(P/3))
    tol = inv(U(2))^P

    value, err, C0 = setprecision(U, Int(floor(1.3 * P))) do
        v, err = quadgk(x -> integrandC(x, τ), cutoff, U(Inf), rtol=tol, order=21)

        c0 = (2 / τ - 3//2 + τ / 6) * cutoff +
             (5//12 - 1 / τ + τ / 12) * cutoff^2 +
             (4 / (9 * τ) - 2//9 + 1//54 * τ - 1//270 * τ^3) * cutoff^3

        return v, err, c0
    end

    return (1 / (2 * τ)) * log(2 * oftype(τ, π)) - value - C0, err
end

function integrandD(x, τ)
    return x*exp((1-τ)*x)/(sinh(x)*sinh(τ*x)) - exp(-2*x)/(τ*x)
end

function modularD(τ::T)::Tuple{T, T} where {T}
    P = precision(real(τ))
    #temporarily increase precision to avoid artificial divergence around zero
    tol = inv(big(2))^P
    value, error = setprecision(BigFloat, Int(floor(1.3 * P))) do
        quadgk(x -> integrandD(x, τ), big(0), big(Inf), rtol=tol, order=21)
    end
    return value, error
end

@memoize function modularcoeff_a(τ)
    modC, _ = modularC(τ)
    return 1/2*τ*log(2*(π*τ)) + 1/2*log(τ) - τ*modC
end

@memoize function modularcoeff_b(τ)
    modD, _ = modularD(τ)
    return -τ*log(τ) - τ^2*modD
end

function log_Barnes_GN(N, z, τ::T)::T where {T}
    # keep only the minimum precision
    #compute the sum
    r = 0
    r -= log(τ) + loggamma(complex(z))
    r += modularcoeff_a(τ) * z / τ + modularcoeff_b(τ) * z^2 / (2 * τ^2)
    for m in 1:N
        mτ = m * τ
        d = digamma(mτ)
        t = trigamma(mτ)
        r += loggamma(complex(mτ)) - loggamma(complex(z + mτ)) +
             z * d + z^2 / 2 * t
    end
    return r
end

@memoize factorial_big(n) = if n < 21 factorial(n) else factorial(big(n)) end

@memoize function polynomial_Pn_coeff(n, j, τ::T)::T where {T}
    n == 1 && j == 0 && return one(τ)/6
    j >= n && return zero(τ)
    j == n-1 && return one(τ)/(factorial_big(n+2))
    return -sum(
        ((1 + τ)^(k + 2) - 1 - τ^(k + 2)) / factorial_big(k + 2) / τ *
        polynomial_Pn_coeff(n - k, j, τ)
        for k in 1:n-1-j
    )
end

polynomial_Pn(n, τ) = [polynomial_Pn_coeff(n, j, τ) for j in 0:n-1]

function rest_RMN(M, N, z, τ::T)::T where {T}
    mτ = -τ
    p(k::Int) = evalpoly(z, polynomial_Pn(k, mτ))
    coeffs_sum = vcat(0, [factorial_big(k-1) * p(k) for k in 1:M])
    return -inv(τ)*evalpoly(inv(-N*τ), coeffs_sum)
end

"""
    log_barnesdoublegamma(z, τ)

Logarithm of Barne's G-function ``\\log(G(z; τ))``.
Can get very expensive for high precision.

# Examples

```jldoctest
julia> z = 1; τ = sqrt(3); log_barnesdoublegamma(z, τ)
4.7036054521762405e-12

julia> z = sqrt(big"2"); τ = sqrt(big"3"); log_barnesdoublegamma(z, τ)
0.293394920968626643183216869255154162603276161914057065216580158409738193128475
```
"""
function _log_barnesdoublegamma(z::Complex, τ::Complex)
    d = precision(real(τ))
    M = floor(Int, 0.5/log(20)*d)
    N = 50*M
    return log_Barnes_GN(N, z, τ) + z^3*rest_RMN(M, N, z, τ)
end

_log_barnesdoublegamma(z::Real, τ::Real) =
    real(_log_barnesdoublegamma(complex(z), complex(τ)))

_log_barnesdoublegamma(z, τ::Complex) =
    _log_barnesdoublegamma(complex(z), complex(τ))

_log_barnesdoublegamma(z::Complex, τ) =
    _log_barnesdoublegamma(complex(z), complex(τ))

log_barnesdoublegamma(z::Number, τ::Number) = _log_barnesdoublegamma(float(z), float(τ))

"""
    barnesdoublegamma(z, τ)

Barne's G-function ``G(z, τ)``.
Can get very expensive for high precision.

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
barnesdoublegamma(z, τ) =
    exp(log_barnesdoublegamma(z, τ))

# function barnesdoublegamma(z, τ)
#     if tol !== missing
#         exp(log_barnesdoublegamma(z, τ, tol))
#     else
#         exp(log_barnesdoublegamma(z, τ))
#     end
# end

"""
    loggamma2(w, β)

Logarithm of the ``Γ_2(w, β)`` function.
"""
function loggamma2(w, β)
    β = real(β-1/β) < 0 ? 1/β : β # change β -> 1/β if needed
    l = log_barnesdoublegamma(w / β, 1 / β^2)
    return w/(2*β)*log(2*oftype(β, π)) + (w/2*(w-β-1/β)+1)*log(β) - l
end

"""
    gamma2(w, β)

``Γ_2(w, β)`` function.
"""
function gamma2(w, β)
    return exp(loggamma2(w, β))
end

"""
        logdoublegamma(w, β)

Compute the logarithm of the double gamma function ``Γ_β(w, β)``.
"""
function logdoublegamma(w, β)
    return loggamma2(w, β) - loggamma2((β+1/β)/2, β)
end

"""
        doublegamma(w, β)

Compute the double gamma function ``Γ_β(w)``.
"""
function doublegamma(w, β)
    exp(logdoublegamma(w, β))
end
