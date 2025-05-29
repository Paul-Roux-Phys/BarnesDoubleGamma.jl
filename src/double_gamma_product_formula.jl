#=
# See arXiv:2208.13876
=#

#====================================================
Cache
=====================================================#
struct BDGCache{T}
    M::Int
    N::Int
    τ::T
    logτ::T
    mτ_pows::Vector{T}
    mτp1_pows::Vector{T}
    loggammas::Vector{T}
    digammas::Vector{T}
    trigammas::Vector{T}
    factorials::Vector{BigInt}
    polynomial_Pns::Matrix{T}
    modularcoeff_a::T
    modularcoeff_b::T
    
end

function BDGCache(τ::T) where {T <: Complex}
    d = precision(real(τ))
    M = floor(Int, 0.5/log(20)*d)
    N = 50*M
    logτ = log(τ)
    mτ_pows = zeros(T, M+2)
    mτp1_pows = zeros(T, M+2)
    mτ_pows[1] = -τ
    mτp1_pows[1] = -τ + 1
    for i in 2:M+2
        mτ_pows[i] = -τ * mτ_pows[i-1]
        mτp1_pows[i] = (-τ+1) * mτp1_pows[i-1]
    end
    loggammas = [loggamma(m*τ) for m in 1:N]
    digammas = [digamma(m*τ) for m in 1:N]
    trigammas = [trigamma(m*τ) for m in 1:N]
    factorials = [factorial(big(n)) for n in 1:M+2]
    poly = polynomial_Pns(M, mτ_pows, mτp1_pows, factorials)
    moda = modularcoeff_a(τ)
    modb = modularcoeff_b(τ)
    return BDGCache{T}(
        M, N,
        τ, logτ, mτ_pows, mτp1_pows,
        loggammas, digammas, trigammas,
        factorials, poly, moda, modb
    )
end

BDGCache(τ::Real) = BDGCache(complex(τ))

#===============================================================================
Barnes Double Gamma G(z, τ)
===============================================================================#
struct LogBDoubleGamma{T}
    cache::BDGCache{T}
end

LogBDoubleGamma(τ::T) where {T} = LogBDoubleGamma(BDGCache(τ))
(f::LogBDoubleGamma)(z::Complex) = _log_barnesdoublegamma(z, f.cache)

struct BDoubleGamma{T}
    logBDG::LogBDoubleGamma{T}
end

BDoubleGamma(τ::T) where {T} = BDoubleGamma(LogBDoubleGamma(τ))
(f::BDoubleGamma)(z::Complex) = exp(f.logBDG(z))

#===============================================================================
Barnes Gamma2 Γ_2(w, β)
===============================================================================#
struct LogGamma2{T}
    logBDG::LogBDoubleGamma{T}
    β::T
end

function LogGamma2(β::T) where {T}
    β = real(β - 1/β) < 0 ? inv(β) : β
    τ = inv(β^2)
    return LogGamma2(LogBDoubleGamma(τ), β)
end

function (f::LogGamma2)(w)
    β = f.β
    l = f.logBDG(w / β) # log_barnesdoublegamma(w / β, 1/β^2)
    return w/(2*β)*log(2*oftype(β, π)) + (w/2*(w-β-1/β)+1)*log(β) - l
end 

struct Gamma2{T}
    loggamma2::LogGamma2{T}
end

Gamma2(β::T) where {T} = Gamma2{T}(LogGamma2(β))

#===============================================================================
Double Gamma Γ_β(w)
===============================================================================#
struct LogDoubleGamma{T}
    loggamma2::LogGamma2{T}
    refval::T # value at (β + 1/β)/2
end

function LogDoubleGamma(β::T) where {T}
    lg = LogGamma2(β)
    refval = lg((β + inv(β)) / 2)
    return LogDoubleGamma{T}(lg, refval)
end

(f::LogDoubleGamma)(w) = f.loggamma2(w) - f.refval

struct DoubleGamma{T}
    inner::LogDoubleGamma{T}
end

DoubleGamma(β::T) where {T} = DoubleGamma{T}(LogDoubleGamma(β))
(f::DoubleGamma)(w::T) where {T} = exp(f.inner(w))

#===============================================================================
Implementation
===============================================================================#
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

function modularcoeff_a(τ)
    modC, _ = modularC(τ)
    return 1/2*τ*log(2*(π*τ)) + 1/2*log(τ) - τ*modC
end

function modularcoeff_b(τ)
    modD, _ = modularD(τ)
    return -τ*log(τ) - τ^2*modD
end

function log_Barnes_GN(N, z, cache::BDGCache{T}) where {T}
    # keep only the minimum precision
    #compute the sum
    τ = cache.τ
    r = zero(T)
    r -= cache.logτ + loggamma(z)
    r += cache.modularcoeff_a * z / cache.τ + cache.modularcoeff_b * z^2 / (2 * τ^2)
    for m in 1:N
        mτ = m * τ
        d = cache.digammas[m]
        t = cache.trigammas[m]
        r += cache.loggammas[m] - loggamma(z + mτ) + z * d + z^2 / 2 * t
    end
    return r
end

function polynomial_Pns(M, τ_pows::Vector{T}, τp1_pows, factorials)::Matrix{T} where {T}
    coeffs = zeros(T, M+1, M+1)
    τfactors = [
        (τp1_pows[k+2] - 1 - τ_pows[k+2]) /
        factorials[k + 2] / τ_pows[1]
        for k in 1:M
    ]
    for n in 1:M+1
        coeffs[n, n] = one(T)/factorials[n+1]
        for j in 1:n-2
            acc = zero(T)
            for k in 1:(n-1-j)
                acc += τfactors[k] * coeffs[n - k, j+1]
            end
            coeffs[n, j+1] = -acc
        end
    end
    return coeffs
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

function rest_RMN(M, N, z, cache::BDGCache{T})::T where {T}
    τ = cache.τ
    coeffs_sum = zeros(typeof(z), M + 1)
    coeffs_sum[2] = evalpoly(z, cache.polynomial_Pns[2, :])
    for k in 2:M
        coeffs_sum[k+1] = cache.factorials[k - 1] * evalpoly(z, cache.polynomial_Pns[k+1, :])
    end
    return -evalpoly(inv(-N * τ), coeffs_sum) / τ
end

function _log_barnesdoublegamma(z::Complex, cache::BDGCache)
    M, N = cache.M, cache.N
    return log_Barnes_GN(N, z, cache) + z^3*rest_RMN(M, N, z, cache)
end

# function _log_barnesdoublegamma(z::Complex, τ::T) where {T<:Complex}
#     cache = BDGCache(τ)
#     return _log_barnesdoublegamma(z, cache)
# end

# _log_barnesdoublegamma(z::Real, τ::Real) =
#     real(_log_barnesdoublegamma(complex(z), complex(τ)))

# _log_barnesdoublegamma(z, τ::Complex) =
#     _log_barnesdoublegamma(complex(z), complex(τ))

# _log_barnesdoublegamma(z::Complex, τ) =
#     _log_barnesdoublegamma(complex(z), complex(τ))

# """
#     log_barnesdoublegamma(z, τ)

# Logarithm of Barne's G-function ``\\log(G(z; τ))``.
# Can get very expensive for high precision.

# # Examples

# ```jldoctest
# julia> z = 1; τ = sqrt(3); log_barnesdoublegamma(z, τ)
# -3.5564013784958066e-9

# julia> z = sqrt(big"2"); τ = sqrt(big"3"); log_barnesdoublegamma(z, τ)
# 0.293394920968626643183216869255154162603276275448888004730390602371621786480874
# ```
# """
# log_barnesdoublegamma(z::Number, τ::Number) = _log_barnesdoublegamma(float(z), float(τ))

# """
#     barnesdoublegamma(z, τ)

# Barne's G-function ``G(z, τ)``.
# Can get very expensive for high precision.

# # Examples

# ```jldoctest
# julia> z = big"1"; τ = sqrt(big"3"); barnesdoublegamma(z, τ)
# 0.9999999999999999999999999999999999999999999999265240691883395187060685710095162

# julia> z = sqrt(big"2"); τ = sqrt(big"3"); barnesdoublegamma(z, τ)
# 1.340972263940081256497568500074283394055091822104168575112011391955491855627026

# julia> z = big(sqrt(3)); τ = sqrt(big"3"); barnesdoublegamma(z, τ)
# 1.488928335365086422942328604671778776079655676458875630387076246079377268516627
# ```
# """
# barnesdoublegamma(z, τ) =
#     exp(log_barnesdoublegamma(z, τ))

# function barnesdoublegamma(z, τ)
#     if tol !== missing
#         exp(log_barnesdoublegamma(z, τ, tol))
#     else
#         exp(log_barnesdoublegamma(z, τ))
#     end
# end

# function loggamma2(w, β::LogGamma2)
#     l = _log_barnesdoublegamma(w / β, β.cache)
#     return w/(2*β)*log(2*oftype(β, π)) + (w/2*(w-β-1/β)+1)*log(β) - l
# end

# """
#     loggamma2(w, β)

# Logarithm of the ``Γ_2(w, β)`` function.
# """
# function loggamma2(w, β)
#     β = real(β-1/β) < 0 ? 1/β : β # change β -> 1/β if needed
#     l = log_barnesdoublegamma(w / β, 1 / β^2)
#     return w/(2*β)*log(2*oftype(β, π)) + (w/2*(w-β-1/β)+1)*log(β) - l
# end

# function gamma2(w, β::Gamma2{T}) where {T}
#     β = LogGamma2{T}(β.cache, β.β)
#     return exp(loggamma2(w, β))
# end

# """
#     gamma2(w, β)

# ``Γ_2(w, β)`` function.
# """
# function gamma2(w, β)
#     return exp(loggamma2(w, β))
# end

# function logdoublegamma(w, β::LogDoubleGamma{T}) where {T}
#     lg2 = LogGamma2{T}(β.cache, β.β)
#     return loggamma2(w, lg2) - loggamma2(
# end

# """
#         logdoublegamma(w, β)

# Compute the logarithm of the double gamma function ``Γ_β(w, β)``.
# """
# function logdoublegamma(w, β)
#     return loggamma2(w, β) - loggamma2((β+1/β)/2, β)
# end

# """
#         doublegamma(w, β)

# Compute the double gamma function ``Γ_β(w)``.
# """
# function doublegamma(w, β)
#     exp(logdoublegamma(w, β))
# end
