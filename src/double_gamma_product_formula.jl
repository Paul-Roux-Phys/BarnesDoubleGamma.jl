#=
# See arXiv:2208.13876
=#

#====================================================
Cache
=====================================================#
struct BDGCache
    M::Int
    N::Int
    τ::Acb
    logτ::Acb
    mτ_pows::Vector{Acb} # (-τ)^n
    mτp1_pows::Vector{Acb} # (-τ+1)^n
    Nmτ_pows::Vector{Acb} # (N - τ)^n
    loggammas::Vector{Acb}
    digammas::Vector{Acb}
    trigammas::Vector{Acb}
    factorials::Vector{BigInt}
    polynomial_Pns::Vector{Vector{Acb}}
    modularcoeff_a::Acb
    modularcoeff_b::Acb
end

function BDGCache(τ::Acb)
    if real(τ) <= 0
        throw(DomainError(τ, "Accepted range is Re(τ) > 0"))
    end
    d = precision(real(τ))
    M = 10
    N = 50*M
    logτ = log(τ)
    mτ_pows = Vector{Acb}(undef, M+2)
    mτp1_pows = Vector{Acb}(undef, M+2)
    Nmτ_pows = Vector{Acb}(undef, M+2)
    mτ_pows[1] = -τ
    mτp1_pows[1] = -τ + 1
    Nmτ_pows[1] = 1
    for i in 2:M+2
        mτ_pows[i] = -τ * mτ_pows[i-1]
        mτp1_pows[i] = (-τ+1) * mτp1_pows[i-1]
        Nmτ_pows[i] = - 1 / (N * τ) * Nmτ_pows[i-1]
    end
    loggammas = [loggamma(m*τ) for m in 1:N]
    digammas = [digamma(m*τ) for m in 1:N]
    trigammas = [trigamma(m*τ) for m in 1:N]
    factorials = [factorial(big(n)) for n in 1:M+2]
    poly = polynomial_Pns(M, mτ_pows, mτp1_pows, factorials)
    moda = modularcoeff_a(τ)
    modb = modularcoeff_b(τ)
    return BDGCache(
        M, N,
        τ, logτ, mτ_pows, mτp1_pows, Nmτ_pows,
        loggammas, digammas, trigammas,
        factorials, poly, moda, modb
    )
end

BDGCache(τ::Real) = BDGCache(complex(τ))

#===============================================================================
Barnes Double Gamma G(z, τ)
===============================================================================#
struct LogBDoubleGamma
    cache::BDGCache
end

LogBDoubleGamma(τ::Acb) = LogBDoubleGamma(BDGCache(τ))
LogBDoubleGamma(τ::Complex) = LogBDoubleGamma(Acb(τ))
(f::LogBDoubleGamma)(z::Union{Complex,Acb}) = _log_barnesdoublegamma(z, f.cache)
(f::LogBDoubleGamma)(x::Real) = f(complex(x))
(f::LogBDoubleGamma)(x::Union{Int, BigInt, Complex{Int}, Complex{BigInt}}) = f(float(x))
function (f::LogBDoubleGamma)(xs::AbstractVector{<:Complex})
    cache = f.cache
    res = Vector{Acb}(undef, length(xs))
    Threads.@threads for i in eachindex(xs)
        @inbounds res[i] = _log_barnesdoublegamma(xs[i], cache)
    end
    return res
end

struct BDoubleGamma
    logBDG::LogBDoubleGamma
end

BDoubleGamma(τ::Number) = BDoubleGamma(LogBDoubleGamma(τ))
(f::BDoubleGamma)(z) = exp(f.logBDG(z))
(f::BDoubleGamma)(z::AbstractVector) = exp.(f.logBDG(z))

#===============================================================================
Barnes Gamma2 Γ_2(w, β)
===============================================================================#
struct LogGamma2
    logBDG::LogBDoubleGamma
    β::Acb
    invβ::Acb
    logβ::Acb
    c1::Acb
    c2::Acb
end

function LogGamma2(β::Acb)
    β = real(β - 1/β) < 0 ? inv(β) : β
    invβ = 1/β
    τ = invβ^2
    logβ = log(β)
    c1 = 1/2*invβ*(log(2*convert(Acb, π)))
    c2 =  -β-invβ
    return LogGamma2(LogBDoubleGamma(τ), β, invβ, logβ, c1, c2)
end

function (f::LogGamma2)(w)
    β = f.β
    l = f.logBDG.(w ./ β) # log_barnesdoublegamma(w / β, 1/β^2)
    logβ = f.logβ
    @. return w * f.c1 + (w / 2 * (w + f.c2) + 1) * logβ - l
end 

struct Gamma2
    loggamma2::LogGamma2
end

Gamma2(β::Number) = Gamma2(LogGamma2(β))
(f::Gamma2)(w) = exp(f.loggamma2(w))
(f::Gamma2)(w::AbstractVector) = exp.(f.loggamma2(w))

#===============================================================================
Double Gamma Γ_β(w)
===============================================================================#
struct LogDoubleGamma
    loggamma2::LogGamma2
    refval::Acb # value at (β + 1/β)/2
end

function LogDoubleGamma(β::Acb)
    lg = LogGamma2(β)
    refval = lg((β + inv(β)) / 2)
    return LogDoubleGamma(lg, refval)
end

(f::LogDoubleGamma)(w) = f.loggamma2(w) - f.refval
(f::LogDoubleGamma)(w::AbstractVector) = f.loggamma2(w) .- f.refval

struct DoubleGamma
    inner::LogDoubleGamma
end

DoubleGamma(β::Number) = DoubleGamma(LogDoubleGamma(β))
(f::DoubleGamma)(w) = exp(f.inner(w))
(f::DoubleGamma)(w::AbstractVector) = exp.(f.inner(w))

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

function modularC(τ::Complex{BigFloat})
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

function modularC(τ::Acb)
    val, err = modularC(Complex{BigFloat}(τ))
    return Acb(val), err
end

function integrandD(x, τ)
    return x*exp((1-τ)*x)/(sinh(x)*sinh(τ*x)) - exp(-2*x)/(τ*x)
end

function modularD(τ::Complex{BigFloat})
    P = precision(real(τ))
    # temporarily increase precision to avoid artificial divergence around zero
    tol = inv(big(2))^P
    value, error = setprecision(BigFloat, Int(floor(1.3 * P))) do
        quadgk(x -> integrandD(x, τ), big(0), big(Inf), rtol=tol, order=21)
    end
    return value, error
end

function modularD(τ::Acb)
    val, err = modularD(Complex{BigFloat}(τ))
    return Acb(val), err
end

function modularcoeff_a(τ)
    modC, _ = modularC(τ)
    return 1/2*τ*log(2*(π*τ)) + 1/2*log(τ) - τ*modC
end

function modularcoeff_b(τ)
    modD, _ = modularD(τ)
    return -τ*log(τ) - τ^2*modD
end

function log_Barnes_GN(z, cache::BDGCache)
    #compute the sum
    N = cache.N
    τ = cache.τ
    res = zero(τ)
    buf = zero(τ)
    Arblib.lgamma!(res, z)
    Arblib.sub!(res, res, cache.logτ)
    res += cache.modularcoeff_a * z / cache.τ + cache.modularcoeff_b * z^2 / (2 * τ^2)
    # the following loop is by far the bottleneck of the program, by one order of magnitude.
    # it is the only part that can't be pre-computed easily (depends on z and τ)
    # the bottleneck is the computation of loggamma(z + m τ)
    # for m in 1:N
    #     Arblib.lgamma!(lgams[m], z + m * τ)
    # end

    for m in 1:N
        # d = cache.digammas[m]
        # t = cache.trigammas[m]
        Arblib.lgamma!(buf, z + m * τ)
        Arblib.add!(res, res, cache.loggammas[m])
        Arblib.sub!(res, res, buf)
        Arblib.addmul!(res, z, cache.digammas[m])
        Arblib.addmul!(res, z^2 / 2, cache.trigammas[m])
        # res += cache.loggammas[m] - lgams[m] + z * d + z^2 / 2 * t
    end
    return res
end

function polynomial_Pns(M, τ_pows::Vector{Acb}, τp1_pows, factorials)::Vector{Vector{Acb}}
    coeffs = [zeros(Acb, n) for n in 1:M]
    τfactors = [
        (τp1_pows[k+2] - 1 - τ_pows[k+2]) /
        factorials[k + 2] / τ_pows[1]
        for k in 1:M
    ]
    for n in 1:M
        coeffs[n][n] = one(Acb)/factorials[n+2]
        for j in 0:n-2
            acc = zero(Acb)
            for k in 1:(n-1-j)
                acc += τfactors[k] * coeffs[n - k][j+1]
            end
            coeffs[n][j+1] = -acc
        end
    end
    return coeffs
end

function rest_RMN(z, cache::BDGCache)::Acb
    M, N = cache.M, cache.N
    mτ = - cache.τ
    coeffs_sum = zeros(typeof(z), M + 1)
    coeffs_sum[2] = evalpoly(z, cache.polynomial_Pns[1])
    for k in 2:M
        coeffs_sum[k+1] = cache.factorials[k - 1] * evalpoly(z, cache.polynomial_Pns[k])
    end
    return sum(
        cache.Nmτ_pows[i] * coeffs_sum[i]
        for i in eachindex(coeffs_sum)
    ) / mτ
end

function _log_barnesdoublegamma(z::Union{Complex,Acb}, cache::BDGCache)
    N = cache.N
    return log_Barnes_GN(z, cache) + z^3*rest_RMN(z, cache)
end

function _log_barnesdoublegamma(zs::AbstractArray{<:Union{Complex,Acb}}, cache::BDGCache)
    N = cache.N
    return [_log_barnesdoublegamma(z, cache) for z in zs]
end
   
"""
    log_barnesdoublegamma(z, τ)

Logarithm of Barne's G-function ``\\log(G(z; τ))``.
Can get very expensive for high precision.

# Examples

```jldoctest
julia> z = 1; τ = sqrt(3); log_barnesdoublegamma(z, τ) ≈ -3.5564013784958066e-9
true

julia> z = sqrt(big"2"); τ = sqrt(big"3"); log_barnesdoublegamma(z, τ) ≈ 0.293394920968626643183216869255154162603276275448888004730390602371621786480874
true
```
"""
log_barnesdoublegamma(z, τ) = LogBDoubleGamma(τ)(z)

"""
    barnesdoublegamma(z, τ)

Barne's G-function ``G(z, τ)``.
Can get very expensive for high precision.

# Examples

```jldoctest
julia> z = big"1"; τ = sqrt(big"3"); barnesdoublegamma(z, τ) ≈ 1
true

julia> z = sqrt(big"2"); τ = sqrt(big"3"); barnesdoublegamma(z, τ) ≈ 1.340972263940081256497568500074283394055091822104168575112011391955491855627026
true

julia> s3 = sqrt(big"3"); z = s3; τ = s3; barnesdoublegamma(z, τ) ≈ (2big(π))^((τ-1)/2)/sqrt(τ)
true
```
"""
barnesdoublegamma(z, τ) = BDoubleGamma(τ)(z)
 
"""
    loggamma2(w, β)

Logarithm of the ``Γ_2(w, β)`` function.
"""
loggamma2(w, β) = LogGamma2(β)(w)

"""
    gamma2(w, β)

``Γ_2(w, β)`` function.
"""
gamma2(w, β) = exp(loggamma2(w, β))

"""
        logdoublegamma(w, β)

Compute the logarithm of the double gamma function ``Γ_β(w, β)``.
"""
logdoublegamma(w, β) = loggamma2(w, β) - loggamma2((β+1/β)/2, β)

"""
        doublegamma(w, β)

Compute the double gamma function ``Γ_β(w)``.
"""
doublegamma(w, β) = exp(logdoublegamma(w, β))

function Base.show(io::IO, DG::DoubleGamma)
    print(io, "Γ(β = $(DG.inner.loggamma2.β))")
end
