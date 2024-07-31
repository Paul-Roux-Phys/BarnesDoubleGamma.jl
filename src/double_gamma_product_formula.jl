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
    res -= log(τ) + loggamma(z)
    res += convert_precision(modularcoeff_a(τ)*z/τ + modularcoeff_b(τ)*z^2/(2*τ^2), prec)
    for m in 1:N
        mτ = m*τ
        d = digamma(mτ)
        t = trigamma(mτ)
        res += loggamma(mτ) - loggamma(z + mτ) + z*d + z^2 / 2 * t
    end
    return res
end

@memoize function factorial_big(n)::BigInt
    return factorial(big(n))
end

"""Polynomial in τ which is the coefficient of z^j in Pn"""
@memoize function polynomial_Pn_coeff(n, j, τ)
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
Numerical approximation of the logarithm of Barne's G-function, up to a given tolerance.
The tolerance is rather pessimistic, on tests I found that a tolerance of 1e-30 gave about
40 correct digits for .
"""
function log_barnes_doublegamma(z, τ, tol)
    d = -log(tol)/log(10)
    M = floor(Int, 0.8*log(10)/log(20)*d)
    N = 30*M
    return log_Barnes_GN(N, z, τ) + z^3*rest_RMN(M, N, z, τ)
end

function barnes_doublegamma(z, τ, tol)
    return exp(log_barnes_doublegamma(z, τ, tol))
end

function log_gamma2(w, β, tol)
    β = real(β-1/β) < 0 ? 1/β : β # change β -> 1/β if needed
    return w/(2*β)*log(2*oftype(β, π)) + (w/2*(w-β-1/β)+1)*log(β) - log_barnes_doublegamma(w/β, 1/β^2, tol)
end

function gamma2(w, β, tol)
    return exp(log_gamma2(w, β, tol))
end

"""
        logdoublegamma(w, β, tol) = Γ_β(w)

Compute the logarithm of the double gamma function Γ_β(w, β) as defined in
[https://en.wikipedia.org/wiki/Multiple_gamma_function](Multiple Gamma function),
with precision tol
"""
function logdoublegamma(w, β, tol)
    return log_gamma2(w, β, tol) - log_gamma2((β+1/β)/2, β, tol)
end

"""
        doublegamma(w, β, tol)

Compute the double gamma function Γ_β(w)as defined in
[https://en.wikipedia.org/wiki/Multiple_gamma_function](Multiple Gamma function),
with precision tol
"""
function doublegamma(w, β, tol)
    exp(logdoublegamma(w, β, tol))
end