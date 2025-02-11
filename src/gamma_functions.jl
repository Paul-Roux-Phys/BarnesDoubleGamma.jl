import SpecialFunctions: gamma, loggamma, digamma, trigamma
# import ArbNumerics.lgamma as arb_loggamma
import ArbNumerics: ArbComplex as ArbC
import ArbNumerics.lgamma
import ArbNumerics: gamma, digamma, polygamma

# Extend gamma, loggamma, trigamma
gamma(z::Complex{BigFloat}) = Complex{BigFloat}(
    gamma(ArbC(z, bits=precision(BigFloat)))
) 
loggamma(z::Complex{BigFloat}) = Complex{BigFloat}(
    lgamma(ArbC(z, bits=precision(BigFloat)))
)

# Extend digamma
digamma(z::Complex{BigFloat}) = Complex{BigFloat}(
    digamma(ArbC(z, bits=precision(BigFloat)))
)

# use ArbNumerics.polygamma to compute trigamma in arbitrary precision.
trigamma(z::Complex{BigFloat}) = Complex{BigFloat}(
    polygamma(ArbC(1, bits=precision(BigFloat)), ArbC(z, bits=precision(BigFloat)))
)
trigamma(z::BigFloat) = real(trigamma(complex(z)))

# Extend polygamma
polygamma(n::Integer, z::Complex{BigFloat}) = Complex{BigFloat}(
    polygamma(ArbC(n, bits=precision(BigFloat)), ArbC(z, bits=precision(BigFloat)))
)

# cotpi(x) = cot(π * x)
cotpi(x) = cospi(x) / sinpi(x)

"""
    digamma_reg(z)

Digamma function ``ψ(z)`` regularised at negative integers thanks to the formula

```math
ψ(1-z) - ψ(z) = π \\operatorname{cot}(πz)
```

# Examples

```jldoctest
julia> digamma_reg(-1)
0.4227843350984672
```
"""
function digamma_reg(z)
    if real(z) > 0
        return digamma(z)
    elseif isreal(z) && real(z) <= 0 && real(z)%1 == 0
        return digamma(1-z)
    else
        return digamma(1-z) - oftype(z, π)*cotpi(z)
    end
end

function precision(z::Complex; base=2)
    return precision(real(z), base=base)
end

function convert_precision(x::Number, precision)
    if precision <= 53 # precision of Float64
        if isreal(x)
            return Float64(Real(x))
        else
            return Complex{Float64}(x)
        end
    else # precision > 53
        if isreal(x)
            return BigFloat(Real(x), precision=precision)
        else
            r = BigFloat(real(x), precision=precision)
            i = BigFloat(imag(x), precision=precision)
            return Complex{BigFloat}(r, i)
        end
    end
end

macro convert_precision!(var, precision)
    quote
        $var = convert_precision($var, $precision)
    end |> esc
end
