import SpecialFunctions: loggamma, gamma, digamma, trigamma, polygamma

for func in (:loggamma, :gamma, :digamma, :trigamma)
    @eval function SpecialFunctions.$func(z::Complex{BigFloat})
        az = Acb(z, prec=precision(real(z)))
        res = SpecialFunctions.$func(az)
        return Complex{BigFloat}(real(res), imag(res))
    end
    @eval SpecialFunctions.$func(z::Complex{BigInt}) = SpecialFunctions(float(z))
end

function SpecialFunctions.polygamma(n::Integer, z::Complex{BigFloat})
    az = Acb(z, prec=precision(real(z)))
    an = Acb(n, prec=precision(real(z)))
    res = SpecialFunctions.polygamma(an, az)
    return Complex{BigFloat}(real(res), imag(res))
end

function SpecialFunctions.polygamma(n::Complex, z::Complex)
    az = Acb(z, prec=precision(real(z)))
    an = Acb(n, prec=precision(real(z)))
    res = SpecialFunctions.polygamma(an, az)
    return Complex{typeof(real(z))}(real(res), imag(res))
end

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
