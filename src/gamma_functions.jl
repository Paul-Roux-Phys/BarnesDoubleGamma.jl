import SpecialFunctions: loggamma, gamma, digamma, trigamma, polygamma

# for func in (:loggamma, :digamma, :trigamma)
#     @eval function SpecialFunctions.$func(z::Double64)
#         p = precision(z)
#         az = Arb(big(z), prec=p)
#         res = SpecialFunctions.$func(az)
#         return Double64(big(res))
#     end
# end

for func in (:loggamma, :gamma, :digamma, :trigamma)
    @eval function SpecialFunctions.$func(z::Complex{BigFloat})
        az = Acb(z, prec=precision(real(z)))
        res = SpecialFunctions.$func(az)
        return Complex{BigFloat}(real(res), imag(res))
    end
    # @eval function SpecialFunctions.$func(z::ComplexDF64)
    #     p = precision(z)
    #     az = Acb(complex(big(real(z)), imag(z)), prec=p)
    #     res = SpecialFunctions.$func(az)
    #     return ComplexDF64(big(real(res)), big(imag(res)))
    # end
end

function SpecialFunctions.polygamma(n::Integer, z::Complex{BigFloat})
    az = Acb(z, prec=precision(real(z)))
    an = Acb(n, prec=precision(real(z)))
    res = SpecialFunctions.polygamma(an, az)
    return Complex{BigFloat}(real(res), imag(res))
end

# function SpecialFunctions.polygamma(n::Integer, z::ComplexDF64)
#     az = Acb(big(z), prec=precision(z))
#     an = Acb(n, prec=precision(z))
#     res = SpecialFunctions.polygamma(an, az)
#     return ComplexDF64(real(res), imag(res))
# end

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
