import SpecialFunctions: loggamma, gamma, digamma, trigamma, polygamma

for func in (:loggamma, :gamma, :digamma, :trigamma)
    @eval function SpecialFunctions.$func(z::Complex{BigFloat})
        az = Acb(z, prec=precision(real(z)))
        res = SpecialFunctions.$func(az)
        return Complex{BigFloat}(real(res), imag(res))
    end
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
