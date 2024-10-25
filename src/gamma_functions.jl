function precision(z::Complex)
    return precision(real(z))
end

"""Convert to a standard precision number"""
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

"""Change the precision of x to precision"""
macro convert_precision!(var, precision)
    quote
        $var = convert_precision($var, $precision)
    end |> esc
end

#=
 We create bindings to gamma, digamma, and polygamma functions in arbitrary precision
=#
for f in (:gamma, :digamma)
    @eval $f(z::Union{Real, Complex{Float64}}) = SF.$f(z)
    @eval $f(z::Complex{BigFloat}) = Complex{BigFloat}(ArbNumerics.$f(ArbComplex(z, bits=precision(BigFloat))))
end

trigamma(z::ComplexOrReal{Float64}) = SF.trigamma(z)
trigamma(z::BigFloat) = BigFloat(ArbNumerics.polygamma(ArbComplex(1, bits=precision(BigFloat)), ArbComplex(z, bits=precision(BigFloat))))
trigamma(z::Complex{BigFloat}) = Complex{BigFloat}(ArbNumerics.polygamma(ArbComplex(1, bits=precision(BigFloat)), ArbComplex(z, bits=precision(BigFloat))))

loggamma(z::Union{ComplexOrReal{Float64}, BigFloat}) = SF.loggamma(z)
loggamma(z::Complex{BigFloat}) = Complex{BigFloat}(ArbNumerics.lgamma(ArbComplex(z, bits=precision(BigFloat))))

polygamma(n, z::Union{Real, Complex{Float64}}) = SF.polygamma(n, z)
polygamma(n, z::Complex{BigFloat}) = Complex{BigFloat}(ArbNumerics.polygamma(ArbComplex(n, bits=precision(BigFloat)), ArbComplex(z, bits=precision(BigFloat))))


"""
    cotpi(x) = cot(π * x)
"""
cotpi(x) = cospi(x) / sinpi(x)

"""Regularised digamma function"""
function digamma_reg(z)
    if real(z) > 0
        return digamma(z)
    elseif isreal(z) && real(z) < 0 && real(z)%1 == 0
        return digamma(1-z)
    else
        return digamma(1-z) - oftype(z, π)*cotpi(z)
    end
end
