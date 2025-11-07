struct DSine{T}
    Nmax::Int
    b::T
    Q::T
    minus_ipiover2::T
    sixthQsqp1::T
    twopii_binv::T
    twopii_b::T
    twopii_m_bsq_inv::Vector{T}
    twopii_m_b2::Vector{T}
end

function DSine(b::T) where {T<:Complex}
    @assert imag(b^2) > 0 "Only valid for Im b^2 > 0"
    precision = -ceil(Int64, log(eps(real(typeof(b)))))
    Nmax = precision * ceil(Int64, max(1 / imag(b^2), - 1 / imag(b^(-2))))
    Q = b + 1 / b
    π = oftype(b, Base.π)
    minus_ipiover2 = -im * π / 2
    sixthQsqp1 = (Q^2 + 1) / 6
    twopii_binv = 2 * π * im / b
    twopii_b = 2 * π * im * b
    twopii_m_b2_inv = [2 * π * im * m / b^2 for m in 1:Nmax]
    twopii_m_b2 = [2 * π * im * m * b^2 for m in 1:Nmax]
    return DSine{T}(
        Nmax, b, Q, minus_ipiover2, sixthQsqp1, twopii_binv, twopii_b, twopii_m_b2_inv, twopii_m_b2
    )
end

function (S::DSine)(z)
    res = exp(S.minus_ipiover2 * (z*z - S.Q * z + S.sixthQsqp1))
    res /= 1 - exp(S.twopii_b * z)
    for m = 1:S.Nmax
        res *= (1 - exp(S.twopii_binv * z - S.twopii_m_bsq_inv[m])) /
            (1 - exp(S.twopii_b * z + S.twopii_m_b2[m]))
    end
    return res
end

function Base.show(io::IO, S::DSine)
    print(io, "S_{b = $(S.b)}")
end
