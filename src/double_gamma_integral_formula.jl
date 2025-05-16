using QuadGK

# assume for now w is in the right range for the integral formula.
integrand_I(t, w, β) =
    w / t - w^2 * exp(-2t) - exp(-w * t) * sinh(w * t) / sinh(β * t) / sinh(t / β)

ϵ_terms(w, β, a) = 2*w2 -2*w3/3 + (β^2+1/β^2)/6*w

function _log_doublegamma(w, β)
    if real(2w - β - inv(β)) <= 0
        return log(2oftype(w, π))/2
    end
end
