module BarnesDoubleGamma

using  Memoization # cache intermediary results
using  QuadGK      # numerical integration
import Base: precision

export loggamma,
       gamma,
       trigamma,
       polygamma,      # generic-typed gamma functions
       digamma_reg

export log_barnesdoublegamma,
       barnesdoublegamma,
       loggamma2,
       gamma2,
       logdoublegamma,
       doublegamma

const ComplexOrReal{T} = Union{T, Complex{T}}
const ComplexOrRealFloat = Union{ComplexOrReal{Float64}, ComplexOrReal{BigFloat}}

include("gamma_functions.jl")
include("double_gamma_product_formula.jl")

end
