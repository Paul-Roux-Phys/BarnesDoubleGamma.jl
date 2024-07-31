module BarnesDoubleGamma

import SpecialFunctions as SF
using  Memoization # cache intermediary results
using  ArbNumerics # the SpecialFunctions package has no arbitrary-precision complex-variable
                   # gamma function, however the ArbNumerics does. We use this,
                   # and convert to a Complex{BigFloat}
using  QuadGK      # numerical integration
import Base: precision

export loggamma,
       gamma,
       trigamma,
       polygamma,      # generic-typed gamma functions
       digamma_reg

export log_barnes_doublegamma,
       barnes_doublegamma,
       log_gamma2,
       gamma2,
       logdoublegamma,
       doublegamma

const ComplexOrReal{T} = Union{T, Complex{T}}
const ComplexOrRealFloat = Union{ComplexOrReal{Float64}, ComplexOrReal{BigFloat}}

include("gamma_functions.jl")
include("double_gamma_product_formula.jl")

end
