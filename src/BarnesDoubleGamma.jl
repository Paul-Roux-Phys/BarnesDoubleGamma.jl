module BarnesDoubleGamma

using  Memoization, # cache intermediary results
    QuadGK, # numerical integration
    DoubleFloats # 30-digit floats, fast.

import Base: precision

export Double64, ComplexDF64,
    loggamma,
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

include("gamma_functions.jl")
include("double_gamma_product_formula.jl")

end
