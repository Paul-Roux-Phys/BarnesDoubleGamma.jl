module BarnesDoubleGamma

using QuadGK, # numerical integration
    SpecialFunctions,
    Arblib
    # DoubleFloats

# we extend these methods to Complex{BigFloat} using Arblib
import SpecialFunctions: loggamma, gamma, digamma, trigamma, polygamma

# generic-typed gamma functions
export loggamma,
    gamma,
    trigamma,
    polygamma,      
    digamma_reg

# Double Gamma functions
export LogBDoubleGamma, log_barnesdoublegamma,
       BDoubleGamma, barnesdoublegamma,
       LogGamma2, loggamma2,
       Gamma2, gamma2,
       LogDoubleGamma, logdoublegamma,
       DoubleGamma, doublegamma

include("gamma_functions.jl")
include("double_gamma_product_formula.jl")

end
