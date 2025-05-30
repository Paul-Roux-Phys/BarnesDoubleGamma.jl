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
export LogBDoubleGamma,
       BDoubleGamma,
       LogGamma2,
       Gamma2,
       LogDoubleGamma,
       DoubleGamma

include("gamma_functions.jl")
include("double_gamma_product_formula.jl")

end
