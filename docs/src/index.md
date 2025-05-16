# `BarnesDoubleGamma.jl` documentation

```@meta
CurrentModule = BarnesDoubleGamma
DocTestSetup = quote
    using BarnesDoubleGamma
end
```

This is the documentation page for the `BarnesDoubleGamma` package.
`BarnesDoubleGamma` is a julia package for computing Barnes' double gamma functions, see [this wikipedia article](https://en.wikipedia.org/wiki/Multiple_gamma_function) for a reference. The package works with real and complex numbers, in standard and multiple precision:

```jldoctest
julia> β = 0.74; w = 1.42; doublegamma(w, β) ≈ 1.330272154435545
true

julia> w = big(w); doublegamma(w, β) ≈ 1.330272154422092394938807374719727965657364472505115527465224747340970511072316
true

julia> w += big"0.1"*im; doublegamma(w, β) ≈ 1.320088253471064630736078956865882843925682644428785275208996772003089457017443 + 0.1236629786396930558290716646529048849169199735177974299246859233300114349665599im
true
```

As a by-product, the package also exports gamma and polygamma functions that work with real and complex numbers in arbitrary precision, see [Gamma and polygamma functions](@ref).

## Contents

* [Installation](installation.md)
* [Gamma and polygamma functions](gamma_functions.md)
* [Barnes double gamma functions](double_gamma.md)
