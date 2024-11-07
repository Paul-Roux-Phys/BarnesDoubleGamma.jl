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
julia> β = 0.74; w = 1.42; v1 = doublegamma(w, β)
1.3302721668409918

julia> w = big(w); v2 = doublegamma(w, β)
1.330272166855746161361194084745097684734792530313816671600279858462974443246171

julia> w += big"0.1"*im; doublegamma(w, β)
1.320088265498793958374893745924799011170556730654186142861596806489898950191236 + 0.1236629830881181358927960770656024793760681821275006476859665026276396033441505im
```

As a by-product, the package also exports gamma and polygamma functions that work with real and complex numbers in arbitrary precision, see [Gamma and polygamma functions](@ref).

## Contents

* [Installation](installation.md)
* [Gamma and polygamma functions](gamma_functions.md)
* [Barnes double gamma functions](double_gamma.md)
