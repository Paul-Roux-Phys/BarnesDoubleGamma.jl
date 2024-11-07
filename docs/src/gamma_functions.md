# Gamma and polygamma functions

Since the [SpecialFunctions](https://specialfunctions.juliamath.org/latest/) module does not expose methods for its `gamma`, `digamma`, `trigamma`, and `polygamma` functions for `Complex{BigFloat}` arguments (or `BigFloat` arguments in some cases), the `BarnesDoubleGamma` package imports these methods from the `ArbNumerics` package instead, and exposes these functions. This means that users of this package can compute gamma functions in arbitrary precision:

```jldoctest
julia> using BarnesDoubleGamma

julia> gamma(0.5)
1.772453850905516

julia> gamma(big"0.5")
1.772453850905516027298167483341145182797549456122387128213807789852911284591025

julia> trigamma(big"0.5" + 0.1im)
4.4780986332966090515538950704365981683158900719789372422121210172045798838399969294524 - 1.561594164192146024275706100493907132241385845177793426420827577586671268403448im

julia> gamma(0.5+0.1im)
1.697617826382886 - 0.33284283907262135im
```

The package also exports a regularised digamma function `digamma_reg`

```@docs
digamma_reg(z)
```
