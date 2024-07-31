# BarnesDoubleGamma

This Julia package exports functions
`log_barnesgamma`, `barnesgamma`
that compute the values of the Barne's double gamma function as defined in [arXiv:2208.13879](https://arxiv.org/abs/2208.13876).
It also exports the related functions
`log_gamma2`, `gamma2`,
`log_doublegamma`, `double_gamma`,
corresponding to the $\Gamma_2$ and $\Gamma_b$ functions of [this](https://en.wikipedia.org/wiki/Multiple_gamma_function) wikipedia article.

All of these functions are available for real or complex argument, in standard or multiple precision (using BigFloat and Complex{BigFloat} arguments)

The package also exports
`loggamma`, `gamma`,
`digamma_reg`,
`trigamma`,
`polygamma`
functions that work for real or complex arguments in standard or multiple precision, by calling the appropriate functions from the `SpecialFunctions` package, or from the `ArbNumerics` package.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Paul-Roux-Phys.github.io/BarnesDoubleGamma.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Paul-Roux-Phys.github.io/BarnesDoubleGamma.jl/dev/)
[![Build Status](https://github.com/Paul-Roux-Phys/BarnesDoubleGamma.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Paul-Roux-Phys/BarnesDoubleGamma.jl/actions/workflows/CI.yml?query=branch%3Amain)
