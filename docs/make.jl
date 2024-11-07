push!(LOAD_PATH, joinpath("..", "src"))

using BarnesDoubleGamma
using Documenter

DocMeta.setdocmeta!(
    BarnesDoubleGamma,
    :DocTestSetup,
    :(using BarnesDoubleGamma),
    recursive = true
)

makedocs(
    sitename = "BarnesDoubleGamma.jl",
    repo = "https://github.com/Paul-Roux-Phys/BarnesDoubleGamma.jl",
    modules = [BarnesDoubleGamma],
    format = Documenter.HTML(
        repolink = "https://github.com/Paul-Roux-Phys/BarnesDoubleGamma.jl",
        edit_link = :commit
    ),
    doctest = true,
    pages = [
        "Home" => "index.md",
        "installation.md",
        "gamma_functions.md",
        "double_gamma.md"
    ]
)

deploydocs(
    repo   = "https://github.com/Paul-Roux-Phys/BarnesDoubleGamma.jl",
    push_preview = true,
    modules = [BarnesDoubleGamma],
    format = Documenter.HTML(
        repolink = "https://github.com/Paul-Roux-Phys/BarnesDoubleGamma.jl",
        edit_link = :commit
    ),
    doctest = true,
    pages = [
        "Home" => "index.md",
        "installation.md",
        "gamma_functions.md",
        "double_gamma.md"
    ]
)
