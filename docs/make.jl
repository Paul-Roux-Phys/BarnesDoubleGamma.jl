using BarnesDoubleGamma
using Documenter

DocMeta.setdocmeta!(BarnesDoubleGamma, :DocTestSetup, :(using BarnesDoubleGamma); recursive=true)

makedocs(;
    modules=[BarnesDoubleGamma],
    authors="Paul Roux",
    repo="https://github.com/Paul-Roux-Phys/BarnesDoubleGamma.jl/blob/{commit}{path}#{line}",
    sitename="BarnesDoubleGamma.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Paul-Roux-Phys.github.io/BarnesDoubleGamma.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Paul-Roux-Phys/BarnesDoubleGamma.jl",
    devbranch="main",
)
