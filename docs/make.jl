using MultivariateInterpolation
using Documenter

DocMeta.setdocmeta!(
    MultivariateInterpolation,
    :DocTestSetup,
    :(using MultivariateInterpolation);
    recursive=true,
)

makedocs(;
    modules=[MultivariateInterpolation],
    authors="Uwe Hernandez Acosta",
    repo="https://github.com/MultivariantToolkit/MultivariateInterpolation.jl/blob/{commit}{path}#{line}",
    sitename="MultivariateInterpolation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MultivariantToolkit.github.io/MultivariateInterpolation.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)

deploydocs(;
    repo="github.com/MultivariantToolkit/MultivariateInterpolation.jl", devbranch="main"
)
