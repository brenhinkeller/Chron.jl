using Chron
using Documenter

DocMeta.setdocmeta!(Chron, :DocTestSetup, :(using Chron); recursive=true)

makedocs(;
    modules=[Chron],
    authors="C. Brenhin Keller",
    repo="https://github.com/brenhinkeller/Chron.jl/blob/{commit}{path}#{line}",
    sitename="Chron.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://brenhinkeller.github.io/Chron.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/brenhinkeller/Chron.jl",
    devbranch = "main",
)
