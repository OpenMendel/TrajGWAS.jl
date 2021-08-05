using TrajGWAS
using Documenter

makedocs(;
    modules=[TrajGWAS],
    authors="Seyoon Ko <kos@ucla.edu> and contributors",
    repo="https://github.com/OpenMendel/TrajGWAS.jl/blob/{commit}{path}#L{line}",
    sitename="TrajGWAS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://OpenMendel.github.io/TrajGWAS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/OpenMendel/TrajGWAS.jl",
    devbranch = "main"
)
