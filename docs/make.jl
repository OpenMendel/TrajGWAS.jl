using vGWAS
using Documenter

makedocs(;
    modules=[vGWAS],
    authors="Seyoon Ko <syko0507@snu.ac.kr> and contributors",
    repo="https://github.com/kose-y/vGWAS.jl/blob/{commit}{path}#L{line}",
    sitename="vGWAS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kose-y.github.io/vGWAS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kose-y/vGWAS.jl",
)
