using Documenter, BlockMaps

makedocs(;
    modules=[BlockMaps],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/jagot/BlockMaps.jl/blob/{commit}{path}#L{line}",
    sitename="BlockMaps.jl",
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/jagot/BlockMaps.jl",
)
