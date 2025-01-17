using PaulStrings2
using Documenter

DocMeta.setdocmeta!(PaulStrings2, :DocTestSetup, :(using PaulStrings2); recursive=true)

makedocs(;
    modules=[PaulStrings2],
    authors="John Lapeyre",
    sitename="PaulStrings2.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
