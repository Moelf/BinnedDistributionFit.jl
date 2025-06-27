using BinnedDistributionFit
using Documenter

DocMeta.setdocmeta!(BinnedDistributionFit, :DocTestSetup, :(using BinnedDistributionFit); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [BinnedDistributionFit],
    authors = "Moelf <proton@jling.dev> and contributors",
    repo = "https://github.com/Moelf/BinnedDistributionFit.jl/blob/{commit}{path}#{line}",
    sitename = "BinnedDistributionFit.jl",
    format = Documenter.HTML(; canonical = "https://Moelf.github.io/BinnedDistributionFit.jl"),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/Moelf/BinnedDistributionFit.jl")
