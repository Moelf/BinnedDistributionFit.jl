using BinnedDistributionFit
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(BinnedDistributionFit, :DocTestSetup, :(using BinnedDistributionFit); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [BinnedDistributionFit],
    repo = Remotes.GitHub("Moelf", "BinnedDistributionFit.jl"),
    authors = "Moelf <proton@jling.dev> and contributors",
    sitename = "BinnedDistributionFit.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/Moelf/BinnedDistributionFit.jl/",
    ),
    pages = ["index.md"; numbered_pages],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/Moelf/BinnedDistributionFit.jl",
    target = "build", # this is where Vitepress stores its output
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)
