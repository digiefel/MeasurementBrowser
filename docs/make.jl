pushfirst!(LOAD_PATH, joinpath(@__DIR__, ".."))

using DataBrowser
using Documenter
using DocumenterMermaid

makedocs(
    modules=[DataBrowser],
    sitename="DataBrowser",
    source="public",
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
        collapselevel=1,
    ),
    pages=[
        "Home" => "index.md",
        "How DataBrowser works" => "pipeline.md",
        "Your first project" => "getting-started.md",
        "Project guide" => [
            "Registration API" => "registration.md",
            "Metadata and collections" => "metadata-and-collections.md",
            "Workspaces" => "workspaces.md",
        ],
        "Type API" => "type-api.md",
        "Examples" => "examples.md",
        "API reference" => "reference.md",
    ],
    warnonly=[:missing_docs],
)
