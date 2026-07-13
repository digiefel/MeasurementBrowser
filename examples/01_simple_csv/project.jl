using CSV
using DataBrowser
using DataFrames: DataFrame

length(ARGS) == 1 || error("Usage: julia --project project.jl DATA_DIRECTORY")

project = define_project("CSV tables")

register_item!(project;
    detect = (file::SourceFile) -> endswith(lowercase(file.filename), ".csv"),
    read = (file::SourceFile) -> CSV.read(file.filepath, DataFrame),
)

workspace = open_workspace(project, only(ARGS))
open_browser(workspace)
