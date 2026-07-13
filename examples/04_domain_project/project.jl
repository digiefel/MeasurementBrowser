using DataBrowser

include("FerroelectricProject.jl")
using .FerroelectricProject: define_ferroelectric_project

length(ARGS) == 1 || error("Usage: julia --project project.jl DATA_DIRECTORY")

project = define_ferroelectric_project()
workspace = open_workspace(project, only(ARGS))
open_browser(workspace)
