using DataBrowser

length(ARGS) == 1 || error("Usage: julia --project project.jl DATA_DIRECTORY")

project = define_project("CSV tables")

# A premade recipe supplies `detect` and `read`; every other `register_item!` callback still
# applies, so this remains a starting point rather than a special case.
register_csv!(project)

workspace = open_workspace(project, only(ARGS))
open_browser(workspace)
