"""
Supertype of every DataBrowser project.

A project answers the declarations in this file and implements the pipeline stages in
`stage_contract.jl`. It carries no package-owned state: caches, jobs, indexes, and browser state
belong to the workspace. The registration dialect's recipe-holding project is one implementation
of this contract, not the contract itself.
"""
abstract type AbstractProject end

"""Return the stable name used to identify a project. Defaults to the project type's name."""
project_name(project::AbstractProject)::String = string(nameof(typeof(project)))

"""Return a short human-readable description of a project."""
project_description(::AbstractProject)::String = ""

"""Return the human-readable label for a project item kind."""
kind_label(::AbstractProject, kind::Symbol)::String = string(kind)

"""Return the human-readable label for one logical item record."""
function display_label end

"""Return the project-specific display label for one collection path."""
function collection_path_label end

"""
Projects a GUI session can offer when the caller did not supply one, and the preferred default.

Registries of `AbstractProject`, not of any one dialect: the GUI picks between projects without
knowing how they were defined.
"""
const PROJECTS = AbstractProject[]
const DEFAULT_PROJECT = Ref{Union{AbstractProject,Nothing}}(nothing)
