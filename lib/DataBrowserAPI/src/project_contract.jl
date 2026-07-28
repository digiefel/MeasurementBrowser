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
Resolve one stored item kind back to the concrete type that produced it, or `nothing`.

The cache stores `kind` as a `Symbol`, but rebuilding a cached item through `reconstruct` needs the
type. Returning `nothing` (the default) is always safe: the engine looks for a loaded leaf subtype
of `AbstractDataItem` whose name matches the kind, then falls back to rerunning
`read` → `entries` → `process`. Override when kinds are not type names.
"""
item_type(::AbstractProject, ::Symbol)::Union{Nothing,Type} = nothing

"""
Resolve one stored collection kind back to the concrete type that produced it, or `nothing`.

The collection counterpart of `item_type`. Collections are always rebuilt from their stored row —
they have no rerun-from-source fallback — so a kind that resolves to nothing is an error at the
point of use, not a slow path. Returning `nothing` (the default) is safe: the engine looks for a
loaded leaf subtype of `AbstractCollection` whose name matches.
"""
collection_type(::AbstractProject, ::Symbol)::Union{Nothing,Type} = nothing

"""
Projects a GUI session can offer when the caller did not supply one, and the preferred default.

Registries of `AbstractProject`, not of any one dialect: the GUI picks between projects without
knowing how they were defined.
"""
const PROJECTS = AbstractProject[]
const DEFAULT_PROJECT = Ref{Union{AbstractProject,Nothing}}(nothing)
