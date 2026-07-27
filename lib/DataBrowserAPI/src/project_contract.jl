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

"""
Classify the item kind a project would produce from a source filename, without reading it.

A cheap description query for status surfaces, not a pipeline stage: routing happens inside `read`.
"""
detect_kind(::AbstractProject, ::String)::Symbol = :unknown

"""
Whether a project runs a collection `process` stage for one item kind.

A scheduling declaration, not a capability: the stage dispatches for any project, but the engine
only materializes a collection's members when its project says the work exists. The registration
dialect answers from its collection recipes; a typed project with collection stages overrides this.
"""
_has_collection_process(::AbstractProject, ::Symbol)::Bool = false

"""Whether a project runs a collection `analyze` stage for one item kind. See above."""
_has_collection_analysis(::AbstractProject, ::Symbol)::Bool = false

"""Return the human-readable label for one logical item record."""
function display_label end

"""Return the project-specific display label for one collection path."""
function collection_path_label end

"""
Resolve one stored item kind back to the concrete type that produced it, or `nothing`.

The cache stores `kind` as a `Symbol`, but rebuilding a cached item through `reconstruct` needs the
type. The project answers, because the project is what knows its own item types — no module
scanning, and no module-qualified type name to keep valid across refactors. Returning `nothing`
(the default) is always safe: the engine falls back to rerunning `read` → `entries` → `process`.
"""
item_type(::AbstractProject, ::Symbol)::Union{Nothing,Type} = nothing

"""
Projects a GUI session can offer when the caller did not supply one, and the preferred default.

Registries of `AbstractProject`, not of any one dialect: the GUI picks between projects without
knowing how they were defined.
"""
const PROJECTS = AbstractProject[]
const DEFAULT_PROJECT = Ref{Union{AbstractProject,Nothing}}(nothing)
