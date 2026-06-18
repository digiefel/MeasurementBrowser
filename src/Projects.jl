module Projects

# ---------------------------------------------------------------------------
# Project type
#
# A project is a value built by registration (define_project + register_*). It is the single project
# type in the package. The struct and its recipe types live here, before ItemIndex/Workspace,
# so those modules can name the concrete type directly. The methods that operate on it (registration,
# interpretation, plotting, serialization) are defined later, once ItemIndex types exist.
# ---------------------------------------------------------------------------

"""
One registered item recipe for a kind.

`entries(file, data)` returns the per-item entries (`Vector{<:AbstractDataItem}`): the package's
`DataItem` for the recipe API, or a project's own subtype for the type API. The engine derives the
internal record from each via the contract.
"""
mutable struct ItemRecipe
    kind::Symbol
    detect::Function
    read::Function
    entries::Function
    process::Union{Nothing,Function}
    stats::Union{Nothing,Function}
    label::Union{Nothing,Function}
end

"""One registered cross-item, collection-scoped stat."""
struct CollectionStatRecipe
    kinds::Vector{Symbol}
    group_by::Function
    compute_stats::Function
end

"""
One registered plot recipe for an item kind.

`setup(workspace, items)` returns the `Figure`; `draw(workspace, items, figure)` fills it. `items` are
the loaded, data-bearing items for the selection (`Vector{<:AbstractDataItem}`); each carries its
processed payload as `item.data`. Neither callback resolves data itself.
"""
struct PlotRecipe
    kind::Symbol
    label::String
    setup::Function
    draw::Function
end

"""Accumulated read/stats timing for one item kind in the current scan."""
mutable struct KindProfile
    source_items::Int
    items::Int
    read_seconds::Float64
    stats_seconds::Float64
end
KindProfile() = KindProfile(0, 0, 0.0, 0.0)

"""
A callback project assembled from registered recipes.

Source interpretation, data processing, and presentation are defined by the registered callbacks.
Package-owned cache, job, and browser state does not belong here.
"""
mutable struct Project
    name::String
    description::String
    recipes::Vector{ItemRecipe}
    collection_stats::Dict{Tuple{Vararg{Symbol}},CollectionStatRecipe}
    plots::Dict{Symbol,Dict{String,PlotRecipe}}
    # Transient per-item analysis failures gathered while interpreting source items, as
    # (source item id, item id, message) tuples. Plain string tuples
    # so this early module needs no ItemIndex types.
    stat_failures::Vector{Tuple{String,String,String}}
    stat_failures_lock::ReentrantLock
    # Transient per-kind timing for the latest scan, surfaced in the performance window. Reset at the
    # start of each scan and replaced wholesale, so it stays bounded to one row per item kind.
    scan_profile::Dict{Symbol,KindProfile}
    profile_lock::ReentrantLock
end

const PROJECTS = Project[]
const DEFAULT_PROJECT = Ref{Union{Project,Nothing}}(nothing)

# ---------------------------------------------------------------------------
# The low-level source contract
#
# These are intentionally internal for now. The engine is source-first, but the exported API remains
# the smaller callback surface until the low-level module settles.
# ---------------------------------------------------------------------------

"""A data origin with lifecycle and discovery of source items."""
abstract type AbstractDataSource end

"""Private source adapter implementing the high-level callback API."""
mutable struct RegisteredProjectSource <: AbstractDataSource
    project::Project
    root_path::String
    collection_metadata::Any
end

RegisteredProjectSource(project::Project, root_path::AbstractString)::RegisteredProjectSource =
    RegisteredProjectSource(project, normpath(abspath(expanduser(String(root_path)))), nothing)

"""One addressable unit discovered inside a data source."""
abstract type AbstractDataSourceItem end

"""Stable source identity used for workspace/cache ownership."""
function source_id end

"""Human-readable source name."""
function source_label end

"""Prepare a source for use. Simple immutable sources return themselves."""
open_source(source::AbstractDataSource)::AbstractDataSource = source

"""Release resources owned by a source."""
close_source!(::AbstractDataSource)::Nothing = nothing

"""Return the current source items discovered by a source."""
function source_items end

"""Future live-update hook. `nothing` means the source is static."""
watch_source(::AbstractDataSource, ::Function) = nothing

"""Optional source-wide invalidation token."""
source_fingerprint(::AbstractDataSource) = nothing

"""Stable source-item identity within a source."""
function source_item_id end

"""Human-readable source-item label."""
function source_item_label end

"""Optional source-item invalidation token."""
source_item_fingerprint(::AbstractDataSourceItem) = nothing

"""Filesystem path for a source item, when one exists."""
source_item_path(::AbstractDataSourceItem)::Union{Nothing,String} = nothing

"""Timestamp for a source item, when one exists."""
source_item_timestamp(::AbstractDataSourceItem) = nothing

"""Interpret one source item into lightweight logical data items."""
function data_items end

"""Reload one logical data item with its payload attached."""
function load_data_item end

# ---------------------------------------------------------------------------
# The AbstractDataItem contract
#
# `AbstractDataItem` is the logical browser object. The index/cache store `ItemRecord`s derived from
# this contract; views receive loaded `AbstractDataItem` values.
# ---------------------------------------------------------------------------

abstract type AbstractDataItem end

"""Nestable container an item is placed in. Future home for collection metadata/behaviour."""
abstract type Collection end

"""Stable string identity of an item."""
function item_id end

"""Human-readable label for an item."""
function item_label end

"""The kind symbol the plot registry is keyed on."""
function kind end

"""Canonical tree placement of an item, as nested collection names (`Vector{String}`)."""
function collection end

"""Parsed acquisition parameters of an item (`Dict{Symbol,Any}`)."""
function parameters end

"""Computed statistics of an item (`Dict{Symbol,Any}`)."""
function stats end

"""The materialized payload carried by an item (also reachable as `item.data`)."""
function item_data end

"""Process raw item data into the payload a view consumes. Optional; default passthrough."""
process(::AbstractDataItem, data) = data

"""Whether an item's payload should be persisted by the data cache. Optional; default `false`."""
cacheable(::AbstractDataItem)::Bool = false

"""Optional item-level invalidation token."""
item_fingerprint(::AbstractDataItem) = nothing

"""Optional fold over one collection node's data-less item handles."""
collection_stats(
    ::AbstractDataSource,
    ::Vector{String},
    ::Vector{<:AbstractDataItem},
)::Dict{Symbol,Any} = Dict{Symbol,Any}()

# ---------------------------------------------------------------------------
# Interface functions
#
# Declared here so the package's submodules can call them; the methods are defined later, after the
# types they depend on (ItemRecord, SourceFile, Workspace) are available.
# ---------------------------------------------------------------------------

"""Return the stable name used to identify a project."""
function project_name end

"""Return a short human-readable description of a project."""
function project_description end

"""Classify the item kind represented by a project source filename."""
function detect_kind end

"""Return the human-readable label for a project item kind."""
function kind_label end

"""Return the human-readable label for one logical item."""
function display_label end

"""Return the project-specific display label for one device path."""
function collection_path_label end

"""Clear any per-scan timing a project accumulates. Called once at the start of every scan."""
function reset_scan_profile! end

"""
Per-item-kind timing for the most recent scan, newest scan replacing the last.

Returns one `(; kind, source_items, items, read_seconds, stats_seconds)` row per kind. Surfaced in
the performance window.
"""
function scan_profile_summary end

end
