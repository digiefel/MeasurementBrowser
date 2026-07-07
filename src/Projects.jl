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
`DataItem` for the recipe API, or a project's own subtype for the type API. Later callbacks receive
items and read `item.data`.
"""
mutable struct ItemRecipe
    kind::Symbol
    detect::Function
    read::Function
    entries::Function
    process::Union{Nothing,Function}
    analyze::Union{Nothing,Function}
    label::Union{Nothing,Function}
end

"""
One registered collection recipe for a kind.

`process(items)` rewrites each member (one output per input, same ids) and may change per-item data
or metadata. `analyze(items)` folds the post-process members into a `Dict{Symbol,Any}` attached to
the collection node only.
"""
struct CollectionRecipe
    kind::Symbol
    process::Union{Nothing,Function}
    analyze::Union{Nothing,Function}
end

"""
One registered plot recipe for an item kind.

`setup(workspace, items)` returns the `Figure`; `draw(workspace, items, figure)` fills it. `items` are
the loaded, data-bearing items for the selection (`Vector{<:AbstractDataItem}`); `process(item)`, if
registered, already returned the item shape these callbacks receive. Neither callback resolves data
itself.
"""
struct PlotRecipe
    kind::Symbol
    label::String
    setup::Function
    draw::Function
end

"""Timing retained for one source item in the current scan."""
mutable struct SourceItemProfile
    source_item_id::String
    source_item_label::String
    source_item_path::Union{Nothing,String}
    kind::Symbol
    item_count::Int
    detect_seconds::Float64
    read_seconds::Float64
    entries_seconds::Float64
    process_seconds::Float64
    analyze_seconds::Float64
    total_seconds::Float64
    thread_ids::Set{Int}
end

SourceItemProfile(source_item_id::AbstractString)::SourceItemProfile = SourceItemProfile(
    String(source_item_id), String(source_item_id), nothing,
    :unmatched, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Set{Int}())

"""One item-kind summary derived from the current source-item timings."""
mutable struct KindProfileRow
    kind::Symbol
    source_items::Int
    items::Int
    detect_seconds::Float64
    read_seconds::Float64
    entries_seconds::Float64
    process_seconds::Float64
    analyze_seconds::Float64
    total_seconds::Float64
end

KindProfileRow(kind::Symbol)::KindProfileRow =
    KindProfileRow(kind, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

"""Read-only source-item row returned to the Performance window."""
struct SourceProfileRow
    source_item_id::String
    source_item_label::String
    source_item_path::Union{Nothing,String}
    kind::Symbol
    items::Int
    detect_seconds::Float64
    read_seconds::Float64
    entries_seconds::Float64
    process_seconds::Float64
    analyze_seconds::Float64
    total_seconds::Float64
    thread_ids::Vector{Int}
end

"""
A callback project assembled from registered recipes.

Source interpretation, data processing, and presentation are defined by the registered callbacks.
Package-owned cache, job, and browser state does not belong here.
"""
mutable struct Project
    name::String
    description::String
    recipes::Vector{ItemRecipe}
    collections::Dict{Symbol,CollectionRecipe}
    plots::Dict{Symbol,Dict{String,PlotRecipe}}
    # Transient timing for the latest scan, surfaced in the performance window. One row per source
    # item keeps the profile useful without retaining every expanded data-item event.
    scan_profile::Dict{String,SourceItemProfile}
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

"""One addressable unit discovered inside a data source."""
abstract type AbstractDataSourceItem end

"""
One source-owned change batch.

Sources report physical source-item replacements/removals and source-provided metadata changes
through the same update contract. An empty item batch with `metadata_changed=true` updates
existing logical items without reinterpreting their source items.
"""
struct SourceChanges{S<:AbstractDataSourceItem}
    upserts::Vector{S}
    removals::Vector{String}
    metadata_changed::Bool
end

SourceChanges(
    upserts::Vector{S},
    removals::Vector{String};
    metadata_changed::Bool=false,
) where {S<:AbstractDataSourceItem} = SourceChanges(upserts, removals, metadata_changed)

"""
One recoverable source failure reported through the watch contract.

Failures are expected states, not watcher deaths: a source that cannot currently produce a
consistent update (e.g. a malformed live metadata file) reports the reason as a value and keeps
watching; a later successful update clears it.
"""
struct SourceError
    message::String
end

"""Stable source identity used for workspace/cache ownership."""
function source_id end

"""Human-readable source name."""
function source_label end

"""Prepare a source for use. Simple immutable sources return themselves."""
open_source(source::AbstractDataSource)::AbstractDataSource = source

"""Release resources owned by a source."""
close_source!(::AbstractDataSource)::Nothing = nothing

"""
Source construction options to replay when reopening an equivalent source, as keyword arguments
accepted by the source's `open_workspace` method. Sources without reopen options return `(;)`.
"""
source_open_options(::AbstractDataSource)::NamedTuple = (;)

"""Return the current source items discovered by a source."""
function source_items end

source_items(source::AbstractDataSource; on_progress::Union{Nothing,Function}=nothing) =
    source_items(source)

"""Human noun for source items, used by status surfaces."""
source_item_noun(::AbstractDataSource)::String = "source items"

"""
Watch a source and call `on_change` with each `SourceChanges` batch or recoverable `SourceError`.
`nothing` means the source is static.
"""
watch_source(::AbstractDataSource, ::Function) = nothing

"""Stable source-item identity within a source."""
function source_item_id end

"""Human-readable source-item label."""
function source_item_label end

"""Optional invalidation token for source items and data items."""
function fingerprint end

fingerprint(::AbstractDataSourceItem) = nothing

"""Filesystem path for a source item, when one exists."""
source_item_path(::AbstractDataSourceItem)::Union{Nothing,String} = nothing

"""Timestamp for a source item, when one exists."""
source_item_timestamp(::AbstractDataSourceItem) = nothing

"""Interpret one source item into lightweight logical data items."""
function data_items end

# ---------------------------------------------------------------------------
# The AbstractDataItem contract
#
# `AbstractDataItem` is the logical browser object. The index/cache store `ItemRecord`s derived from
# this contract; views receive loaded `AbstractDataItem` values.
# ---------------------------------------------------------------------------

abstract type AbstractDataItem end

"""Nestable container an item is placed in. Future home for collection behaviour."""
abstract type Collection end

"""Stable string identity of an item."""
function id end

"""Human-readable label for an item."""
function item_label end

"""The kind symbol the plot registry is keyed on."""
function kind end

"""Canonical tree placement of an item, as nested collection names (`Vector{String}`)."""
function collection end

"""Metadata of an item (`Dict{Symbol,Any}`): parsed parameters, computed values, provenance merged."""
function metadata end

"""The materialized data carried by an item (also reachable as `item.data`)."""
function item_data end

"""Process an item into the item a view consumes. Optional; default identity."""
process(item::AbstractDataItem) = item

"""Whether an item's data should be persisted by the data cache. Optional; default `false`."""
cacheable(::AbstractDataItem)::Bool = false

fingerprint(::AbstractDataItem) = nothing

"""Optional rewrite of a collection's members (one output per input). Internal workspace hook."""
function process_collection end

"""Optional fold over a collection's members into collection-node metadata. Internal workspace hook."""
function analyze_collection end

"""Collection metadata contributed by a source for one collection path."""
collection_metadata(::AbstractDataSource, ::AbstractVector{<:AbstractString})::Dict{Symbol,Any} =
    Dict{Symbol,Any}()

"""Whether a source supplied collection metadata for this scan."""
has_collection_metadata(::AbstractDataSource)::Bool = false

"""Per-item analysis metadata computed after indexing. Internal workspace hook."""
function analyze_item end

"""Whether one item kind has a registered collection `process` stage."""
function has_collection_process end

"""Whether one item kind has a registered collection `analyze` stage."""
function has_collection_analysis end

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

"""Return the project-specific display label for one collection path."""
function collection_path_label end

"""Clear any per-scan timing a project accumulates. Called once at the start of every scan."""
function reset_scan_profile! end

"""Record one callback phase for a source item in the current scan."""
function record_scan_phase! end

"""Finish one source-item timing after all expanded items have been processed."""
function finish_source_profile! end

"""
Per-item-kind timing for the most recent scan, newest scan replacing the last.

Returns one aggregate row per item kind. Surfaced in the performance window.
"""
function scan_profile_summary end

"""Return source-item timing rows for the latest scan, slowest first."""
function scan_source_profile end

end
