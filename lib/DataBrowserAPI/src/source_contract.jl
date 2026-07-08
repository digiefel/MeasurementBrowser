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
