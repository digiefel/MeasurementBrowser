using CancellationTokens: CancellationToken

"""A data origin with lifecycle and discovery of source items."""
abstract type AbstractDataSource end

"""One addressable unit discovered inside a data source."""
abstract type AbstractDataSourceItem end

"""
    metadata(value) -> Dict

Return metadata supplied directly by `value`. It does not include metadata supplied by source
items or collections around that value.
"""
function metadata end

"""Return metadata supplied directly by a source item. The default is an empty `Dict`."""
metadata(::AbstractDataSourceItem)::Dict = Dict()

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
    source_items(source; cancel_token, on_progress, on_item) -> Vector{<:AbstractDataSourceItem}

Return the current source items discovered by a source.

The engine always passes the three keywords, so every method must accept them — a source that
discovers everything at once absorbs them with `; kwargs...` and ignores them. A source slow enough
to be worth streaming calls `on_item` per item as it finds them (the engine then queues each one
immediately rather than waiting for the return), reports counts through `on_progress`, and honors
`cancel_token`. The returned vector is queued only when `on_item` was never called, so a streaming
source may return its items as well without them being processed twice.
"""
function source_items end

"""Human noun for source items, used by status surfaces."""
source_item_noun(::AbstractDataSource)::String = "source items"

"""
Watch a source and call `on_change` with each `SourceChanges` batch or recoverable `SourceError`.
`nothing` means the source is static.
"""
watch_source(::AbstractDataSource, ::Function; cancel_token::CancellationToken) = nothing

"""
    id(value)::String

Stable identity supplied by `value`, used exactly as returned. Sources, source items, data items,
and collections all implement it; none has a default, and none has its answer wrapped or
namespaced. An id has to stay stable across scans and reopenings, and two values answering the same
id collide.
"""
function id end

"""
    label(value) -> String

Human-readable label supplied by `value`. Source items must implement this; data items and
collections have defaults (see the item contract).
"""
function label end

"""Optional invalidation token for source items and data items."""
function fingerprint end

fingerprint(::AbstractDataSourceItem) = nothing

"""Filesystem path for a source item, when one exists."""
source_item_path(::AbstractDataSourceItem)::Union{Nothing,String} = nothing

"""Timestamp for a source item, when one exists."""
source_item_timestamp(::AbstractDataSourceItem) = nothing

"""
Where a source places an item that declares no collection path of its own.

Applied by the engine after `entries`, where the source is legitimately in hand — the pipeline
stages themselves stay pure functions of values. The default leaves such items at the root; a
directory source places them under their directory relative to its root.
"""
default_collection_path(::AbstractDataSource, ::AbstractDataSourceItem) = AbstractCollection[]

"""
Let a source attach its own metadata to the levels of one item's collection path.

The directory source attaches its `metadata.txt` entries to each named level. The default returns
the path unchanged, so a source with nothing to add costs nothing.
"""
annotate_collection_path(::AbstractDataSource, path::AbstractVector) = path
