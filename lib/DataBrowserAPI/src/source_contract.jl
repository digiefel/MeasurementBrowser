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
Source construction options to replay when reopening an equivalent source, as keyword arguments
accepted by the source's `open_workspace` method. Sources without reopen options return `(;)`.
"""
source_open_options(::AbstractDataSource)::NamedTuple = (;)

"""Return the current source items discovered by a source."""
function source_items end

"""Human noun for source items, used by status surfaces."""
source_item_noun(::AbstractDataSource)::String = "source items"

"""
Watch a source and call `on_change` with each `SourceChanges` batch or recoverable `SourceError`.
`nothing` means the source is static.
"""
watch_source(::AbstractDataSource, ::Function; cancel_token::CancellationToken) = nothing

"""
    id(value)

Stable identity supplied by `value`. Source items must implement this: the id has to be stable
within their source across scans and reopenings. Data items and collections have defaults (see
the item contract).
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

"""Interpret one source item into lightweight logical data items."""
function data_items end
