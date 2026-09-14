"""
The typed pipeline stage contract.

Source discovery belongs to the source (`source_items`). Everything after it is these stages:

```text
read → entries → process → analyze → collection process → collection analyze
```

Each stage has a project-aware form, which the engine calls, and a context-free form, which
projects usually implement. The project-aware defaults forward to the context-free ones, so a
project that owns its own source and item types never mentions the project, while a project reusing
a shared source (`DirectorySource`) dispatches on its own project type instead.

`read` obtains the source content; `entries` turns that content into items. The workspace runs
these as separate stages. Read results are cached in memory by default and may be evicted. A
scheduled interpretation retains its input until it finishes. After eviction or reopening, reading
runs again only if a consumer needs the missing input; valid downstream results remain usable.
Source changes invalidate both stages. A failure belongs to the stage that raised it and remains
cached until that stage is invalidated.

`entries` still receives the source item because identity is not payload: the loaded value stays
purely the expensive data, while ids, labels, and collection placement derive from the source item —
which the engine can supply without touching the origin.

`read`'s return value and `entries`' `loaded` argument deliberately share no supertype. The loaded
value is a private handoff between two stages of the same project; the type discipline is enforced
by the concrete signature on the receiving end.
"""

"""
    PipelineStage

The stage identified by a scheduled job, cached payload, or cached completion record.
Values follow pipeline order. `SOURCE_READ` runs `read`; `SOURCE_INTERPRET` runs `entries`
and assigns item identity and collection placement. Analysis stages produce metadata.
"""
@enum PipelineStage::Int8 begin
    SOURCE_READ = 0
    SOURCE_INTERPRET = 1
    ITEM_PROCESS = 2
    ITEM_ANALYZE = 3
    COLLECTION_PROCESS = 4
    COLLECTION_ANALYZE = 5
end

# ---------------------------------------------------------------------------
# read
# ---------------------------------------------------------------------------

"""
    read(project, source, item) -> loaded
    read(source, item) -> loaded

Perform the one expensive source operation for a source item and return the project's loaded value.

This is the only stage that sees the source. Implement the context-free form when the project owns
its source item type, or the project-aware form when it reuses a shared source.
"""
read(project::AbstractProject, source::AbstractDataSource, item::AbstractDataSourceItem) =
    read(source, item)

read(source::AbstractDataSource, item::AbstractDataSourceItem) = error(
    "No read stage for $(typeof(item)) from $(typeof(source)). Implement " *
    "`DataBrowser.read(::$(typeof(source)), ::$(typeof(item)))`, or the project-aware " *
    "`DataBrowser.read(::YourProject, source, item)` when reusing a shared source.",
)

# ---------------------------------------------------------------------------
# entries
# ---------------------------------------------------------------------------

"""
    entries(project, item, loaded) -> Vector{<:AbstractDataItem}
    entries(item, loaded) -> Vector{<:AbstractDataItem}

Expand one loaded value into zero, one, or many concrete data items without rereading the source.

The default treats the loaded value as a single item. The workspace may call `entries` again with
the same cached `loaded` value. Preserve that value for reuse; copy data that needs modification.
"""
entries(project::AbstractProject, item::AbstractDataSourceItem, loaded) = entries(item, loaded)

entries(::AbstractDataSourceItem, loaded) = [loaded]

# ---------------------------------------------------------------------------
# process and analyze
# ---------------------------------------------------------------------------

"""
    process(project, item) -> AbstractDataItem

Process one interpreted item. The default forwards to the context-free `process(item)`.

Metadata returned by `metadata(processed_item)` is merged over the interpreted metadata. The merged
metadata is supplied to later reconstruction and analysis.
"""
process(project::AbstractProject, item::AbstractDataItem) = process(item)

"""
    analyze(project, item) -> Dict

Analyze one processed item into additional metadata. The returned metadata is merged over the
processed item's metadata and supplied to later reconstruction. The default forwards to
`analyze(item)`.
"""
analyze(project::AbstractProject, item::AbstractDataItem) = analyze(item)

"""
    process(project, collection, items) -> Vector{<:AbstractDataItem}
    process(collection, items) -> Vector{<:AbstractDataItem}

Rewrite one collection's members, returning exactly one output per input with matching ids. The
default passes the members through unchanged.
"""
process(project::AbstractProject, collection::AbstractCollection, items::AbstractVector) =
    process(collection, items)

process(::AbstractCollection, items::AbstractVector) = items

"""
    analyze(project, collection, items) -> Dict
    analyze(collection, items) -> Dict

Fold one collection's processed members into metadata attached to the collection. The default is
empty.
"""
analyze(project::AbstractProject, collection::AbstractCollection, items::AbstractVector) =
    analyze(collection, items)

analyze(::AbstractCollection, items::AbstractVector)::Dict = Dict()

# ---------------------------------------------------------------------------
# reconstruct
# ---------------------------------------------------------------------------

"""
    reconstruct(::Type{T}, id, data, metadata::Dict) -> Union{Nothing,T}

Rebuild one concrete item from its stored id, cached payload, and the metadata available at the
calling stage. The engine can call this method before item analysis, after item analysis, or while
preparing collection input. The item id, payload, and metadata must therefore contain everything
the method needs at that point.

The default returns the payload when it is already a `T`, preserving the identity of custom items
whose `item_data(item)` is the item itself. Otherwise it returns `nothing`. While preparing an
interpreted input, the engine reruns `entries`, reusing its cached read input or calling `read` if
that input is absent. While preparing a processed input, it also reruns `process`. The fallback restores those stage outputs, but it cannot add later metadata
to a custom item. A type that needs later metadata must retain it in its `reconstruct` result.
Cached payloads remain available for views either way; reconstruction is needed to run further
project dispatch on a cached item.

Rehydration must be a pure function of cached content. Anything a type needs to rebuild itself
belongs in its id, data, or metadata, never in live workspace state.
"""
reconstruct(::Type{T}, id::AbstractString, data, metadata::Dict) where {T} =
    data isa T ? data : nothing

"""
    reconstruct(::Type{T}, id, metadata::Dict) -> T

Rebuild one collection value from its stored id and its own metadata. Collections hold no payload,
so those two are everything the row keeps about the value.

Unlike the item method there is no default and no way to opt out: an item that cannot be rebuilt is
recreated by rerunning `read` → `entries`, but a collection has no such path, so a project that
defines collection types must make them rebuildable. The engine rebuilds collections from the very
first interpretation, not only after reopening, so a missing method fails immediately.
"""
reconstruct(::Type{T}, id::AbstractString, metadata::Dict) where {T<:AbstractCollection} =
    error("Collection type $T must implement reconstruct(::Type{$T}, id, metadata)")
