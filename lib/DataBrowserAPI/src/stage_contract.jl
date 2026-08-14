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

**The source appears exactly once, at `read`.** Every later stage is a pure function of values: the
source item (an address — fingerprintable, recordable), the loaded payload, the items. That rule is
what makes a cached stage result sufficient on its own: rerunning `entries` without rereading,
warm reopen, eviction recovery, and cross-machine rehydration are sound by construction rather than
by discipline. A `read` returning a live handle (an HDF5 group, a database cursor) is legitimate
but uncacheable; that project's durable boundary is `process` instead.

`entries` still receives the source item because identity is not payload: the loaded value stays
purely the expensive data, while ids, labels, and collection placement derive from the source item —
which the engine can supply without touching the origin.

`read`'s return value and `entries`' `loaded` argument deliberately share no supertype. The loaded
value is a private handoff between two stages of the same project; the type discipline is enforced
by the concrete signature on the receiving end.
"""

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

The default treats the loaded value as a single item.
"""
entries(project::AbstractProject, item::AbstractDataSourceItem, loaded) = entries(item, loaded)

entries(::AbstractDataSourceItem, loaded) = [loaded]

# ---------------------------------------------------------------------------
# process and analyze
# ---------------------------------------------------------------------------

"""
    process(project, item) -> AbstractDataItem

Process one interpreted item. The default forwards to the context-free `process(item)`.
"""
process(project::AbstractProject, item::AbstractDataItem) = process(item)

"""
    analyze(project, item) -> Dict

Analyze one processed item into additional metadata. The default forwards to `analyze(item)`.
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

Rebuild one concrete item from its stored id, cached payload, and metadata. The default returns
`nothing`; the engine then reruns `read` → `entries` → `process`, which is always correct and only
slower. Cached payloads are still delivered to views either way; this is needed only to run further
project dispatch on a cached item.

Rehydration must be a pure function of cached content. Anything a type needs to rebuild itself
belongs in its id or metadata, never in live workspace state.
"""
reconstruct(::Type, id::AbstractString, data, metadata::Dict) = nothing

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
