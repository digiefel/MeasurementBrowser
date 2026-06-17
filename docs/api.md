# Public API

MeasurementBrowser has **two layers**:

- a **low-level, type-based source API** — the foundation the engine is written against. You model a
  data origin as an `AbstractDataSource`, the units it discovers as `AbstractDataSourceItem`s, and the
  logical browsable objects as `AbstractDataItem`s. No registration, no callbacks: you implement an
  interface and hand the engine a configured source value.
- a **high-level callback API** — a convenience built privately on top of the low level. You
  `define_project` and `register_*` callbacks; the engine wraps that in a built-in source adapter.

The high-level callback API is the **exported** surface today and is what most projects use. The
low-level types are currently engine-internal — reachable as `MeasurementBrowser.name` and staged for
a future dedicated submodule — but they are first-party and the callback API is implemented entirely
through them. The package ships **no domain types**: a source/project supplies its kinds, readers,
and plots.

## Run the app

```julia
ws = open_workspace(mysource)                       # low level: any AbstractDataSource
ws = open_workspace("/path/to/data"; project)       # high level: a project on a data root
open_browser(ws)                                    # launch the CImGui browser on the workspace
close_workspace!(ws)                                # stop background work and release the source
select_items!(ws, items)                            # set the selection (Vector{<:AbstractDataItem})
```

`open_workspace`'s first argument is reinterpreted: a `String` is a data root (the high-level path —
the `project` keyword's recipes interpret it, via the built-in `RegisteredProjectSource` adapter), an
`AbstractDataSource` is opened directly (no `project` needed). `project` stays a keyword.

---

## The low-level source API

Three abstract types, each a small documented interface. A value "is a source / source item / data
item" iff it answers the interface — there is no required field or supertype field to declare.

### `AbstractDataSource` — lifecycle and discovery

A source owns the lifecycle (open/close, optional watch) and discovery of source items. It may be a
dataset root, a database query, an instrument session, a live stream, or a remote service.

| Method | Required | Meaning |
|---|---|---|
| `source_id(source)::String` | yes | stable identity for cache / workspace ownership |
| `source_label(source)::String` | yes | user-facing source name |
| `open_source(source)::AbstractDataSource` | yes | prepare resources; return the opened source (simple sources return themselves) |
| `close_source!(source)::Nothing` | yes | release files, sockets, tasks, sessions, streams |
| `source_items(source)::Vector{<:AbstractDataSourceItem}` | yes | scan and return the current discovered units |
| `collection_stats(source, collection, items)::Dict{Symbol,Any}` | no (→ `Dict()`) | cross-item fold over one collection node (see below) |
| `source_fingerprint(source)` | no (→ `nothing`) | invalidation token for the whole source |
| `watch_source(source, on_change)` | no (→ `nothing`) | future live-update hook; `nothing` means static |

`open_source` returns the opened source rather than mutating in place, so it carries no bang;
`close_source!` mutates the source (releases its resources) and does. `source_items` and
`collection_stats` are getters — they return values and never mutate their arguments.

### `AbstractDataSourceItem` — one discovered unit

A source item is one addressable unit a source discovers — a file, a database row, a run id, a stream
channel, an image-stack member. It is the unit of scan progress, failure reporting, and source-level
invalidation. **It is not necessarily a browsable item**: one source item produces zero, one, or many
data items.

| Method | Required | Meaning |
|---|---|---|
| `source_item_id(item)::String` | yes | stable within `source_id(source)` |
| `source_item_label(item)::String` | yes | user-facing label for progress / errors |
| `data_items(source, source_item)::Vector{<:AbstractDataItem}` | yes | interpret one unit into lightweight (data-less) logical items for indexing |
| `load_data_item(source, source_item_id, item_id)::AbstractDataItem` | yes | reload one logical item later, with data available via `item_data` |
| `source_item_fingerprint(item)` | no (→ `nothing`) | invalidates records/payloads from this source item |
| `source_item_path(item)` | no (→ `nothing`) | filesystem path, when the unit has one |
| `source_item_timestamp(item)` | no (→ `nothing`) | acquisition/modification time, when known |

`data_items` enumerates and returns lightweight handles; `load_data_item` rebuilds one item with its
payload attached. Both are getters (no bang). The package's built-in file-backed source item is
`SourceFile` (below); a bare `SourceFile` has **no** universal `data_items` — the source decides what
a file means.

### `AbstractDataItem` — the logical browser item

The object the app indexes, selects, loads, inspects, and plots. Shared by both APIs (a callback's
`entries` returns `AbstractDataItem`s too), so these getters **are exported**.

| Method | Required | Meaning |
|---|---|---|
| `item_id(item)::String` | yes | stable within its source item |
| `item_label(item)::String` | yes | display title |
| `kind(item)::Symbol` | yes | coarse UI/plot key — *not* the source of semantic truth |
| `collection(item)::Vector{String}` | yes | canonical hierarchy path |
| `parameters(item)::Dict{Symbol,Any}` | yes | metadata known at interpretation time |
| `stats(item)::Dict{Symbol,Any}` | yes | computed values |
| `item_data(item)` | yes | loaded payload consumed by views/plots (`nothing` on a data-less handle) |
| `process(item, data)` | no (→ `data`) | optional transform applied to loaded data |
| `cacheable(item)::Bool` | no (→ `false`) | opt into persistent payload caching |
| `item_fingerprint(item)` | no (→ `nothing`) | invalidates this item's cached payload |

### A complete low-level source

No `register!` anywhere — a configured value passed to `open_workspace`:

```julia
struct PhotoDataset <: AbstractDataSource
    root::String
end
struct PhotoFile <: AbstractDataSourceItem
    path::String
    fingerprint::FileFingerprint
end
struct Photo <: AbstractDataItem
    id::String; label::String; collection::Vector{String}
    exposure::Float64; path::String; data::Union{Nothing,Matrix{Float64}}
end

# source
source_id(ds::PhotoDataset)    = abspath(ds.root)
source_label(::PhotoDataset)   = "Photo dataset"
open_source(ds::PhotoDataset)  = ds
close_source!(::PhotoDataset)  = nothing
source_items(ds::PhotoDataset) = [PhotoFile(p, file_fingerprint(p)) for p in photo_paths(ds.root)]

# source item
source_item_id(f::PhotoFile)          = f.path
source_item_label(f::PhotoFile)       = basename(f.path)
source_item_fingerprint(f::PhotoFile) = f.fingerprint
function data_items(::PhotoDataset, f::PhotoFile)
    m = read_photo_header(f.path)                                   # cheap, data-less
    [Photo(stable_id(f.path, m), m.label, m.collection, m.exposure, f.path, nothing)]
end
function load_data_item(ds::PhotoDataset, source_item_id, item_id)
    m = read_photo_header(source_item_id)
    Photo(item_id, m.label, m.collection, m.exposure, source_item_id, read_photo_matrix(source_item_id))
end

# data item
item_id(p::Photo)    = p.id
item_label(p::Photo) = p.label
kind(::Photo)        = :photo
collection(p::Photo) = p.collection
parameters(p::Photo) = Dict{Symbol,Any}(:exposure => p.exposure)
stats(p::Photo)      = Dict{Symbol,Any}()
item_data(p::Photo)  = p.data

ws = open_workspace(PhotoDataset("/path/to/photos"))
open_browser(ws)
```

### Low-level plotting

Plots receive loaded data items, never internal records. The source declares and draws them:

```julia
plot_kinds(source, items)::Vector{Type{<:PlotKind}}            # available plots for these items
plot_kind_label(source, plot_kind)::String                    # menu label
setup_plot(source, plot_kind, items)::Figure                  # build the figure
plot_data!(source, plot_kind, items, figure)::Nothing         # fill it; reads item_data
```

`PlotKind` and these methods are engine-internal for now (staged with the rest of the low level).

---

## The high-level callback API (convenience)

The exported, batteries-included path. A project is a value built by mutation; each `register_*`
points at plain callbacks. Re-calling with the same key replaces that recipe in place, so editing and
re-running a line updates a live project. Internally a project becomes a `RegisteredProjectSource`.

```julia
project = define_project("MyProject"; description="…")

register_item!(project, :iv;
    detect  = file -> endswith(file.filename, ".csv"),   # Bool; first matching recipe wins
    read    = file -> DataFrame(...),                    # whole file, parsed once
    entries = (file, data) -> [DataItem(kind=:iv, collection=[dev], parameters=p) for …],
    process = (item, data) -> data,                      # optional; default passthrough
    stats   = (item, processed) -> Dict{Symbol,Any}(…),  # optional; merged into the item's stats
    label   = item -> "…")                               # optional; display label

register_collection_stat!(project; kinds=[:iv],          # cross-item fold over each collection node
    compute_stats = group -> …)                          # returns a Dict stored on the node

register_plot!(project, :iv; label="I–V", setup=…, draw=…)   # one or more plots per kind
```

- **`entries(file, data)`** enumerates the items in one file, returning `Vector{<:AbstractDataItem}`
  — the package's `DataItem` (the **recipe API**) or your own subtype (the **type API**). The adapter
  fills each item's source-item identity (`filepath`/`timestamp`) from the `SourceFile` and mints
  `item_id` from filepath + kind + `parameters`, then derives the internal record via the contract.
  **Project code never constructs or names the internal record.**
- **`register_collection_stat!`** is the high-level form of the source's `collection_stats` hook: its
  `compute_stats` fold runs over each collection node's items and the result is stored on the node.
- **`register_plot!`** `setup(workspace, items)` returns the `Figure`; `draw(workspace, items, figure)`
  fills it, each item carrying its payload as `item.data`. See
  [plans/plotting-api-design.md](plans/plotting-api-design.md).

## Types you name

- **`Project`** — the value `define_project` returns and `register_*` mutate.
- **`SourceFile`** — the built-in file-backed `AbstractDataSourceItem` (`.filename`, `.filepath`,
  `.timestamp`, `.fingerprint`), passed to callbacks.
- **`DataItem`** — the normal concrete `AbstractDataItem`; subtype `AbstractDataItem` for the type API,
  or let `entries` build a `DataItem`.

See [ARCHITECTURE.md](ARCHITECTURE.md) for how these fit the engine, and
[plans/data-model-generalization.md](plans/data-model-generalization.md) for the design rationale.
</content>
</invoke>
