# Public API

MeasurementBrowser has **two layers**:

- a **low-level, type-based source API** — the foundation the engine is written against. You model a
  data origin as an `AbstractDataSource`, the units it discovers as `AbstractDataSourceItem`s, and the
  logical browsable objects as `AbstractDataItem`s. A workspace is still opened for a `Project`; the
  source supplies discovery while project/source methods interpret and load items.
- a **high-level callback API** — a convenience built on the same source engine. You `define_project`
  and `register_*` callbacks; the built-in `DirectorySource` discovers files and the project recipes
  interpret them.

The high-level callback API is the **exported** surface today and is what most projects use. The
low-level types are currently engine-internal — reachable as `MeasurementBrowser.name` and staged for
a future dedicated submodule — but they are first-party and the callback API is implemented entirely
through them. The package ships **no domain types**: a source/project supplies its kinds, readers,
and plots.

## Run the app

```julia
ws = open_workspace(project, mysource)              # project-owned low-level source
ws = open_workspace(project, "/path/to/data")       # high level: a project on a data root
open_browser(ws)                                    # launch the CImGui browser on the workspace
close_workspace!(ws)                                # stop background work and release the source
select_items!(ws, items_or_records)                 # programmatic selection
```

Every workspace is associated with a `Project`; the first argument to `open_workspace` is always that
project. Passing `root_path` is the exported directory convenience path: it constructs a
`DirectorySource(root_path)` and opens the workspace with that source. Passing an
`AbstractDataSource` keeps the source explicit while still associating the workspace with the project
that owns labels, registered plots, and browser state.

`select_items!` accepts indexed records, exact item ids, or loaded `AbstractDataItem` values whose
ids are present in the current workspace.

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
| `collection_stats(project, source, collection, items)::Dict{Symbol,Any}` | no (→ `Dict()`) | background cross-item fold over one collection node |
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
| `data_items(project, source, source_item)::Vector{<:AbstractDataItem}` | yes | read/interpret one unit into loaded logical items for processing, stats, and indexing |
| `load_data_item(project, source, source_item_id, id)::AbstractDataItem` | yes | reload one logical item later, with data available via `item_data` |
| `fingerprint(item::AbstractDataSourceItem)` | no (→ `nothing`) | invalidates records and cached item data from this source item |
| `source_item_path(item)` | no (→ `nothing`) | filesystem path, when the unit has one |
| `source_item_timestamp(item)` | no (→ `nothing`) | acquisition/modification time, when known |

`data_items` reads/interprets one source item and returns its logical data items. Their loaded data is
processed, analyzed, and optionally cached before that source-item pass ends. `load_data_item`
rebuilds one item later when neither memory nor cache can serve it. Both are getters (no bang). The
package's built-in file-backed source item is
`SourceFile` (below); a bare `SourceFile` has **no** universal `data_items` — the source decides what
a file means.

### `AbstractDataItem` — the logical browser item

The object the app indexes, selects, loads, inspects, and plots. Shared by both APIs (a callback's
`entries` returns `AbstractDataItem`s too), so these getters **are exported**.

| Method | Required | Meaning |
|---|---|---|
| `id(item)::String` | yes | stable within its source |
| `item_label(item)::String` | yes | display title |
| `kind(item)::Symbol` | yes | coarse UI/plot key — *not* the source of semantic truth |
| `collection(item)::Vector{String}` | yes | canonical hierarchy path |
| `parameters(item)::Dict{Symbol,Any}` | yes | item parameters; entry callbacks provide local params, indexed/loaded `DataItem`s carry effective params |
| `stats(item)::Dict{Symbol,Any}` | yes | computed values filled during the source-item pass |
| `item_data(item)` | yes | loaded data consumed by processing/views; internal index handles carry no data |
| `process(item)` | no (→ `item`) | optional transform returning the item views receive |
| `cacheable(item)::Bool` | no (→ `false`) | opt into persistent item-data caching |
| `fingerprint(item::AbstractDataItem)` | no (→ `nothing`) | invalidates this item's cached data |

### A complete low-level source

No item registration is involved — the project still owns the workspace, while the configured source
value supplies discovery and loading:

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
fingerprint(f::PhotoFile) = f.fingerprint
function data_items(::Project, ::PhotoDataset, f::PhotoFile)
    m = read_photo_header(f.path)                                   # cheap, data-less
    [Photo(stable_id(f.path, m), m.label, m.collection, m.exposure, f.path, nothing)]
end
function load_data_item(::Project, ds::PhotoDataset, source_item_id, id)
    m = read_photo_header(source_item_id)
    Photo(id, m.label, m.collection, m.exposure, source_item_id, read_photo_matrix(source_item_id))
end

# data item
id(p::Photo)    = p.id
item_label(p::Photo) = p.label
kind(::Photo)        = :photo
collection(p::Photo) = p.collection
parameters(p::Photo) = Dict{Symbol,Any}(:exposure => p.exposure)
stats(p::Photo)      = Dict{Symbol,Any}()
item_data(p::Photo)  = p.data

project = define_project("Photos")
ws = open_workspace(project, PhotoDataset("/path/to/photos"))
open_browser(ws)
```

### Plotting

The external plot API is `register_plot!(project, kind; label, setup, draw)`. The workspace
materializes selected records into loaded data items before invoking the callbacks, so
`setup(workspace, items)` and `draw(workspace, items, figure)` receive loaded `AbstractDataItem`s and
never internal records.

`PlotKind`, `RegisteredPlot`, `setup_plot`, and `plot_data!` are engine internals used by the GUI to
dispatch registered plot callbacks. A source-level plot hook belongs with the future low-level export,
not the exported API today.

---

## The high-level callback API (convenience)

The exported, batteries-included path. A project is a value built by mutation; each `register_*`
points at plain callbacks. Re-calling with the same key replaces that recipe in place, so editing and
re-running a line updates a live project. `open_workspace(project, root_path)` uses `DirectorySource`
for directory discovery; project recipes interpret each discovered `SourceFile`.

```julia
project = define_project("MyProject"; description="…")

register_item!(project, :iv;
    detect  = file -> endswith(file.filename, ".csv"),   # Bool; first matching recipe wins
    read    = file -> DataFrame(...),                    # whole file, parsed once
    entries = (file, data) -> [DataItem(kind=:iv, collection=[dev], parameters=p, data=slice) for …],
    process = item -> DataItem(item, clean(item.data)),  # optional; default passthrough
    stats   = item -> Dict{Symbol,Any}(…),               # optional; computed after indexing
    label   = item -> "…")                               # optional; display label

register_collection_stat!(project; kinds=[:iv],          # cross-item fold over each collection node
    compute_stats = group -> …)                          # returns a Dict stored on the node

register_plot!(project, :iv; label="I–V", setup=…, draw=…)   # one or more plots per kind
```

- **`entries(file, data)`** enumerates the items in one file, returning `Vector{<:AbstractDataItem}`
  — the package's `DataItem` (the **recipe API**) or your own subtype (the **type API**). It attaches
  the raw per-item data as `item.data`; optional `process(item)` returns the item that stats and views
  receive. The source worker runs `process` and `stats` before releasing the file data, while records
  are streamed to the workspace as each source item finishes. The adapter derives each internal record from the returned item and the
  `SourceFile`. When a recipe entry does not provide an id, the adapter mints one from the
  source-item path, kind, and `parameters`.
  **Project code never constructs or names the internal record.**
- **`register_collection_stat!`** is the high-level form of the source's `collection_stats` hook: its
  `compute_stats` fold runs in background analysis over each collection node's items and the result
  is stored on the node.
- **`register_plot!`** `setup(workspace, items)` returns the `Figure`; `draw(workspace, items, figure)`
  fills it, each item carrying its payload as `item.data`. See
  [plans/plotting-api-design.md](plans/plotting-api-design.md).

## Types you name

- **`Project`** — the value `define_project` returns and `register_*` mutate.
- **`SourceFile`** — the built-in file-backed `AbstractDataSourceItem` (`.filename`, `.filepath`,
  `.timestamp`, `.fingerprint`) discovered by `DirectorySource` and passed to callbacks.
- **`DataItem`** — the normal concrete `AbstractDataItem`; subtype `AbstractDataItem` for the type API,
  or let `entries` build a `DataItem`.

See [ARCHITECTURE.md](ARCHITECTURE.md) for how these fit the engine, and
[plans/data-model-generalization.md](plans/data-model-generalization.md) for the design rationale.
</content>
</invoke>
