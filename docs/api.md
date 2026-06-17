# Public API

This is the complete public surface of `MeasurementBrowser` — the names it exports. Everything else
is engine-internal and reachable only as `MeasurementBrowser.name` when you genuinely need it (and
may change without notice). The package ships **no domain types**: a project supplies its kinds,
readers, and plots through this API.

## Run the app

```julia
ws = open_workspace(project, "/path/to/data")   # open a data root; scans + caches in the background
open_browser(ws)                                # launch the CImGui browser on the workspace
close_workspace!(ws)                            # stop background work and release the root
select_items!(ws, items)                        # set the workspace selection (Vector{<:AbstractDataItem})
```

## Build a project

A project is a value built by mutation. Each `register_*` points at plain callbacks; re-calling with
the same key replaces that recipe in place, so editing and re-running a line updates a live project.

```julia
project = define_project("MyProject"; description="…")

register_item!(project, :iv;
    detect  = file -> endswith(file.filename, ".csv"),   # Bool; first matching recipe wins
    read    = file -> DataFrame(...),                    # whole file, parsed once
    entries = (file, data) -> [DataItem(kind=:iv, collection=[dev], parameters=p) for …],
    process = (item, data) -> data,                      # optional; default passthrough
    stats   = (item, processed) -> Dict{Symbol,Any}(…),  # optional; merged into the item's stats
    label   = item -> "…")                               # optional; display label

register_collection_stat!(project; kinds=[:iv],          # cross-item fold over each collection
    compute_stats = group -> …)                          # mutates each record's `stats` in place

register_plot!(project, :iv; label="I–V", setup=…, draw=…)   # one or more plots per kind
```

- **`entries(file, data)`** enumerates the items in one file. It returns `Vector{<:AbstractDataItem}`
  — the package's `DataItem` (the **recipe API**) or your own subtype (the **type API**, below). The
  engine fills each item's `filepath`/`timestamp` from the `SourceFile` and mints its `unique_id` from
  filepath + kind + `parameters`, then derives the internal record via the contract. **Project code
  never constructs or names the internal record.**
- **`register_plot!`** `setup(workspace, items)` returns the `Figure`; `draw(workspace, items, figure)`
  fills it. `items` are the loaded items for the selection, each carrying its payload as `item.data`.
  See [plans/plotting-api-design.md](plans/plotting-api-design.md).

## Implement a custom item — the `AbstractDataItem` contract

To go beyond the generic `DataItem`, subtype `AbstractDataItem`, carry metadata as typed fields and
the payload directly, and have your recipe's `entries` return instances of it. Implement:

```julia
struct Photo <: AbstractDataItem
    id::String; collection::Vector{String}; exposure::Float64; data::Matrix{RGB}
end
item_id(p::Photo)      = p.id
item_label(p::Photo)   = "exp=$(p.exposure)"
kind(::Photo)          = :photo
collection(p::Photo)   = p.collection
parameters(p::Photo)   = Dict{Symbol,Any}(:exposure => p.exposure)
stats(p::Photo)        = Dict{Symbol,Any}()
item_data(p::Photo)    = p.data
# optional, with defaults: process(::AbstractDataItem, data) = data ; cacheable(::AbstractDataItem) = false
# read_data is the engine-internal loader for the type API.

register_item!(project, :photo; detect=…, read=…, entries=(file,data)->[Photo(…) for …])
```

The engine indexes the subtype through the same contract the package's `DataItem` answers, and at
view time re-runs `read`/`entries` to hand your items (with `item.data`) straight to the plot.

## Types you name

- **`Project`** — the value `define_project` returns and `register_*` mutate.
- **`SourceFile`** — one physical file (`.filename`, `.filepath`, `.timestamp`), passed to callbacks.
- **`DataItem`** — the normal concrete item; subtype it (`AbstractDataItem`) or let `entries` build it.

See [ARCHITECTURE.md](ARCHITECTURE.md) for how these fit the engine, and
[plans/data-model-generalization.md](plans/data-model-generalization.md) for the design rationale.
