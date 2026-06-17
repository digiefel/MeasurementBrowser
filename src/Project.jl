"""
Project.jl - methods of the registration-based project API.

The `Project` type and its recipe types are defined in the `Projects` module (early, so other
submodules can name the concrete type). This file defines the methods: `define_project`, the
`register_*` functions, and the project's implementation of the engine interface (interpretation,
data access, plotting, serialization).

A project is built by mutation with `register_item!` / `register_collection_stat!` /
`register_plot!`, each pointing at plain callbacks, which the engine then drives:

Pipeline per item kind (fixed order, each fed the previous output):
    detect(file)          -> Bool
    read(file)            -> data               # whole file, parsed once
    entries(file, data)   -> Vector{<:AbstractDataItem}  # one entry per item (a single one is a vector of one)
    process(item, data)   -> data               # optional; default passthrough
    stats(item, processed)-> Dict{Symbol,Any}   # optional
    label(item)           -> String             # optional

`entries` returns the package's `DataItem` (recipe API) or a project subtype (type API); the engine
derives the internal record from each via the `AbstractDataItem` contract. `file` is a `SourceFile`
(has `.filename`, `.filepath`, `.timestamp`).
"""

using DataFrames: DataFrame
import Serialization

# The project is reachable from a cached SourceScan (twice: source.project and hierarchy.project), so
# it is serialized into the HDF5 cache. Only the registered recipes are persisted; transient scan
# state (read cache, profiling, locks) is rebuilt empty on load. This keeps the cache format stable
# when transient fields change, so adding scan instrumentation never silently invalidates a cache.
# serialize_cycle/deserialize_cycle preserve object identity, so the two references still resolve to
# one project after load.
function Serialization.serialize(s::Serialization.AbstractSerializer, project::Project)
    Serialization.serialize_cycle(s, project) && return nothing
    Serialization.serialize_type(s, Project, true)
    Serialization.serialize(s, project.name)
    Serialization.serialize(s, project.description)
    Serialization.serialize(s, project.recipes)
    Serialization.serialize(s, project.collection_stats)
    Serialization.serialize(s, project.plots)
    return nothing
end

function Serialization.deserialize(s::Serialization.AbstractSerializer, ::Type{Project})
    project = Project(
        "", "", ItemRecipe[],
        Dict{Tuple{Vararg{Symbol}},CollectionStatRecipe}(),
        Dict{Symbol,Dict{String,PlotRecipe}}(),
        Tuple{String,String,String}[], ReentrantLock(),
        Dict{Symbol,KindProfile}(), ReentrantLock(),
    )
    Serialization.deserialize_cycle(s, project)
    project.name = Serialization.deserialize(s)
    project.description = Serialization.deserialize(s)
    project.recipes = Serialization.deserialize(s)
    project.collection_stats = Serialization.deserialize(s)
    project.plots = Serialization.deserialize(s)
    return project
end

# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------

"""Name of a callback, or "λ" for an anonymous function."""
_callback_name(f)::String = (n = string(nameof(f)); startswith(n, "#") ? "λ" : n)

"""Count with a pluralized noun, e.g. `1 plot`, `3 plots`."""
_plural(n::Integer, word::AbstractString)::String = "$n $word" * (n == 1 ? "" : "s")

Base.show(io::IO, project::Project) =
    print(io, "Project(\"$(project.name)\", ", _plural(length(project.recipes), "item kind"), ")")

function Base.show(io::IO, ::MIME"text/plain", project::Project)
    print(io, "Project \"$(project.name)\"")
    isempty(project.description) || print(io, " · ", project.description)

    print(io, "\n  ", _plural(length(project.recipes), "item kind"))
    isempty(project.recipes) || print(io, " (detection order):")
    pad = isempty(project.recipes) ? 0 : maximum(length(string(r.kind)) for r in project.recipes)
    for recipe in project.recipes
        optional = [name for (name, f) in
            (("process", recipe.process), ("stats", recipe.stats), ("label", recipe.label))
            if f !== nothing]
        steps = isempty(optional) ? "" : " · " * join(optional, " · ")
        print(io, "\n    ", rstrip(rpad(string(recipe.kind), pad) * steps))
    end

    print(io, "\n  ", _plural(length(project.collection_stats), "collection stat"))
    isempty(project.collection_stats) || print(io, ":")
    for recipe in values(project.collection_stats)
        print(io, "\n    (", join(recipe.kinds, ", "), ") → ", _callback_name(recipe.compute_stats))
    end

    plot_count = sum(length, values(project.plots); init=0)
    print(io, "\n  ", _plural(plot_count, "plot"))
    isempty(project.plots) || print(io, ":")
    for kind in sort!(collect(keys(project.plots)); by=String)
        for label in sort!(collect(keys(project.plots[kind])))
            print(io, "\n    ", kind, " → \"", label, "\"")
        end
    end
end

"""Create an empty project. Register items/stats/plots into it afterwards."""
function define_project(name::AbstractString; description::AbstractString="")::Project
    return Project(
        String(name),
        String(description),
        ItemRecipe[],
        Dict{Tuple{Vararg{Symbol}},CollectionStatRecipe}(),
        Dict{Symbol,Dict{String,PlotRecipe}}(),
        Tuple{String,String,String}[],
        ReentrantLock(),
        Dict{Symbol,KindProfile}(),
        ReentrantLock(),
    )
end

"""Accumulate one timed file read (and the measurements it expanded to) into the scan profile."""
function _record_read!(
    project::Project,
    kind::Symbol,
    seconds::Float64,
    measurements::Int,
)::Nothing
    lock(project.profile_lock) do
        entry = get!(KindProfile, project.scan_profile, kind)
        entry.files += 1
        entry.read_seconds += seconds
        entry.measurements += measurements
    end
    return nothing
end

"""Accumulate one timed process+stats computation into the per-kind scan profile."""
function _record_stats!(project::Project, kind::Symbol, seconds::Float64)::Nothing
    lock(project.profile_lock) do
        entry = get!(KindProfile, project.scan_profile, kind)
        entry.stats_seconds += seconds
    end
    return nothing
end

function reset_scan_profile!(project::Project)::Nothing
    lock(project.profile_lock) do
        empty!(project.scan_profile)
    end
    lock(project.stat_failures_lock) do
        empty!(project.stat_failures)
    end
    return nothing
end

"""Log and record one per-measurement analysis failure raised while interpreting a file."""
function _record_stat_failure!(
    project::Project,
    filepath::String,
    measurement_id::String,
    step::Symbol,
    err,
)::Nothing
    @error("Measurement analysis failed", file=filepath, measurement=measurement_id, step,
        exception=(err, catch_backtrace()))
    lock(project.stat_failures_lock) do
        push!(project.stat_failures, (filepath, measurement_id, "step=$step: " * sprint(showerror, err)))
    end
    return nothing
end

function scan_profile_summary(project::Project)::Vector{NamedTuple}
    rows = lock(project.profile_lock) do
        [(; kind, files=p.files, measurements=p.measurements,
            read_seconds=p.read_seconds, stats_seconds=p.stats_seconds)
         for (kind, p) in project.scan_profile]
    end
    sort!(rows; by=row -> row.read_seconds + row.stats_seconds, rev=true)
    return rows
end

"""
Register (or replace) the recipe for one item kind.

`detect`, `read`, and `entries` are required. `entries(file, data)` returns the per-item entries as
`AbstractDataItem`s — the package's `DataItem` (recipe API) or your own subtype (type API); a
single-item file just returns a vector of one. Optional `process`/`stats`/`label` refine the recipe
API. Re-calling with the same `kind` replaces the recipe in place, so editing and re-running the line
updates a live project.
"""
function register_item!(
    project::Project,
    kind::Symbol;
    detect::Function,
    read::Function,
    entries::Function,
    process::Union{Nothing,Function}=nothing,
    stats::Union{Nothing,Function}=nothing,
    label::Union{Nothing,Function}=nothing,
)::Project
    recipe = ItemRecipe(kind, detect, read, entries, process, stats, label)
    index = findfirst(r -> r.kind === kind, project.recipes)
    index === nothing ? push!(project.recipes, recipe) : (project.recipes[index] = recipe)
    return project
end

"""
Register (or replace) a cross-item stat folded over each collection's items.

`compute_stats` receives the group `Vector{ItemRecord}` and fills each item's `stats` in place.
The engine groups by `group_by` (default: collection path) and runs this after the full scan.
Re-calling with the same `kinds` replaces it.
"""
function register_collection_stat!(
    project::Project;
    kinds::Vector{Symbol},
    compute_stats::Function,
    group_by::Function=_default_collection_group,
)::Project
    key = Tuple(sort(kinds))
    project.collection_stats[key] = CollectionStatRecipe(kinds, group_by, compute_stats)
    return project
end

"""
Register (or replace) one plot recipe for an item kind.

`setup(workspace, items)` builds and returns the `Figure`; `draw(workspace, items, figure)` fills it
in. `items` are the loaded, data-bearing items for the selection (`Vector{<:AbstractDataItem}`); each
carries its processed payload as `item.data`, materialized by the package through its processed-data
cache, so neither callback calls `process_item_data` itself. A kind may have multiple plots;
re-registering the same `label` for the same kind replaces that plot, which keeps REPL iteration
stable.
"""
function register_plot!(
    project::Project,
    kind::Symbol;
    label::AbstractString,
    setup::Function,
    draw::Function,
)::Project
    label_string = String(label)
    recipes = get!(project.plots, kind) do
        Dict{String,PlotRecipe}()
    end
    recipes[label_string] = PlotRecipe(kind, label_string, setup, draw)
    return project
end

_default_collection_group(record::ItemRecord)::String = collection_path_key(record.collection)

# ---------------------------------------------------------------------------
# Recipe lookup helpers
# ---------------------------------------------------------------------------

_recipe(project::Project, kind::Symbol)::Union{Nothing,ItemRecipe} =
    (i = findfirst(r -> r.kind === kind, project.recipes); i === nothing ? nothing : project.recipes[i])

function _detect_recipe(project::Project, file::SourceFile)::Union{Nothing,ItemRecipe}
    for recipe in project.recipes
        recipe.detect(file)::Bool && return recipe
    end
    return nothing
end

"""Engine-generated stable identity: file + kind + the params that distinguish siblings."""
function _mint_id(filepath::AbstractString, kind::Symbol, params::Dict{Symbol,Any})::String
    isempty(params) && return "$(filepath)#kind=$(kind)"
    ordered = sort!(collect(keys(params)))
    return "$(filepath)#kind=$(kind)," * join(("$(k)=$(params[k])" for k in ordered), ",")
end

# ---------------------------------------------------------------------------
# Engine interface implementation
# ---------------------------------------------------------------------------

project_name(project::Project)::String = project.name
project_description(project::Project)::String = project.description

function detect_kind(project::Project, filename::String)::Symbol
    recipe = _detect_recipe(project, index_source_file(filename))
    return recipe === nothing ? :unknown : recipe.kind
end

function parse_collection_info(::Project, file::SourceFile)::Vector{String}
    error("Project sets collection placement inside `measurements`, not via parse_collection_info ($(file.filepath))")
end

function kind_label(::Project, kind::Symbol)::String
    return string(kind)
end

function display_label(project::Project, measurement::ItemRecord)::String
    recipe = _recipe(project, measurement.kind)
    (recipe === nothing || recipe.label === nothing) && return _default_label(measurement)
    return recipe.label(measurement)::String
end

function _default_label(measurement::ItemRecord)::String
    parts = Any[measurement.timestamp, string(measurement.kind)]
    return join(filter(!isnothing, parts), " ")
end

"""
Interpret one file: read it, enumerate its items, and compute their per-item stats while the parse is
still in hand.

`entries` returns `AbstractDataItem`s (the package's `DataItem` or a project subtype); the engine
derives the internal `ItemRecord` from each via the contract, stamping the file's path/timestamp and
an engine-minted id. Folding stats into this single pass frees the (often large) parse as soon as the
file's items are analyzed, so scan memory stays bounded to the files in flight. Per-item failures are
recorded on the project and drained by `compute_and_add_item_stats!`.
"""
function interpret_file(project::Project, file::SourceFile)::Vector{ItemRecord}
    recipe = _detect_recipe(project, file)
    recipe === nothing && return ItemRecord[]
    read_started = time_ns()
    data = recipe.read(file)
    read_seconds = (time_ns() - read_started) / 1e9
    items = recipe.entries(file, data)::Vector{<:AbstractDataItem}
    # Note which API this recipe uses, so the view path knows how to re-materialize data: a type-API
    # `entries` returns data-bearing items, a recipe-API one returns metadata-only `DataItem`s.
    isempty(items) || (recipe.entries_carry_data = item_data(first(items)) !== nothing)
    records = ItemRecord[_record_from_item(recipe, file, item) for item in items]
    _record_read!(project, recipe.kind, read_seconds, length(records))
    if recipe.stats !== nothing
        for (item, record) in zip(items, records)
            try
                stats_started = time_ns()
                processed = recipe.process === nothing ? data : recipe.process(item, data)
                merge!(record.stats, recipe.stats(item, processed))
                _record_stats!(project, recipe.kind, (time_ns() - stats_started) / 1e9)
            catch err
                is_job_cancelled(err) && rethrow()
                _record_stat_failure!(project, record.filepath, record.unique_id, :stats, err)
            end
        end
    end
    # `data` falls out of scope here, so the parse becomes collectable once this file is done.
    return records
end

"""Derive the internal record from one `entries` item, stamping file context and a minted id."""
_record_from_item(recipe::ItemRecipe, file::SourceFile, item::AbstractDataItem)::ItemRecord =
    ItemRecord(
        item;
        filepath=file.filepath,
        filename=file.filename,
        timestamp=file.timestamp,
        kind=recipe.kind,
        unique_id=_mint_id(file.filepath, recipe.kind, parameters(item)),
    )

"""Direct data for an item is the whole-file read; `process` slices it later (recipe API)."""
function load_source_data(
    project::Project,
    source_file::SourceFile;
    record::Union{Nothing,ItemRecord}=nothing,
)
    kind = record === nothing ?
        detect_kind(project, source_file.filename) : record.kind
    recipe = _recipe(project, kind)
    recipe === nothing && error("No registered item recipe for kind $kind")
    return recipe.read(source_file)
end

"""Processed data: run the recipe's `process` over the cached whole-file read (recipe API)."""
function process_item_data(
    workspace::Workspace.Workspace{Project},
    record::ItemRecord,
)
    recipe = _recipe(workspace.project, record.kind)
    recipe === nothing && error("No registered item recipe for kind $(record.kind)")
    raw = only(read_item_data(workspace, [record]))
    recipe.process === nothing && return raw
    return recipe.process(DataItem(record, nothing), raw)
end

"""
Materialize the loaded, data-bearing items for a selection of records.

Recipe-API records resolve their payload through the cached `read`/`process` path and become a
`DataItem` carrying it. Type-API records (whose `entries` returns data-bearing items) re-run
`read`/`entries` for their file and return the rebuilt items directly, matched to records by id — the
same re-read the recipe path already pays, with no parallel data array.
"""
function materialize_items(
    workspace::Workspace.Workspace{Project},
    records::Vector{ItemRecord},
)::Vector{AbstractDataItem}
    result = Vector{AbstractDataItem}(undef, length(records))
    positions_by_file = Dict{String,Vector{Int}}()
    for (position, record) in pairs(records)
        push!(get!(positions_by_file, record.filepath, Int[]), position)
    end
    for (filepath, positions) in positions_by_file
        recipe = _recipe(workspace.project, records[positions[1]].kind)
        if recipe !== nothing && recipe.entries_carry_data
            file = index_source_file(filepath)
            raw = recipe.read(file)
            built = Dict(
                _mint_id(filepath, recipe.kind, parameters(item)) => item
                for item in recipe.entries(file, raw)::Vector{<:AbstractDataItem}
            )
            for position in positions
                result[position] = built[records[position].unique_id]
            end
        else
            selected = records[positions]
            processed = process_item_data(workspace, selected)
            for (offset, position) in pairs(positions)
                result[position] = DataItem(records[position], processed[offset])
            end
        end
    end
    return result
end

"""
Run the cross-measurement device stats, then return every analysis failure for the scan.

Per-measurement process/stats already ran inside `interpret_file` (one pass over each file's parse),
so this only does the disjoint device-scoped folds, which read the measurements' computed `stats`
rather than any raw data. Per-measurement failures collected during interpretation are drained and
returned alongside any device-fold failures.
"""
function compute_and_add_item_stats!(
    project::Project,
    measurements::Vector{ItemRecord},
    files::Vector{SourceFile};
    on_progress::Union{Nothing,Function}=nothing,
)::Vector{ItemFailure}
    # The cancel callback lives in task-local storage, which spawned tasks do not inherit, so capture
    # it here and re-establish it inside each task via `with_cancel`.
    failures = ItemFailure[]
    failures_lock = ReentrantLock()
    add_failure! = (filepath, id, step, err) -> lock(failures_lock) do
        push!(failures, _step_failure(filepath, id, step, err))
    end
    cancel_requested = get(task_local_storage(), CANCEL_CALLBACK_KEY, nothing)

    # Cross-item collection-stat folds. Groups within one recipe are disjoint, so they fold in
    # parallel; separate stat recipes stay sequential.
    recipes = collect(values(project.collection_stats))
    total = length(recipes)
    done = Base.Threads.Atomic{Int}(0)
    for recipe in recipes
        check_cancel()
        @sync for group in values(_group_for_collection_stat(measurements, recipe))
            isempty(group) && continue
            Base.Threads.@spawn with_cancel(cancel_requested) do
                check_cancel()
                try
                    recipe.compute_stats(group)
                catch err
                    is_job_cancelled(err) && rethrow()
                    add_failure!(group[1].filepath, group[1].unique_id, :compute_stats, err)
                end
            end
        end
        if on_progress !== nothing
            processed = Base.Threads.atomic_add!(done, 1) + 1
            on_progress((total=total, processed=processed, current_path=""))
        end
    end

    # Drain the per-measurement failures gathered while interpreting files.
    lock(project.stat_failures_lock) do
        for (filepath, measurement_id, message) in project.stat_failures
            push!(failures, ItemFailure(filepath, measurement_id, message))
        end
        empty!(project.stat_failures)
    end
    return failures
end

function _group_for_collection_stat(
    measurements::Vector{ItemRecord},
    recipe::CollectionStatRecipe,
)::Dict{Any,Vector{ItemRecord}}
    kinds = Set(recipe.kinds)
    groups = Dict{Any,Vector{ItemRecord}}()
    for measurement in measurements
        measurement.kind in kinds || continue
        push!(get!(groups, recipe.group_by(measurement), ItemRecord[]), measurement)
    end
    return groups
end

function _step_failure(filepath::String, measurement_id::String, step::Symbol, err)::ItemFailure
    @error("Measurement analysis failed", file=filepath, measurement=measurement_id, step, exception=(err, catch_backtrace()))
    return ItemFailure(filepath, measurement_id, "step=$step: " * sprint(showerror, err))
end
