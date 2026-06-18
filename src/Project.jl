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

# RegisteredProjectSource serializes its project into the HDF5 cache. Only registered recipes are
# persisted; transient scan state (read cache, profiling, locks) is rebuilt empty on load. This keeps
# the cache format stable when transient fields change, so adding scan instrumentation never silently
# invalidates a cache.
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

"""Accumulate one timed source-item read and the items it expanded to into the scan profile."""
function _record_read!(
    project::Project,
    kind::Symbol,
    seconds::Float64,
    item_count::Int,
)::Nothing
    lock(project.profile_lock) do
        entry = get!(KindProfile, project.scan_profile, kind)
        entry.source_items += 1
        entry.read_seconds += seconds
        entry.items += item_count
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

"""Log and record one per-item analysis failure raised while interpreting a source item."""
function _record_stat_failure!(
    project::Project,
    filepath::String,
    item_id::String,
    step::Symbol,
    err,
)::Nothing
    @error("Item analysis failed", source_item=filepath, item=item_id, step,
        exception=(err, catch_backtrace()))
    lock(project.stat_failures_lock) do
        push!(project.stat_failures, (filepath, item_id, "step=$step: " * sprint(showerror, err)))
    end
    return nothing
end

function scan_profile_summary(project::Project)::Vector{NamedTuple}
    rows = lock(project.profile_lock) do
        [(; kind, source_items=p.source_items, items=p.items,
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

`compute_stats` receives the group `Vector{<:AbstractDataItem}` for one collection and returns a
`Dict{Symbol,Any}` stored on that collection node. The callback adapter implements this through the
low-level `collection_stats(source, collection, items)` hook. Re-calling with the same `kinds`
replaces the recipe.
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
cache, so neither callback resolves item data itself. A kind may have multiple plots;
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

_default_collection_group(item::AbstractDataItem)::String = collection_path_key(collection(item))

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

"""Engine-generated item identity within one source item."""
function _mint_item_id(kind::Symbol, params::Dict{Symbol,Any})::String
    isempty(params) && return "kind=$(kind)"
    ordered = sort!(collect(keys(params)))
    return "kind=$(kind)," * join(("$(k)=$(params[k])" for k in ordered), ",")
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

function kind_label(::Project, kind::Symbol)::String
    return string(kind)
end

function display_label(project::Project, item::ItemRecord)::String
    recipe = _recipe(project, item.kind)
    (recipe === nothing || recipe.label === nothing) && return item.item_label
    return recipe.label(DataItem(item, nothing))::String
end

Projects.source_id(source::RegisteredProjectSource)::String = source.root_path
Projects.source_label(source::RegisteredProjectSource)::String = project_name(source.project)
Projects.source_fingerprint(::RegisteredProjectSource) = nothing

function Projects.source_items(source::RegisteredProjectSource)::Vector{SourceFile}
    reset_scan_profile!(source.project)
    source.collection_metadata = load_scan_metadata(source.root_path)
    return collect_source_files(source.root_path)
end

ItemIndex.source_collection_metadata(source::RegisteredProjectSource) = source.collection_metadata

function _normalize_entry(recipe::ItemRecipe, item::AbstractDataItem)::AbstractDataItem
    id = item_id(item)
    isempty(id) || return item
    return DataItem(;
        kind=recipe.kind,
        collection=collection(item),
        label=item_label(item),
        parameters=parameters(item),
        stats=stats(item),
        data=item_data(item),
        unique_id=_mint_item_id(recipe.kind, parameters(item)),
    )
end

function Projects.data_items(
    project::Project,
    source::RegisteredProjectSource,
    file::SourceFile,
)::Vector{<:AbstractDataItem}
    recipe = _detect_recipe(project, file)
    recipe === nothing && return DataItem[]
    read_started = time_ns()
    data = recipe.read(file)
    read_seconds = (time_ns() - read_started) / 1e9
    items = AbstractDataItem[
        _normalize_entry(recipe, item)
        for item in recipe.entries(file, data)::Vector{<:AbstractDataItem}
    ]
    _record_read!(project, recipe.kind, read_seconds, length(items))
    if recipe.stats !== nothing
        for item in items
            try
                stats_started = time_ns()
                processed = recipe.process === nothing ? data : recipe.process(item, data)
                merge!(stats(item), recipe.stats(item, processed))
                _record_stats!(project, recipe.kind, (time_ns() - stats_started) / 1e9)
            catch err
                is_job_cancelled(err) && rethrow()
                _record_stat_failure!(project, file.filepath, item_id(item), :stats, err)
            end
        end
    end
    return items
end

function Projects.load_data_item(
    project::Project,
    source::RegisteredProjectSource,
    source_item_id::AbstractString,
    id::AbstractString,
)::AbstractDataItem
    file = index_source_file(source_item_id)
    recipe = _detect_recipe(project, file)
    recipe === nothing && error("No registered item recipe for source item $(source_item_id)")
    raw = recipe.read(file)
    items = AbstractDataItem[
        _normalize_entry(recipe, item)
        for item in recipe.entries(file, raw)::Vector{<:AbstractDataItem}
    ]
    index = findfirst(item -> item_id(item) == id, items)
    index === nothing && error("Source item $(source_item_id) did not produce item id $(id)")
    item = items[index]
    processed = recipe.process === nothing ? raw : recipe.process(item, raw)
    if recipe.process === nothing && !(item isa DataItem)
        return item
    end
    return DataItem(;
        kind=kind(item),
        collection=collection(item),
        label=item_label(item),
        parameters=parameters(item),
        stats=stats(item),
        data=processed,
        unique_id=item_id(item),
    )
end

"""Interpret one physical file through the high-level callback adapter."""
function items_for_file(
    project::Project,
    filepath::AbstractString;
    meta::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}}=nothing,
)::Vector{ItemRecord}
    source = RegisteredProjectSource(project, dirname(filepath))
    source.collection_metadata = meta
    records, _ = ItemIndex.interpret_source_item(
        project,
        source,
        index_source_file(filepath),
        source.collection_metadata,
    )
    return records
end

function Projects.collection_stats(
    project::Project,
    source::RegisteredProjectSource,
    ::Vector{String},
    items::Vector{<:AbstractDataItem},
)::Dict{Symbol,Any}
    merged = Dict{Symbol,Any}()
    for recipe in values(project.collection_stats)
        selected = AbstractDataItem[item for item in items if kind(item) in recipe.kinds]
        isempty(selected) && continue
        stats = recipe.compute_stats(selected)::Dict{Symbol,Any}
        merge!(merged, stats)
    end
    return merged
end

function _item_key(source_item_id::AbstractString, item_id::AbstractString)::String
    return "$(ncodeunits(source_item_id)):$(source_item_id)$(item_id)"
end

function ItemIndex.source_analysis_failures(
    project::Project,
    source::RegisteredProjectSource,
)::Vector{ItemFailure}
    failures = ItemFailure[]
    lock(project.stat_failures_lock) do
        for (source_item_id, item_id, message) in project.stat_failures
            push!(
                failures,
                ItemFailure(
                    source_item_id,
                    _item_key(source_item_id, item_id),
                    message,
                ),
            )
        end
        empty!(project.stat_failures)
    end
    return failures
end
