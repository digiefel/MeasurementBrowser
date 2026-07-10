"""
Project.jl - engine methods for the registration-based project API.

Construction (`define_project`, `register_*`) and the payload-agnostic contracts live in
`DataBrowserAPI`. This file keeps serialization, scan profiling, and callback-driving methods that
touch engine types (`ItemRecord`, `SourceFile`, `DataFrame`, …).
"""

using DataBrowserAPI
using DataBrowserAPI: AbstractDataItem, AbstractDataSource, KindProfileRow, SourceProfileRow, collection, id, item_data, item_label, kind, metadata
using DataBrowserSources: DirectorySource, SourceFile, index_source_file
using DataFrames: DataFrame
import Serialization
import DataBrowserAPI:
    _analyze_collection,
    _analyze_item,
    collection_path_label,
    data_items,
    detect_kind,
    display_label,
    finish_source_profile!,
    _has_collection_analysis,
    _has_collection_process,
    kind_label,
    process,
    _process_collection,
    project_description,
    project_name,
    record_scan_phase!,
    reset_scan_profile!,
    scan_profile_summary,
    scan_source_profile
import DataBrowserAPI.ItemIndex: DataItem, ItemRecord, interpret_source_item

# Only registered recipes are persisted; transient scan state (read cache, profiling, locks) is
# rebuilt empty on load. This keeps the cache format stable when transient fields change, so adding
# scan instrumentation never silently invalidates a cache.
function Serialization.serialize(s::Serialization.AbstractSerializer, project::Project)
    Serialization.serialize_cycle(s, project) && return nothing
    Serialization.serialize_type(s, Project, true)
    Serialization.serialize(s, project.name)
    Serialization.serialize(s, project.description)
    Serialization.serialize(s, project.recipes)
    Serialization.serialize(s, project.collections)
    Serialization.serialize(s, project.plots)
    return nothing
end

function Serialization.deserialize(s::Serialization.AbstractSerializer, ::Type{Project})
    project = Project(
        "", "", ItemRecipe[],
        Dict{Symbol,CollectionRecipe}(),
        Dict{Symbol,Dict{String,PlotRecipe}}(),
        Dict{String,SourceItemProfile}(), ReentrantLock(),
    )
    Serialization.deserialize_cycle(s, project)
    project.name = Serialization.deserialize(s)
    project.description = Serialization.deserialize(s)
    project.recipes = Serialization.deserialize(s)
    project.collections = Serialization.deserialize(s)
    project.plots = Serialization.deserialize(s)
    return project
end

function record_scan_phase!(
    project::Project,
    source_item_id::AbstractString,
    kind::Symbol,
    phase::Symbol,
    seconds::Float64,
    thread_id::Int,
)::Nothing
    lock(project.profile_lock) do
        entry = get!(() -> SourceItemProfile(source_item_id), project.scan_profile, source_item_id)
        kind !== :unmatched && (entry.kind = kind)
        phase === :detect ? (entry.detect_seconds += seconds) :
        phase === :read ? (entry.read_seconds += seconds) :
        phase === :entries ? (entry.entries_seconds += seconds) :
        phase === :process ? (entry.process_seconds += seconds) :
        phase === :analyze ? (entry.analyze_seconds += seconds) :
        error("Unknown scan timing phase: $phase")
        push!(entry.thread_ids, thread_id)
    end
    return nothing
end

function finish_source_profile!(
    project::Project,
    source_item_id::AbstractString,
    source_item_label::AbstractString,
    source_item_path::Union{Nothing,String},
    kind::Symbol,
    item_count::Int,
    total_seconds::Float64,
    thread_ids::Set{Int},
)::Nothing
    lock(project.profile_lock) do
        entry = get!(() -> SourceItemProfile(source_item_id), project.scan_profile, source_item_id)
        entry.source_item_label = String(source_item_label)
        entry.source_item_path = source_item_path
        entry.kind = kind
        entry.item_count = item_count
        # process_seconds and analyze_seconds are accumulated separately during analysis via
        # record_scan_phase!, so finish (which runs at scan time) must not overwrite them.
        entry.total_seconds = total_seconds
        union!(entry.thread_ids, thread_ids)
    end
    return nothing
end

function reset_scan_profile!(project::Project)::Nothing
    lock(project.profile_lock) do
        empty!(project.scan_profile)
    end
    return nothing
end

function scan_profile_summary(project::Project)::Vector{KindProfileRow}
    rows = lock(project.profile_lock) do
        by_kind = Dict{Symbol,KindProfileRow}()
        for profile in values(project.scan_profile)
            row = get!(() -> KindProfileRow(profile.kind), by_kind, profile.kind)
            row.source_items += 1
            row.items += profile.item_count
            row.detect_seconds += profile.detect_seconds
            row.read_seconds += profile.read_seconds
            row.entries_seconds += profile.entries_seconds
            row.process_seconds += profile.process_seconds
            row.analyze_seconds += profile.analyze_seconds
            row.total_seconds += profile.total_seconds
        end
        collect(values(by_kind))
    end
    sort!(rows; by=row -> row.total_seconds, rev=true)
    return rows
end

function scan_source_profile(project::Project)::Vector{SourceProfileRow}
    rows = lock(project.profile_lock) do
        [SourceProfileRow(
            profile.source_item_id,
            profile.source_item_label,
            profile.source_item_path,
            profile.kind,
            profile.item_count,
            profile.detect_seconds,
            profile.read_seconds,
            profile.entries_seconds,
            profile.process_seconds,
            profile.analyze_seconds,
            profile.total_seconds,
            sort!(collect(profile.thread_ids)),
        ) for profile in values(project.scan_profile)]
    end
    sort!(rows; by=row -> row.total_seconds, rev=true)
    return rows
end

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

"""Default item identity for recipe entries that do not provide one."""
function _mint_id(
    source_item_id::AbstractString,
    kind::Symbol,
    metadata::Dict{Symbol,Any},
)::String
    source_id = String(source_item_id)
    isempty(metadata) && return source_id
    ordered = sort!(collect(keys(metadata)))
    return "$source_id#kind=$(kind)," * join(("$(k)=$(metadata[k])" for k in ordered), ",")
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

function _normalize_entry(
    recipe::ItemRecipe,
    source_item_id::AbstractString,
    item::AbstractDataItem,
)::AbstractDataItem
    entry_id = id(item)
    isempty(entry_id) || return item
    return DataItem(;
        kind=recipe.kind,
        collection=collection(item),
        label=item_label(item),
        metadata=metadata(item),
        data=item_data(item),
        id=_mint_id(source_item_id, recipe.kind, metadata(item)),
    )
end

function _processed_item(recipe::ItemRecipe, item::AbstractDataItem)::AbstractDataItem
    recipe.process === nothing && return item
    result = recipe.process(item)
    result isa AbstractDataItem && return result
    error(
        "process callback for kind $(recipe.kind) must return an AbstractDataItem; " *
        "got $(typeof(result)). Return DataItem(item, processed) or a custom AbstractDataItem instead."
    )
end

"""
Process one loaded item through its registered callback or the low-level `process(item)` hook.
"""
function process(
    project::Project,
    ::AbstractDataSource,
    item::AbstractDataItem,
)::AbstractDataItem
    recipe = _recipe(project, kind(item))
    (recipe === nothing || recipe.process === nothing) && return process(item)
    return _processed_item(recipe, item)
end

function data_items(
    project::Project,
    source::AbstractDataSource,
    file::SourceFile,
)::Vector{<:AbstractDataItem}
    recipe = nothing
    started = time_ns()
    profiler = Profiling.current_session()
    try
        recipe = Profiling.@profile_span profiler :project :detect Profiling.ProfileAttributes(
            source_id=file.filepath,
        ) begin
            _detect_recipe(project, file)
        end
    finally
        record_scan_phase!(
            project, file.filepath,
            recipe === nothing ? :unmatched : recipe.kind,
            :detect, (time_ns() - started) / 1e9, Base.Threads.threadid())
    end
    recipe === nothing && return DataItem[]
    started = time_ns()
    data = try
        Profiling.@profile_span profiler :project :read Profiling.ProfileAttributes(
            kind=recipe.kind,
            source_id=file.filepath,
        ) begin
            recipe.read(file)
        end
    finally
        record_scan_phase!(
            project, file.filepath, recipe.kind, :read,
            (time_ns() - started) / 1e9, Base.Threads.threadid())
    end
    started = time_ns()
    items = try
        Profiling.@profile_span profiler :project :entries Profiling.ProfileAttributes(
            kind=recipe.kind,
            source_id=file.filepath,
        ) begin
            AbstractDataItem[
                _normalize_entry(recipe, file.filepath, item)
                for item in recipe.entries(file, data)::Vector{<:AbstractDataItem}
            ]
        end
    finally
        record_scan_phase!(
            project, file.filepath, recipe.kind, :entries,
            (time_ns() - started) / 1e9, Base.Threads.threadid())
    end
    return items
end

"""Interpret one physical file through the high-level callback adapter."""
function items_for_file(
    project::Project,
    filepath::AbstractString;
    meta::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}}=nothing,
)::Vector{ItemRecord}
    source = DirectorySource(dirname(filepath); metadata_file=nothing)
    if meta !== nothing
        source.collection_metadata_entries = meta
        source.has_metadata = true
    end
    interpretation = interpret_source_item(
        project,
        source,
        index_source_file(filepath),
    )
    return interpretation.records
end

"""
Rewrite one collection's members through the kind's registered collection `process`.

Members are grouped by kind; a kind without a registered `process` passes its members through
unchanged. The callback returns one output per input; the adapter validates ids and count.
"""
function _process_collection(
    project::Project,
    ::AbstractDataSource,
    ::Vector{String},
    items::Vector{<:AbstractDataItem},
)::Vector{<:AbstractDataItem}
    rewritten = AbstractDataItem[]
    for group in _group_by_kind(items)
        recipe = get(project.collections, kind(first(group)), nothing)
        if recipe === nothing || recipe.process === nothing
            append!(rewritten, group)
            continue
        end
        outputs = recipe.process(group)::Vector{<:AbstractDataItem}
        length(outputs) == length(group) || error(
            "collection process for kind $(kind(first(group))) must return one item per input; " *
            "got $(length(outputs)) for $(length(group)) members",
        )
        input_ids = Set(id(item) for item in group)
        for output in outputs
            id(output) in input_ids || error(
                "collection process for kind $(kind(first(group))) returned unknown item id " *
                "'$(id(output))'",
            )
        end
        append!(rewritten, outputs)
    end
    return rewritten
end

"""Fold one collection's post-process members into collection-node metadata."""
function _analyze_collection(
    project::Project,
    ::AbstractDataSource,
    ::Vector{String},
    items::Vector{<:AbstractDataItem},
)::Dict{Symbol,Any}
    merged = Dict{Symbol,Any}()
    for group in _group_by_kind(items)
        recipe = get(project.collections, kind(first(group)), nothing)
        (recipe === nothing || recipe.analyze === nothing) && continue
        merge!(merged, recipe.analyze(group)::Dict{Symbol,Any})
    end
    return merged
end

"""Group items by kind, preserving first-seen kind order."""
function _group_by_kind(items::Vector{<:AbstractDataItem})::Vector{Vector{AbstractDataItem}}
    groups = Dict{Symbol,Vector{AbstractDataItem}}()
    order = Symbol[]
    for item in items
        item_kind = kind(item)
        haskey(groups, item_kind) || push!(order, item_kind)
        push!(get!(() -> AbstractDataItem[], groups, item_kind), item)
    end
    return Vector{AbstractDataItem}[groups[item_kind] for item_kind in order]
end

function _analyze_item(
    project::Project,
    ::AbstractDataSource,
    item::AbstractDataItem,
)::Dict{Symbol,Any}
    recipe = _recipe(project, kind(item))
    (recipe === nothing || recipe.analyze === nothing) && return Dict{Symbol,Any}()
    return recipe.analyze(item)::Dict{Symbol,Any}
end

function _has_collection_process(project::Project, item_kind::Symbol)::Bool
    recipe = get(project.collections, item_kind, nothing)
    return recipe !== nothing && recipe.process !== nothing
end

function _has_collection_analysis(project::Project, item_kind::Symbol)::Bool
    recipe = get(project.collections, item_kind, nothing)
    return recipe !== nothing && recipe.analyze !== nothing
end
