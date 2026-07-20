"""
Project.jl - engine methods for the registration-based project API.

Construction (`define_project`, `register_item!`, `register_collection_analysis!`) and the
payload-agnostic contracts live in `DataBrowserAPI`. This file keeps scan profiling and
callback-driving methods that touch engine types (`ItemRecord`, `SourceFile`, `DataFrame`, …).
"""

using DataBrowserAPI
using DataBrowserAPI: AbstractCollection, AbstractDataItem, AbstractDataSource, AbstractDataSourceItem, KindProfileRow, SourceProfileRow, analyze, fingerprint, source_id, source_item_id, source_item_label, source_item_path, source_item_timestamp
using DataBrowserSources: DirectorySource, SourceFile, index_source_file
using DataFrames: DataFrame
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
    collection,
    id,
    item_data,
    item_label,
    kind,
    metadata,
    process,
    _process_collection,
    project_description,
    project_name,
    record_scan_phase!,
    reset_scan_profile!,
    scan_profile_summary,
    scan_source_profile
import DataBrowserAPI.ItemIndex:
    CollectionIndex,
    CollectionInput,
    ItemFailure,
    ItemRecord,
    MetadataDict,
    RegisteredDataItem,
    collection_inputs,
    effective_record,
    resolve_collection_path!,
    metadata_dict

function _with_data(item::RegisteredDataItem, data)::RegisteredDataItem
    return RegisteredDataItem(item, data)
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

"""Stable private identity for one value returned by a registration."""
function _mint_id(
    source_item_id::AbstractString,
    kind::Symbol,
    position::Integer,
    supplied_key=nothing,
)::String
    suffix = supplied_key === nothing ? string(position) : string(supplied_key)
    return "$(source_item_id)#$(kind):$(suffix)"
end

"""Split the documented `(data=..., metadata=Dict(...))` result form."""
function _data_and_metadata(value)::Tuple{Any,MetadataDict}
    if value isa NamedTuple && keys(value) == (:data, :metadata)
        return value.data, metadata_dict(value.metadata)
    end
    return value, MetadataDict()
end

function _default_collection(
    source::AbstractDataSource,
    source_item::AbstractDataSourceItem,
)::Vector{String}
    path = source_item_path(source_item)
    path === nothing && return String[]
    root = source_id(source)
    isabspath(root) || return String[]
    relative_directory = dirname(relpath(path, root))
    relative_directory == "." && return String[]
    return splitpath(relative_directory)
end

function _registration_collection_path(value)::Vector{String}
    value isa AbstractVector || throw(ArgumentError(
        "a register_item! collection callback must return a vector of strings; got $(typeof(value))",
    ))
    all(segment -> segment isa AbstractString, value) || throw(ArgumentError(
        "a register_item! collection callback must return only strings; got $(repr(value))",
    ))
    return String[segment for segment in value]
end

function _typed_collection_path(item::AbstractDataItem)::Vector{AbstractCollection}
    value = collection(item)
    value isa AbstractVector || throw(ArgumentError(
        "collection(::$(typeof(item))) must return a vector of AbstractCollection values; " *
        "got $(typeof(value))",
    ))
    all(segment -> segment isa AbstractCollection, value) || throw(ArgumentError(
        "collection(::$(typeof(item))) must return only AbstractCollection values; got $(repr(value))",
    ))
    return AbstractCollection[segment for segment in value]
end

function _registered_item(
    recipe::ItemRecipe,
    source::AbstractDataSource,
    source_item::AbstractDataSourceItem,
    position::Integer,
    value,
    inherited_metadata::MetadataDict,
)::RegisteredDataItem
    data, entry_metadata = _data_and_metadata(value)
    local_metadata = merge(copy(inherited_metadata), entry_metadata)
    collection_value = recipe.collection === nothing ?
        _default_collection(source, source_item) : recipe.collection(data, local_metadata)
    collection_path = _registration_collection_path(collection_value)
    supplied_key = recipe.id === nothing ? nothing : recipe.id(data, local_metadata)
    label = recipe.label === nothing ? "" : String(recipe.label(data, local_metadata))
    return RegisteredDataItem(
        _mint_id(source_item_id(source_item), recipe.kind, position, supplied_key),
        label,
        recipe.kind,
        collection_path,
        data,
        local_metadata,
    )
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
    return item.item_label
end

function _processed_item(
    recipe::ItemRecipe,
    item::RegisteredDataItem,
)::RegisteredDataItem
    recipe.process === nothing && return item
    return _with_data(item, recipe.process(item.data, item.metadata))
end

"""
Process one loaded item through its registered callback or the low-level `process(item)` hook.
"""
function process(
    project::Project,
    ::AbstractDataSource,
    item::AbstractDataItem,
)::AbstractDataItem
    if !(item isa RegisteredDataItem)
        result = process(item)
        result isa AbstractDataItem || error(
            "process(::$(typeof(item))) must return an AbstractDataItem; got $(typeof(result))",
        )
        return result
    end
    recipe = _recipe(project, item.registration)
    recipe === nothing && error("Missing registration $(item.registration)")
    return _processed_item(recipe, item)
end

function data_items(
    project::Project,
    source::AbstractDataSource,
    file::SourceFile,
)::Vector{<:AbstractDataItem}
    recipe = nothing
    started = time_ns()
    try
        recipe = @time_dbg _detect_recipe(project, file)
    finally
        record_scan_phase!(
            project, file.filepath,
            recipe === nothing ? :unmatched : recipe.kind,
            :detect, (time_ns() - started) / 1e9, Base.Threads.threadid())
    end
    recipe === nothing && return AbstractDataItem[]
    started = time_ns()
    read_result = try
        @time_dbg "read" recipe.read(file)
    finally
        record_scan_phase!(
            project, file.filepath, recipe.kind, :read,
            (time_ns() - started) / 1e9, Base.Threads.threadid())
    end
    loaded_data, read_metadata = _data_and_metadata(read_result)
    inherited_metadata = merge(metadata_dict(metadata(file)), read_metadata)
    started = time_ns()
    items = try
        @time_dbg "entries" begin
            values = recipe.entries === nothing ? Any[loaded_data] :
                recipe.entries(loaded_data, inherited_metadata)
            values isa AbstractVector || error(
                "entries callback for registration $(recipe.kind) must return a vector; " *
                "got $(typeof(values))",
            )
            AbstractDataItem[
                _registered_item(
                    recipe, source, file, position, value, inherited_metadata)
                for (position, value) in pairs(values)
            ]
        end
    finally
        record_scan_phase!(
            project, file.filepath, recipe.kind, :entries,
            (time_ns() - started) / 1e9, Base.Threads.threadid())
    end
    return items
end

"""
One completed source-item pass.

`records` are retained by the index. `interpreted_items` carry effective item data only on direct
interpretation paths; workspace workers put that data in the memory cache before publishing a
lightweight completion.
"""
struct SourceItemInterpretation
    records::Vector{ItemRecord}
    collection_paths::Vector{Vector{CollectionInput}}
    interpreted_items::Vector{AbstractDataItem}
    failures::Vector{ItemFailure}
end

"""
Interpret every logical data item produced by one source item.

The source item is read only by `data_items`. Concrete typed items remain unchanged; registered data
uses a private carrier. Processing and analysis belong to the workspace work graph.
"""
function interpret_source_item(
    project::Project,
    source::AbstractDataSource,
    source_item::AbstractDataSourceItem,
)::SourceItemInterpretation
    source_item_id_value = source_item_id(source_item)
    source_item_path_value = source_item_path(source_item)
    source_item_label_value = if source_item_path_value !== nothing &&
                                 isabspath(source_id(source))
        relpath(source_item_path_value, source_id(source))
    else
        source_item_label(source_item)
    end
    source_started = time_ns()
    handles = try
        @time_dbg data_items(project, source, source_item)
    catch
        finish_source_profile!(
            project, source_item_id_value, source_item_label_value, source_item_path_value,
            :unmatched, 0,
            (time_ns() - source_started) / 1e9, Set([Base.Threads.threadid()]))
        rethrow()
    end
    item_count = length(handles)
    item_kinds = unique(kind(handle) for handle in handles)
    source_kind = length(item_kinds) == 1 ? only(item_kinds) :
        isempty(item_kinds) ? :unmatched : :mixed
    records = Vector{ItemRecord}(undef, item_count)
    collection_paths = Vector{Vector{CollectionInput}}(undef, item_count)
    interpreted_items = Vector{AbstractDataItem}(undef, item_count)
    source_metadata = metadata_dict(metadata(source_item))
    for (index, handle) in pairs(handles)
        item_metadata = handle isa RegisteredDataItem ? handle.metadata : metadata_dict(metadata(handle))
        record_metadata = merge(copy(source_metadata), item_metadata)
        item_id = isempty(id(handle)) ?
            _mint_id(source_item_id_value, kind(handle), index) : id(handle)
        label = item_label(handle)
        normalized_collection = handle isa RegisteredDataItem ?
            collection_inputs(source, collection(handle)) :
            collection_inputs(_typed_collection_path(handle))
        record = ItemRecord(;
            source_item_id=source_item_id_value,
            source_item_path=source_item_path_value,
            source_item_timestamp=source_item_timestamp(source_item),
            id=item_id,
            item_label=isempty(label) ? source_item_label(source_item) : label,
            kind=kind(handle),
            collection_key=nothing,
            metadata=record_metadata,
            item_fingerprint=fingerprint(handle),
        )
        records[index] = record
        collection_paths[index] = normalized_collection
        interpreted_items[index] = handle isa RegisteredDataItem ?
            RegisteredDataItem(record, item_data(handle), collection(handle)) : handle
    end
    finish_source_profile!(
        project, source_item_id_value, source_item_label_value, source_item_path_value,
        source_kind, item_count,
        (time_ns() - source_started) / 1e9, Set([Base.Threads.threadid()]))
    return SourceItemInterpretation(records, collection_paths, interpreted_items, ItemFailure[])
end

"""
Interpret one physical file into the data items produced by project code.

Typed items retain their concrete collection values. Registered items retain the string path
returned by their registration callback.
"""
function items_for_file(
    project::Project,
    filepath::AbstractString;
    meta::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}}=nothing,
)::Vector{AbstractDataItem}
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
    collections = CollectionIndex(source_id(source))
    records = ItemRecord[
        ItemRecord(record; collection_key=resolve_collection_path!(collections, path))
        for (record, path) in zip(
            interpretation.records,
            interpretation.collection_paths,
        )
    ]
    return AbstractDataItem[
        item isa RegisteredDataItem ?
            RegisteredDataItem(
                effective_record(collections, record),
                item_data(item),
                collection(item),
            ) :
            item
        for (record, item) in zip(records, interpretation.interpreted_items)
    ]
end

"""
Rewrite one collection's members through the kind's registered collection `process`.

Members are grouped by kind; a kind without a registered `process` passes its members through
unchanged. The callback returns one output per input; the adapter validates ids and count.
"""
function _process_collection(
    project::Project,
    ::AbstractDataSource,
    items::Vector{<:AbstractDataItem},
)::Vector{<:AbstractDataItem}
    rewritten = AbstractDataItem[]
    for positions in _group_positions_by_kind(items)
        group = items[positions]
        first(group) isa RegisteredDataItem || (append!(rewritten, group); continue)
        recipe = get(project.collections, first(group).registration, nothing)
        if recipe === nothing || recipe.process === nothing
            append!(rewritten, group)
            continue
        end
        output_data = recipe.process(item_data.(group), metadata.(group))
        output_data isa AbstractVector || error(
            "collection process for kind $(kind(first(group))) must return a vector; " *
            "got $(typeof(output_data))",
        )
        outputs = AbstractDataItem[
            _with_data(input::RegisteredDataItem, data)
            for (input, data) in zip(group, output_data)
        ]
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
    items::Vector{<:AbstractDataItem},
)::Dict{Symbol,Any}
    merged = Dict{Symbol,Any}()
    for positions in _group_positions_by_kind(items)
        group = items[positions]
        first(group) isa RegisteredDataItem || continue
        recipe = get(project.collections, first(group).registration, nothing)
        (recipe === nothing || recipe.analyze === nothing) && continue
        merge!(merged, metadata_dict(recipe.analyze(
            item_data.(group), metadata.(group))))
    end
    return merged
end

"""Group item positions by kind, preserving first-seen kind order."""
function _group_positions_by_kind(items::Vector{<:AbstractDataItem})::Vector{Vector{Int}}
    groups = Dict{Tuple{Bool,Symbol},Vector{Int}}()
    order = Tuple{Bool,Symbol}[]
    for (position, item) in pairs(items)
        key = (item isa RegisteredDataItem, kind(item))
        haskey(groups, key) || push!(order, key)
        push!(get!(() -> Int[], groups, key), position)
    end
    return Vector{Int}[groups[key] for key in order]
end

function _analyze_item(
    project::Project,
    ::AbstractDataSource,
    item::AbstractDataItem,
)::Dict{Symbol,Any}
    if !(item isa RegisteredDataItem)
        return metadata_dict(analyze(item))
    end
    recipe = _recipe(project, item.registration)
    recipe === nothing && error("Missing registration $(item.registration)")
    recipe.analyze === nothing && return Dict{Symbol,Any}()
    return metadata_dict(recipe.analyze(item_data(item), metadata(item)))
end

function _has_collection_process(project::Project, item_kind::Symbol)::Bool
    recipe = get(project.collections, item_kind, nothing)
    return recipe !== nothing && recipe.process !== nothing
end

function _has_collection_analysis(project::Project, item_kind::Symbol)::Bool
    recipe = get(project.collections, item_kind, nothing)
    return recipe !== nothing && recipe.analyze !== nothing
end
