"""
Project.jl - engine methods for the registration-based project API.

Construction (`define_project`, `register_item!`, `register_collection_analysis!`) and the
payload-agnostic contracts live in `DataBrowserAPI`. This file keeps the callback-driving methods
that touch engine types (`ItemRecord`, `SourceFile`, `DataFrame`, …).
"""

using DataBrowserAPI
using DataBrowserAPI: AbstractCollection, AbstractDataItem, AbstractDataSource, AbstractDataSourceItem, AbstractProject, annotate_collection_path, default_collection_path, source_id, source_item_path, source_item_timestamp
using DataBrowserSources: DirectorySource, SourceFile, index_source_file
using DataFrames: DataFrame
import DataBrowserAPI:
    analyze,
    attach_record,
    detect_kind,
    entries,
    item_type,
    _has_collection_analysis,
    _has_collection_process,
    collection,
    id,
    item_data,
    kind,
    label,
    metadata,
    process,
    project_description,
    project_name,
    read
import DataBrowserAPI.ItemIndex:
    CollectionIndex,
    CollectionInput,
    ItemFailure,
    ItemRecord,
    MetadataDict,
    NoMatch,
    RegisteredDataItem,
    RegisteredReadResult,
    collection_inputs,
    effective_record,
    named_collection_path,
    resolve_collection_path!,
    metadata_dict

function _with_data(item::RegisteredDataItem, data)::RegisteredDataItem
    return RegisteredDataItem(item, data)
end

# ---------------------------------------------------------------------------
# Recipe lookup helpers
# ---------------------------------------------------------------------------

_recipe(project::Project, kind::Symbol)::Union{Nothing,ItemRecipe} =
    (i = findfirst(r -> r.kind === kind, project.recipes); i === nothing ? nothing : project.recipes[i])

function _detect_recipe(
    project::Project,
    source_item::AbstractDataSourceItem,
)::Union{Nothing,ItemRecipe}
    for recipe in project.recipes
        recipe.detect(source_item)::Bool && return recipe
    end
    return nothing
end

"""
Mint one item's final id from its source item, kind, position, and optional sibling key.

This is the single identity rule for every interpreted item: registered and typed items both
supply at most a sibling key (`recipe.id` callback or `id(item)`), and the engine namespaces it
under the source item and kind. An absent or empty key falls back to the item's returned position.
"""
function _mint_id(
    source_item_id::AbstractString,
    kind::Symbol,
    position::Integer,
    supplied_key=nothing,
)::String
    suffix = supplied_key === nothing || supplied_key == "" ?
        string(position) : string(supplied_key)
    return "$(source_item_id)#$(kind):$(suffix)"
end

"""Split the documented `(data=..., metadata=Dict(...))` result form."""
function _data_and_metadata(value)::Tuple{Any,MetadataDict}
    if value isa NamedTuple && keys(value) == (:data, :metadata)
        return value.data, metadata_dict(value.metadata)
    end
    return value, MetadataDict()
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

"""
Resolve where one interpreted item is placed, with the source in hand.

An item declaring no path of its own takes the source's default placement; every path is then
offered to the source, which may attach its own metadata to each level. This is the one point at
which placement sees the source, so the stages before it stay pure functions of values.
"""
function _placed_collection_path(
    source::AbstractDataSource,
    source_item::AbstractDataSourceItem,
    item::AbstractDataItem,
)::Vector{AbstractCollection}
    path = _collection_path(item)
    isempty(path) && (path = default_collection_path(source, source_item))
    isempty(path) && return path
    return AbstractCollection[segment for segment in annotate_collection_path(source, path)]
end

"""Validate one item's `collection(item)` contract result and return its concrete segments."""
function _collection_path(item::AbstractDataItem)::Vector{AbstractCollection}
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

"""
Adapt one registration callback result into a `RegisteredDataItem`.

This translates the registration dialect (entry tuples, callback-supplied keys, labels, and
collection strings) into the item contract; it mints nothing and never sees the source. A recipe
without a `collection` callback leaves the path empty, and the engine applies the source's default
placement afterwards — the same rule typed items get. The carrier's `id` holds only the
callback-supplied sibling key (or `""`); final ids are minted once, in `interpret_source_item`,
through the same path typed items take.
"""
function _registered_item(
    recipe::ItemRecipe,
    value,
    inherited_metadata::MetadataDict,
)::RegisteredDataItem
    data, entry_metadata = _data_and_metadata(value)
    local_metadata = merge(copy(inherited_metadata), entry_metadata)
    collection_path = recipe.collection === nothing ? AbstractCollection[] :
        named_collection_path(
            _registration_collection_path(recipe.collection(data, local_metadata)))
    supplied_key = recipe.id === nothing ? nothing : recipe.id(data, local_metadata)
    label = recipe.label === nothing ? "" : String(recipe.label(data, local_metadata))
    return RegisteredDataItem(
        supplied_key === nothing ? "" : string(supplied_key),
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

"""Every registered kind is carried by the same package-owned type."""
item_type(project::Project, kind::Symbol)::Union{Nothing,Type} =
    _recipe(project, kind) === nothing ? nothing : RegisteredDataItem

project_name(project::Project)::String = project.name
project_description(project::Project)::String = project.description

function detect_kind(project::Project, filename::String)::Symbol
    recipe = _detect_recipe(project, index_source_file(filename))
    return recipe === nothing ? :unknown : recipe.kind
end

function _processed_item(
    recipe::ItemRecipe,
    item::RegisteredDataItem,
)::RegisteredDataItem
    recipe.process === nothing && return item
    return _with_data(item, recipe.process(item.data, item.metadata))
end

"""Process one registered item through its `process` callback."""
function process(project::Project, item::RegisteredDataItem)::AbstractDataItem
    recipe = _recipe(project, item.registration)
    recipe === nothing && error("Missing registration $(item.registration)")
    return _processed_item(recipe, item)
end

"""
Run the matching registration's `read` callback, tagging its result with the registration name.

Detection happens here rather than in a pipeline stage of its own: `read` runs for every changed
source item anyway, and `detect` costs filename-predicate time. A source item matching no
registration yields `NoMatch`, whose `entries` produces no items.
"""
function read(project::Project, ::AbstractDataSource, item::AbstractDataSourceItem)
    recipe = @timed_dbg _detect_recipe(project, item)
    recipe === nothing && return NoMatch()
    return @timed_dbg "read" _invoke_read(recipe, item)
end

"""Function barrier: the concrete recipe specializes the callback call and result construction."""
function _invoke_read(recipe::ItemRecipe, item::AbstractDataSourceItem)
    data, read_metadata = _data_and_metadata(recipe.read(item))
    return RegisteredReadResult{recipe.kind}(data, read_metadata)
end

entries(::Project, ::AbstractDataSourceItem, ::NoMatch)::Vector{AbstractDataItem} =
    AbstractDataItem[]

"""
Expand one registration's read result into its data items.

The source item is still in hand because identity is not payload: collection defaults and inherited
metadata derive from it, while the loaded value stays purely the expensive data.
"""
function entries(
    project::Project,
    source_item::AbstractDataSourceItem,
    loaded::RegisteredReadResult{K},
)::Vector{AbstractDataItem} where {K}
    recipe = _recipe(project, K)
    recipe === nothing && error("Missing registration $K")
    inherited_metadata = merge(metadata_dict(metadata(source_item)), loaded.metadata)
    return @timed_dbg "entries" begin
        values = recipe.entries === nothing ? Any[loaded.data] :
            recipe.entries(loaded.data, inherited_metadata)
        values isa AbstractVector || error(
            "entries callback for registration $K must return a vector; got $(typeof(values))",
        )
        AbstractDataItem[
            _registered_item(recipe, value, inherited_metadata) for value in values
        ]
    end
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
    # The source item's display label, resolved once here so the UI and cache never rerun
    # project label code.
    source_item_label::String
end

"""
Interpret every logical data item produced by one source item.

The source is touched only by `read`; `entries` expands its result into items without going back to
the origin. Collection placement is applied here, where the source is legitimately in hand: an item
declaring no path of its own takes the source's default, and every path is offered to the source to
annotate. Concrete typed items remain unchanged; registered data uses a private carrier. Processing
and analysis belong to the workspace work graph. `source_item_key` is the workspace-minted surrogate
stamped on every record; transient interpretations outside a workspace (such as `items_for_file`)
leave it at 0.
"""
function interpret_source_item(
    project::AbstractProject,
    source::AbstractDataSource,
    source_item::AbstractDataSourceItem;
    source_item_key::Int64=Int64(0),
)::SourceItemInterpretation
    source_item_id_value = id(source_item)
    source_item_path_value = source_item_path(source_item)
    source_item_label_value = if source_item_path_value !== nothing &&
                                 isabspath(source_id(source))
        relpath(source_item_path_value, source_id(source))
    else
        label(source_item)
    end
    loaded = read(project, source, source_item)
    handles = entries(project, source_item, loaded)
    handles isa AbstractVector || error(
        "entries(::$(typeof(project)), ::$(typeof(source_item)), ::$(typeof(loaded))) must " *
        "return a vector of items; got $(typeof(handles))",
    )
    item_count = length(handles)
    records = Vector{ItemRecord}(undef, item_count)
    collection_paths = Vector{Vector{CollectionInput}}(undef, item_count)
    interpreted_items = Vector{AbstractDataItem}(undef, item_count)
    source_metadata = metadata_dict(metadata(source_item))
    for (index, handle) in pairs(handles)
        record_metadata = merge(copy(source_metadata), metadata_dict(metadata(handle)))
        item_label_value = label(handle)
        record = ItemRecord(;
            source_item_key,
            source_item_path=source_item_path_value,
            source_item_timestamp=source_item_timestamp(source_item),
            id=_mint_id(source_item_id_value, kind(handle), index, id(handle)),
            label=isempty(item_label_value) ? label(source_item) : item_label_value,
            kind=kind(handle),
            collection_key=nothing,
            metadata=record_metadata,
        )
        records[index] = record
        collection_paths[index] = collection_inputs(
            _placed_collection_path(source, source_item, handle))
        interpreted_items[index] = attach_record(handle, record)
    end
    return SourceItemInterpretation(
        records, collection_paths, interpreted_items, ItemFailure[],
        String(source_item_label_value))
end

"""
Interpret one physical file into the data items produced by project code.

Every item answers `collection(item)` with its concrete collection segments; registered items
carry the normalized segments adapted from their registration callback's string path.
"""
function items_for_file(
    project::AbstractProject,
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
        attach_record(item, effective_record(collections, record))
        for (record, item) in zip(records, interpretation.interpreted_items)
    ]
end

"""
Rewrite one collection's members through each kind's registered collection `process`.

Members are grouped by kind; a kind without a registered `process` falls through to the typed
`process(collection, items)`, whose default passes them along unchanged. The callback returns one
output per input; the adapter validates ids and count.
"""
function process(
    project::Project,
    collection_value::AbstractCollection,
    items::AbstractVector,
)::Vector{<:AbstractDataItem}
    rewritten = AbstractDataItem[]
    for positions in _group_positions_by_kind(items)
        group = items[positions]
        recipe = first(group) isa RegisteredDataItem ?
            get(project.collections, first(group).registration, nothing) : nothing
        if recipe === nothing || recipe.process === nothing
            append!(rewritten, process(collection_value, group))
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
function analyze(
    project::Project,
    collection_value::AbstractCollection,
    items::AbstractVector,
)::Dict{Symbol,Any}
    merged = Dict{Symbol,Any}()
    for positions in _group_positions_by_kind(items)
        group = items[positions]
        recipe = first(group) isa RegisteredDataItem ?
            get(project.collections, first(group).registration, nothing) : nothing
        if recipe === nothing || recipe.analyze === nothing
            merge!(merged, metadata_dict(analyze(collection_value, group)))
            continue
        end
        merge!(merged, metadata_dict(recipe.analyze(
            item_data.(group), metadata.(group))))
    end
    return merged
end

"""Group item positions by kind, preserving first-seen kind order."""
function _group_positions_by_kind(items::AbstractVector)::Vector{Vector{Int}}
    groups = Dict{Tuple{Bool,Symbol},Vector{Int}}()
    order = Tuple{Bool,Symbol}[]
    for (position, item) in pairs(items)
        key = (item isa RegisteredDataItem, kind(item))
        haskey(groups, key) || push!(order, key)
        push!(get!(() -> Int[], groups, key), position)
    end
    return Vector{Int}[groups[key] for key in order]
end

"""Analyze one registered item through its `analyze` callback."""
function analyze(project::Project, item::RegisteredDataItem)::Dict{Symbol,Any}
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
