"""
The dialect's implementation of the typed stage contract.

Every registration callback runs through these adapter methods over the same stage runners and work
keys as a typed project's methods. Nothing here is privileged: the engine calls `read`, `entries`,
`process`, and `analyze`, and this file answers them.
"""

"""Split the documented `(data=..., metadata=Dict(...))` result form."""
function _data_and_metadata(value)::Tuple{Any,MetadataDict}
    if value isa NamedTuple && keys(value) == (:data, :metadata)
        return value.data, metadata_dict(value.metadata)
    end
    return value, MetadataDict()
end

function _with_data(item::RegisteredDataItem, data)::RegisteredDataItem
    return RegisteredDataItem(item, data)
end


"""
Select the recipe registered under one name.

The one dynamic step in the dialect. Everything downstream takes the returned recipe by value, so
its concrete parameters specialize the rest of the call.
"""
_recipe(project::Project, kind::Symbol)::Union{Nothing,ItemRecipe} =
    (i = findfirst(r -> r.kind === kind, project.recipes); i === nothing ? nothing : project.recipes[i])

"""Select the first recipe whose `detect` accepts a source item; registration order decides ties."""
function _detect_recipe(
    project::Project,
    source_item::AbstractDataSourceItem,
)::Union{Nothing,ItemRecipe}
    for recipe in project.recipes
        _detects(recipe, source_item) && return recipe
    end
    return nothing
end

_detects(recipe::ItemRecipe, source_item::AbstractDataSourceItem)::Bool =
    recipe.detect(source_item)::Bool


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
Adapt one registration callback result into a `RegisteredDataItem`.

This translates the registration dialect (entry tuples, callback-supplied keys, labels, and
collection strings) into the item contract. A recipe without a `collection` callback leaves the path
empty, and the engine applies the source's default placement afterwards — the same rule typed items
get.

Identity is the dialect's job, not the engine's: `id(item)` is verbatim for every item type, so this
builds one from the source item, the registration, and either the `id` callback's value or the
entry's position. A `register_item!` user never supplies or sees it.
"""
function _registered_item(
    recipe::ItemRecipe{<:Any,<:Any,<:Any,<:Any,<:Any,Label,Collection,Id},
    value,
    inherited_metadata::MetadataDict,
    source_item_id::AbstractString,
    position::Integer,
)::RegisteredDataItem where {Label,Collection,Id}
    data, entry_metadata = _data_and_metadata(value)
    local_metadata = merge(copy(inherited_metadata), entry_metadata)
    collection_path = recipe.collection === nothing ? AbstractCollection[] :
        named_collection_path(
            _registration_collection_path(recipe.collection(data, local_metadata)))
    supplied_key = recipe.id === nothing ? nothing : recipe.id(data, local_metadata)
    key = supplied_key === nothing || supplied_key == "" ?
        string(position) : string(supplied_key)
    label = recipe.label === nothing ? "" : String(recipe.label(data, local_metadata))
    return RegisteredDataItem{recipe.kind}(
        "$(source_item_id)#$(recipe.kind):$(key)",
        label,
        collection_path,
        data,
        local_metadata,
    )
end


project_name(project::Project)::String = project.name
project_description(project::Project)::String = project.description

function _processed_item(
    recipe::ItemRecipe,
    item::RegisteredDataItem,
)::RegisteredDataItem
    recipe.process === nothing && return item
    return _with_data(item, recipe.process(item.data, item.metadata))
end

"""Process one registered item through its `process` callback."""
function process(project::Project, item::RegisteredDataItem{K})::AbstractDataItem where {K}
    recipe = _recipe(project, K)
    recipe === nothing && error("Missing registration $K")
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
            _registered_item(recipe, value, inherited_metadata, id(source_item), position)
            for (position, value) in pairs(values)
        ]
    end
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
            get(project.collections, label(typeof(first(group))), nothing) : nothing
        if recipe === nothing || recipe.process === nothing
            append!(rewritten, process(collection_value, group))
            continue
        end
        output_data = recipe.process(item_data.(group), metadata.(group))
        output_data isa AbstractVector || error(
            "collection process for kind $(label(typeof(first(group)))) must return a vector; " *
            "got $(typeof(output_data))",
        )
        outputs = AbstractDataItem[
            _with_data(input::RegisteredDataItem, data)
            for (input, data) in zip(group, output_data)
        ]
        length(outputs) == length(group) || error(
            "collection process for kind $(label(typeof(first(group)))) must return one item per input; " *
            "got $(length(outputs)) for $(length(group)) members",
        )
        input_ids = Set(id(item) for item in group)
        for output in outputs
            id(output) in input_ids || error(
                "collection process for kind $(label(typeof(first(group)))) returned unknown item id " *
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
            get(project.collections, label(typeof(first(group))), nothing) : nothing
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
        key = (item isa RegisteredDataItem, label(typeof(item)))
        haskey(groups, key) || push!(order, key)
        push!(get!(() -> Int[], groups, key), position)
    end
    return Vector{Int}[groups[key] for key in order]
end

"""Analyze one registered item through its `analyze` callback."""
function analyze(project::Project, item::RegisteredDataItem{K})::Dict{Symbol,Any} where {K}
    recipe = _recipe(project, K)
    recipe === nothing && error("Missing registration $K")
    recipe.analyze === nothing && return Dict{Symbol,Any}()
    return metadata_dict(recipe.analyze(item_data(item), metadata(item)))
end


