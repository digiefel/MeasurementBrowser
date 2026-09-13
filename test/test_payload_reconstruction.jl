using DataBrowser
using DataBrowserAPI
using DataBrowserAPI: ITEM_PROCESS, COLLECTION_PROCESS
using DataBrowserCache
using Test

const PAYLOAD_API = DataBrowserAPI

mutable struct PayloadCounters
    reads::Threads.Atomic{Int}
    processes::Threads.Atomic{Int}
    collection_processes::Threads.Atomic{Int}
    self_reconstructs_with_analysis::Threads.Atomic{Int}
end

PayloadCounters() = PayloadCounters(
    Threads.Atomic{Int}(0), Threads.Atomic{Int}(0), Threads.Atomic{Int}(0),
    Threads.Atomic{Int}(0))

const PAYLOAD_COUNTERS = Dict{String,PayloadCounters}()

struct PayloadContractProject <: PAYLOAD_API.AbstractProject
    name::String
end

struct PayloadContractSource <: PAYLOAD_API.AbstractDataSource
    name::String
    mode::Symbol
end

struct PayloadContractSourceItem <: PAYLOAD_API.AbstractDataSourceItem
    source_name::String
    mode::Symbol
end

struct PayloadContractCollection <: PAYLOAD_API.AbstractCollection
    name::String
end

struct PayloadContractItem{K,D} <: PAYLOAD_API.AbstractDataItem
    item_id::String
    payload::D
    item_metadata::Dict{Symbol,Any}
end

PayloadContractItem{K}(
    item_id::AbstractString,
    payload::D,
    item_metadata::Dict{Symbol,Any},
) where {K,D} = PayloadContractItem{K,D}(String(item_id), payload, item_metadata)

struct PayloadSelfItem <: PAYLOAD_API.AbstractDataItem
    item_id::String
    value::Int
end

const PAYLOAD_SELF_ITEMS = Dict{String,PayloadSelfItem}()

struct SourceFallbackItem <: PAYLOAD_API.AbstractDataItem
    item_id::String
    payload::NamedTuple{(:x,),Tuple{Vector{Int}}}
    item_metadata::Dict{Symbol,Any}
end

PAYLOAD_API.project_name(project::PayloadContractProject) = project.name
PAYLOAD_API.source_id(source::PayloadContractSource) = source.name
PAYLOAD_API.source_label(source::PayloadContractSource) = source.name
PAYLOAD_API.source_items(source::PayloadContractSource; kwargs...) =
    [PayloadContractSourceItem(source.name, source.mode)]
PAYLOAD_API.id(item::PayloadContractSourceItem) = "$(item.source_name)-source"
PAYLOAD_API.label(item::PayloadContractSourceItem) = String(item.mode)
PAYLOAD_API.fingerprint(item::PayloadContractSourceItem) = item.mode

PAYLOAD_API.id(collection::PayloadContractCollection) = collection.name
PAYLOAD_API.label(collection::PayloadContractCollection) = collection.name
PAYLOAD_API.reconstruct(
    ::Type{PayloadContractCollection}, identity::AbstractString, ::Dict,
) = PayloadContractCollection(String(identity))

PAYLOAD_API.id(item::PayloadContractItem) = item.item_id
PAYLOAD_API.item_data(item::PayloadContractItem) = item.payload
PAYLOAD_API.metadata(item::PayloadContractItem) = item.item_metadata
PAYLOAD_API.collection(::PayloadContractItem) =
    PAYLOAD_API.AbstractCollection[PayloadContractCollection("members")]
function PAYLOAD_API.reconstruct(
    ::Type{<:PayloadContractItem{K}},
    item_id::AbstractString,
    payload,
    metadata::Dict,
) where {K}
    return PayloadContractItem{K}(String(item_id), payload, Dict{Symbol,Any}(metadata))
end

PAYLOAD_API.id(item::PayloadSelfItem) = item.item_id
function PAYLOAD_API.reconstruct(
    ::Type{PayloadSelfItem},
    item_id::AbstractString,
    data::PayloadSelfItem,
    metadata::Dict,
)
    if get(metadata, :analysis_only, false)
        name = first(split(String(item_id), "-item"; limit=2))
        Threads.atomic_add!(PAYLOAD_COUNTERS[name].self_reconstructs_with_analysis, 1)
    end
    return data
end

PAYLOAD_API.id(item::SourceFallbackItem) = item.item_id
PAYLOAD_API.item_data(item::SourceFallbackItem) = item.payload
PAYLOAD_API.metadata(item::SourceFallbackItem) = item.item_metadata
PAYLOAD_API.reconstruct(
    ::Type{SourceFallbackItem}, ::AbstractString, data, ::Dict,
) = nothing

function PAYLOAD_API.read(
    project::PayloadContractProject,
    ::PayloadContractSource,
    source_item::PayloadContractSourceItem,
)
    Threads.atomic_add!(PAYLOAD_COUNTERS[project.name].reads, 1)
    return source_item.mode
end

function PAYLOAD_API.entries(
    ::PayloadContractProject,
    source_item::PayloadContractSourceItem,
    mode::Symbol,
)
    item_id = "$(source_item.source_name)-item"
    mode === :array && return [PayloadContractItem{:array}(
        item_id, [1, 2], Dict{Symbol,Any}(:entry => "array"))]
    mode === :table && return [PayloadContractItem{:table}(
        item_id, (x=[1, 2],), Dict{Symbol,Any}(:entry => "table"))]
    mode === :nothing && return [PayloadContractItem{:nothing}(
        item_id, nothing, Dict{Symbol,Any}(:entry => "nothing"))]
    mode === :self && return [PayloadSelfItem(item_id, 7)]
    mode === :fallback && return [SourceFallbackItem(
        item_id, (x=[1, 2],), Dict{Symbol,Any}(:entry => "fallback"))]
    error("unknown payload test mode $mode")
end

function PAYLOAD_API.process(
    project::PayloadContractProject,
    item::PayloadContractItem{K},
) where {K}
    Threads.atomic_add!(PAYLOAD_COUNTERS[project.name].processes, 1)
    haskey(item.item_metadata, :analysis_only) &&
        error("item process received downstream analysis metadata")
    payload = K === :array ? item.payload .* 2 :
        K === :table ? (x=item.payload.x .* 2,) : nothing
    metadata = merge(copy(item.item_metadata), Dict{Symbol,Any}(:processed => String(K)))
    return PayloadContractItem{K}(item.item_id, payload, metadata)
end

function PAYLOAD_API.process(project::PayloadContractProject, item::PayloadSelfItem)
    Threads.atomic_add!(PAYLOAD_COUNTERS[project.name].processes, 1)
    PAYLOAD_SELF_ITEMS[project.name] = item
    return item
end

function PAYLOAD_API.process(project::PayloadContractProject, item::SourceFallbackItem)
    Threads.atomic_add!(PAYLOAD_COUNTERS[project.name].processes, 1)
    haskey(item.item_metadata, :analysis_only) &&
        error("fallback process received downstream analysis metadata")
    return SourceFallbackItem(
        item.item_id,
        (x=item.payload.x .* 2,),
        merge(copy(item.item_metadata), Dict{Symbol,Any}(:processed => "fallback")),
    )
end

function PAYLOAD_API.analyze(
    ::PayloadContractProject,
    item::PayloadContractItem,
)
    @assert haskey(item.item_metadata, :processed)
    return Dict{Symbol,Any}(:analysis_only => true)
end

PAYLOAD_API.analyze(::PayloadContractProject, ::PayloadSelfItem) =
    Dict{Symbol,Any}(:analysis_only => true)
PAYLOAD_API.analyze(::PayloadContractProject, ::SourceFallbackItem) =
    Dict{Symbol,Any}(:analysis_only => true)

function PAYLOAD_API.process(
    project::PayloadContractProject,
    ::PayloadContractCollection,
    items::AbstractVector,
)
    all(item -> haskey(PAYLOAD_API.metadata(item), :processed), items) ||
        error("collection process did not receive item process metadata")
    all(item -> haskey(PAYLOAD_API.metadata(item), :analysis_only), items) ||
        error("collection process did not receive item analysis metadata")
    Threads.atomic_add!(PAYLOAD_COUNTERS[project.name].collection_processes, 1)
    return items
end

function _payload_value(item::PayloadContractItem{:table})
    return item.payload.x
end

_payload_value(item::PayloadContractItem) = item.payload

@testset "payload reconstruction is independent of cache storage" begin
    for (mode, disk_cache) in ((:array, false), (:nothing, false), (:table, false), (:table, true))
        name = "PayloadContract_$(mode)_$(time_ns())"
        counters = PayloadCounters()
        PAYLOAD_COUNTERS[name] = counters
        project = PayloadContractProject(name)
        source = PayloadContractSource(name, mode)
        item_id = "$(name)-item"
        workspace = DataBrowser.open_workspace(
            project, source; cache=disk_cache, background_processing=true)
        try
            DataBrowser.wait_workspace_idle!(workspace; timeout=30)
            @test isempty(DataBrowser.workspace_status(workspace).errors)
            item = only(DataBrowser.materialize_items(workspace, [item_id]))
            @test PAYLOAD_API.id(item) == item_id
            @test _payload_value(item) == (mode === :nothing ? nothing : [2, 4])
            @test PAYLOAD_API.metadata(item)[:processed] == String(mode)
            @test PAYLOAD_API.metadata(item)[:analysis_only] === true
            @test counters.reads[] == 1
            @test counters.processes[] == 1
            @test counters.collection_processes[] == 1
            record = workspace.index.items[item_id]
            @test DataBrowserCache.has_payload(
                workspace.cache.db, item_id; stage=ITEM_PROCESS)
            @test !DataBrowserCache.has_payload(
                workspace.cache.db, "missing-item"; stage=ITEM_PROCESS)
            hit = only(DataBrowserCache.read_payload(
                workspace.cache.db, [record]; stage=ITEM_PROCESS))
            @test hit isa Some
            @test something(hit) == PAYLOAD_API.item_data(item)
            @test only(DataBrowserCache.read_payload(
                workspace.cache.db, [record]; stage=COLLECTION_PROCESS)) === nothing

            if mode === :array && !disk_cache
                processes_before = counters.processes[]
                DataBrowserCache.clear_cached_result_state!(
                    workspace.cache.db, ITEM_PROCESS, item_id)
                @test DataBrowserCache.has_payload(
                    workspace.cache.db, item_id; stage=ITEM_PROCESS)
                @test DataBrowserCache.cached_result_state(
                    workspace.cache.db,
                    ITEM_PROCESS,
                    item_id,
                ) === nothing
                only(DataBrowser.materialize_items(workspace, [item_id]))
                @test counters.processes[] == processes_before + 1

                processes_before = counters.processes[]
                delete!(workspace.cache.db.processed_memory, item_id)
                @test !DataBrowserCache.has_payload(
                    workspace.cache.db, item_id; stage=ITEM_PROCESS)
                only(DataBrowser.materialize_items(workspace, [item_id]))
                @test counters.processes[] == processes_before + 1
            end
        finally
            DataBrowser.close_workspace!(workspace)
        end

        if disk_cache && mode === :table
            reopened = DataBrowser.open_workspace(
                project, source; cache=true, background_processing=true)
            try
                DataBrowser.wait_workspace_idle!(reopened; timeout=30)
                item = only(DataBrowser.materialize_items(reopened, [item_id]))
                @test _payload_value(item) == [2, 4]
                @test PAYLOAD_API.metadata(item)[:processed] == "table"
                @test PAYLOAD_API.metadata(item)[:analysis_only] === true
                @test counters.reads[] == 1
                @test counters.processes[] == 1
            finally
                DataBrowser.close_workspace!(reopened)
            end
        end
        delete!(PAYLOAD_COUNTERS, name)
    end
end

@testset "an item can be its own cached payload" begin
    name = "PayloadSelf_$(time_ns())"
    counters = PayloadCounters()
    PAYLOAD_COUNTERS[name] = counters
    project = PayloadContractProject(name)
    workspace = DataBrowser.open_workspace(
        project,
        PayloadContractSource(name, :self);
        cache=false,
        background_processing=true,
    )
    try
        DataBrowser.wait_workspace_idle!(workspace; timeout=30)
        item = only(DataBrowser.materialize_items(workspace, ["$(name)-item"]))
        @test item isa PayloadSelfItem
        @test item === PAYLOAD_SELF_ITEMS[name]
        @test item.value == 7
        @test counters.self_reconstructs_with_analysis[] >= 1
        @test counters.reads[] == 1
        @test counters.processes[] == 1
    finally
        DataBrowser.close_workspace!(workspace)
        delete!(PAYLOAD_COUNTERS, name)
        delete!(PAYLOAD_SELF_ITEMS, name)
    end
end

@testset "source fallback runs process once per reconstruction" begin
    for disk_cache in (false, true)
        name = "PayloadFallback_$(disk_cache)_$(time_ns())"
        counters = PayloadCounters()
        PAYLOAD_COUNTERS[name] = counters
        project = PayloadContractProject(name)
        source = PayloadContractSource(name, :fallback)
        item_id = "$(name)-item"
        workspace = DataBrowser.open_workspace(
            project, source; cache=disk_cache, background_processing=true)
        try
            DataBrowser.wait_workspace_idle!(workspace; timeout=30)
            @test isempty(DataBrowser.workspace_status(workspace).errors)
            processes_before = counters.processes[]
            reads_before = counters.reads[]
            item = only(DataBrowser.materialize_items(workspace, [item_id]))
            @test item.payload.x == [2, 4]
            @test counters.processes[] == processes_before + 1
            @test counters.reads[] == reads_before + 1
        finally
            DataBrowser.close_workspace!(workspace)
        end

        if disk_cache
            processes_before = counters.processes[]
            reads_before = counters.reads[]
            reopened = DataBrowser.open_workspace(
                project, source; cache=true, background_processing=true)
            try
                DataBrowser.wait_workspace_idle!(reopened; timeout=30)
                @test counters.processes[] == processes_before
                @test counters.reads[] == reads_before
                item = only(DataBrowser.materialize_items(reopened, [item_id]))
                @test item.payload.x == [2, 4]
                @test counters.processes[] == processes_before + 1
                @test counters.reads[] == reads_before + 1
                @test isempty(DataBrowser.workspace_status(reopened).errors)
            finally
                DataBrowser.close_workspace!(reopened)
            end
        end
        delete!(PAYLOAD_COUNTERS, name)
    end
end
