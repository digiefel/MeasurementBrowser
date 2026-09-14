using Test
using DataBrowser
import DataBrowserAPI, DataBrowserCore, DataBrowserCache
using DataBrowserAPI: SOURCE_READ, SOURCE_INTERPRET

const SOURCE_API = DataBrowserAPI
const SOURCE_WORK = DataBrowserCore.Workspace
const SOURCE_CACHE = DataBrowserCache

mutable struct SplitSourceProject <: SOURCE_API.AbstractProject
    name::String
    reads::Threads.Atomic{Int}
    interpretations::Threads.Atomic{Int}
    fail_read::Bool
    fail_interpret::Bool
    empty::Bool
    nothing_input::Bool
    blocked_stage::Symbol
    entered::Channel{Nothing}
    release::Channel{Nothing}
end
SplitSourceProject() = SplitSourceProject("SourceStages_$(time_ns())",
    Threads.Atomic{Int}(0), Threads.Atomic{Int}(0), false, false, false, false,
    :none, Channel{Nothing}(1), Channel{Nothing}(1))

struct SplitSource <: SOURCE_API.AbstractDataSource
    name::String
    version::Int
end
struct SplitSourceItem <: SOURCE_API.AbstractDataSourceItem
    version::Int
end
struct SplitItem <: SOURCE_API.AbstractDataItem
    data::NamedTuple{(:x,),Tuple{Vector{Int}}}
end
SOURCE_API.project_name(project::SplitSourceProject) = project.name
SOURCE_API.source_id(source::SplitSource) = source.name
SOURCE_API.source_label(source::SplitSource) = source.name
SOURCE_API.source_items(source::SplitSource; kwargs...) = [SplitSourceItem(source.version)]
SOURCE_API.id(::SplitSourceItem) = "source"
SOURCE_API.label(::SplitSourceItem) = "source"
SOURCE_API.fingerprint(item::SplitSourceItem) = item.version
SOURCE_API.id(::SplitItem) = "item"
SOURCE_API.label(::SplitItem) = "item"
SOURCE_API.item_data(item::SplitItem) = item.data
SOURCE_API.reconstruct(::Type{SplitItem}, ::AbstractString, data, ::Dict) = SplitItem(data)

function SOURCE_API.read(project::SplitSourceProject, ::SplitSource, item::SplitSourceItem)
    Threads.atomic_add!(project.reads, 1)
    if project.blocked_stage === :read && item.version == 1
        put!(project.entered, nothing)
        take!(project.release)
    end
    project.fail_read && error("read fixture failed")
    return project.nothing_input ? nothing : item.version
end
function SOURCE_API.entries(project::SplitSourceProject, item::SplitSourceItem, loaded)
    Threads.atomic_add!(project.interpretations, 1)
    if project.blocked_stage === :interpret && item.version == 1
        put!(project.entered, nothing)
        take!(project.release)
    end
    project.fail_interpret && error("interpret fixture failed")
    project.empty && return SplitItem[]
    return [SplitItem((x=[something(loaded, 0)],))]
end

function source_stage_state(workspace, stage)
    key = SOURCE_CACHE.source_item_key(workspace.cache.db, "source")
    return SOURCE_CACHE.cached_result_state(workspace.cache.db, stage, key)
end
function finish_source_work(workspace)
    DataBrowser.wait_workspace_idle!(workspace; timeout=20)
    @test isempty(workspace.work.nodes)
end
function reinterpret_source(workspace)
    key = SOURCE_CACHE.source_item_key(workspace.cache.db, "source")
    SOURCE_CACHE.clear_cached_result_state!(workspace.cache.db, SOURCE_INTERPRET, key)
    result = SOURCE_WORK.ensure_uptodate!(workspace, SOURCE_WORK.WorkKey(SOURCE_INTERPRET, key))
    @test result.failure === nothing
    finish_source_work(workspace)
end

@testset "source read is independently reusable" begin
    for disk in (false, true), nothing_input in (false, true)
        project = SplitSourceProject()
        project.nothing_input = nothing_input
        workspace = DataBrowser.open_workspace(project, SplitSource(project.name, 1); cache=disk)
        try
            finish_source_work(workspace)
            @test project.reads[] == project.interpretations[] == 1
            @test source_stage_state(workspace, SOURCE_READ).status == Int8(SOURCE_CACHE.RESULT_READY)
            @test source_stage_state(workspace, SOURCE_INTERPRET).status == Int8(SOURCE_CACHE.RESULT_READY)
            key = SOURCE_CACHE.source_item_key(workspace.cache.db, "source")
            @test SOURCE_CACHE.read_payload(workspace.cache.db, key) == Some(nothing_input ? nothing : 1)
            reinterpret_source(workspace)
            @test project.reads[] == 1
            @test project.interpretations[] == 2
            # Eviction is independent of completion. Interpretation must pull a missing read input.
            SOURCE_CACHE.clear!(workspace.cache.db.source_reads)
            reinterpret_source(workspace)
            @test project.reads[] == 2
            @test project.interpretations[] == 3
            @test isempty(DataBrowser.workspace_status(workspace).errors)
        finally
            DataBrowser.close_workspace!(workspace)
        end
    end
end

@testset "source failures retain their stage and fingerprint across reopen" begin
    for failed_stage in (:read, :interpret)
        project = SplitSourceProject()
        project.fail_read = failed_stage === :read
        project.fail_interpret = failed_stage === :interpret
        source = SplitSource(project.name, 1)
        for opening in 1:2
            workspace = DataBrowser.open_workspace(project, source)
            try
                finish_source_work(workspace)
                @test project.reads[] == 1
                @test project.interpretations[] == (failed_stage === :interpret ? 1 : 0)
                stage = failed_stage === :read ? SOURCE_READ : SOURCE_INTERPRET
                @test source_stage_state(workspace, stage).status == Int8(SOURCE_CACHE.RESULT_FAILED)
                @test isempty(DataBrowser.query_items(workspace))
                @test DataBrowser.workspace_status(workspace).counts.sources_pending == 0
                counts = SOURCE_CACHE.cache_stage_summary(workspace.cache.db)
                @test counts.failed_read == (failed_stage === :read ? 1 : 0)
                @test counts.failed_interpret == (failed_stage === :interpret ? 1 : 0)
                failed_stage === :read && @test source_stage_state(workspace, SOURCE_INTERPRET) === nothing
            finally
                DataBrowser.close_workspace!(workspace)
            end
        end
        project.fail_read = project.fail_interpret = false
        workspace = DataBrowser.open_workspace(project, SplitSource(project.name, 2))
        try
            finish_source_work(workspace)
            @test project.reads[] == 2
            @test DataBrowser.query_items(workspace) == ["item"]
            @test isempty(DataBrowser.workspace_status(workspace).errors)
        finally
            DataBrowser.close_workspace!(workspace)
        end
    end
end

@testset "interpretation failure can retry with its successful read" begin
    project = SplitSourceProject()
    project.fail_interpret = true
    workspace = DataBrowser.open_workspace(project, SplitSource(project.name, 1); cache=false)
    try
        finish_source_work(workspace)
        @test project.reads[] == project.interpretations[] == 1
        project.fail_interpret = false
        reinterpret_source(workspace)
        @test project.reads[] == 1
        @test project.interpretations[] == 2
        @test DataBrowser.query_items(workspace) == ["item"]
    finally
        DataBrowser.close_workspace!(workspace)
    end
end

@testset "empty interpretations and warm processed results need no source replay" begin
    for empty_result in (false, true)
        project = SplitSourceProject()
        project.empty = empty_result
        source = SplitSource(project.name, 1)
        for opening in 1:2
            workspace = DataBrowser.open_workspace(project, source; background_processing=true)
            try
                finish_source_work(workspace)
                @test project.reads[] == project.interpretations[] == 1
                if !empty_result
                    @test only(DataBrowser.materialize_items(workspace, ["item"])).data.x == [1]
                end
                if opening == 2
                    key = SOURCE_CACHE.source_item_key(workspace.cache.db, "source")
                    @test SOURCE_CACHE.read_payload(workspace.cache.db, key) === nothing
                end
            finally
                DataBrowser.close_workspace!(workspace)
            end
        end
    end
end

@testset "source changes supersede both source stages" begin
    for stage in (:read, :interpret)
        project = SplitSourceProject()
        project.blocked_stage = stage
        workspace = DataBrowser.open_workspace(project, SplitSource(project.name, 1); cache=false)
        try
            @test timedwait(() -> isready(project.entered), 20) === :ok
            changes = SOURCE_API.SourceChanges([SplitSourceItem(2)], String[])
            SOURCE_WORK.publish_source_event!(workspace, changes)
            put!(project.release, nothing)
            finish_source_work(workspace)
            @test DataBrowser.query_items(workspace) == ["item"]
            @test only(DataBrowser.materialize_items(workspace, ["item"])).data.x == [2]
            @test isempty(DataBrowser.workspace_status(workspace).errors)
        finally
            isready(project.release) || put!(project.release, nothing)
            DataBrowser.close_workspace!(workspace)
        end
    end
end

@testset "interpretation replacement retires old items on disk" begin
    project = SplitSourceProject()
    source = SplitSource(project.name, 1)
    workspace = DataBrowser.open_workspace(project, source; background_processing=true)
    try
        finish_source_work(workspace)
        @test DataBrowser.query_items(workspace) == ["item"]
        project.empty = true
        reinterpret_source(workspace)
        @test project.reads[] == 1
        @test isempty(DataBrowser.query_items(workspace))
    finally
        DataBrowser.close_workspace!(workspace)
    end
    workspace = DataBrowser.open_workspace(project, source)
    try
        finish_source_work(workspace)
        @test isempty(DataBrowser.query_items(workspace))
        @test project.reads[] == 1
    finally
        DataBrowser.close_workspace!(workspace)
    end
end

@testset "an unfinished source resumes on reopen" begin
    project = SplitSourceProject()
    source = SplitSource(project.name, 1)
    workspace = DataBrowser.open_workspace(project, source)
    try
        finish_source_work(workspace)
        key = SOURCE_CACHE.source_item_key(workspace.cache.db, "source")
        SOURCE_CACHE.delete_source_output!(workspace.cache.db, key, collect(values(workspace.index.items)))
    finally
        DataBrowser.close_workspace!(workspace)
    end
    workspace = DataBrowser.open_workspace(project, source)
    try
        finish_source_work(workspace)
        @test project.reads[] == project.interpretations[] == 2
        @test DataBrowser.query_items(workspace) == ["item"]
    finally
        DataBrowser.close_workspace!(workspace)
    end
end

@testset "memory cache weights values and preserves acquired references" begin
    cache = SOURCE_CACHE.MemoryStore{Int,Any}(; capacity=100, weight=Base.summarysize)
    payload = [1, 2, 3]
    append!(cache, 1, payload)
    acquired = SOURCE_CACHE.read_hit(cache, 1)
    append!(cache, 2, zeros(100))
    @test !haskey(cache, 1)
    @test something(acquired) === payload
    @test haskey(cache, 2) # One oversized entry is usable.
    append!(cache, 3, nothing)
    @test !haskey(cache, 2)
    @test SOURCE_CACHE.read_hit(cache, 3) === Some(nothing)
    append!(cache, 3, payload)
    push!(payload, 4)
    delete!(cache, 3)
    @test cache.size == 0 # Mutation did not change the recorded insertion weight.
    SOURCE_CACHE.close!(cache)
end
