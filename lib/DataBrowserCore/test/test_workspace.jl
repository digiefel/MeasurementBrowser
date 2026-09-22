using DataBrowserAPI: AbstractProject, AbstractDataSource, AbstractDataSourceItem, AbstractDataItem, id, item_data, metadata
using Test, DataBrowserAPI, DataBrowserCore
import DataBrowserCore.Workspace as W

mutable struct EngineProject <: AbstractProject
    name::String
    reads::Int
end
DataBrowserAPI.project_name(p::EngineProject) = p.name
struct EngineSource <: AbstractDataSource
    name::String
    version::Int
    fail::Bool
end
struct EngineSourceItem <: AbstractDataSourceItem
    version::Int
    fail::Bool
end
struct EngineItem <: AbstractDataItem
    key::String
    data::Any
    meta::Dict
end
DataBrowserAPI.source_id(s::EngineSource) = s.name
DataBrowserAPI.source_label(s::EngineSource) = s.name
DataBrowserAPI.source_items(s::EngineSource; kwargs...) =
    s.version == 0 ? EngineSourceItem[] : [EngineSourceItem(s.version, s.fail)]
DataBrowserAPI.id(::EngineSourceItem) = "source"
DataBrowserAPI.label(::EngineSourceItem) = "source"
DataBrowserAPI.fingerprint(s::EngineSourceItem) = (s.version, s.fail)
DataBrowserAPI.id(i::EngineItem) = i.key
DataBrowserAPI.item_data(i::EngineItem) = i.data
DataBrowserAPI.metadata(i::EngineItem) = i.meta
DataBrowserAPI.reconstruct(::Type{EngineItem}, key::AbstractString, data, meta::Dict) = EngineItem(key, data, meta)
function DataBrowserAPI.read(p::EngineProject, ::EngineSource, s::EngineSourceItem)
    p.reads += 1
    s.fail && error("fixture read failed")
    return s.version
end
DataBrowserAPI.entries(::EngineProject, ::EngineSourceItem, version) =
    [EngineItem(string(n), (x=[version],), Dict(:input => version, :total => -1)) for n in 1:version]
DataBrowserAPI.process(::EngineProject, i::EngineItem) = EngineItem(i.key, (x=i.data.x .* 2,), merge(i.meta, Dict(:total => -2)))
DataBrowserAPI.analyze(::EngineProject, i::EngineItem) = Dict(:total => sum(i.data.x))

function settled(ws)
    W.wait_workspace_idle!(ws; timeout=20)
    @test !W.workspace_status(ws).busy
end

@testset "reconstruction, reuse and source replacement (cache=$cache)" for cache in (true, false)
    mktempdir() do depot
        pushfirst!(DEPOT_PATH, depot)
        try
            project = EngineProject("engine", 0)
            ws = W.open_workspace(project, EngineSource("source", 2, false); cache, background_processing=true)
            try
                settled(ws)
                items = W.materialize_items(ws, sort(W.query_items(ws)))
                @test id.(items) == ["1", "2"]
                @test all(i -> item_data(i).x == [4] && metadata(i)[:total] == 4, items)
            finally
                W.close_workspace!(ws)
            end
            ws = W.open_workspace(project, EngineSource("source", 2, false); cache, background_processing=true)
            try
                settled(ws)
                @test metadata(only(W.materialize_items(ws, ["1"])))[:total] == 4
                @test project.reads == (cache ? 1 : 2)
                W.modify_workspace!(ws; source=EngineSource("source", 1, false))
                settled(ws)
                @test W.query_items(ws) == ["1"]
                @test item_data(only(W.materialize_items(ws, ["1"]))).x == [2]
                W.modify_workspace!(ws; source=EngineSource("source", 0, false))
                settled(ws)
                @test isempty(W.query_items(ws))
                @test W.workspace_status(ws).level !== :error
                @test isempty(W.workspace_status(ws).errors)
            finally
                W.close_workspace!(ws)
            end
        finally
            popfirst!(DEPOT_PATH)
        end
    end
end

@testset "failure reuse and recovery (cache=$cache)" for cache in (true, false)
    mktempdir() do depot
        pushfirst!(DEPOT_PATH, depot)
        try
            project = EngineProject("failure", 0)
            for _ in 1:2
                ws = W.open_workspace(project, EngineSource("source", 1, true); cache)
                try
                    settled(ws)
                    @test !isempty(W.workspace_status(ws).errors)
                finally
                    W.close_workspace!(ws)
                end
            end
            @test project.reads == (cache ? 1 : 2)
            ws = W.open_workspace(project, EngineSource("source", 1, false); cache)
            try
                settled(ws)
                @test isempty(W.workspace_status(ws).errors)
                @test item_data(only(W.materialize_items(ws, ["1"]))).x == [2]
            finally
                W.close_workspace!(ws)
            end
        finally
            popfirst!(DEPOT_PATH)
        end
    end
end
