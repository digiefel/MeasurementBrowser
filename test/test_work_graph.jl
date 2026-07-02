using MeasurementBrowser
using Test

const WORK = MeasurementBrowser.Workspace
const CACHE = MeasurementBrowser.Cache

@testset "work graph dependencies" begin
    graph = WORK.WorkDependencyGraph()
    process_key = WORK.WorkKey(WORK.PROCESS_ITEM, "item-1")
    item_a = WORK.WorkKey(WORK.ITEM_STATS, "item-1")
    item_b = WORK.WorkKey(WORK.ITEM_STATS, "item-2")
    collection = WORK.WorkKey(WORK.COLLECTION_STATS, "device-A")

    lock(graph.lock) do
        WORK.replace_work_node!(
            graph,
            WORK.WorkNode(process_key, UInt64(1), :failed, 0, WORK.WorkKey[], Channel{Any}[], 0),
        )
        blocked_item_stats = WORK.WorkNode(
            item_a, UInt64(1), :waiting, 1, WORK.WorkKey[process_key], Channel{Any}[], 0)
        WORK.replace_work_node!(graph, blocked_item_stats)
        @test !WORK.dependencies_ready(graph, blocked_item_stats)

        graph.nodes[process_key].state = :ready
        @test WORK.dependencies_ready(graph, blocked_item_stats)

        graph.nodes[item_a].state = :ready
        WORK.replace_work_node!(
            graph,
            WORK.WorkNode(item_b, UInt64(1), :missing, 0, WORK.WorkKey[], Channel{Any}[], 0),
        )
        collection_node = WORK.WorkNode(
            collection, UInt64(1), :waiting, 1, WORK.WorkKey[item_a, item_b],
            Channel{Any}[], 0)
        WORK.replace_work_node!(graph, collection_node)

        WORK.wake_ready_dependents!(graph, item_a)
        @test collection_node.state === :waiting
        @test isempty(graph.queue)

        graph.nodes[item_b].state = :failed
        WORK.wake_ready_dependents!(graph, item_b)
        @test collection_node.state === :queued
        @test graph.queue == [(collection, UInt64(1))]
    end
end

@testset "semantic collection-stat cache delete" begin
    mktempdir() do dir
        source = MeasurementBrowser.DirectorySource(dir)
        identity = CACHE.project_cache_identity("WorkGraphDelete_$(basename(dir))", source)
        cache = CACHE.open_cache_db(identity)
        try
            CACHE.store_collection_stats!(cache, "device-A", Dict(:mean => 1.0, :count => 2))
            CACHE.delete_collection_stats!(cache, ["device-A"])

            @test !haskey(
                read(cache.metadata),
                (Int8(CACHE.SCOPE_NODE_STATS), "device-A", "mean"),
            )
            @test !haskey(
                read(cache.metadata),
                (Int8(CACHE.SCOPE_NODE_STATS), "device-A", "count"),
            )
            @test !haskey(
                read(cache.result_states),
                (Int8(CACHE.COLLECTION_STATS_RESULT), "device-A"),
            )
        finally
            CACHE.close_cache_db!(cache)
            rm(dirname(identity.cache_path); force=true, recursive=true)
        end
    end
end
