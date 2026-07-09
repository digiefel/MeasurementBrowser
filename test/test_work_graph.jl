using DataBrowser
using Test
using DBInterface
using DuckDB

const WORK = DataBrowserCore.Workspace
const CACHE = DataBrowserCore.Cache

@testset "work graph dependencies" begin
    graph = WORK.WorkDependencyGraph()
    process_key = WORK.WorkKey(WORK.ITEM_PROCESS, "item-1")
    item_a = WORK.WorkKey(WORK.ITEM_ANALYZE, "item-1")
    item_b = WORK.WorkKey(WORK.ITEM_ANALYZE, "item-2")
    collection = WORK.WorkKey(WORK.COLLECTION_ANALYZE, "device-A")

    lock(graph.lock) do
        process_node = WORK.WorkNode(
            process_key, UInt16(1), :running, 0, Set{WORK.WorkKey}(), UInt64(0),
            Channel{Any}[], 0)
        graph.nodes[process_key] = process_node

        blocked_item_analyze = WORK.WorkNode(
            item_a, UInt16(1), :waiting, 1, Set{WORK.WorkKey}(), UInt64(0),
            Channel{Any}[], 0)
        WORK.seed_node_dependencies!(graph, blocked_item_analyze, WORK.WorkKey[process_key])
        graph.nodes[item_a] = blocked_item_analyze
        @test blocked_item_analyze.pending == 1
        @test item_a in process_node.dependents
        @test !WORK.dependencies_ready(blocked_item_analyze)

        WORK.wake_ready_dependents!(graph, process_node)
        delete!(graph.nodes, process_key)
        @test blocked_item_analyze.pending == 0
        @test WORK.dependencies_ready(blocked_item_analyze)
        WORK.queue_ready_node!(graph, blocked_item_analyze)
        @test blocked_item_analyze.state === :queued
        blocked_item_analyze.state = :waiting
        blocked_item_analyze.pending = 0
        empty!(graph.queue)

        collection_node = WORK.WorkNode(
            collection, UInt16(1), :waiting, 1, Set{WORK.WorkKey}(), UInt64(0),
            Channel{Any}[], 0)
        WORK.seed_node_dependencies!(graph, collection_node, WORK.WorkKey[item_a, item_b])
        graph.nodes[collection] = collection_node
        @test collection_node.pending == 1

        WORK.wake_ready_dependents!(graph, blocked_item_analyze)
        delete!(graph.nodes, item_a)
        @test collection_node.pending == 0
        @test collection_node.state === :queued
        @test graph.queue == Dict(1 => [(collection, UInt16(1))])
        @test WORK.pop_queued_node!(graph) === collection_node
        @test collection_node.state === :running
        @test WORK.pop_queued_node!(graph) === nothing
    end
end

@testset "stale cache schema requires explicit rebuild" begin
    mktempdir() do dir
        source = DataBrowser.DirectorySource(dir)
        identity = CACHE.project_cache_identity("StaleSchema_$(basename(dir))", source)
        mkpath(dirname(identity.cache_path))
        db = DBInterface.connect(DuckDB.DB, identity.cache_path)
        connection = DBInterface.connect(db)
        try
            DBInterface.execute(connection, "CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT)")
            DBInterface.execute(connection, "INSERT INTO meta VALUES ('schema_version', '0')")
        finally
            DBInterface.close!(connection)
            DBInterface.close!(db)
        end

        @test_throws CACHE.ProjectCacheSchemaError CACHE.open_cache_db(identity)
        @test isfile(identity.cache_path)

        db = DBInterface.connect(DuckDB.DB, identity.cache_path)
        connection = DBInterface.connect(db)
        try
            tables = String[
                String(row.table_name)
                for row in DBInterface.execute(
                    connection,
                    "SELECT table_name FROM information_schema.tables " *
                    "WHERE table_schema = 'main' ORDER BY table_name",
                )
            ]
            @test tables == ["meta"]
        finally
            DBInterface.close!(connection)
            DBInterface.close!(db)
        end

        cache = CACHE.open_cache_db(identity; rebuild=true)
        try
            @test isfile(identity.cache_path)
        finally
            CACHE.close_cache_db!(cache)
            rm(dirname(identity.cache_path); force=true, recursive=true)
        end
    end
end

@testset "workspace falls back to memory cache on stale disk schema" begin
    mktempdir() do dir
        project_name = "StaleWorkspace_$(basename(dir))"
        project = DataBrowser.define_project(project_name)
        source = DataBrowser.DirectorySource(dir)
        identity = CACHE.project_cache_identity(project_name, source)
        mkpath(dirname(identity.cache_path))
        db = DBInterface.connect(DuckDB.DB, identity.cache_path)
        connection = DBInterface.connect(db)
        try
            DBInterface.execute(connection, "CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT)")
            DBInterface.execute(connection, "INSERT INTO meta VALUES ('schema_version', '0')")
        finally
            DBInterface.close!(connection)
            DBInterface.close!(db)
        end

        workspace = DataBrowser.open_workspace(project, source)
        try
            @test workspace.cache.db isa CACHE.MemoryCacheDB
            @test workspace.cache.disk_error isa CACHE.ProjectCacheSchemaError
            @test !workspace.background_processing
        finally
            DataBrowser.close_workspace!(workspace)
            rm(dirname(identity.cache_path); force=true, recursive=true)
        end
    end
end

@testset "semantic collection-metadata cache delete" begin
    mktempdir() do dir
        source = DataBrowser.DirectorySource(dir)
        identity = CACHE.project_cache_identity("WorkGraphDelete_$(basename(dir))", source)
        cache = CACHE.open_cache_db(identity)
        try
            CACHE.store_collection_metadata!(
                cache,
                "device-A",
                Dict(:mean => 1.0, :count => 2),
            )
            @test read(cache.analyzed_collection_metadata)["device-A"][:mean] == 1.0
            @test haskey(
                read(cache.result_states),
                (Int8(CACHE.COLLECTION_ANALYSIS_RESULT), "device-A"),
            )

            CACHE.delete_collection_metadata!(cache, ["device-A"])

            @test !haskey(read(cache.analyzed_collection_metadata), "device-A")
            @test !haskey(
                read(cache.result_states),
                (Int8(CACHE.COLLECTION_ANALYSIS_RESULT), "device-A"),
            )
        finally
            CACHE.close_cache_db!(cache)
            rm(dirname(identity.cache_path); force=true, recursive=true)
        end
    end
end

@testset "cache stage summary follows persisted stage rows" begin
    mktempdir() do dir
        write_test_source(joinpath(dir, "a.csv"))
        workspace = DataBrowser.open_workspace(
            TEST_PROJECT,
            DataBrowser.DirectorySource(dir);
            background_processing=true,
        )
        try
            DataBrowserCore.Workspace.wait_workspace_idle!(workspace)
            summary = CACHE.cache_stage_summary(workspace.cache.db)
            @test summary.cached_sources == 1
            @test summary.interpreted_items == 1
            @test summary.processed == 1
            @test summary.analyzed == 1

            source_item_id = only(keys(workspace.index.items_by_source))
            records = DataBrowserCore.Workspace.source_item_records(
                workspace.index, source_item_id)
            CACHE.delete_source_item!(workspace.cache.db, source_item_id, records)
            summary = CACHE.cache_stage_summary(workspace.cache.db)
            @test summary.cached_sources == 0
            @test summary.interpreted_items == 0
            @test summary.processed == 0
            @test summary.analyzed == 0
        finally
            DataBrowser.close_workspace!(workspace)
        end
    end
end
