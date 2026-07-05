using MeasurementBrowser
using Test
using DBInterface
using DuckDB

const WORK = MeasurementBrowser.Workspace
const CACHE = MeasurementBrowser.Cache

@testset "work graph dependencies" begin
    graph = WORK.WorkDependencyGraph()
    process_key = WORK.WorkKey(WORK.ITEM_PROCESS, "item-1")
    item_a = WORK.WorkKey(WORK.ITEM_ANALYZE, "item-1")
    item_b = WORK.WorkKey(WORK.ITEM_ANALYZE, "item-2")
    collection = WORK.WorkKey(WORK.COLLECTION_ANALYZE, "device-A")

    lock(graph.lock) do
        WORK.replace_work_node!(
            graph,
            WORK.WorkNode(process_key, UInt64(1), :failed, 0, WORK.WorkKey[], Channel{Any}[], 0),
        )
        blocked_item_analyze = WORK.WorkNode(
            item_a, UInt64(1), :waiting, 1, WORK.WorkKey[process_key], Channel{Any}[], 0)
        WORK.replace_work_node!(graph, blocked_item_analyze)
        @test !WORK.dependencies_ready(graph, blocked_item_analyze)

        graph.nodes[process_key].state = :ready
        @test WORK.dependencies_ready(graph, blocked_item_analyze)

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
        @test graph.queue == Dict(1 => [(collection, UInt64(1))])
        @test WORK.pop_queued_node!(graph) === collection_node
        @test collection_node.state === :running
        @test WORK.pop_queued_node!(graph) === nothing
    end
end

@testset "stale cache schema requires explicit rebuild" begin
    mktempdir() do dir
        source = MeasurementBrowser.DirectorySource(dir)
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
        project = MeasurementBrowser.define_project(project_name)
        source = MeasurementBrowser.DirectorySource(dir)
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

        workspace = MeasurementBrowser.open_workspace(project, source)
        try
            @test workspace.cache.db isa CACHE.MemoryCacheDB
            @test workspace.cache.disk_error isa CACHE.ProjectCacheSchemaError
            @test !workspace.background_processing
        finally
            MeasurementBrowser.close_workspace!(workspace)
            rm(dirname(identity.cache_path); force=true, recursive=true)
        end
    end
end

@testset "semantic collection-metadata cache delete" begin
    mktempdir() do dir
        source = MeasurementBrowser.DirectorySource(dir)
        identity = CACHE.project_cache_identity("WorkGraphDelete_$(basename(dir))", source)
        cache = CACHE.open_cache_db(identity)
        try
            CACHE.store_collection_metadata!(
                cache,
                "device-A",
                Dict(:mean => 1.0, :count => 2),
            )
            @test read(cache.collection_metadata)["device-A"][:mean] == 1.0
            @test haskey(
                read(cache.result_states),
                (Int8(CACHE.COLLECTION_ANALYSIS_RESULT), "device-A"),
            )

            CACHE.delete_collection_metadata!(cache, ["device-A"])

            @test !haskey(read(cache.collection_metadata), "device-A")
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
