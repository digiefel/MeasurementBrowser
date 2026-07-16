using DataBrowser
using Test
using DBInterface
using DuckDB

const WORK = DataBrowserCore.Workspace
const CACHE = DataBrowserCache

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

@testset "workspace status does not block on publication" begin
    mktempdir() do dir
        project = DataBrowser.define_project("StatusSnapshot_$(basename(dir))")
        workspace = DataBrowser.open_workspace(project, DataBrowser.DirectorySource(dir); cache=false)
        try
            DataBrowser.wait_workspace_idle!(workspace)
            lock(workspace.publish_lock) do
                workspace.index.analysis_errors["item"] = "failed"
            end
            DataBrowserCore.Workspace.refresh_status!(workspace)
            baseline = workspace.status

            locked = Channel{Nothing}(1)
            release = Channel{Nothing}(1)
            holder = Threads.@spawn lock(workspace.publish_lock) do
                put!(locked, nothing)
                take!(release)
            end
            take!(locked)
            started = time()
            try
                snapshot = DataBrowser.workspace_status(workspace)
                @test time() - started < 1
                @test snapshot.errors == baseline.errors
            finally
                put!(release, nothing)
                wait(holder)
            end
        finally
            DataBrowser.close_workspace!(workspace)
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
                Int64(1),
                Dict(:mean => 1.0, :count => 2),
            )
            @test read(cache.analyzed_collection_metadata)[Int64(1)][:mean] == 1.0
            @test haskey(
                read(cache.collection_result_states),
                (Int8(CACHE.COLLECTION_ANALYSIS_RESULT), Int64(1)),
            )

            CACHE.delete_collection_metadata!(cache, Int64[1])

            @test !haskey(read(cache.analyzed_collection_metadata), Int64(1))
            @test !haskey(
                read(cache.collection_result_states),
                (Int8(CACHE.COLLECTION_ANALYSIS_RESULT), Int64(1)),
            )
        finally
            CACHE.close_cache_db!(cache)
            rm(dirname(identity.cache_path); force=true, recursive=true)
        end
    end
end
