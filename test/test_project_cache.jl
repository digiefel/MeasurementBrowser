using MeasurementBrowser
using MeasurementBrowser: ItemFailure
using Test

const ProjectCache = MeasurementBrowser.Cache

@testset "project HDF5 cache" begin
    mktempdir() do dir
        source = test_source(TEST_PROJECT, dir)
        @test basename(ProjectCache.project_cache_identity(
            "20260430_120001",
            source,
        ).cache_path) == "20260430_120001.h5"
    end

    mktempdir() do first_root
        mktempdir() do second_root
            @test ProjectCache.project_cache_id(test_source(TEST_PROJECT, first_root)) !=
                  ProjectCache.project_cache_id(test_source(TEST_PROJECT, second_root))
        end
    end

    mktempdir() do dir
        write_test_source(joinpath(dir, "first.csv"))
        write_test_source(joinpath(dir, "second.csv"), 10)
        write(joinpath(dir, "metadata.txt"), "collection_path,wafer\ntest,A\n")
        source = scan_test_source(TEST_PROJECT, dir)
        adapter = test_source(TEST_PROJECT, dir)
        identity = ProjectCache.project_cache_identity(
            "20260430_120005",
            adapter,
        )
        rm(identity.cache_path; force=true)

        try
            progress = NamedTuple[]
            snapshot = ProjectCache.write_project_cache!(
                identity,
                source;
                replace=true,
                on_progress=event -> push!(progress, event),
            )
            @test length(snapshot.source.hierarchy.all_items) == 2
            @test ProjectCache.cache_status(snapshot, source).fresh_source_items == 2
            @test any(event -> event.phase == :cache_finalize, progress)

            loaded = ProjectCache.load_project_cache(identity)
            @test loaded.source.hierarchy.index[("test", "first")].parameters[:wafer] == "A"
            @test Set(
                item.id
                for item in loaded.source.hierarchy.all_items
            ) == Set(
                item.id
                for item in source.hierarchy.all_items
            )
            @test all(
                item -> item.parameters[:wafer] == "A",
                loaded.source.hierarchy.all_items,
            )
            @test ProjectCache.cached_item_data(
                loaded,
                loaded.source.hierarchy.all_items,
            ) == Any[nothing, nothing]

            sleep(0.05)
            touch(first(source.hierarchy.all_items).source_item_path)
            rm(joinpath(dir, "second.csv"))
            write_test_source(joinpath(dir, "third.csv"), 20)
            updated_source = scan_test_source(TEST_PROJECT, dir)
            status = ProjectCache.cache_status(loaded, updated_source)
            @test status.stale_source_items == 1
            @test status.new_source_items == 1
            @test status.deleted_source_items == 1

            updated = ProjectCache.write_project_cache!(loaded.identity, updated_source)
            updated_status = ProjectCache.cache_status(updated, updated_source)
            @test updated_status.fresh_source_items == 2
            @test updated_status.stale_source_items == 0
            @test updated_status.new_source_items == 0
            @test updated_status.deleted_source_items == 0

            failed_source = MeasurementBrowser.SourceScan(
                updated_source.source_id,
                updated_source.source_label,
                updated_source.source_item_fingerprints,
                updated_source.hierarchy,
                [ItemFailure(
                    first(updated_source.hierarchy.all_items).source_item_id,
                    first(updated_source.hierarchy.all_items).id,
                    "synthetic analysis failure",
                )],
            )
            failed = ProjectCache.write_project_cache!(loaded.identity, failed_source)
            @test only(values(failed.analysis_errors)) == "synthetic analysis failure"
        finally
            rm(identity.cache_path; force=true)
        end
    end

    @testset "opening a project updates a stale cache" begin
        mktempdir() do dir
            write_test_source(joinpath(dir, "first.csv"))
            source = scan_test_source(TEST_PROJECT, dir)
            adapter = test_source(TEST_PROJECT, dir)
            identity = ProjectCache.project_cache_identity(
                ProjectCache.project_cache_id(adapter),
                adapter,
            )
            rm(identity.cache_path; force=true)
            try
                ProjectCache.write_project_cache!(identity, source; replace=true)
                write_test_source(joinpath(dir, "second.csv"), 10)

                state = MeasurementBrowser.Browser.BrowserState(project_locked=true)
                MeasurementBrowser.Browser._open_project_path!(
                    state,
                    dir;
                    project=TEST_PROJECT,
                    persist=false,
                )
                saw_writing = false
                deadline = time() + 10
                while time() < deadline
                    workspace = state.workspace
                    MeasurementBrowser.Workspace.poll_workspace!(workspace)
                    saw_writing |= workspace.cache_job.state == :writing
                    status = workspace.cache.status
                    if saw_writing &&
                       workspace.cache_job.state == :ready &&
                       status isa ProjectCache.ProjectCacheStatus &&
                       status.cached_source_items == 2 &&
                       status.new_source_items == 0
                        break
                    end
                    sleep(0.02)
                end

                workspace = state.workspace
                status = workspace.cache.status
                @test saw_writing
                @test workspace.scan.state == :done
                @test workspace.cache_job.state == :ready
                @test status isa ProjectCache.ProjectCacheStatus
                @test status.cached_source_items == 2
                @test status.new_source_items == 0
                close_workspace!(workspace)
            finally
                rm(identity.cache_path; force=true)
            end
        end
    end
end
