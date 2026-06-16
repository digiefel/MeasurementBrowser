using DataFrames: AbstractDataFrame, nrow
using MeasurementBrowser
using Test

const ProjectCache = MeasurementBrowser.Cache

"""Return whether two dataframes have the same columns and values."""
function same_dataframe(left::AbstractDataFrame, right::AbstractDataFrame)::Bool
    names(left) == names(right) || return false
    nrow(left) == nrow(right) || return false
    return all(name -> left[!, name] == right[!, name], names(left))
end

@testset "project HDF5 cache" begin
    @test basename(ProjectCache.project_cache_identity(
        "20260430_120001",
        TEST_PROJECT,
        pwd(),
    ).cache_path) == "20260430_120001.h5"

    mktempdir() do first_root
        mktempdir() do second_root
            @test ProjectCache.project_cache_id(first_root) !=
                  ProjectCache.project_cache_id(second_root)
        end
    end

    mktempdir() do dir
        write_test_source(joinpath(dir, "first.csv"))
        write_test_source(joinpath(dir, "second.csv"), 10)
        source = MeasurementBrowser.scan_source(dir; project=TEST_PROJECT)
        identity = ProjectCache.project_cache_identity(
            "20260430_120005",
            TEST_PROJECT,
            source.root_path,
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
            @test length(snapshot.source.hierarchy.all_measurements) == 2
            @test ProjectCache.cache_status(snapshot, source).fresh_files == 2

            loaded = ProjectCache.load_project_cache(identity)
            @test Set(
                measurement.unique_id
                for measurement in loaded.source.hierarchy.all_measurements
            ) == Set(
                measurement.unique_id
                for measurement in source.hierarchy.all_measurements
            )

            measurement = first(loaded.source.hierarchy.all_measurements)
            source_workspace = MeasurementBrowser.Workspace.Workspace(TEST_PROJECT, dir)
            source_data = only(read_item_data(source_workspace, [measurement]))

            cached_workspace = MeasurementBrowser.Workspace.Workspace(TEST_PROJECT, dir)
            cached_workspace.cache.index = loaded
            @test only(ProjectCache.cached_measurement_data(
                loaded,
                [measurement],
            )) === nothing

            cached_data = only(read_item_data(cached_workspace, [measurement]))
            @test same_dataframe(cached_data, source_data)
            @test same_dataframe(
                only(ProjectCache.cached_measurement_data(loaded, [measurement])),
                source_data,
            )

            @test only(ProjectCache.cached_measurement_data(
                loaded,
                [measurement];
                processed=true,
            )) === nothing
            processed = only(process_item_data(cached_workspace, [measurement]))
            @test names(processed) == ["x", "y", "processed"]
            @test same_dataframe(
                only(ProjectCache.cached_measurement_data(
                    loaded,
                    [measurement];
                    processed=true,
                )),
                processed,
            )

            sleep(0.05)
            touch(measurement.filepath)
            @test only(ProjectCache.cached_measurement_data(
                loaded,
                [measurement],
            )) === nothing

            rm(joinpath(dir, "second.csv"))
            write_test_source(joinpath(dir, "third.csv"), 20)
            updated_source = MeasurementBrowser.scan_source(dir; project=TEST_PROJECT)
            status = ProjectCache.cache_status(loaded, updated_source)
            @test status.stale_files == 1
            @test status.new_files == 1
            @test status.deleted_files == 1

            updated = ProjectCache.write_project_cache!(loaded.identity, updated_source)
            updated_status = ProjectCache.cache_status(updated, updated_source)
            @test updated_status.fresh_files == 2
            @test updated_status.stale_files == 0
            @test updated_status.new_files == 0
            @test updated_status.deleted_files == 0

            failed_source = MeasurementBrowser.SourceScan(
                updated_source.root_path,
                updated_source.project,
                updated_source.files,
                updated_source.hierarchy,
                [ItemFailure(
                    first(updated_source.files).filepath,
                    first(updated_source.hierarchy.all_measurements).unique_id,
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
            source = MeasurementBrowser.scan_source(dir; project=TEST_PROJECT)
            identity = ProjectCache.project_cache_identity(
                ProjectCache.project_cache_id(dir),
                TEST_PROJECT,
                dir,
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
                       status.cached_files == 2 &&
                       status.new_files == 0
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
                @test status.cached_files == 2
                @test status.new_files == 0
                close_workspace!(workspace)
            finally
                rm(identity.cache_path; force=true)
            end
        end
    end
end
