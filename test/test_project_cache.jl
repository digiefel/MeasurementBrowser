using MeasurementBrowser
using DataFrames: nrow
using Test

const _CACHE_FIXTURES = (
    breakdown="Break of oxide 15V [RuO2test_A9_VI_D1D2(1) ; 2025-10-01 16_27_53].csv",
    pund="3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv",
    tlm="TLM_4P [RuO2test_A9_VI_TLML100W2(12) ; 2025-10-01 16_21_45].csv",
)

function _copy_cache_fixture!(dir::AbstractString, name::AbstractString)
    src = joinpath(@__DIR__, "fixtures", "RuO2", name)
    dst = joinpath(dir, name)
    cp(src, dst)
    return dst
end

function _same_dataframe(left, right)
    names(left) == names(right) || return false
    nrow(left) == nrow(right) || return false
    return all(name -> left[!, name] == right[!, name], names(left))
end

function _cache_test_ui_state()
    ui_state = Dict{Symbol,Any}()
    MeasurementBrowser._init_tag_state!(ui_state)
    MeasurementBrowser._init_figure_script_state!(ui_state)
    MeasurementBrowser._init_plot_state!(ui_state)
    ui_state[:project_preference] = "RuO2"
    ui_state[:recent_projects] = Dict{String,String}[]
    return ui_state
end

@testset "project HDF5 cache" begin
    @test basename(MeasurementBrowser.project_cache_identity(
        "20260430_120001",
        MeasurementBrowser.RUO2_PROJECT,
        pwd(),
    ).cache_path) == "20260430_120001.h5"
    mktempdir() do dir_a
        mktempdir() do dir_b
            id_a = MeasurementBrowser.project_cache_id(dir_a)
            id_b = MeasurementBrowser.project_cache_id(dir_b)
            @test id_a != id_b
            @test occursin(r"^[a-z0-9_.-]+-[0-9a-f]{12}$", id_a)
        end
    end

    mktempdir() do dir
        cache_id = "20260430_120005"
        cache_identity = MeasurementBrowser.project_cache_identity(
            cache_id,
            MeasurementBrowser.RUO2_PROJECT,
            dir,
        )
        try
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.breakdown)
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.pund)

            source = MeasurementBrowser.scan_source(
                dir;
                project=MeasurementBrowser.RUO2_PROJECT,
            )
            cache_identity = MeasurementBrowser.project_cache_identity(
                cache_id,
                MeasurementBrowser.RUO2_PROJECT,
                source.root_path,
            )
            rm(cache_identity.cache_path; force=true)
            progress_events = NamedTuple[]
            snapshot = MeasurementBrowser.write_project_cache!(
                cache_identity,
                source;
                replace=true,
                on_progress=progress -> push!(progress_events, progress),
            )
            @test snapshot.identity.cache_id == cache_id
            @test snapshot.source.hierarchy.has_device_metadata == false
            @test length(snapshot.source.hierarchy.all_measurements) == 3
            @test snapshot.source.hierarchy.skipped_count == 0
            built_status = MeasurementBrowser.cache_status(snapshot, source)
            @test built_status.total_files == 2
            @test built_status.cached_files == 2
            @test built_status.fresh_files == 2

            load_progress = NamedTuple[]
            loaded = MeasurementBrowser.load_project_cache(
                cache_identity;
                on_progress=progress -> push!(load_progress, progress),
            )
            @test length(load_progress) == 1
            @test only(load_progress).phase == :cache_load
            @test only(load_progress).total_csv == 1
            @test only(load_progress).loaded_measurements == 3
            @test sort([m.unique_id for m in loaded.source.hierarchy.all_measurements]) ==
                sort([m.unique_id for m in snapshot.source.hierarchy.all_measurements])

            pund_measurement = only(filter(
                m -> m.measurement_kind === :pund,
                loaded.source.hierarchy.all_measurements,
            ))
            source_workspace = MeasurementBrowser.Workspace.Workspace(
                MeasurementBrowser.RUO2_PROJECT,
                dir,
            )
            source_data = only(read_measurement_data(
                source_workspace,
                [pund_measurement],
            ))
            cached_workspace = MeasurementBrowser.Workspace.Workspace(
                MeasurementBrowser.RUO2_PROJECT,
                dir,
            )
            cached_workspace.cache.index = loaded
            @test only(MeasurementBrowser.cached_measurement_data(
                loaded,
                [pund_measurement],
            )) === nothing
            loaded_data = only(read_measurement_data(
                cached_workspace,
                [pund_measurement],
            ))
            @test _same_dataframe(loaded_data, source_data)
            cached_data = only(MeasurementBrowser.cached_measurement_data(
                loaded,
                [pund_measurement],
            ))
            @test cached_data !== nothing
            @test _same_dataframe(cached_data, source_data)
            @test only(MeasurementBrowser.cached_measurement_data(
                loaded,
                [pund_measurement];
                processed=true,
            )) === nothing
            processed_data = only(process_measurement_data(
                cached_workspace,
                [pund_measurement],
            ))
            @test all(name in names(processed_data) for name in [
                "time_us",
                "current_uA",
                "i_fe_uA",
                "q_fe_pC",
                "q_centered_pC",
                "pulse_group_size",
            ])
            @test only(MeasurementBrowser.cached_measurement_data(
                loaded,
                [pund_measurement];
                processed=true,
            )) !== nothing
            @test _same_dataframe(
                only(MeasurementBrowser.cached_measurement_data(
                    loaded,
                    [pund_measurement],
                )),
                source_data,
            )
            @test (
                pund_measurement.stats[:V_base],
                pund_measurement.stats[:V_min],
                pund_measurement.stats[:V_max],
                pund_measurement.stats[:V_amp],
            ) == (0.0, -3.0, 3.0, 3.0)

            sleep(0.05)
            touch(pund_measurement.filepath)
            @test only(MeasurementBrowser.cached_measurement_data(
                loaded,
                [pund_measurement],
            )) === nothing
            rm(joinpath(dir, _CACHE_FIXTURES.breakdown))
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.tlm)

            updated_source = MeasurementBrowser.scan_source(
                dir;
                project=MeasurementBrowser.RUO2_PROJECT,
            )
            status = MeasurementBrowser.cache_status(loaded, updated_source)
            @test status.total_files == 2
            @test status.cached_files == 2
            @test status.fresh_files == 0
            @test status.stale_files == 1
            @test status.new_files == 1
            @test status.deleted_files == 1

            update_progress = NamedTuple[]
            updated = MeasurementBrowser.write_project_cache!(
                loaded.identity,
                updated_source;
                on_progress=progress -> push!(update_progress, progress),
            )
            updated_status = MeasurementBrowser.cache_status(updated, updated_source)
            @test updated_status.total_files == 2
            @test updated_status.cached_files == 2
            @test updated_status.fresh_files == 2
            @test updated_status.stale_files == 0
            @test updated_status.new_files == 0
            @test updated_status.deleted_files == 0

            failed_source = MeasurementBrowser.SourceScan(
                updated_source.root_path,
                updated_source.project,
                updated_source.files,
                updated_source.hierarchy,
                [MeasurementBrowser.MeasurementAnalysisFailure(
                    first(updated_source.files).filepath,
                    first(updated_source.hierarchy.all_measurements).unique_id,
                    "synthetic analysis failure",
                )],
            )
            analysis_error_snapshot = MeasurementBrowser.write_project_cache!(
                loaded.identity,
                failed_source;
            )
            @test length(analysis_error_snapshot.analysis_errors) == 1
            @test only(values(analysis_error_snapshot.analysis_errors)) ==
                "synthetic analysis failure"
        finally
            rm(cache_identity.cache_path; force=true)
        end
    end

    @testset "project open updates stale cache" begin
        mktempdir() do dir
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.tlm)
            source = MeasurementBrowser.scan_source(
                dir;
                project=MeasurementBrowser.RUO2_PROJECT,
            )
            identity = MeasurementBrowser.project_cache_identity(
                MeasurementBrowser.project_cache_id(dir),
                MeasurementBrowser.RUO2_PROJECT,
                source.root_path,
            )
            try
                MeasurementBrowser.write_project_cache!(identity, source; replace=true)
                _copy_cache_fixture!(dir, _CACHE_FIXTURES.pund)

                ui_state = _cache_test_ui_state()
                MeasurementBrowser._open_project_path!(ui_state, dir; persist=false)
                saw_writing = false
                deadline = time() + 10
                while time() < deadline
                    workspace = ui_state[:workspace]
                    MeasurementBrowser.Workspace.poll_workspace!(workspace)
                    saw_writing |= workspace.cache_job.state == :writing
                    status = workspace.cache.status
                    if saw_writing &&
                       workspace.cache_job.state == :ready &&
                       status isa MeasurementBrowser.ProjectCacheStatus &&
                       status.cached_files == 2 &&
                       status.new_files == 0
                        break
                    end
                    sleep(0.02)
                end

                workspace = ui_state[:workspace]
                status = workspace.cache.status
                @test saw_writing
                @test workspace.scan.state == :done
                @test workspace.cache_job.state == :ready
                @test status isa MeasurementBrowser.ProjectCacheStatus
                @test status.cached_files == 2
                @test status.new_files == 0
                MeasurementBrowser.close_workspace!(workspace)
            finally
                rm(identity.cache_path; force=true)
            end
        end
    end

end
