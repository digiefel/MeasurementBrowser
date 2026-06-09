using MeasurementBrowser
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

function _remove_cache_file(cache_id::AbstractString)
    path = MeasurementBrowser.project_cache_path(cache_id, MeasurementBrowser.RUO2_PROJECT)
    rm(path; force=true)
end

function _same_dataframe(left, right)
    names(left) == names(right) || return false
    nrow(left) == nrow(right) || return false
    return all(name -> left[!, name] == right[!, name], names(left))
end

function _cache_test_ui_state()
    ui_state = Dict{Symbol,Any}()
    MeasurementBrowser._init_scan_state!(ui_state)
    MeasurementBrowser._init_cache_state!(ui_state)
    MeasurementBrowser._init_tag_state!(ui_state)
    MeasurementBrowser._init_figure_script_state!(ui_state)
    MeasurementBrowser._init_plot_state!(ui_state)
    ui_state[:project_preference] = "RuO2"
    ui_state[:recent_projects] = Dict{String,String}[]
    return ui_state
end

@testset "project HDF5 cache" begin
    cache_id = MeasurementBrowser.new_project_cache_id()
    @test occursin(r"^\d{8}_\d{6}$", cache_id)
    @test basename(MeasurementBrowser.project_cache_path(
        "20260430_120001",
        MeasurementBrowser.RUO2_PROJECT,
    )) == "20260430_120001.h5"
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
        _remove_cache_file(cache_id)
        try
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.breakdown)
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.pund)

            source = scan_source(dir; project=MeasurementBrowser.RUO2_PROJECT)
            identity = project_cache_identity(cache_id, MeasurementBrowser.RUO2_PROJECT, source.root_path)
            progress_events = NamedTuple[]
            snapshot = write_project_cache!(
                identity,
                source;
                mode=:build,
                on_progress=progress -> push!(progress_events, progress),
            )
            @test snapshot.identity.cache_id == cache_id
            @test snapshot.hierarchy.has_device_metadata == false
            @test length(snapshot.hierarchy.all_measurements) == 3
            @test snapshot.hierarchy.skipped_count == 0
            @test snapshot.status.total_files == 2
            @test snapshot.status.cached_files == 2
            @test snapshot.status.fresh_files == 2
            @test any(event -> event.phase == :cache_update, progress_events)
            @test last(filter(event -> event.phase == :cache_update, progress_events)).processed_csv == 2

            load_progress = NamedTuple[]
            loaded = load_project_cache(identity; on_progress=progress -> push!(load_progress, progress))
            @test length(load_progress) == 1
            @test only(load_progress).phase == :cache_load
            @test only(load_progress).total_csv == 1
            @test only(load_progress).loaded_measurements == 3
            @test sort([m.unique_id for m in loaded.hierarchy.all_measurements]) ==
                sort([m.unique_id for m in snapshot.hierarchy.all_measurements])

            pund_measurement = only(filter(
                m -> m.measurement_kind === :pund,
                loaded.hierarchy.all_measurements,
            ))
            source_data = only(read_measurement_data(
                MeasurementBrowser.RUO2_PROJECT,
                [pund_measurement],
            ))
            MeasurementBrowser._set_active_project_cache!(loaded.identity)
            @test only(MeasurementBrowser._cached_measurements_data(
                MeasurementBrowser.RUO2_PROJECT,
                [pund_measurement],
            )) === nothing
            loaded_data = only(read_measurement_data(
                MeasurementBrowser.RUO2_PROJECT,
                [pund_measurement],
            ))
            @test _same_dataframe(loaded_data, source_data)
            cached_data = only(MeasurementBrowser._cached_measurements_data(
                MeasurementBrowser.RUO2_PROJECT,
                [pund_measurement],
            ))
            @test cached_data !== nothing
            @test _same_dataframe(cached_data, source_data)
            @test only(MeasurementBrowser._cached_measurements_data(
                MeasurementBrowser.RUO2_PROJECT,
                [pund_measurement];
                processed=true,
            )) === nothing
            processed_data = only(process_measurement_data(
                MeasurementBrowser.RUO2_PROJECT,
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
            @test only(MeasurementBrowser._cached_measurements_data(
                MeasurementBrowser.RUO2_PROJECT,
                [pund_measurement];
                processed=true,
            )) !== nothing
            @test _same_dataframe(
                only(MeasurementBrowser._cached_measurements_data(
                    MeasurementBrowser.RUO2_PROJECT,
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
            @test only(MeasurementBrowser._cached_measurements_data(
                MeasurementBrowser.RUO2_PROJECT,
                [pund_measurement],
            )) === nothing
            rm(joinpath(dir, _CACHE_FIXTURES.breakdown))
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.tlm)

            updated_source = scan_source(dir; project=MeasurementBrowser.RUO2_PROJECT)
            status = cache_status(loaded.identity, updated_source)
            @test status.total_files == 2
            @test status.cached_files == 2
            @test status.fresh_files == 0
            @test status.stale_files == 1
            @test status.new_files == 1
            @test status.deleted_files == 1

            update_progress = NamedTuple[]
            updated = write_project_cache!(
                loaded.identity,
                updated_source;
                mode=:update,
                on_progress=progress -> push!(update_progress, progress),
            )
            @test updated.status.total_files == 2
            @test updated.status.cached_files == 2
            @test updated.status.fresh_files == 2
            @test updated.status.stale_files == 0
            @test updated.status.new_files == 0
            @test updated.status.deleted_files == 0
            cache_update_events = filter(event -> event.phase == :cache_update, update_progress)
            @test last(cache_update_events).processed_csv == 3
            @test last(cache_update_events).total_csv == 3

            failed_source = SourceScan(
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
            analysis_error_snapshot = write_project_cache!(
                loaded.identity,
                failed_source;
                mode=:update,
            )
            @test analysis_error_snapshot.status.error_files == 1
            @test only(analysis_error_snapshot.errors).message == "synthetic analysis failure"
        finally
            MeasurementBrowser._set_active_project_cache!(nothing)
            _remove_cache_file(cache_id)
        end
    end

    @testset "project open updates stale cache" begin
        mktempdir() do dir
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.tlm)
            source = scan_source(dir; project=MeasurementBrowser.RUO2_PROJECT)
            identity = project_cache_identity(
                project_cache_id(dir),
                MeasurementBrowser.RUO2_PROJECT,
                source.root_path,
            )
            try
                write_project_cache!(identity, source; mode=:build)
                _copy_cache_fixture!(dir, _CACHE_FIXTURES.pund)

                ui_state = _cache_test_ui_state()
                MeasurementBrowser._open_project_path!(ui_state, dir; persist=false)
                saw_writing = false
                deadline = time() + 10
                while time() < deadline
                    MeasurementBrowser._poll_cache_events!(ui_state)
                    MeasurementBrowser._poll_source_scan_events!(ui_state)
                    saw_writing |= get(ui_state, :cache_state, :idle) == :writing
                    status = get(ui_state, :cache_status, nothing)
                    if saw_writing &&
                       get(ui_state, :cache_state, :idle) == :ready &&
                       status isa ProjectCacheStatus &&
                       status.cached_files == 2 &&
                       status.new_files == 0
                        break
                    end
                    sleep(0.02)
                end

                status = get(ui_state, :cache_status, nothing)
                @test saw_writing
                @test get(ui_state, :source_scan_state, :idle) == :done
                @test get(ui_state, :cache_state, :idle) == :ready
                @test status isa ProjectCacheStatus
                @test status.cached_files == 2
                @test status.new_files == 0
            finally
                rm(identity.cache_path; force=true)
            end
        end
    end

end
