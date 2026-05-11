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

function _remove_cache_file(cache_id::AbstractString)
    path = MeasurementBrowser.project_cache_path(cache_id, MeasurementBrowser.RUO2_PROJECT)
    rm(path; force=true)
end

function _write_bad_pund_fixture!(dir::AbstractString)
    dst = joinpath(dir, _CACHE_FIXTURES.pund)
    open(dst, "w") do io
        write(io, "This file has a PUND filename but no FE PUND data header.\n")
    end
    return dst
end

function _init_cache_test_ui(cache_id::AbstractString, dir::AbstractString)
    ui_state = Dict{Symbol,Any}()
    MeasurementBrowser._init_scan_state!(ui_state)
    MeasurementBrowser._init_cache_state!(ui_state)
    MeasurementBrowser._init_tag_state!(ui_state)
    MeasurementBrowser._init_figure_script_state!(ui_state)
    MeasurementBrowser._init_plot_state!(ui_state)
    ui_state[:project_preference] = "RuO2"
    ui_state[:recent_projects] = [Dict{String,String}(
        "path" => realpath(dir),
        "project_preference" => "RuO2",
        "figure_script_output_dir" => "",
        "cache_id" => String(cache_id),
    )]
    return ui_state
end

function _drain_scan_job!(ui_state)
    for _ in 1:500
        MeasurementBrowser._poll_cache_events!(ui_state)
        get(ui_state, :cache_events, nothing) === nothing && return
        sleep(0.01)
    end
    error("Timed out waiting for cache job to finish")
end

function _drain_background_jobs!(ui_state)
    for _ in 1:500
        MeasurementBrowser._poll_cache_events!(ui_state)
        MeasurementBrowser._poll_source_scan_events!(ui_state)
        get(ui_state, :cache_events, nothing) === nothing &&
            get(ui_state, :source_scan_events, nothing) === nothing &&
            return
        sleep(0.01)
    end
    error("Timed out waiting for background jobs to finish")
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

            loaded = load_project_cache(identity)
            @test sort([m.id for m in loaded.hierarchy.all_measurements]) ==
                sort([m.id for m in snapshot.hierarchy.all_measurements])
            @test loaded.semantic_fields[:measurement] == [
                :measurement_id,
                :measurement_kind,
                :timestamp,
                :source_file,
            ]

            pund_measurement = only(filter(
                m -> m.measurement_kind === :pund,
                loaded.hierarchy.all_measurements,
            ))
            @test get(pund_measurement.parameters, :global_pund_pulse_count, nothing) == 1
            @test :global_pund_pulse_count in loaded.semantic_fields[:summary]
            cached_plot = MeasurementBrowser._measurement_group_for_cached_plot(
                loaded.identity,
                pund_measurement,
            )
            @test hasproperty(cached_plot, :df)
            @test hasproperty(cached_plot, :pulse_groups)
            @test cached_plot.debug == false
            @test nrow(cached_plot.df) > 0
            cached_stats = MeasurementBrowser.compute_pund_stats_from_analyzed_plot(
                cached_plot,
                pund_measurement.device_info.parameters,
            )
            @test haskey(cached_stats, :voltage_max_V)
            @test haskey(cached_stats, :frequency_kHz)

            debug_ui = Dict{Symbol,Any}(
                :cache_identity => loaded.identity,
                :cache_status => loaded.status,
                :debug_plot_mode => true,
            )
            debug_job = MeasurementBrowser._plot_job(
                debug_ui,
                MeasurementBrowser.RUO2_PROJECT,
                [pund_measurement],
                :pund,
                :main;
                target_id="main",
            )
            cached_debug_plot = MeasurementBrowser._run_plot_job(debug_job, () -> false)
            @test cached_debug_plot.debug == true
            @test hasproperty(cached_debug_plot, :debug_boundaries)
            @test hasproperty(cached_debug_plot, :debug_labels)
            @test nrow(cached_debug_plot.df) == nrow(cached_plot.df)

            sleep(0.05)
            touch(pund_measurement.filepath)
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
        finally
            _remove_cache_file(cache_id)
        end
    end

    mktempdir() do dir
        cache_id = "20260430_120006"
        _remove_cache_file(cache_id)
        try
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.breakdown)
            bad_path = _write_bad_pund_fixture!(dir)

            source = scan_source(dir; project=MeasurementBrowser.RUO2_PROJECT)
            identity = project_cache_identity(cache_id, MeasurementBrowser.RUO2_PROJECT, source.root_path)
            snapshot = write_project_cache!(
                identity,
                source;
                mode=:build,
            )

            # breakdown expands to D1+D2 (2 ok), bad_pund is error (1 included in hierarchy)
            @test length(snapshot.hierarchy.all_measurements) == 3
            @test snapshot.hierarchy.skipped_count == 0
            @test snapshot.status.total_files == 2
            @test snapshot.status.cached_files == 2
            @test snapshot.status.fresh_files == 1
            @test snapshot.status.error_files == 1
            @test length(snapshot.errors) == 1
            @test snapshot.errors[1].path == normpath(abspath(bad_path))
            @test occursin("Could not find FE PUND data header", snapshot.errors[1].message)

            loaded = load_project_cache(identity)
            @test length(loaded.hierarchy.all_measurements) == 3
            @test loaded.status.error_files == 1
            @test length(loaded.errors) == 1

            load_progress = NamedTuple[]
            cached_only = MeasurementBrowser._load_project_cache_contents(
                identity;
                on_progress=progress -> push!(load_progress, progress),
            )
            @test length(cached_only.hierarchy.all_measurements) == 3
            @test any(event -> event.phase == :cache_load, load_progress)
        finally
            _remove_cache_file(cache_id)
        end
    end

    mktempdir() do dir
        cache_id = MeasurementBrowser.project_cache_id(dir)
        _remove_cache_file(cache_id)
        try
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.breakdown)
            source = scan_source(dir; project=MeasurementBrowser.RUO2_PROJECT)
            identity = project_cache_identity(cache_id, MeasurementBrowser.RUO2_PROJECT, source.root_path)
            write_project_cache!(
                identity,
                source;
                mode=:build,
            )

            ui_state = _init_cache_test_ui(cache_id, dir)
            write(joinpath(dir, "bad_measurements"), "device A\n")
            MeasurementBrowser._open_project_path!(ui_state, dir; persist=false)
            @test "bad" in ui_state[:tag_state].assignments["A"]
            _drain_background_jobs!(ui_state)
            @test ui_state[:cache_state] == :ready
            MeasurementBrowser._launch_source_scan_job!(
                ui_state,
                dir,
                MeasurementBrowser.RUO2_PROJECT;
                persist=false,
            )
            _drain_background_jobs!(ui_state)
            @test ui_state[:cache_state] == :ready
            @test ui_state[:cache_status].fresh_files == 1
            @test length(ui_state[:all_measurements]) == 2
            @test get(ui_state, :source_scan_state, :idle) == :done
        finally
            _remove_cache_file(cache_id)
        end
    end

    mktempdir() do dir
        cache_id = MeasurementBrowser.project_cache_id(dir)
        _remove_cache_file(cache_id)
        try
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.breakdown)
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.pund)
            source = scan_source(dir; project=MeasurementBrowser.RUO2_PROJECT)
            identity = project_cache_identity(cache_id, MeasurementBrowser.RUO2_PROJECT, source.root_path)
            write_project_cache!(
                identity,
                source;
                mode=:build,
            )

            source_path = String(dir)
            ui_state = _init_cache_test_ui(cache_id, source_path)
            rm(source_path; recursive=true)
            MeasurementBrowser._open_project_path!(ui_state, source_path; persist=false)
            _drain_scan_job!(ui_state)

            @test ui_state[:cache_state] == :ready
            @test ui_state[:cache_source_checked] == false
            @test get(ui_state, :source_scan_state, :idle) == :idle
            @test length(ui_state[:all_measurements]) == 3
            @test ui_state[:cache_status].cached_files == 2

            pund_meas = only(filter(
                m -> m.measurement_kind === :pund,
                ui_state[:all_measurements],
            ))
            job = MeasurementBrowser._plot_job(
                ui_state,
                MeasurementBrowser.RUO2_PROJECT,
                [pund_meas],
                :pund,
                :main;
                target_id="main",
            )
            cached_pund = MeasurementBrowser._run_plot_job(job, () -> false)
            @test hasproperty(cached_pund, :df)
        finally
            _remove_cache_file(cache_id)
        end
    end
end
