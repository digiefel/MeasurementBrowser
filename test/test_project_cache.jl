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
            @test isempty(loaded.semantic_fields)

            pund_measurement = only(filter(
                m -> m.measurement_kind === :pund,
                loaded.hierarchy.all_measurements,
            ))
            @test (
                pund_measurement.stats[:V_base],
                pund_measurement.stats[:V_min],
                pund_measurement.stats[:V_max],
                pund_measurement.stats[:V_amp],
            ) == (0.0, -3.0, 3.0, 3.0)

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
end
