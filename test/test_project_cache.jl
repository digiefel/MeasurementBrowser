using MeasurementBrowser
using Test

const _CACHE_FIXTURES = (
    wakeup="Wakeup 3V [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_10_48].csv",
    pund="3V FE PUND [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv",
    tlm="TLM_4P [RuO2test_A9_VI_TLML100W2(12) ; 2025-10-01 16_21_45].csv",
)

function _copy_cache_fixture!(dir::AbstractString, name::AbstractString)
    src = joinpath(@__DIR__, name)
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

@testset "project HDF5 cache" begin
    cache_id = MeasurementBrowser.new_project_cache_id()
    @test occursin(r"^\d{8}_\d{6}$", cache_id)
    @test basename(MeasurementBrowser.project_cache_path(
        "20260430_120001",
        MeasurementBrowser.RUO2_PROJECT,
    )) == "20260430_120001.h5"

    mktempdir() do dir
        cache_id = "20260430_120005"
        _remove_cache_file(cache_id)
        try
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.wakeup)
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.pund)

            snapshot = build_project_cache!(
                dir,
                MeasurementBrowser.RUO2_PROJECT,
                cache_id;
                full_rebuild=true,
            )
            @test snapshot.identity.cache_id == cache_id
            @test snapshot.hierarchy.has_device_metadata == false
            @test length(snapshot.hierarchy.all_measurements) == 2
            @test snapshot.hierarchy.skipped_count == 0
            @test snapshot.status.total_files == 2
            @test snapshot.status.cached_files == 2
            @test snapshot.status.fresh_files == 2

            loaded = load_project_cache(dir, MeasurementBrowser.RUO2_PROJECT, cache_id)
            @test [m.id for m in loaded.hierarchy.all_measurements] ==
                [m.id for m in snapshot.hierarchy.all_measurements]
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
            cached_plot = MeasurementBrowser._measurement_group_for_cached_plot(
                loaded.identity,
                pund_measurement,
            )
            @test hasproperty(cached_plot, :df)
            @test hasproperty(cached_plot, :pulse_groups)
            @test nrow(cached_plot.df) > 0

            sleep(0.05)
            touch(pund_measurement.filepath)
            rm(joinpath(dir, _CACHE_FIXTURES.wakeup))
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.tlm)

            status = cache_status(loaded.identity)
            @test status.total_files == 2
            @test status.cached_files == 2
            @test status.fresh_files == 0
            @test status.stale_files == 1
            @test status.new_files == 1
            @test status.deleted_files == 1
        finally
            _remove_cache_file(cache_id)
        end
    end

    mktempdir() do dir
        cache_id = "20260430_120006"
        _remove_cache_file(cache_id)
        try
            _copy_cache_fixture!(dir, _CACHE_FIXTURES.wakeup)
            bad_path = _write_bad_pund_fixture!(dir)

            snapshot = build_project_cache!(
                dir,
                MeasurementBrowser.RUO2_PROJECT,
                cache_id;
                full_rebuild=true,
            )

            @test length(snapshot.hierarchy.all_measurements) == 1
            @test snapshot.hierarchy.skipped_count == 1
            @test snapshot.status.total_files == 2
            @test snapshot.status.cached_files == 2
            @test snapshot.status.fresh_files == 1
            @test snapshot.status.error_files == 1
            @test length(snapshot.errors) == 1
            @test snapshot.errors[1].path == bad_path
            @test occursin("Could not find FE PUND data header", snapshot.errors[1].message)

            loaded = load_project_cache(dir, MeasurementBrowser.RUO2_PROJECT, cache_id)
            @test length(loaded.hierarchy.all_measurements) == 1
            @test loaded.status.error_files == 1
            @test length(loaded.errors) == 1
        finally
            _remove_cache_file(cache_id)
        end
    end
end
