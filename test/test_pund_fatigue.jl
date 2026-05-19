using Test
using DataFrames

using DataLoader
using MeasurementBrowser

@testset "PUND fatigue regression" begin
    fixture = joinpath(@__DIR__, "fixtures", "RuO2", "RuO2test_A12_XI_FeCap_D6_20260511_222203_PUND_Fatigue_V2.csv")

    @test isfile(fixture)

    @testset "RuO2 fatigue file helpers" begin
        fatigue_df = MeasurementBrowser._load_ruo2_pund_fatigue_file(fixture)
        @test fatigue_df isa DataFrame
        @test names(fatigue_df) == ["cycle", "time", "voltage", "current"]
        cycles = sort(unique(fatigue_df.cycle))
        @test length(cycles) > 1

        cycle = cycles[2]
        cycle_df = MeasurementBrowser._select_pund_fatigue_cycle(fatigue_df, cycle)
        @test cycle_df isa DataFrame
        @test names(cycle_df) == ["time", "current", "voltage"]
        @test nrow(cycle_df) == count(==(cycle), fatigue_df.cycle)
        @test cycle_df.voltage == fatigue_df.voltage[fatigue_df.cycle .== cycle]
        @test cycle_df.current == fatigue_df.current[fatigue_df.cycle .== cycle]
    end

    @testset "RuO2 expansion and staged plot loading" begin
        meas = MeasurementInfo(fixture, RUO2_PROJECT)
        @test meas.measurement_kind == :pund_fatigue

        expanded = measurements_for_file(RUO2_PROJECT, fixture)
        cycles = sort(unique(MeasurementBrowser._load_ruo2_pund_fatigue_file(fixture).cycle))
        @test length(expanded) == length(cycles)
        @test all(m -> m.measurement_kind == :pund, expanded)
        @test all(m -> haskey(m.parameters, :wakeup_count), expanded)
        @test all(m -> haskey(m.parameters, :wakeup_f), expanded)
        @test all(m -> haskey(m.parameters, :wakeup_V), expanded)
        @test all(m -> haskey(m.parameters, :fatigue_count), expanded)
        @test all(m -> haskey(m.parameters, :fatigue_f), expanded)
        @test all(m -> haskey(m.parameters, :fatigue_V), expanded)
        @test all(m -> get(m.parameters, :wakeup_count, NaN) == 0, expanded)
        @test all(m -> isnan(get(m.parameters, :wakeup_f, NaN)), expanded)
        @test all(m -> isnan(get(m.parameters, :wakeup_V, NaN)), expanded)
        @test Set(get(m.parameters, :fatigue_f, NaN) for m in expanded) == Set([100000.0])
        @test Set(get(m.parameters, :fatigue_V, NaN) for m in expanded) == Set([2.9])
        @test [get(m.parameters, :fatigue_count, NaN) for m in expanded] == cycles
        @test all(
            m -> (
                m.stats[:V_base],
                m.stats[:V_min],
                m.stats[:V_max],
                m.stats[:V_amp],
            ) == (0.6, -2.3, 3.5, 2.9),
            expanded,
        )

        selected = expanded[2]
        selected_count = get(selected.parameters, :fatigue_count, NaN)
        @test selected_count == cycles[2]
        selected_cycle = isfinite(selected_count) ? Int(selected_count) : cycles[2]
        expected_df = MeasurementBrowser._select_pund_fatigue_cycle(
            MeasurementBrowser._load_ruo2_pund_fatigue_file(fixture),
            selected_cycle,
        )
        device_params = merge(selected.device_info.parameters, selected.parameters)
        loaded = MeasurementBrowser.load_plot_for_file(
            RUO2_PROJECT,
            fixture,
            :pund;
            device_params,
        )

        @test loaded !== nothing
        @test loaded.area_um2 == get(device_params, :area_um2, nothing)
        @test loaded.df == expected_df
    end
end
