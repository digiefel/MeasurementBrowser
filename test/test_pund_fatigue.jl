using Test
using DataFrames

using DataLoader
using MeasurementBrowser

const RUO2_PROJECT = MeasurementBrowser.RUO2_PROJECT

@testset "PUND fatigue regression" begin
    fixture = joinpath(@__DIR__, "fixtures", "RuO2", "RuO2test_A12_XI_FeCap_D6_20260511_222203_PUND_Fatigue_V2.csv")
    headerless_fixture = joinpath(@__DIR__, "fixtures", "RuO2", "RuO2test_A11_VII_FeCap_A2_20260130_155835_273K_PUND_Fatigue.csv")

    @test isfile(fixture)
    @test isfile(headerless_fixture)

    @testset "multi-loop Pr stats" begin
        df = DataFrame(
            pulse_idx=[2, 2, 4, 4, 7, 7, 9, 9, 10],
            Q_FE=[0.0, 2e-6, 4e-6, 6e-6, 0.0, 4e-6, 8e-6, 12e-6, NaN],
        )
        @test MeasurementBrowser.pund_pr_value(df, 1e8) ≈ 4.5
    end

    @testset "RuO2 fatigue file helpers" begin
        fatigue_df = MeasurementBrowser.read_pund_fatigue_file(fixture)
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
        @test detect_kind(RUO2_PROJECT, basename(fixture)) == :pund_fatigue

        expanded = measurements_for_file(RUO2_PROJECT, fixture)
        cycles = sort(unique(MeasurementBrowser.read_pund_fatigue_file(fixture).cycle))
        @test length(expanded) == length(cycles)
        @test all(m -> m.measurement_kind == :pund, expanded)
        @test all(m -> m.stats[:wakeup_count] == 0, expanded)
        @test all(m -> !haskey(m.stats, :wakeup_f), expanded)
        @test all(m -> !haskey(m.stats, :wakeup_V), expanded)
        @test Set(m.parameters[:fatigue_f] for m in expanded) == Set([100000.0])
        @test Set(m.parameters[:fatigue_Vamp] for m in expanded) == Set([2.9])
        @test Set(m.parameters[:fatigue_Vbase] for m in expanded) == Set([0.6])
        @test all(m -> m.stats[:fatigue_f] == m.stats[:frequency_kHz] * 1000, expanded)
        @test Set(m.stats[:fatigue_V] for m in expanded) == Set([2.9])
        @test [m.stats[:fatigue_count] for m in expanded] == cycles
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
        selected_count = selected.parameters[:fatigue_idx]
        @test selected_count == cycles[2]
        expected_df = MeasurementBrowser._select_pund_fatigue_cycle(
            MeasurementBrowser.read_pund_fatigue_file(fixture),
            Int(selected_count),
        )
        workspace = MeasurementBrowser.Workspace.Workspace(RUO2_PROJECT, dirname(fixture))
        data = MeasurementBrowser.read_measurement_data(workspace, [selected])

        @test length(data) == 1
        @test only(data) == expected_df
    end

    @testset "RuO2 headerless fatigue file expansion" begin
        expanded = measurements_for_file(RUO2_PROJECT, headerless_fixture)
        cycles = sort(unique(MeasurementBrowser.read_pund_fatigue_file(headerless_fixture).cycle))
        @test length(expanded) == length(cycles)
        @test all(m -> m.measurement_kind == :pund, expanded)
        @test [m.parameters[:fatigue_idx] for m in expanded] == cycles
        @test all(m -> isnan(m.parameters[:fatigue_f]), expanded)
        @test all(m -> isnan(m.parameters[:fatigue_Vamp]), expanded)
        @test all(m -> isnan(m.parameters[:fatigue_Vbase]), expanded)
        @test [m.stats[:fatigue_count] for m in expanded] == cycles
        @test all(m -> m.stats[:wakeup_count] == 0, expanded)
        @test all(m -> m.stats[:fatigue_f] == m.stats[:frequency_kHz] * 1000, expanded)
        @test all(m -> m.stats[:fatigue_V] == m.stats[:V_amp], expanded)
    end
end
