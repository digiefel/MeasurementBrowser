using Test
using DataFrames

using DataLoader
using MeasurementBrowser

@testset "PUND fatigue regression" begin
    fixture = joinpath(@__DIR__, "3V PUND Fatigue [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv")

    @test isfile(fixture)

    @testset "DataLoader cycle readers" begin
        cycles = read_pund_fatigue_cycles(basename(fixture), @__DIR__)
        @test cycles == [1, 2]

        cycle_df = read_pund_fatigue_cycle(basename(fixture), @__DIR__, 2)
        @test cycle_df isa DataFrame
        @test names(cycle_df) == ["time", "current", "voltage"]
        @test nrow(cycle_df) == 3
        @test cycle_df.voltage == [0.0, 3.0, -3.0]
        @test cycle_df.current == [1.5e-9, 2.5e-6, -2.5e-6]
    end

    @testset "RuO2 expansion and plot loading" begin
        meas = MeasurementInfo(fixture)
        @test meas.measurement_kind == :pund_fatigue

        expanded = expand_measurement(RUO2_PROJECT, meas)
        @test length(expanded) == 2
        @test all(m -> m.measurement_kind == :pund, expanded)
        @test [m.parameters[:fatigue_cycle] for m in expanded] == [1, 2]
        @test all(m -> m.parameters[:voltage_V] == 3.0, expanded)
        @test display_label(RUO2_PROJECT, expanded[1]) == "2025-10-01T17:12:33 FE PUND 3.0V cycle 1 (fatigue)"

        device_params = merge(expanded[2].device_info.parameters, expanded[2].parameters)
        loaded = MeasurementBrowser.load_plot_input_for_file(
            RUO2_PROJECT,
            fixture,
            :pund;
            device_params,
        )

        @test loaded !== nothing
        @test loaded.title == "3V PUND Fatigue [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33] cycle 2 (fatigue)"
        @test loaded.area_um2 === nothing
        @test nrow(loaded.df) == 3
        @test loaded.df.voltage == [0.0, 3.0, -3.0]
        @test loaded.df.current == [1.5e-9, 2.5e-6, -2.5e-6]
    end
end
