using Test
using DataFrames

using DataLoader
using MeasurementBrowser

@testset "PUND fatigue regression" begin
    fixture = joinpath(@__DIR__, "fixtures", "RuO2", "3V PUND Fatigue [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv")

    @test isfile(fixture)

    @testset "DataLoader cycle readers" begin
        cycles = read_pund_fatigue_cycles(basename(fixture), dirname(fixture))
        @test cycles == [1, 2]

        cycle_df = read_pund_fatigue_cycle(basename(fixture), dirname(fixture), 2)
        @test cycle_df isa DataFrame
        @test names(cycle_df) == ["time", "current", "voltage"]
        @test nrow(cycle_df) == 3
        @test cycle_df.voltage == [0.0, 3.0, -3.0]
        @test cycle_df.current == [1.5e-9, 2.5e-6, -2.5e-6]
    end

    @testset "RuO2 fatigue scan requires monotonic cycle blocks" begin
        mktempdir() do dir
            bad_path = joinpath(dir, "bad_fatigue.csv")
            write(bad_path, join([
                "Cycle,Time_s,Voltage_V,Current_A",
                "1,0.0,1.0,0.0",
                "2,1.0,2.0,0.0",
                "1,2.0,3.0,0.0",
                "",
            ], "\n"))
            err = try
                MeasurementBrowser._ruo2_scan_fatigue_file(bad_path)
                nothing
            catch caught
                caught
            end
            @test err isa ErrorException
            @test occursin("non-monotonic cycle blocks", sprint(showerror, err))
        end
    end

    @testset "FE PUND loader errors explicitly on malformed files" begin
        mktempdir() do dir
            bad_path = joinpath(dir, "bad_pund.csv")
            write(bad_path, "not,a,real,pund,file\n1,2,3,4\n")
            err = try
                read_fe_pund(basename(bad_path), dir)
                nothing
            catch caught
                caught
            end
            @test err isa ErrorException
            @test occursin("Could not find FE PUND data header", sprint(showerror, err))
        end
    end

    @testset "FE PUND loader accepts compact source CSVs" begin
        mktempdir() do dir
            path = joinpath(dir, "compact_pund.csv")
            write(path, join([
                "Time_s,Current_A,Voltage_V",
                "0.0,1.0e-9,0.0",
                "1.0e-6,2.0e-6,3.0",
                "2.0e-6,-2.0e-6,-3.0",
                "",
            ], "\n"))

            df = read_fe_pund(basename(path), dir)
            @test names(df) == ["time", "current", "voltage", "current_time", "voltage_time"]
            @test df.time == [0.0, 1.0e-6, 2.0e-6]
            @test df.current == [1.0e-9, 2.0e-6, -2.0e-6]
            @test df.voltage == [0.0, 3.0, -3.0]
            @test df.current_time == df.time
            @test df.voltage_time == df.time
        end
    end

    @testset "RuO2 expansion and staged plot loading" begin
        meas = MeasurementInfo(fixture, RUO2_PROJECT)
        @test meas.measurement_kind == :pund_fatigue

        expanded = measurements_for_file(RUO2_PROJECT, fixture)
        @test length(expanded) == 2
        @test all(m -> m.measurement_kind == :pund, expanded)
        @test [m.parameters[:fatigue_cycle] for m in expanded] == [1, 2]
        @test all(m -> m.parameters[:voltage_V] == 3.0, expanded)
        @test display_label(RUO2_PROJECT, expanded[1]) == "2025-10-01T17:12:33 FE PUND 3.0V cycle 1 (fatigue)"

        device_params = merge(expanded[2].device_info.parameters, expanded[2].parameters)
        loaded = MeasurementBrowser.load_plot_for_file(
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
