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
        cycles = sort(unique(MeasurementBrowser._load_ruo2_pund_fatigue_file(fixture).cycle))
        @test length(expanded) == length(cycles)
        @test all(m -> m.measurement_kind == :pund, expanded)
        @test [m.parameters[:fatigue_cycle] for m in expanded] == cycles
        @test all(m -> m.parameters[:voltage_V] > 0, expanded)
        @test occursin("cycle $(cycles[1]) (fatigue)", display_label(RUO2_PROJECT, expanded[1]))

        selected = expanded[2]
        selected_cycle = Int(selected.parameters[:fatigue_cycle])
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
        @test occursin("cycle $selected_cycle (fatigue)", loaded.title)
        @test loaded.area_um2 == get(device_params, :area_um2, nothing)
        @test loaded.df == expected_df
    end
end
