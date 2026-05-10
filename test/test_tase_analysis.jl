using MeasurementBrowser
using Test

@testset "TASE Analysis" begin
    fixture1 = joinpath(@__DIR__, "TASESNS1c1f_A_2TSNJunction_11_20260224_111623_298K_FourTerminalIV.csv")
    fixture2 = joinpath(@__DIR__, "TASESNS1c1f_A_2TSNJunction_31_20260224_111700_298K_FourTerminalIV.csv")
    source1 = index_source_file(fixture1)
    measurements = [
        MeasurementInfo(fixture1, TASE_PROJECT),
        MeasurementInfo(fixture2, TASE_PROJECT),
    ]

    @test parse_device_info(TASE_PROJECT, source1).location == measurements[1].device_info.location

    @testset "available analyses" begin
        options = available_analyses(TASE_PROJECT, measurements)
        @test length(options) == 2
        @test Set(opt.key for opt in options) == Set([:iv_trace_table, :iv_device_summary])
    end

    @testset "iv trace table" begin
        result = run_analysis(TASE_PROJECT, :iv_trace_table, measurements)
        @test result isa AnalysisResult
        @test result.key == :iv_trace_table
        @test result.row_kind == :iv_point
        @test nrow(result.table) == 6
        @test Set(result.table.num_wires) == Set([1, 4])
        @test Set(result.table.wire_diameter_nm) == Set([32.0])
        @test result.table.current_uA == [1.0, 2.0, 3.0, 1.0, 2.0, 3.0]
        @test result.table.voltage_mV == [1.0, 2.0, 3.0, 2.0, 4.0, 6.0]
    end

    @testset "iv device summary" begin
        result = run_analysis(TASE_PROJECT, :iv_device_summary, measurements)
        @test result isa AnalysisResult
        @test result.key == :iv_device_summary
        @test result.row_kind == :device_measurement
        @test nrow(result.table) == 2

        sort!(result.table, :R_ohm)
        g0 = 77.480_917_30e-6
        @test result.table.device_id == ["11", "31"]
        @test result.table.num_wires == [1, 4]
        @test result.table.R_ohm ≈ [1000.0, 2000.0]
        @test result.table.R_ci_low_ohm ≈ [1000.0, 2000.0]
        @test result.table.R_ci_high_ohm ≈ [1000.0, 2000.0]
        @test result.table.G_S ≈ [0.001, 0.0005]
        @test result.table.G_over_G0 ≈ [0.001 / g0, 0.0005 / g0]
        @test result.table.G_per_wire_G0 ≈ [0.001 / g0, 0.0005 / (4 * g0)]
    end
end
