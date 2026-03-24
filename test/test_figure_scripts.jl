using MeasurementBrowser
using Test

function _copy_fixture(temp_root::AbstractString, fixture_name::AbstractString)
    source = joinpath(@__DIR__, String(fixture_name))
    target = joinpath(temp_root, String(fixture_name))
    cp(source, target; force=true)
    return target
end

@testset "figure script resolution" begin
    mktempdir() do temp_root
        fixture_path = _copy_fixture(
            temp_root,
            "3V PUND Fatigue [RuO2test_A9_VI_D1(2) ; 2025-10-01 17_12_33].csv",
        )
        expanded = expand_measurement(
            MeasurementBrowser.RUO2_PROJECT,
            MeasurementInfo(fixture_path, MeasurementBrowser.RUO2_PROJECT),
        )
        cycle_one = only(filter(measurement -> measurement.parameters[:fatigue_cycle] == 1, expanded))
        cycle_two = only(filter(measurement -> measurement.parameters[:fatigue_cycle] == 2, expanded))

        resolved = MeasurementBrowser._resolve_measurement_lookup(
            temp_root,
            MeasurementBrowser.RUO2_PROJECT,
            [cycle_two.id, cycle_one.id],
        )

        @test resolved[cycle_two.id].id == cycle_two.id
        @test resolved[cycle_two.id].measurement_kind == :pund
        @test resolved[cycle_two.id].parameters[:fatigue_cycle] == 2
        @test resolved[cycle_one.id].parameters[:fatigue_cycle] == 1
    end
end

@testset "figure script preparation" begin
    mktempdir() do temp_root
        fixture_path = _copy_fixture(
            temp_root,
            "TASESNS1c1f_A_2TSNJunction_11_20260224_111623_298K_FourTerminalIV.csv",
        )
        measurement = MeasurementInfo(fixture_path, MeasurementBrowser.TASE_PROJECT)

        groups = [
            NamedMeasurementGroup("single", [measurement.id]),
            NamedMeasurementGroup("repeat", [measurement.id]),
        ]
        data = prepare_measurement_groups(temp_root, MeasurementBrowser.TASE_PROJECT, groups)

        @test [entry.measurement.id for entry in data["single"]] == [measurement.id]
        @test data["single"][1].measurement.measurement_kind == :four_terminal_iv
        @test !isnothing(data["single"][1].loaded)
        @test !isnothing(data["single"][1].analyzed)
        @test [entry.measurement.id for entry in data["repeat"]] == [measurement.id]
    end
end

@testset "figure script validation" begin
    mktempdir() do temp_root
        fixture_path = _copy_fixture(
            temp_root,
            "TASESNS1c1f_A_2TSNJunction_11_20260224_111623_298K_FourTerminalIV.csv",
        )
        measurement = MeasurementInfo(fixture_path, MeasurementBrowser.TASE_PROJECT)

        duplicate_names = [
            NamedMeasurementGroup("group", [measurement.id]),
            NamedMeasurementGroup("group", [measurement.id]),
        ]
        duplicate_ids = [
            NamedMeasurementGroup("group", [measurement.id, measurement.id]),
        ]
        missing_ids = [
            NamedMeasurementGroup("missing", [joinpath(temp_root, "missing.csv")]),
        ]
        measurement_lookup = Dict(measurement.id => measurement)

        @test_throws MeasurementBrowser.FigureScriptValidationError prepare_measurement_groups(
            temp_root,
            MeasurementBrowser.TASE_PROJECT,
            duplicate_names,
        )
        @test_throws MeasurementBrowser.FigureScriptValidationError prepare_measurement_groups(
            temp_root,
            MeasurementBrowser.TASE_PROJECT,
            duplicate_ids,
        )
        @test_throws MeasurementBrowser.FigureScriptResolutionError prepare_measurement_groups(
            temp_root,
            MeasurementBrowser.TASE_PROJECT,
            missing_ids,
        )
        @test_throws MeasurementBrowser.FigureScriptValidationError MeasurementBrowser.write_figure_script(
            "",
            temp_root,
            MeasurementBrowser.TASE_PROJECT,
            "figure_1",
            [NamedMeasurementGroup("group", [measurement.id])],
            measurement_lookup,
        )
        @test_throws MeasurementBrowser.FigureScriptValidationError MeasurementBrowser.write_figure_script(
            "relative/output",
            temp_root,
            MeasurementBrowser.TASE_PROJECT,
            "figure_1",
            [NamedMeasurementGroup("group", [measurement.id])],
            measurement_lookup,
        )
    end
end

@testset "figure script writing" begin
    mktempdir() do temp_root
        _copy_fixture(
            temp_root,
            "TASESNS1c1f_A_2TSNJunction_11_20260224_111623_298K_FourTerminalIV.csv",
        )
        hierarchy = scan_directory(temp_root; project=MeasurementBrowser.TASE_PROJECT)
        measurement_lookup = Dict(measurement.id => measurement for measurement in hierarchy.all_measurements)
        measurement = only(hierarchy.all_measurements)
        groups = [NamedMeasurementGroup("reference", [measurement.id])]
        output_dir = joinpath(temp_root, "scripts_out")

        script_path = MeasurementBrowser.write_figure_script(
            output_dir,
            temp_root,
            MeasurementBrowser.TASE_PROJECT,
            "figure_1",
            groups,
            measurement_lookup,
        )
        @test script_path == joinpath(output_dir, "figure_1.jl")
        @test isfile(script_path)

        contents = read(script_path, String)
        @test occursin("/home/dgfl/code/julia/MeasurementBrowser/docs/figure_scripts.md", contents)
        @test occursin("MeasurementBrowser.NamedMeasurementGroup", contents)
        @test occursin(measurement.id, contents)
        @test occursin("data = MeasurementBrowser.prepare_measurement_groups", contents)

        @test_throws MeasurementBrowser.FigureScriptExistsError MeasurementBrowser.write_figure_script(
            output_dir,
            temp_root,
            MeasurementBrowser.TASE_PROJECT,
            "figure_1",
            groups,
            measurement_lookup,
        )

        module_name = gensym(:FigureScriptInclude)
        module_ref = Main.eval(:(module $module_name
            using MeasurementBrowser
        end))
        Base.include(module_ref, script_path)
        data = Core.eval(module_ref, :data)

        @test [entry.measurement.id for entry in data["reference"]] == [measurement.id]
    end
end
