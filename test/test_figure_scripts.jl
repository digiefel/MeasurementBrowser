using MeasurementBrowser
using Test

function _copy_fixture(temp_root::AbstractString, fixture_name::AbstractString)
    source = joinpath(@__DIR__, String(fixture_name))
    target = joinpath(temp_root, String(fixture_name))
    cp(source, target; force=true)
    return target
end

function _test_measurement(
    id::AbstractString,
    filepath::AbstractString,
    kind::Symbol,
    device_path::Vector{String};
    device_parameters::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    parameters::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)
    return MeasurementInfo(
        String(id),
        basename(String(filepath)),
        String(filepath),
        String(id),
        kind,
        nothing,
        DeviceInfo(copy(device_path), deepcopy(device_parameters)),
        deepcopy(parameters),
        nothing,
    )
end

@testset "figure script group inference" begin
    cycle_1000_a = _test_measurement(
        "cycle_1000_a",
        "/tmp/f1.csv",
        :pund,
        ["RuO2test", "A9", "VI", "D1"];
        parameters=Dict{Symbol,Any}(:fatigue_cycle => 1000),
    )
    cycle_2000_a = _test_measurement(
        "cycle_2000_a",
        "/tmp/f1.csv",
        :pund,
        ["RuO2test", "A9", "VI", "D1"];
        parameters=Dict{Symbol,Any}(:fatigue_cycle => 2000),
    )
    cycle_1000_b = _test_measurement(
        "cycle_1000_b",
        "/tmp/f2.csv",
        :pund,
        ["RuO2test", "A10", "VI", "D2"];
        parameters=Dict{Symbol,Any}(:fatigue_cycle => 1000),
    )
    wakeup_1000 = _test_measurement(
        "wakeup_1000",
        "/tmp/f3.csv",
        :wakeup,
        ["RuO2test", "A11", "VI", "D3"];
        parameters=Dict{Symbol,Any}(:fatigue_cycle => 1000),
    )

    all_measurements = [cycle_1000_a, cycle_2000_a, cycle_1000_b, wakeup_1000]
    group = infer_measurement_group("cycle_1000", [cycle_1000_a, cycle_1000_b], all_measurements)

    matched = MeasurementBrowser._matching_measurements(all_measurements, group)
    @test [measurement.id for measurement in matched] == ["cycle_1000_a", "cycle_1000_b"]
end

@testset "figure script path inference" begin
    selected_a = _test_measurement("selected_a", "/tmp/a.csv", :iv, ["chip", "A", "VI", "D1"])
    selected_b = _test_measurement("selected_b", "/tmp/b.csv", :iv, ["chip", "A", "FE", "D2"])
    outside = _test_measurement("outside", "/tmp/c.csv", :iv, ["chip", "B", "VI", "D3"])

    group = infer_measurement_group("site_a", [selected_a, selected_b], [selected_a, selected_b, outside])

    matched = MeasurementBrowser._matching_measurements([selected_a, selected_b, outside], group)
    @test [measurement.id for measurement in matched] == ["selected_a", "selected_b"]
end

@testset "figure script exact selector inference" begin
    split_a1 = _test_measurement("split_a1", "/tmp/shared.csv", :breakdown, ["chip", "VI", "FecapBD", "A1"])
    split_a2 = _test_measurement("split_a2", "/tmp/shared.csv", :breakdown, ["chip", "VI", "FecapBD", "A2"])
    unrelated = _test_measurement("unrelated", "/tmp/other.csv", :breakdown, ["chip", "VI", "FecapBD", "A3"])

    group = infer_measurement_group("only_a1", [split_a1], [split_a1, split_a2, unrelated])

    @test MeasurementBrowser._matching_measurements([split_a1, split_a2, unrelated], group) == [split_a1]
end

@testset "figure script fact cover inference" begin
    cycle_1000_a = _test_measurement(
        "cycle_1000_a",
        "/tmp/f1.csv",
        :pund,
        ["RuO2test", "A9", "VI", "D1"];
        parameters=Dict{Symbol,Any}(:fatigue_cycle => 1000),
    )
    cycle_2000_a = _test_measurement(
        "cycle_2000_a",
        "/tmp/f1.csv",
        :pund,
        ["RuO2test", "A9", "VI", "D1"];
        parameters=Dict{Symbol,Any}(:fatigue_cycle => 2000),
    )
    cycle_1000_b = _test_measurement(
        "cycle_1000_b",
        "/tmp/f2.csv",
        :pund,
        ["RuO2test", "A10", "VI", "D2"];
        parameters=Dict{Symbol,Any}(:fatigue_cycle => 1000),
    )
    wakeup_1000 = _test_measurement(
        "wakeup_1000",
        "/tmp/f3.csv",
        :wakeup,
        ["RuO2test", "A11", "VI", "D3"];
        parameters=Dict{Symbol,Any}(:fatigue_cycle => 1000),
    )
    all_measurements = [cycle_1000_a, cycle_2000_a, cycle_1000_b, wakeup_1000]

    group = infer_measurement_group(
        "cycle_1000",
        [cycle_1000_a, cycle_1000_b],
        all_measurements,
    )

    matched = MeasurementBrowser._matching_measurements(all_measurements, group)
    @test [measurement.id for measurement in matched] == ["cycle_1000_a", "cycle_1000_b"]
    @test length(group.filter.clauses) == 1
    clause = only(group.filter.clauses)
    @test clause.source_file === nothing
    @test clause.measurement_kind == :pund
    @test isempty(clause.device_path)
    @test clause.parameter_conditions == Pair{Symbol,Any}[:fatigue_cycle => 1000]
end

@testset "figure script fact cover prefers shared path clauses" begin
    selected_a = _test_measurement("selected_a", "/tmp/a.csv", :iv, ["chip", "A", "VI", "D1"])
    selected_b = _test_measurement("selected_b", "/tmp/b.csv", :iv, ["chip", "A", "FE", "D2"])
    outside = _test_measurement("outside", "/tmp/c.csv", :iv, ["chip", "B", "VI", "D3"])

    group, _ = MeasurementBrowser._infer_measurement_group_profiled(
        "site_a",
        [selected_a, selected_b],
        [selected_a, selected_b, outside],
    )

    @test length(group.filter.clauses) == 1
    clause = only(group.filter.clauses)
    @test clause.source_file === nothing
    @test clause.measurement_kind === nothing
    @test clause.device_path == ["chip", "A"]
    @test isempty(clause.parameter_conditions)
end

@testset "figure script inference profiling" begin
    cycle_1000_a = _test_measurement(
        "cycle_1000_a",
        "/tmp/f1.csv",
        :pund,
        ["RuO2test", "A9", "VI", "D1"];
        parameters=Dict{Symbol,Any}(:fatigue_cycle => 1000),
    )
    cycle_2000_a = _test_measurement(
        "cycle_2000_a",
        "/tmp/f1.csv",
        :pund,
        ["RuO2test", "A9", "VI", "D1"];
        parameters=Dict{Symbol,Any}(:fatigue_cycle => 2000),
    )
    cycle_1000_b = _test_measurement(
        "cycle_1000_b",
        "/tmp/f2.csv",
        :pund,
        ["RuO2test", "A10", "VI", "D2"];
        parameters=Dict{Symbol,Any}(:fatigue_cycle => 1000),
    )
    all_measurements = [cycle_1000_a, cycle_2000_a, cycle_1000_b]

    group, profile = MeasurementBrowser._infer_measurement_group_profiled(
        "cycle_1000",
        [cycle_1000_a, cycle_1000_b],
        all_measurements,
    )

    @test group.name == "cycle_1000"
    @test profile.group_name == "cycle_1000"
    @test profile.selected_count == 2
    @test profile.measurement_count == 3
    @test profile.total_ms >= 0.0
    @test any(section -> section.key == :collect_candidates, profile.sections)
    @test any(section -> section.key == :candidate_match_mask, profile.sections)
    @test get(profile.counters, :candidate_count, 0) >= 1
    @test get(profile.counters, :final_clause_count, 0) >= 1
end

@testset "figure script inference profiling on error" begin
    measurement = _test_measurement("m1", "/tmp/f.csv", :kind, ["chip", "A"])

    err = try
        MeasurementBrowser._infer_measurement_group_profiled("empty", MeasurementInfo[], [measurement])
        nothing
    catch caught
        caught
    end
    @test err isa MeasurementBrowser.FigureScriptProfiledError
    @test err.profile.group_name == "empty"
    @test err.profile.selected_count == 0
    @test err.profile.measurement_count == 1
    @test err.profile.total_ms >= 0.0
end

@testset "figure script data preparation" begin
    mktempdir() do temp_root
        fixture_path = _copy_fixture(
            temp_root,
            "TASESNS1c1f_A_2TSNJunction_11_20260224_111623_298K_FourTerminalIV.csv",
        )
        measurement = MeasurementInfo(fixture_path, MeasurementBrowser.TASE_PROJECT)
        group = NamedMeasurementGroup(
            "reference",
            MeasurementGroupFilter([
                MeasurementFilterClause(
                    source_file=measurement.filepath,
                    measurement_kind=measurement.measurement_kind,
                    device_path=measurement.device_info.location,
                ),
            ]),
        )

        data = prepare_figure_script_data(
            temp_root,
            MeasurementBrowser.TASE_PROJECT,
            [measurement.filepath],
            [group],
        )

        @test haskey(data, "reference")
        @test keys(data) == ["reference"]
        @test data.source_files == [measurement.filepath]
        @test length(data.measurements) == 1
        @test data["reference"][1] isa MeasurementBrowser.FigureMeasurement
        @test data["reference"][1].measurement.id == measurement.id
        @test data["reference"][1].parameters == merge(measurement.device_info.parameters, measurement.parameters)
        @test !isnothing(data["reference"][1].loaded)
        @test !isnothing(data["reference"][1].analyzed)
        @test_throws MeasurementBrowser.FigureScriptResolutionError data["missing"]
    end
end

@testset "figure script validation" begin
    measurement = _test_measurement("m1", "/tmp/f.csv", :kind, ["chip", "A"])
    valid_group = NamedMeasurementGroup(
        "group",
        MeasurementGroupFilter([
            MeasurementFilterClause(device_path=["chip"]),
        ]),
    )
    duplicate_groups = [valid_group, NamedMeasurementGroup("group", valid_group.filter)]
    missing_group = NamedMeasurementGroup(
        "missing",
        MeasurementGroupFilter([
            MeasurementFilterClause(device_path=["other"]),
        ]),
    )

    @test_throws MeasurementBrowser.FigureScriptValidationError MeasurementBrowser._validate_named_measurement_groups(duplicate_groups)
    @test_throws MeasurementBrowser.FigureScriptValidationError prepare_figure_script_data(
        "/tmp",
        MeasurementBrowser.TASE_PROJECT,
        String[],
        [valid_group],
    )
    @test_throws MeasurementBrowser.FigureScriptResolutionError MeasurementBrowser._group_matches(
        [missing_group],
        [measurement],
    )
end

@testset "figure script writing" begin
    mktempdir() do temp_root
        _copy_fixture(
            temp_root,
            "TASESNS1c1f_A_2TSNJunction_11_20260224_111623_298K_FourTerminalIV.csv",
        )
        hierarchy = scan_directory(temp_root; project=MeasurementBrowser.TASE_PROJECT)
        measurement = only(hierarchy.all_measurements)
        group = NamedMeasurementGroup(
            "reference",
            MeasurementGroupFilter([
                MeasurementFilterClause(
                    source_file=measurement.filepath,
                    measurement_kind=measurement.measurement_kind,
                    device_path=measurement.device_info.location,
                ),
            ]),
        )
        output_dir = joinpath(temp_root, "scripts_out")

        script_path = MeasurementBrowser.write_figure_script(
            output_dir,
            temp_root,
            MeasurementBrowser.TASE_PROJECT,
            "figure_1",
            [group],
            hierarchy.all_measurements,
        )
        @test script_path == joinpath(output_dir, "figure_1.jl")
        @test isfile(script_path)

        contents = read(script_path, String)
        docs_path = abspath(joinpath(@__DIR__, "..", "docs", "figure_scripts.md"))
        @test occursin(docs_path, contents)
        @test occursin("MeasurementBrowser.MeasurementFilterClause", contents)
        @test occursin("source_files = [", contents)
        @test occursin("prepare_figure_script_data", contents)
        @test !occursin("#cycle=", contents)
        @test !occursin("#split=", contents)

        @test_throws MeasurementBrowser.FigureScriptExistsError MeasurementBrowser.write_figure_script(
            output_dir,
            temp_root,
            MeasurementBrowser.TASE_PROJECT,
            "figure_1",
            [group],
            hierarchy.all_measurements,
        )

        module_name = gensym(:FigureScriptInclude)
        module_ref = Main.eval(:(module $module_name
            using MeasurementBrowser
        end))
        Base.include(module_ref, script_path)
        data = Core.eval(module_ref, :data)

        @test data isa MeasurementBrowser.FigureScriptData
        @test [entry.measurement.id for entry in data["reference"]] == [measurement.id]
    end
end
