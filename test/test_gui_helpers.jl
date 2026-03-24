using MeasurementBrowser
using Test

@testset "gui helper utilities" begin
    @test MeasurementBrowser._sanitize_project_preference("auto") == "auto"
    @test MeasurementBrowser._sanitize_project_preference("invalid") == "auto"
    @test MeasurementBrowser._sanitize_project_preference("RuO2") == "RuO2"

    path_a = pwd()
    path_b = joinpath(pwd(), "missing_folder")
    prefs = Dict{String,Any}(
        "recent_projects" => Any[
            Dict(
                "path" => path_a,
                "project_preference" => "RuO2",
                "figure_script_output_dir" => "/tmp/figures_a",
            ),
            Dict("path" => path_b, "project_preference" => 99),
            Dict("path" => 42, "project_preference" => "TASE"),
            "bad_entry",
        ],
    )
    recents = MeasurementBrowser._parse_recent_projects(prefs)
    @test length(recents) == 2
    @test recents[1]["path"] == MeasurementBrowser._normalize_project_path(path_a)
    @test recents[1]["project_preference"] == "RuO2"
    @test recents[1]["figure_script_output_dir"] == "/tmp/figures_a"
    @test recents[2]["path"] == MeasurementBrowser._normalize_project_path(path_b)
    @test recents[2]["project_preference"] == "auto"
    @test recents[2]["figure_script_output_dir"] == ""

    items = [1, 2, 3, 4]
    selected = Int[]
    MeasurementBrowser._update_multi_selection!(selected, 2, items, false, false)
    @test selected == [2]
    MeasurementBrowser._update_multi_selection!(selected, 4, items, false, true)
    @test selected == [2, 4]
    MeasurementBrowser._update_multi_selection!(selected, 2, items, false, true)
    @test selected == [4]
    MeasurementBrowser._update_multi_selection!(selected, 1, items, true, false)
    @test selected == [4, 1, 2, 3]

    m1 = MeasurementInfo(
        "m1",
        "m1.csv", "/tmp/m1.csv", "m1", :kind_a, nothing,
        DeviceInfo(["A"], Dict{Symbol,Any}()), Dict{Symbol,Any}(), nothing
    )
    m2 = MeasurementInfo(
        "m2",
        "m2.csv", "/tmp/m2.csv", "m2", :kind_b, nothing,
        DeviceInfo(["B"], Dict{Symbol,Any}()), Dict{Symbol,Any}(), nothing
    )
    leaf_a = HierarchyNode("A", :leaf, HierarchyNode[], MeasurementInfo[m1])
    leaf_b = HierarchyNode("B", :leaf, HierarchyNode[], MeasurementInfo[m2])
    branch = HierarchyNode("branch", :level, HierarchyNode[leaf_a, leaf_b], MeasurementInfo[])
    root = HierarchyNode("/", :root, HierarchyNode[branch], MeasurementInfo[])

    ui = Dict{Symbol,Any}(
        :hierarchy_root => root,
        :selected_devices => HierarchyNode[leaf_a, leaf_b],
    )
    devices = MeasurementBrowser._all_devices(ui)
    @test length(devices) == 2
    @test devices[1].name == "A"
    @test devices[2].name == "B"

    selected_meas = MeasurementBrowser._selected_measurements(ui)
    @test length(selected_meas) == 2
    @test selected_meas[1].filename == "m1.csv"
    @test selected_meas[2].filename == "m2.csv"

    hierarchy = MeasurementHierarchy([m1, m2], "/tmp", false, MeasurementBrowser.RUO2_PROJECT, 0)
    ui2 = Dict{Symbol,Any}(
        :scan_hierarchy => hierarchy,
        :measurement_index => Dict("m1" => m1, "m2" => m2),
        :selected_device_paths => ["A", "B"],
        :selected_measurement_ids => ["m1", "m2"],
        :bad_registry => MeasurementBrowser.BadRegistry(Set(["B"]), Set{String}()),
        :bad_registry_error => "",
        :show_bad => true,
    )
    MeasurementBrowser._apply_visible_selection!(ui2)
    @test [node.name for node in ui2[:selected_devices]] == ["A", "B"]
    @test [m.id for m in ui2[:selected_measurements]] == ["m1", "m2"]
    @test ui2[:selected_device_paths] == ["A", "B"]
    @test ui2[:selected_measurement_ids] == ["m1", "m2"]

    ui2[:show_bad] = false
    MeasurementBrowser._apply_visible_selection!(ui2)
    @test [node.name for node in ui2[:selected_devices]] == ["A"]
    @test [m.id for m in ui2[:selected_measurements]] == ["m1"]
    @test ui2[:selected_device_paths] == ["A", "B"]
    @test ui2[:selected_measurement_ids] == ["m1", "m2"]

    ui2[:show_bad] = true
    MeasurementBrowser._apply_visible_selection!(ui2)
    @test [node.name for node in ui2[:selected_devices]] == ["A", "B"]
    @test [m.id for m in ui2[:selected_measurements]] == ["m1", "m2"]

    @test MeasurementBrowser._selection_targets([1, 2, 3], 2) == [1, 2, 3]
    @test MeasurementBrowser._selection_targets([1, 2, 3], 4) == [4]

    mktempdir() do dir
        ui3 = Dict{Symbol,Any}(
            :root_path => dir,
            :bad_registry => MeasurementBrowser.BadRegistry(),
            :bad_registry_error => "",
            :show_bad => true,
            :scan_hierarchy => hierarchy,
            :measurement_index => Dict("m1" => m1, "m2" => m2),
            :selected_device_paths => ["A", "B"],
            :selected_measurement_ids => ["m1", "m2"],
        )
        MeasurementBrowser._apply_visible_selection!(ui3)

        @test MeasurementBrowser._set_devices_bad!(ui3, ["B"], true)
        @test ui3[:bad_registry].devices == Set(["B"])
        @test isfile(joinpath(dir, "bad_measurements"))

        @test MeasurementBrowser._set_measurements_bad!(ui3, ["m1"], true)
        @test ui3[:bad_registry].measurements == Set(["m1"])

        @test MeasurementBrowser._set_devices_bad!(ui3, ["B"], false)
        @test isempty(ui3[:bad_registry].devices)
    end

    ui4 = Dict{Symbol,Any}(
        :bad_registry => nothing,
        :bad_registry_error => "bad registry unavailable",
        :show_bad => false,
    )
    @test !MeasurementBrowser._bad_registry_ready(ui4)
    @test_throws ErrorException MeasurementBrowser._device_is_visible(ui4, "A")
    @test_throws ErrorException MeasurementBrowser._measurement_is_visible(ui4, m1)

    mktempdir() do dir
        write(joinpath(dir, "bad_measurements"), "bogus\n")
        ui5 = Dict{Symbol,Any}(:show_bad => false)
        MeasurementBrowser._load_bad_registry_for_root!(ui5, dir)
        @test ui5[:show_bad] == true
        @test ui5[:bad_registry] === nothing
        @test !isempty(ui5[:bad_registry_error])
    end

    m3 = MeasurementInfo(
        "m3",
        "m3.csv", "/tmp/m3.csv", "m3", :kind_c, nothing,
        DeviceInfo(["A"], Dict{Symbol,Any}()), Dict{Symbol,Any}(), nothing
    )
    device_leaf = HierarchyNode("A", :leaf, HierarchyNode[], MeasurementInfo[m1, m2, m3])
    ui6 = Dict{Symbol,Any}(
        :selected_devices => HierarchyNode[device_leaf],
        :selected_measurement_ids => ["m3", "m1"],
        :all_measurements => [m1, m2, m3],
    )
    MeasurementBrowser._init_figure_script_state!(ui6)

    @test [measurement.id for measurement in MeasurementBrowser._selected_measurements_in_panel_order(ui6)] == ["m1", "m3"]
    @test MeasurementBrowser._figure_script_output_path(ui6) === nothing

    MeasurementBrowser._set_buffer_string!(ui6[:figure_script_output_dir_buffer], "/tmp/figure_scripts")
    MeasurementBrowser._set_buffer_string!(ui6[:figure_script_name_buffer], "figure_1")
    @test MeasurementBrowser._current_figure_script_output_dir(ui6) == "/tmp/figure_scripts"
    @test MeasurementBrowser._figure_script_output_path(ui6) == "/tmp/figure_scripts/figure_1.jl"

    MeasurementBrowser._set_buffer_string!(ui6[:figure_script_group_name_buffer], "primary")
    MeasurementBrowser._create_figure_script_group_from_selection!(ui6)
    @test length(MeasurementBrowser._figure_script_groups(ui6)) == 1
    @test MeasurementBrowser._figure_script_groups(ui6)[1].name == "primary"
    @test [measurement.id for measurement in MeasurementBrowser._group_measurements_in_current_scan(
        ui6,
        MeasurementBrowser._figure_script_groups(ui6)[1],
    )] == ["m1", "m3"]

    ui6[:selected_measurement_ids] = ["m2", "m1"]
    MeasurementBrowser._add_selection_to_figure_script_group!(ui6)
    @test [measurement.id for measurement in MeasurementBrowser._group_measurements_in_current_scan(
        ui6,
        MeasurementBrowser._figure_script_groups(ui6)[1],
    )] == ["m1", "m2", "m3"]

    MeasurementBrowser._set_buffer_string!(ui6[:figure_script_group_name_buffer], "renamed")
    MeasurementBrowser._rename_selected_figure_script_group!(ui6)
    @test MeasurementBrowser._figure_script_groups(ui6)[1].name == "renamed"

    ui6[:selected_measurement_ids] = ["m3"]
    MeasurementBrowser._remove_selection_from_figure_script_group!(ui6)
    @test [measurement.id for measurement in MeasurementBrowser._group_measurements_in_current_scan(
        ui6,
        MeasurementBrowser._figure_script_groups(ui6)[1],
    )] == ["m1", "m2"]

    MeasurementBrowser._delete_selected_figure_script_group!(ui6)
    @test isempty(MeasurementBrowser._figure_script_groups(ui6))
    @test ui6[:figure_script_selected_group] == 0

    ui7 = Dict{Symbol,Any}(
        :recent_projects => Dict{String,String}[
            Dict(
                "path" => "/tmp/project_a",
                "project_preference" => "RuO2",
                "figure_script_output_dir" => "/tmp/project_a/figures",
            ),
        ],
    )
    MeasurementBrowser._init_figure_script_state!(ui7)
    MeasurementBrowser._reset_figure_script_state!(ui7, "/tmp/project_a")
    @test MeasurementBrowser._figure_script_output_dir_for_path(ui7, "/tmp/project_a") == "/tmp/project_a/figures"
    @test MeasurementBrowser._buffer_string(ui7[:figure_script_output_dir_buffer]) == "/tmp/project_a/figures"
end
