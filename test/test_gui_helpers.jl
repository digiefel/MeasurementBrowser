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
            Dict("path" => path_a, "project_preference" => "RuO2"),
            Dict("path" => path_b, "project_preference" => 99),
            Dict("path" => 42, "project_preference" => "TASE"),
            "bad_entry",
        ],
    )
    recents = MeasurementBrowser._parse_recent_projects(prefs)
    @test length(recents) == 2
    @test recents[1]["path"] == MeasurementBrowser._normalize_project_path(path_a)
    @test recents[1]["project_preference"] == "RuO2"
    @test recents[2]["path"] == MeasurementBrowser._normalize_project_path(path_b)
    @test recents[2]["project_preference"] == "auto"

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
    MeasurementBrowser._sync_visible_selection!(ui2)
    @test [node.name for node in ui2[:selected_devices]] == ["A", "B"]
    @test [m.id for m in ui2[:selected_measurements]] == ["m1", "m2"]

    ui2[:show_bad] = false
    MeasurementBrowser._sync_visible_selection!(ui2)
    @test [node.name for node in ui2[:selected_devices]] == ["A"]
    @test [m.id for m in ui2[:selected_measurements]] == ["m1"]

    ui2[:show_bad] = true
    MeasurementBrowser._sync_visible_selection!(ui2)
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
        MeasurementBrowser._sync_visible_selection!(ui3)

        @test MeasurementBrowser._set_devices_bad!(ui3, ["B"], true)
        @test ui3[:bad_registry].devices == Set(["B"])
        @test isfile(joinpath(dir, "bad_measurements"))

        @test MeasurementBrowser._set_measurements_bad!(ui3, ["m1"], true)
        @test ui3[:bad_registry].measurements == Set(["m1"])

        @test MeasurementBrowser._set_devices_bad!(ui3, ["B"], false)
        @test isempty(ui3[:bad_registry].devices)
    end
end
