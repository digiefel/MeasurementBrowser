using MeasurementBrowser
using Annotations
using Test

function _drain_figure_script_job!(ui_state)
    for _ in 1:200
        MeasurementBrowser._poll_figure_script_job_events!(ui_state)
        MeasurementBrowser._figure_script_job_running(ui_state) || return
        sleep(0.01)
    end
    error("Timed out waiting for figure-script job to finish")
end

@testset "gui helper utilities" begin
    @test MeasurementBrowser._sanitize_project_preference("auto") == "auto"
    @test MeasurementBrowser._sanitize_project_preference("invalid") == "auto"
    @test MeasurementBrowser._sanitize_project_preference("RuO2") == "RuO2"
    @test MeasurementBrowser._sanitize_figure_script_output_dir("  /tmp/out  ") == "/tmp/out"
    @test MeasurementBrowser._sanitize_figure_script_output_dir(SubString("  /tmp/sub  ", 3, 10)) == "/tmp/sub"

    path_a = pwd()
    path_b = joinpath(pwd(), "missing_folder")
    prefs = Dict{String,Any}(
        "recent_projects" => Any[
            Dict(
                "path" => path_a,
                "project_preference" => "RuO2",
                "figure_script_output_dir" => "/tmp/figures_a",
                "cache_id" => "20260430_120001",
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
    @test recents[1]["cache_id"] == "20260430_120001"
    @test recents[2]["path"] == MeasurementBrowser._normalize_project_path(path_b)
    @test recents[2]["project_preference"] == "auto"
    @test recents[2]["figure_script_output_dir"] == ""
    @test recents[2]["cache_id"] == ""

    recents2 = copy(recents)
    MeasurementBrowser._update_recent_projects(
        recents2,
        path_a,
        SubString(" RuO2 ", 2, 5),
        SubString(" /tmp/out2 ", 2, 10),
        "20260430_120002",
    )
    @test recents2[1]["project_preference"] == "RuO2"
    @test recents2[1]["figure_script_output_dir"] == "/tmp/out2"
    @test recents2[1]["cache_id"] == "20260430_120002"

    @test MeasurementBrowser._sanitize_cache_id("  folder-abc123  ") == "folder-abc123"
    @test_throws ErrorException MeasurementBrowser._sanitize_cache_id("bad/cache-id")

    no_cache_model = MeasurementBrowser._cache_toolbar_model(Dict{Symbol,Any}())
    @test no_cache_model.label == "Cache: No Project"

    cache_identity = MeasurementBrowser.project_cache_identity(
        "20260430_120004",
        MeasurementBrowser.RUO2_PROJECT,
        pwd(),
    )
    missing_model = MeasurementBrowser._cache_toolbar_model(Dict{Symbol,Any}(
        :cache_identity => cache_identity,
        :cache_state => :missing,
    ))
    @test missing_model.label == "Cache: Missing"

    fresh_model = MeasurementBrowser._cache_toolbar_model(Dict{Symbol,Any}(
        :cache_identity => cache_identity,
        :cache_state => :ready,
        :cache_status => MeasurementBrowser.ProjectCacheStatus(2, 2, 2, 0, 0, 0, 0),
    ))
    @test fresh_model.label == "Cache: Fresh"

    unchecked_model = MeasurementBrowser._cache_toolbar_model(Dict{Symbol,Any}(
        :cache_identity => cache_identity,
        :cache_state => :ready,
        :cache_source_checked => false,
        :cache_status => MeasurementBrowser.ProjectCacheStatus(2, 2, 2, 0, 0, 0, 0),
    ))
    @test unchecked_model.label == "Cache: Loaded"
    @test occursin("source not checked", unchecked_model.detail)

    stale_model = MeasurementBrowser._cache_toolbar_model(Dict{Symbol,Any}(
        :cache_identity => cache_identity,
        :cache_state => :ready,
        :cache_status => MeasurementBrowser.ProjectCacheStatus(3, 2, 1, 1, 1, 0, 0),
    ))
    @test stale_model.label == "Cache: Stale"

    errors_model = MeasurementBrowser._cache_toolbar_model(Dict{Symbol,Any}(
        :cache_identity => cache_identity,
        :cache_state => :ready,
        :cache_status => MeasurementBrowser.ProjectCacheStatus(3, 2, 1, 0, 0, 0, 1),
    ))
    @test errors_model.label == "Cache: Errors"
    discovery_model = MeasurementBrowser._cache_activity_model(Dict{Symbol,Any}(
        :scan_state => :cache_discovery,
        :scan_progress => Dict{Symbol,Any}(
            :phase => :cache_discovery,
            :total_csv => 0,
            :processed_csv => 12,
            :loaded_measurements => 0,
            :skipped_csv => 0,
            :current_path => "",
        ),
    ))
    @test discovery_model.title == "Cache: Preparing Build"
    @test discovery_model.progress == "Found 12 source CSV files"

    load_model = MeasurementBrowser._cache_activity_model(Dict{Symbol,Any}(
        :scan_state => :cache_load,
        :scan_progress => Dict{Symbol,Any}(
            :phase => :cache_load,
            :total_csv => 5,
            :processed_csv => 2,
            :loaded_measurements => 14,
            :skipped_csv => 0,
            :current_path => "",
        ),
    ))
    @test load_model.title == "Cache: Loading"
    @test load_model.progress == "Read 2/5 cached files, loaded 14 measurements"

    check_model = MeasurementBrowser._cache_activity_model(Dict{Symbol,Any}(
        :scan_state => :cache_check,
        :scan_progress => Dict{Symbol,Any}(
            :phase => :cache_check,
            :total_csv => 5,
            :processed_csv => 3,
            :loaded_measurements => 0,
            :skipped_csv => 0,
            :current_path => "",
        ),
    ))
    @test check_model.title == "Source: Checking"
    @test check_model.progress == "Checked 3/5 source CSV files"
    source_models = MeasurementBrowser._source_progress_models(Dict{Symbol,Any}(
        :scan_state => :cache_check,
        :source_check_progress => Dict{Symbol,Any}(
            :phase => :cache_check,
            :total_csv => 5,
            :processed_csv => 3,
            :loaded_measurements => 0,
            :skipped_csv => 0,
            :current_path => "",
        ),
    ))
    @test length(source_models) == 1
    @test source_models[1].title == "Source: Checking"
    rescan_models = MeasurementBrowser._source_progress_models(Dict{Symbol,Any}(
        :source_scan_state => :cache_check,
        :source_scan_progress => Dict{Symbol,Any}(
            :phase => :cache_check,
            :total_csv => 8,
            :processed_csv => 2,
            :loaded_measurements => 0,
            :skipped_csv => 0,
            :current_path => "",
        ),
    ))
    @test length(rescan_models) == 1
    @test rescan_models[1].title == "Rescan: Checking Source"
    @test rescan_models[1].progress == "Checking 2/8 source CSV files"
    @test rescan_models[1].show_bar == true
    finding_models = MeasurementBrowser._source_progress_models(Dict{Symbol,Any}(
        :source_scan_state => :cache_check,
        :source_scan_progress => Dict{Symbol,Any}(
            :phase => :cache_check,
            :total_csv => 0,
            :processed_csv => 12,
            :loaded_measurements => 0,
            :skipped_csv => 0,
            :current_path => "",
        ),
    ))
    @test finding_models[1].progress == "Finding source CSV files: 12 found"
    @test finding_models[1].show_bar == false
    done_rescan_models = MeasurementBrowser._source_progress_models(Dict{Symbol,Any}(
        :source_scan_state => :done,
        :source_scan_progress => Dict{Symbol,Any}(
            :phase => :cache_check,
            :total_csv => 8,
            :processed_csv => 8,
            :loaded_measurements => 0,
            :skipped_csv => 0,
            :current_path => "",
        ),
    ))
    @test done_rescan_models[1].title == "Rescan: Complete"
    @test done_rescan_models[1].fraction == 1.0f0
    @test done_rescan_models[1].show_bar == true
    @test isempty(MeasurementBrowser._cache_progress_models(Dict{Symbol,Any}(
        :scan_state => :cache_check,
        :source_check_progress => Dict{Symbol,Any}(),
    )))
    @test isempty(MeasurementBrowser._cache_progress_models(Dict{Symbol,Any}(
        :cache_state => :checking,
        :scan_state => :cache_check,
        :cache_load_progress => Dict{Symbol,Any}(
            :phase => :cache_load,
            :total_csv => 8,
            :processed_csv => 8,
            :loaded_measurements => 32,
            :skipped_csv => 0,
            :current_path => "",
        ),
    )))
    @test !MeasurementBrowser._cache_action_blocked(Dict{Symbol,Any}(
        :cache_state => :checking,
        :scan_state => :cache_check,
    ))
    @test MeasurementBrowser._cache_action_blocked(Dict{Symbol,Any}(
        :cache_state => :loading,
        :scan_state => :cache_load,
    ))
    checking_model = MeasurementBrowser._cache_toolbar_model(Dict{Symbol,Any}(
        :cache_identity => cache_identity,
        :cache_state => :checking,
        :scan_state => :cache_check,
    ))
    @test checking_model.label == "Cache: Loaded"
    @test MeasurementBrowser.MakieImguiIntegration._texture_display_size((1600, 1000), 2.0f0) == (800.0, 500.0)
    @test MeasurementBrowser.MakieImguiIntegration._texture_display_size((800, 500), 1.0f0) == (800.0, 500.0)
    @test MeasurementBrowser.MakieImguiIntegration._imgui_mouse_to_makie(
        (x=150.0, y=125.0),
        (x=100.0, y=100.0),
        (800.0, 500.0),
    ) == (50.0, 475.0)
    makie_mouse = MeasurementBrowser.MakieImguiIntegration.Makie.Mouse
    @test MeasurementBrowser.MakieImguiIntegration._makie_button_for_imgui_right(false) == makie_mouse.right
    @test MeasurementBrowser.MakieImguiIntegration._makie_button_for_imgui_right(true) == makie_mouse.left
    @test MeasurementBrowser.MakieImguiIntegration._makie_button_for_imgui_right(false, true) == makie_mouse.left
    @test MeasurementBrowser.MakieImguiIntegration._should_open_axis_popup(false, false, false)
    @test !MeasurementBrowser.MakieImguiIntegration._should_open_axis_popup(true, false, false)
    @test !MeasurementBrowser.MakieImguiIntegration._should_open_axis_popup(false, true, false)
    @test !MeasurementBrowser.MakieImguiIntegration._should_open_axis_popup(false, false, true)

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
        :selected_all_measurements => MeasurementInfo[m1, m2],
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
    bad_tag_state = Annotations.Tags.TagState(
        [Annotations.Tags.TagDef("bad", (0xff, 0x30, 0x30), 100)],
        Dict{String,Set{String}}("B" => Set(["bad"])),
    )
    ui2 = Dict{Symbol,Any}(
        :scan_hierarchy => hierarchy,
        :measurement_index => Dict("m1" => m1, "m2" => m2),
        :selected_device_paths => ["A", "B"],
        :selected_measurement_ids => ["m1", "m2"],
        :tag_state => bad_tag_state,
        :tag_state_error => "",
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

    mktempdir() do dir
        ui_reload = Dict{Symbol,Any}(
            :scan_hierarchy => hierarchy,
            :measurement_index => Dict("m1" => m1, "m2" => m2),
            :selected_device_paths => ["A", "B"],
            :selected_measurement_ids => ["m1", "m2"],
            :tag_state => Annotations.Tags.TagState(),
            :tag_state_error => "",
            :show_bad => false,
        )
        MeasurementBrowser._apply_visible_selection!(ui_reload)
        @test [node.name for node in ui_reload[:selected_devices]] == ["A", "B"]

        write(joinpath(dir, "bad_measurements"), "device B\n")
        MeasurementBrowser._load_tag_state_for_root!(ui_reload, dir)
        @test "bad" in ui_reload[:tag_state].assignments["B"]
        @test [node.name for node in ui_reload[:selected_devices]] == ["A"]
        @test [m.id for m in ui_reload[:selected_measurements]] == ["m1"]
    end

    @test MeasurementBrowser._selection_targets([1, 2, 3], 2) == [1, 2, 3]
    @test MeasurementBrowser._selection_targets([1, 2, 3], 4) == [4]

    mktempdir() do dir
        ui3 = Dict{Symbol,Any}(
            :root_path => dir,
            :tag_state => Annotations.Tags.TagState(),
            :tag_state_error => "",
            :show_bad => true,
            :scan_hierarchy => hierarchy,
            :measurement_index => Dict("m1" => m1, "m2" => m2),
            :selected_device_paths => ["A", "B"],
            :selected_measurement_ids => ["m1", "m2"],
        )
        MeasurementBrowser._apply_visible_selection!(ui3)

        @test MeasurementBrowser._set_devices_bad!(ui3, ["B"], true)
        @test "bad" in ui3[:tag_state].assignments["B"]
        @test isfile(joinpath(dir, "bad_measurements"))

        @test MeasurementBrowser._set_measurements_bad!(ui3, ["m1"], true)
        @test "bad" in ui3[:tag_state].assignments["m1"]

        @test MeasurementBrowser._set_devices_bad!(ui3, ["B"], false)
        @test !haskey(ui3[:tag_state].assignments, "B")
    end

    ui4 = Dict{Symbol,Any}(
        :tag_state => nothing,
        :tag_state_error => "tag state unavailable",
        :show_bad => false,
    )
    @test !MeasurementBrowser._tag_state_ready(ui4)
    @test MeasurementBrowser._show_bad_effective(ui4)
    @test MeasurementBrowser._device_is_visible(ui4, "A")
    @test MeasurementBrowser._measurement_is_visible(ui4, m1)

    mktempdir() do dir
        write(joinpath(dir, "tags.txt"),
            "[catalog]\nbad\tff3030\t100\n\n[assignments]\nonlyone\n")
        ui5 = Dict{Symbol,Any}(
            :scan_hierarchy => hierarchy,
            :measurement_index => Dict("m1" => m1, "m2" => m2),
            :selected_device_paths => ["A", "B"],
            :selected_measurement_ids => ["m1", "m2"],
            :show_bad => false,
        )
        MeasurementBrowser._load_tag_state_for_root!(ui5, dir)
        @test ui5[:show_bad] == false
        @test MeasurementBrowser._show_bad_effective(ui5)
        @test ui5[:tag_state] === nothing
        @test !isempty(ui5[:tag_state_error])
        @test [node.name for node in ui5[:selected_devices]] == ["A", "B"]
        @test [m.id for m in ui5[:selected_measurements]] == ["m1", "m2"]

        write(joinpath(dir, "tags.txt"),
            "[catalog]\nbad\tff3030\t100\n\n[assignments]\nB\tbad\n")
        MeasurementBrowser._load_tag_state_for_root!(ui5, dir)
        @test MeasurementBrowser._tag_state_ready(ui5)
        @test !MeasurementBrowser._show_bad_effective(ui5)
        @test [node.name for node in ui5[:selected_devices]] == ["A"]
        @test [m.id for m in ui5[:selected_measurements]] == ["m1"]
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
        :selected_all_measurements => MeasurementInfo[m1, m2, m3],
        :selected_measurement_id_set => Set(["m3", "m1"]),
        :all_measurements => [m1, m2, m3],
        :active_scan_id => 1,
        :root_path => "/tmp/project",
    )
    MeasurementBrowser._init_figure_script_state!(ui6)

    @test [measurement.id for measurement in MeasurementBrowser._selected_measurements_in_panel_order(ui6)] == ["m1", "m3"]
    @test MeasurementBrowser._figure_script_output_path(ui6) === nothing
    @test ui6[:figure_script_group_matches_valid] == false
    @test ui6[:figure_script_fact_index_valid] == false

    MeasurementBrowser._set_buffer_string!(ui6[:figure_script_output_dir_buffer], "/tmp/figure_scripts")
    MeasurementBrowser._set_buffer_string!(ui6[:figure_script_name_buffer], "figure_1")
    @test MeasurementBrowser._current_figure_script_output_dir(ui6) == "/tmp/figure_scripts"
    @test MeasurementBrowser._figure_script_output_path(ui6) == "/tmp/figure_scripts/figure_1.jl"

    MeasurementBrowser._set_buffer_string!(ui6[:figure_script_group_name_buffer], "primary")
    MeasurementBrowser._create_figure_script_group_from_selection!(ui6)
    @test MeasurementBrowser._figure_script_job_running(ui6)
    _drain_figure_script_job!(ui6)
    @test length(ui6[:figure_script_job_profiles]) == 1
    @test ui6[:figure_script_job_profiles][1].profile.group_name == "primary"
    @test any(section -> section.key == :collect_candidates, ui6[:figure_script_job_profiles][1].profile.sections)
    @test length(MeasurementBrowser._figure_script_groups(ui6)) == 1
    @test MeasurementBrowser._figure_script_groups(ui6)[1].name == "primary"
    @test ui6[:figure_script_group_matches_valid] == false
    @test [measurement.id for measurement in MeasurementBrowser._group_measurements_in_current_scan(
        ui6,
        MeasurementBrowser._figure_script_groups(ui6)[1],
    )] == ["m1", "m3"]
    @test ui6[:figure_script_group_matches_valid] == true
    @test ui6[:figure_script_fact_index_valid] == true

    ui6[:selected_measurement_ids] = ["m2", "m1"]
    ui6[:selected_measurement_id_set] = Set(["m2", "m1"])
    MeasurementBrowser._add_selection_to_figure_script_group!(ui6)
    @test MeasurementBrowser._figure_script_job_running(ui6)
    _drain_figure_script_job!(ui6)
    @test length(ui6[:figure_script_job_profiles]) == 2
    @test ui6[:figure_script_group_matches_valid] == false
    @test [measurement.id for measurement in MeasurementBrowser._group_measurements_in_current_scan(
        ui6,
        MeasurementBrowser._figure_script_groups(ui6)[1],
    )] == ["m1", "m2", "m3"]
    @test ui6[:figure_script_group_matches_valid] == true
    @test ui6[:figure_script_fact_index_valid] == true

    MeasurementBrowser._set_buffer_string!(ui6[:figure_script_group_name_buffer], "renamed")
    MeasurementBrowser._rename_selected_figure_script_group!(ui6)
    @test MeasurementBrowser._figure_script_groups(ui6)[1].name == "renamed"
    @test ui6[:figure_script_group_matches_valid] == false
    @test [measurement.id for measurement in MeasurementBrowser._group_measurements_in_current_scan(
        ui6,
        MeasurementBrowser._figure_script_groups(ui6)[1],
    )] == ["m1", "m2", "m3"]

    ui6[:selected_measurement_ids] = ["m3"]
    ui6[:selected_measurement_id_set] = Set(["m3"])
    MeasurementBrowser._remove_selection_from_figure_script_group!(ui6)
    @test MeasurementBrowser._figure_script_job_running(ui6)
    _drain_figure_script_job!(ui6)
    @test length(ui6[:figure_script_job_profiles]) == 3
    @test ui6[:figure_script_group_matches_valid] == false
    @test [measurement.id for measurement in MeasurementBrowser._group_measurements_in_current_scan(
        ui6,
        MeasurementBrowser._figure_script_groups(ui6)[1],
    )] == ["m1", "m2"]

    MeasurementBrowser._delete_selected_figure_script_group!(ui6)
    @test isempty(MeasurementBrowser._figure_script_groups(ui6))
    @test ui6[:figure_script_group_matches_valid] == false
    @test ui6[:figure_script_selected_group] == 0

    ui8 = Dict{Symbol,Any}(
        :selected_devices => HierarchyNode[device_leaf],
        :selected_measurement_ids => ["m1"],
        :selected_all_measurements => MeasurementInfo[m1, m2, m3],
        :selected_measurement_id_set => Set(["m1"]),
        :all_measurements => [m1, m2, m3],
        :active_scan_id => 10,
        :root_path => "/tmp/project",
    )
    MeasurementBrowser._init_figure_script_state!(ui8)
    MeasurementBrowser._set_buffer_string!(ui8[:figure_script_group_name_buffer], "stale")
    MeasurementBrowser._create_figure_script_group_from_selection!(ui8)
    @test MeasurementBrowser._figure_script_job_running(ui8)
    @test !isempty(ui8[:figure_script_status])
    MeasurementBrowser._ensure_figure_script_fact_index!(ui8)
    @test ui8[:figure_script_fact_index_valid] == true
    MeasurementBrowser._cancel_figure_script_job!(ui8)
    @test !MeasurementBrowser._figure_script_job_running(ui8)
    @test isempty(ui8[:figure_script_status])
    @test isempty(ui8[:figure_script_error])
    ui8[:active_scan_id] = 11
    MeasurementBrowser._invalidate_figure_script_scan_cache!(ui8)
    @test ui8[:figure_script_fact_index_valid] == false
    _drain_figure_script_job!(ui8)
    @test isempty(MeasurementBrowser._figure_script_groups(ui8))

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
