using MeasurementBrowser
using Test

@testset "project view state" begin
    root_path = mktempdir()
    @test basename(MeasurementBrowser._project_view_file_path(root_path)) == "measurementbrowser.toml"
    default_view = MeasurementBrowser._load_project_view(root_path)
    @test default_view.project == ""
    @test default_view.main_plot.id == "main"
    @test default_view.main_plot.live == true

    view = MeasurementBrowser.PersistedProjectView(
        project="RuO2",
        tree=MeasurementBrowser.PersistedTreeView(
            expanded=["chip/device"],
            selected=["chip/device"],
            filter="tlm",
        ),
        measurements=MeasurementBrowser.PersistedMeasurementsView(
            selected=["measurement-1", "measurement-2"],
            filter="298K",
        ),
        plot_kinds=Dict("iv_sweep" => "RuO2IVSweepPlot"),
        main_plot=MeasurementBrowser.PersistedPlotView(
            id="main",
            title="Plot Area",
            plot_kind="RuO2TLM4PointPlot",
            live=true,
            measurements=["measurement-1"],
        ),
        plot_windows=[
            MeasurementBrowser.PersistedPlotView(
                id="plot_1",
                title="Detached",
                plot_kind="RuO2IVSweepPlot",
                live=false,
                measurements=["measurement-2"],
            ),
        ],
    )

    MeasurementBrowser._save_project_view(root_path, view)
    @test isfile(joinpath(root_path, "measurementbrowser.toml"))
    loaded = MeasurementBrowser._load_project_view(root_path)
    @test loaded.project == view.project
    @test loaded.tree.expanded == view.tree.expanded
    @test loaded.tree.selected == view.tree.selected
    @test loaded.tree.filter == view.tree.filter
    @test loaded.measurements.selected == view.measurements.selected
    @test loaded.measurements.filter == view.measurements.filter
    @test loaded.main_plot.plot_kind == view.main_plot.plot_kind
    @test loaded.main_plot.measurements == view.main_plot.measurements
    @test only(loaded.plot_windows).plot_kind == only(view.plot_windows).plot_kind
    @test only(loaded.plot_windows).measurements == only(view.plot_windows).measurements

    toml_data = MeasurementBrowser._project_view_to_toml(view)
    @test Set(keys(toml_data)) == Set([
        "project",
        "tree",
        "measurements",
        "plot_kinds",
        "main_plot",
        "plot_windows",
    ])

    parsed = MeasurementBrowser._project_view_from_toml(MeasurementBrowser.PersistedProjectView, Dict{String,Any}(
        "project" => "RuO2",
        "tree" => Dict{String,Any}(
            "expanded" => ["chip/device"],
            "selected" => ["chip/device"],
            "filter" => "tlm",
        ),
        "measurements" => Dict{String,Any}(
            "selected" => ["measurement-1"],
            "filter" => "298K",
        ),
        "plot_kinds" => Dict{String,Any}("iv_sweep" => "RuO2IVSweepPlot"),
        "main_plot" => Dict{String,Any}(
            "id" => "main",
            "title" => "Plot Area",
            "plot_kind" => "RuO2TLM4PointPlot",
            "live" => true,
            "measurements" => ["measurement-1"],
        ),
        "plot_windows" => Any[
            Dict{String,Any}(
                "id" => "plot_1",
                "title" => "Detached",
                "plot_kind" => "RuO2IVSweepPlot",
                "live" => false,
                "measurements" => ["measurement-2"],
            ),
        ],
    ))

    @test_throws ErrorException MeasurementBrowser._project_view_from_toml(MeasurementBrowser.PersistedProjectView, Dict{String,Any}(
        "project" => "RuO2",
        "tree" => Dict{String,Any}(
            "expanded" => Any["chip/device", 1],
            "selected" => ["chip/device"],
            "filter" => "tlm",
        ),
        "measurements" => Dict{String,Any}(
            "selected" => ["measurement-1"],
            "filter" => "298K",
        ),
        "plot_kinds" => Dict{String,Any}("iv_sweep" => "RuO2IVSweepPlot"),
        "main_plot" => Dict{String,Any}(
            "id" => "main",
            "title" => "Plot Area",
            "plot_kind" => "RuO2TLM4PointPlot",
            "live" => true,
            "measurements" => ["measurement-1"],
        ),
        "plot_windows" => Any[],
    ))

    @test parsed.project == "RuO2"
    @test parsed.tree.expanded == ["chip/device"]
    @test parsed.measurements.selected == ["measurement-1"]
    @test parsed.main_plot.id == "main"
    @test parsed.main_plot.title == "Plot Area"
    @test parsed.main_plot.live == true
    @test only(parsed.plot_windows).measurements == ["measurement-2"]

    project = MeasurementBrowser.RuO2Project()
    measurement_1 = MeasurementBrowser.MeasurementInfo(;
        unique_id="measurement-1",
        filepath=joinpath(root_path, "measurement-1.csv"),
        clean_title="Measurement 1",
        measurement_kind=:iv_sweep,
        device_info=MeasurementBrowser.DeviceInfo(["chip", "device-1"]),
    )
    measurement_2 = MeasurementBrowser.MeasurementInfo(;
        unique_id="measurement-2",
        filepath=joinpath(root_path, "measurement-2.csv"),
        clean_title="Measurement 2",
        measurement_kind=:iv_sweep,
        device_info=MeasurementBrowser.DeviceInfo(["chip", "device-2"]),
    )
    hierarchy = MeasurementBrowser.MeasurementHierarchy(
        [measurement_1, measurement_2],
        root_path,
        true,
        project,
    )
    ui_state = Dict{Symbol,Any}(
        :root_path => root_path,
        :project => project,
        :scan_hierarchy => hierarchy,
        :hierarchy_root => hierarchy.root,
        :all_measurements => hierarchy.all_measurements,
        :measurement_index => MeasurementBrowser._build_measurement_index(hierarchy.all_measurements),
        :selected_device_paths => String[],
        :selected_measurement_ids => String[],
        :selected_devices => MeasurementBrowser.HierarchyNode[],
        :selected_measurements => MeasurementBrowser.MeasurementInfo[],
        :selected_all_measurements => MeasurementBrowser.MeasurementInfo[],
        :selected_measurement_id_set => Set{String}(),
        :selected_path => String[],
        :_plot_window_counter => 0,
    )

    MeasurementBrowser._apply_project_view!(ui_state, view)
    @test ui_state[:selected_device_paths] == ["chip/device"]
    @test ui_state[:selected_measurement_ids] == ["measurement-1", "measurement-2"]
    @test ui_state[:tree_filter] == "tlm"
    @test ui_state[:measurement_filter] == "298K"
    @test ui_state[:plot_kind_by_measurement_kind] == Dict("iv_sweep" => "RuO2IVSweepPlot")
    @test ui_state[:main_plot_live] == true
    @test ui_state[:main_plot_kind] == MeasurementBrowser.RuO2TLM4PointPlot
    @test only(ui_state[:main_plot_measurements]).unique_id == "measurement-1"
    @test only(ui_state[:open_plot_windows])[:target_id] == "plot_1"
    @test only(only(ui_state[:open_plot_windows])[:measurements]).unique_id == "measurement-2"

    saved_view = MeasurementBrowser._project_view_from_ui_state(ui_state)
    @test saved_view.project == "RuO2"
    @test saved_view.tree.selected == ["chip/device"]
    @test saved_view.tree.expanded == ["chip/device"]
    @test saved_view.measurements.selected == ["measurement-1", "measurement-2"]
    @test saved_view.plot_kinds == Dict("iv_sweep" => "RuO2IVSweepPlot")
    @test saved_view.main_plot.plot_kind == "RuO2TLM4PointPlot"
    @test only(saved_view.plot_windows).plot_kind == "RuO2IVSweepPlot"

    MeasurementBrowser._save_project_view_if_changed!(ui_state)
    loaded_after_ui_save = MeasurementBrowser._load_project_view(root_path)
    @test MeasurementBrowser._same_project_view(loaded_after_ui_save, saved_view)

    ui_state[:project_view_loaded] = view
    MeasurementBrowser._begin_scan!(ui_state, root_path, project, true)
    @test ui_state[:selected_device_paths] == ["chip/device"]
    @test ui_state[:selected_measurement_ids] == ["measurement-1", "measurement-2"]
    @test haskey(ui_state, :main_plot_kind)
end
