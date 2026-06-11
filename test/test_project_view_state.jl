using MeasurementBrowser
using Test

const Browser = MeasurementBrowser.Browser

@testset "project view state" begin
    root_path = mktempdir()
    @test basename(Browser._project_view_file_path(root_path)) == "measurementbrowser.toml"
    default_view = Browser._load_project_view(root_path)
    @test default_view.project == ""
    @test default_view.main_plot.id == "main"
    @test default_view.main_plot.live == true

    view = Browser.PersistedProjectView(
        project="RuO2",
        tree=Browser.PersistedTreeView(
            expanded=["chip/device"],
            selected=["chip/device"],
            filter="tlm",
        ),
        measurements=Browser.PersistedMeasurementsView(
            selected=["measurement-1", "measurement-2"],
            filter="298K",
        ),
        plot_kinds=Dict("iv_sweep" => "RuO2IVSweepPlot"),
        main_plot=Browser.PersistedPlotView(
            id="main",
            title="Plot Area",
            plot_kind="RuO2TLM4PointPlot",
            live=true,
            measurements=["measurement-1"],
        ),
        plot_windows=[
            Browser.PersistedPlotView(
                id="plot_1",
                title="Detached",
                plot_kind="RuO2IVSweepPlot",
                live=false,
                measurements=["measurement-2"],
            ),
        ],
    )

    Browser._save_project_view(root_path, view)
    @test isfile(joinpath(root_path, "measurementbrowser.toml"))
    loaded = Browser._load_project_view(root_path)
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

    toml_data = Browser._project_view_to_toml(view)
    @test Set(keys(toml_data)) == Set([
        "project",
        "tree",
        "measurements",
        "plot_kinds",
        "main_plot",
        "plot_windows",
    ])

    parsed = Browser._project_view_from_toml(Browser.PersistedProjectView, Dict{String,Any}(
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

    @test_throws ErrorException Browser._project_view_from_toml(Browser.PersistedProjectView, Dict{String,Any}(
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
    workspace = MeasurementBrowser.Workspace.Workspace(project, root_path)
    MeasurementBrowser.Workspace.replace_measurement_index!(workspace, hierarchy)
    state = Browser.BrowserState(workspace=workspace)

    Browser._apply_project_view!(state, view)
    @test workspace.selection.device_paths == ["chip/device"]
    @test workspace.selection.measurement_ids == ["measurement-1", "measurement-2"]
    @test state.tree_filter == "tlm"
    @test state.measurement_filter == "298K"
    plots = state.plots
    @test plots.kind_by_measurement ==
          Dict(:iv_sweep => MeasurementBrowser.RuO2IVSweepPlot)
    @test plots.main.live == true
    @test plots.main.plot_kind == MeasurementBrowser.RuO2TLM4PointPlot
    @test plots.main.measurement_ids == ["measurement-1"]
    @test only(plots.windows).id == "plot_1"
    @test only(plots.windows).measurement_ids == ["measurement-2"]

    saved_view = Browser._project_view_from_browser(state)
    @test saved_view.project == "RuO2"
    @test saved_view.tree.selected == ["chip/device"]
    @test saved_view.tree.expanded == ["chip/device"]
    @test saved_view.measurements.selected == ["measurement-1", "measurement-2"]
    @test saved_view.plot_kinds == Dict("iv_sweep" => "RuO2IVSweepPlot")
    @test saved_view.main_plot.plot_kind == "RuO2TLM4PointPlot"
    @test only(saved_view.plot_windows).plot_kind == "RuO2IVSweepPlot"

    Browser._save_project_view_if_changed!(state)
    loaded_after_ui_save = Browser._load_project_view(root_path)
    @test Browser._project_view_to_toml(loaded_after_ui_save) ==
          Browser._project_view_to_toml(saved_view)

end
