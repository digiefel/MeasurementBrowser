using MeasurementBrowser
using Test

const Browser = MeasurementBrowser.Browser

# Persisted plot choices include item-kind and label names; on load they resolve to the
# internal RegisteredPlot identity for that plot.
const ProjectViewIVPlot = MeasurementBrowser.RegisteredPlot{:iv_sweep,:ProjectViewIVPlot}
const ProjectViewTLMPlot = MeasurementBrowser.RegisteredPlot{:iv_sweep,:ProjectViewTLMPlot}

@testset "project view state" begin
    root_path = mktempdir()
    @test basename(Browser._project_view_file_path(root_path)) == "measurementbrowser.toml"
    default_view = Browser._load_project_view(root_path)
    @test default_view.project == ""
    @test default_view.main_plot.id == "main"
    @test default_view.main_plot.live == true

    view = Browser.PersistedProjectView(
        project="ProjectViewTest",
        tree=Browser.PersistedTreeView(
            expanded=["chip/device"],
            selected=["chip/device"],
            filter="tlm",
        ),
        items=Browser.PersistedItemsView(
            selected=["item-1", "item-2"],
            filter="298K",
        ),
        plot_kinds=Dict("iv_sweep" => "iv_sweep::ProjectViewIVPlot"),
        main_plot=Browser.PersistedPlotView(
            id="main",
            title="Plot Area",
            plot_kind="iv_sweep::ProjectViewTLMPlot",
            live=true,
            items=["item-1"],
        ),
        plot_windows=[
            Browser.PersistedPlotView(
                id="plot_1",
                title="Detached",
                plot_kind="iv_sweep::ProjectViewIVPlot",
                live=false,
                items=["item-2"],
            ),
        ],
    )

    Browser._save_project_view(root_path, view)
    @test isfile(joinpath(root_path, "measurementbrowser.toml"))
    loaded = Browser._load_project_view(root_path)
    @test Browser._project_view_to_toml(loaded) == Browser._project_view_to_toml(view)

    legacy_root_path = mktempdir()
    write(joinpath(legacy_root_path, "measurementbrowser.toml"), """
    project = "ProjectViewTest"

    [measurements]
    selected = ["item-1"]
    filter = "298K"

    [main_plot]
    id = "main"
    title = "Plot Area"
    live = true
    measurements = ["item-1"]

    [[plot_windows]]
    id = "plot_1"
    title = "Detached"
    live = false
    measurements = ["item-2"]
    """)
    legacy_loaded = Browser._load_project_view(legacy_root_path)
    @test isempty(legacy_loaded.items.selected)
    @test legacy_loaded.items.filter == ""
    @test isempty(legacy_loaded.main_plot.items)
    @test isempty(only(legacy_loaded.plot_windows).items)

    toml_data = Browser._project_view_to_toml(view)
    @test Set(keys(toml_data)) == Set([
        "project",
        "tree",
        "items",
        "plot_kinds",
        "main_plot",
        "plot_windows",
    ])

    parsed = Browser._project_view_from_toml(Browser.PersistedProjectView, Dict{String,Any}(
        "project" => "ProjectViewTest",
        "tree" => Dict{String,Any}(
            "expanded" => ["chip/device"],
            "selected" => ["chip/device"],
            "filter" => "tlm",
        ),
        "items" => Dict{String,Any}(
            "selected" => ["item-1"],
            "filter" => "298K",
        ),
        "plot_kinds" => Dict{String,Any}("iv_sweep" => "iv_sweep::ProjectViewIVPlot"),
        "main_plot" => Dict{String,Any}(
            "id" => "main",
            "title" => "Plot Area",
            "plot_kind" => "iv_sweep::ProjectViewTLMPlot",
            "live" => true,
            "items" => ["item-1"],
        ),
        "plot_windows" => Any[
            Dict{String,Any}(
                "id" => "plot_1",
                "title" => "Detached",
                "plot_kind" => "iv_sweep::ProjectViewIVPlot",
                "live" => false,
                "items" => ["item-2"],
            ),
        ],
    ))

    @test_throws ErrorException Browser._project_view_from_toml(Browser.PersistedProjectView, Dict{String,Any}(
        "project" => "ProjectViewTest",
        "tree" => Dict{String,Any}(
            "expanded" => Any["chip/device", 1],
            "selected" => ["chip/device"],
            "filter" => "tlm",
        ),
        "items" => Dict{String,Any}(
            "selected" => ["item-1"],
            "filter" => "298K",
        ),
        "plot_kinds" => Dict{String,Any}("iv_sweep" => "iv_sweep::ProjectViewIVPlot"),
        "main_plot" => Dict{String,Any}(
            "id" => "main",
            "title" => "Plot Area",
            "plot_kind" => "iv_sweep::ProjectViewTLMPlot",
            "live" => true,
            "items" => ["item-1"],
        ),
        "plot_windows" => Any[],
    ))

    @test parsed.project == "ProjectViewTest"
    @test parsed.main_plot.plot_kind == "iv_sweep::ProjectViewTLMPlot"
    @test only(parsed.plot_windows).items == ["item-2"]

    project = MeasurementBrowser.define_project("ProjectViewTest")
    source = test_source(project, root_path)
    workspace = MeasurementBrowser.Workspace.Workspace(project, source)
    item_1 = MeasurementBrowser.ItemRecord(;
        source_item_id="file-1",
        source_item_path=joinpath(root_path, "item-1.csv"),
        id="item-1",
        item_label="Item 1",
        kind=:iv_sweep,
        collection=["chip", "device-1"],
    )
    item_2 = MeasurementBrowser.ItemRecord(;
        source_item_id="file-2",
        source_item_path=joinpath(root_path, "item-2.csv"),
        id="item-2",
        item_label="Item 2",
        kind=:iv_sweep,
        collection=["chip", "device-2"],
    )
    hierarchy = MeasurementBrowser.Hierarchy(root_path, true)
    MeasurementBrowser.insert_item!(hierarchy, item_1)
    MeasurementBrowser.insert_item!(hierarchy, item_2)
    MeasurementBrowser.Workspace.replace_item_index!(workspace, hierarchy)
    state = Browser.BrowserState(workspace=workspace)

    bad_plot_view = Browser.PersistedProjectView(
        project="ProjectViewTest",
        plot_kinds=Dict("iv_sweep" => "tlm4p"),
        main_plot=Browser.PersistedPlotView(
            id="main",
            title="Plot Area",
            plot_kind="tlm4p",
            live=true,
        ),
        plot_windows=[
            Browser.PersistedPlotView(
                id="plot_1",
                title="Detached",
                plot_kind="tlm4p",
                live=false,
            ),
        ],
    )
    Browser._apply_project_view!(state, bad_plot_view)
    @test isempty(state.plots.kind_by_item)
    @test state.plots.main.plot_kind === nothing
    @test only(state.plots.windows).plot_kind === nothing

    Browser._apply_project_view!(state, view)
    @test workspace.selection.collection_paths == ["chip/device"]
    @test workspace.selection.item_ids == ["item-1", "item-2"]
    plots = state.plots
    @test plots.kind_by_item == Dict(:iv_sweep => ProjectViewIVPlot)
    @test plots.main.plot_kind == ProjectViewTLMPlot
    @test plots.main.item_ids == ["item-1"]
    @test only(plots.windows).item_ids == ["item-2"]

    saved_view = Browser._project_view_from_browser(state)
    @test saved_view.project == "ProjectViewTest"
    @test saved_view.items.selected == ["item-1", "item-2"]
    @test saved_view.main_plot.plot_kind == "iv_sweep::ProjectViewTLMPlot"

    Browser._save_project_view_if_changed!(state)
    loaded_after_ui_save = Browser._load_project_view(root_path)
    @test Browser._project_view_to_toml(loaded_after_ui_save) ==
          Browser._project_view_to_toml(saved_view)

end
