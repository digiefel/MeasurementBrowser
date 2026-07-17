using DataBrowser
using DataBrowserGUI
using DataBrowserPlots
using Test

const Browser = DataBrowserGUI.Browser

const ProjectViewIVPlot = RegisteredPlot{:iv_sweep,:ProjectViewIVPlot}
const ProjectViewTLMPlot = RegisteredPlot{:iv_sweep,:ProjectViewTLMPlot}

const _PLOTS_EXTENSION_VIEW = Dict{String,Any}(
    "plot_kinds" => Dict("iv_sweep" => "iv_sweep::ProjectViewIVPlot"),
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
)

@testset "project view state" begin
    root_path = mktempdir()
    @test basename(Browser._project_view_file_path(root_path)) == "databrowser.toml"
    default_view = Browser._load_project_view(root_path)
    @test default_view.project == ""
    @test isempty(default_view.extensions)

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
        extensions=Dict("PlotsExtension" => copy(_PLOTS_EXTENSION_VIEW)),
    )

    Browser._save_project_view(root_path, view)
    @test isfile(joinpath(root_path, "databrowser.toml"))
    loaded = Browser._load_project_view(root_path)
    @test Browser._project_view_to_toml(loaded) == Browser._project_view_to_toml(view)

    toml_data = Browser._project_view_to_toml(view)
    @test Set(keys(toml_data)) == Set([
        "project",
        "tree",
        "items",
        "extensions",
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
        "extensions" => Dict{String,Any}(
            "PlotsExtension" => copy(_PLOTS_EXTENSION_VIEW),
        ),
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
        "extensions" => Dict{String,Any}(
            "PlotsExtension" => copy(_PLOTS_EXTENSION_VIEW),
        ),
    ))

    @test parsed.project == "ProjectViewTest"
    @test parsed.extensions["PlotsExtension"]["main_plot"]["plot_kind"] ==
        "iv_sweep::ProjectViewTLMPlot"

    project = DataBrowser.define_project("ProjectViewTest")
    source = test_source(project, root_path)
    workspace = DataBrowserCore.Workspace.Workspace(project, source)
    collections = DataBrowserAPI.ItemIndex.CollectionIndex(root_path)
    path_1 = DataBrowserAPI.ItemIndex.RegisteredCollection.(["chip", "device-1"])
    path_2 = DataBrowserAPI.ItemIndex.RegisteredCollection.(["chip", "device-2"])
    key_1 = DataBrowserAPI.ItemIndex.resolve_collection_path!(
        collections, DataBrowserAPI.ItemIndex.collection_inputs(path_1))
    key_2 = DataBrowserAPI.ItemIndex.resolve_collection_path!(
        collections, DataBrowserAPI.ItemIndex.collection_inputs(path_2))
    item_1 = DataBrowserAPI.ItemIndex.ItemRecord(;
        source_item_key=DataBrowserCache.source_item_key!(
            workspace.cache.db, "file-1"; mint=true),
        source_item_path=joinpath(root_path, "item-1.csv"),
        id="item-1",
        label="Item 1",
        kind=:iv_sweep,
        collection_key=key_1,
    )
    item_2 = DataBrowserAPI.ItemIndex.ItemRecord(;
        source_item_key=DataBrowserCache.source_item_key!(
            workspace.cache.db, "file-2"; mint=true),
        source_item_path=joinpath(root_path, "item-2.csv"),
        id="item-2",
        label="Item 2",
        kind=:iv_sweep,
        collection_key=key_2,
    )
    DataBrowserAPI.ItemIndex.insert_item!(collections, item_1.id, key_1)
    DataBrowserAPI.ItemIndex.insert_item!(collections, item_2.id, key_2)
    DataBrowserCore.Workspace.replace_item_index!(workspace, collections, [item_1, item_2])
    state = Browser.BrowserState(workspace=workspace)
    state.extensions = Browser._instantiate_extensions()

    bad_plot_view = Browser.PersistedProjectView(
        project="ProjectViewTest",
        extensions=Dict(
            "PlotsExtension" => Dict{String,Any}(
                "plot_kinds" => Dict("iv_sweep" => "tlm4p"),
                "main_plot" => Dict{String,Any}(
                    "id" => "main",
                    "title" => "Plot Area",
                    "plot_kind" => "tlm4p",
                    "live" => true,
                    "items" => String[],
                ),
                "plot_windows" => Any[
                    Dict{String,Any}(
                        "id" => "plot_1",
                        "title" => "Detached",
                        "plot_kind" => "tlm4p",
                        "live" => false,
                        "items" => String[],
                    ),
                ],
            ),
        ),
    )
    Browser._apply_project_view!(state, bad_plot_view)
    plots = plots_extension(state).plots
    @test isempty(plots.kind_by_item)
    @test plots.main.plot_kind === nothing
    @test only(plots.windows).plot_kind === nothing

    Browser._apply_project_view!(state, view)
    @test workspace.selection.collection_ids == ["chip/device"]
    @test workspace.selection.item_ids == ["item-1", "item-2"]
    plots = plots_extension(state).plots
    @test plots.kind_by_item == Dict(:iv_sweep => ProjectViewIVPlot)
    @test plots.main.plot_kind == ProjectViewTLMPlot
    @test plots.main.item_ids == ["item-1"]
    @test only(plots.windows).item_ids == ["item-2"]

    saved_view = Browser._project_view_from_browser(state)
    @test saved_view.project == "ProjectViewTest"
    @test saved_view.items.selected == ["item-1", "item-2"]
    @test saved_view.extensions["PlotsExtension"]["main_plot"]["plot_kind"] ==
        "iv_sweep::ProjectViewTLMPlot"

    Browser._save_project_view_if_changed!(state)
    loaded_after_ui_save = Browser._load_project_view(root_path)
    @test Browser._project_view_to_toml(loaded_after_ui_save) ==
          Browser._project_view_to_toml(saved_view)
end

@testset "root items participate in browser selection" begin
    mktempdir() do root_path
        project = DataBrowser.define_project("RootItemSelection_$(basename(root_path))")
        workspace = DataBrowserCore.Workspace.Workspace(
            project,
            test_source(project, root_path),
        )
        try
            collections = DataBrowserAPI.ItemIndex.CollectionIndex(root_path)
            root_item = DataBrowserAPI.ItemIndex.ItemRecord(;
                source_item_key=DataBrowserCache.source_item_key!(
                    workspace.cache.db, "root-file"; mint=true),
                source_item_path=joinpath(root_path, "root.csv"),
                id="root-item",
                label="Root Item",
                kind=:table,
                collection_key=nothing,
            )
            DataBrowserAPI.ItemIndex.insert_item!(collections, root_item.id, nothing)
            DataBrowserCore.Workspace.replace_item_index!(
                workspace,
                collections,
                [root_item],
            )
            state = Browser.BrowserState(workspace=workspace)

            workspace.selection.collection_ids = [Browser.ROOT_COLLECTION_SELECTION_ID]
            workspace.selection.item_ids = [root_item.id]
            selected_collections, selected_items, selected_path =
                Browser._project_visible_selection(state)
            @test isempty(selected_collections)
            @test selected_items == [root_item]
            @test isempty(selected_path)
            @test Browser._items_of_selected_collections(state) == [root_item]

            empty!(workspace.selection.collection_ids)
            empty!(workspace.selection.item_ids)
            @test Browser.select_source_item!(state, "root-file")
            @test workspace.selection.collection_ids ==
                [Browser.ROOT_COLLECTION_SELECTION_ID]
            @test workspace.selection.item_ids == [root_item.id]
        finally
            DataBrowser.close_workspace!(workspace)
        end
    end
end
