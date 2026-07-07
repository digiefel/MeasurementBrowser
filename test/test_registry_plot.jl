using MeasurementBrowser
using MeasurementBrowser: items_for_file, setup_plot, plot_data!
using CSV
using DataFrames: DataFrame, nrow
using GLMakie: Figure, Axis
using Test

const MB = MeasurementBrowser

@testset "registry plot bridge" begin
    dir = mktempdir()
    write(joinpath(dir, "m.csv"), "x,y\n1,2\n3,4\n")

    project = MB.define_project("PlotBridge")
    MB.register_item!(
        project,
        :table;
        detect=file -> endswith(file.filename, ".csv"),
        read=file -> DataFrame(CSV.File(file.filepath)),
        entries=(file, data) -> [DataItem(
            kind=:table,
            collection=["dev"],
            label=file.filename,
            data=data,
        )],
    )
    drew = Ref(0)
    drew_rows = Ref(0)
    MB.register_plot!(
        project,
        :table;
        label="Table Plot",
        setup=(workspace, items) -> Figure(),
        draw=function (workspace, items, figure)
            drew[] = length(items)
            # Each item is an AbstractDataItem carrying its own materialized data as `item.data`,
            # and answering the contract accessors.
            @test all(it -> it isa MB.AbstractDataItem, items)
            @test all(it -> MB.item_data(it) === it.data, items)
            @test all(it -> MB.kind(it) === :table, items)
            @test all(it -> !isempty(MB.collection(it)), items)
            drew_rows[] = sum(nrow(item.data) for item in items)
            Axis(figure[1, 1])
            nothing
        end,
    )
    MB.register_plot!(
        project,
        :table;
        label="Table Summary",
        setup=(workspace, items) -> Figure(),
        draw=(workspace, items, figure) -> nothing,
    )
    MB.register_plot!(
        project,
        :table;
        label="Table Summary",
        setup=(workspace, items) -> Figure(),
        draw=(workspace, items, figure) -> (drew[] = -1; nothing),
    )

    # Identity helpers used by the GUI and persistence.
    table_plot = MB.RegisteredPlot{:table,Symbol("Table Plot")}
    table_summary = MB.RegisteredPlot{:table,Symbol("Table Summary")}
    @test MB.registered_plot_kinds(project, :table) == [table_plot, table_summary]
    @test MB.registered_plot_kinds(project, :missing) == Type{<:MB.PlotKind}[]
    @test MB.plot_kind_label(project, table_plot) == "Table Plot"
    @test MB.plot_kind_name(table_plot) == "table::Table Plot"
    @test MB.plot_kind_from_name("table::Table Plot") === table_plot
    @test MB.plot_kind_from_name("table") === nothing

    # The bridge runs the registered setup/draw callbacks via the engine's plot dispatch.
    items = items_for_file(project, joinpath(dir, "m.csv"))
    workspace = MB.Workspace.Workspace(project, test_source(project, dir))
    for item in items
        workspace.index.items[item.id] = item
    end
    ids = [item.id for item in items]
    @test MB.select_items!(workspace, items) === nothing
    @test workspace.selection.item_ids == ids
    @test MB.select_items!(workspace, ids[1:1]) === nothing
    @test workspace.selection.item_ids == ids[1:1]
    loaded_items = MB.Workspace.materialize_items(workspace, items)
    @test MB.select_items!(workspace, loaded_items) === nothing
    @test workspace.selection.item_ids == ids
    @test MB.select_items!(workspace, []) === nothing
    @test isempty(workspace.selection.item_ids)

    figure = setup_plot(workspace, table_plot, items)
    @test figure isa Figure
    @test plot_data!(workspace, table_plot, items, figure) === nothing
    @test drew[] == length(items)
    @test drew_rows[] == 2   # the 2-row fixture flowed through as item.data
    figure = setup_plot(workspace, table_summary, items)
    @test plot_data!(workspace, table_summary, items, figure) === nothing
    @test drew[] == -1

    # An unregistered kind errors clearly rather than dispatching nowhere.
    @test_throws KeyError setup_plot(
        workspace,
        MB.RegisteredPlot{:missing,Symbol("Missing")},
        items,
    )
end
