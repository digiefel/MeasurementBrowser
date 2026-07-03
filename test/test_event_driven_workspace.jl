using Test
using CSV
using DataFrames: DataFrame, nrow
using GLMakie: Axis, Figure
using MeasurementBrowser

# These tests intentionally never poll: the engine must scan, interpret, process, and publish on
# completion events alone, with blocking callers woken by publications.

"""A project whose process and stats callbacks leave observable marks, with a tunable delay."""
function _event_driven_project(name::AbstractString; process_delay::Real=0.0)
    project = define_project("EventDriven_$name")
    register_item!(
        project,
        :table;
        detect=file -> endswith(lowercase(file.filename), ".csv"),
        read=file -> CSV.read(file.filepath, DataFrame),
        entries=(file, data) -> [DataItem(
            kind=:table,
            collection=["run"],
            label=file.filename,
            data=data,
        )],
        process=function (item)
            process_delay > 0 && sleep(process_delay)
            data = DataFrame(item.data)
            data.sum = data.x .+ data.y
            return DataItem(item, data)
        end,
        stats=item -> Dict{Symbol,Any}(:rows => nrow(item.data)),
    )
    register_collection_stat!(project;
        kinds=[:table],
        compute_stats=items -> Dict{Symbol,Any}(:members => length(items)),
    )
    register_plot!(project, :table; label="Sum",
        setup=(_ws, _items) -> (figure = Figure(); Axis(figure[1, 1]); figure),
        draw=(_ws, _items, _figure) -> nothing,
    )
    return project
end

function _event_driven_dir(count::Int)::String
    dir = mktempdir()
    for index in 1:count
        write_test_source(joinpath(dir, "sample$index.csv"), index)
    end
    return dir
end

@testset "headless scan publishes without polling" begin
    dir = _event_driven_dir(3)
    project = _event_driven_project(basename(dir))
    workspace = MeasurementBrowser.open_workspace(
        project, test_source(project, dir); background_processing=true)
    try
        MeasurementBrowser.Workspace.wait_workspace_idle!(workspace; timeout=30)
        @test !MeasurementBrowser.Workspace.engine_work_running(workspace)
        @test workspace.scan.state in (:done, :unchanged)
        @test length(workspace.index.items) == 3
        @test isempty(workspace.index.analysis_errors)
        # Background processing ran item stats and collection stats purely on completion events.
        @test all(stats -> stats[:rows] == 2, values(workspace.index.item_stats))
        @test length(workspace.index.item_stats) == 3
        node = workspace.index.hierarchy.index[("run",)]
        @test node.stats[:members] == 3
    finally
        MeasurementBrowser.close_workspace!(workspace)
    end
end

@testset "blocking materialize resolves from worker publishes" begin
    dir = _event_driven_dir(2)
    project = _event_driven_project(basename(dir))
    workspace = MeasurementBrowser.open_workspace(
        project, test_source(project, dir); background_processing=false)
    try
        MeasurementBrowser.Workspace.wait_workspace_idle!(workspace; timeout=30)
        records = sort!(collect(values(workspace.index.items)); by=record -> record.id)
        items = MeasurementBrowser.Workspace.materialize_items(workspace, records)
        @test length(items) == 2
        @test all(item -> hasproperty(MeasurementBrowser.item_data(item), :sum), items)
    finally
        MeasurementBrowser.close_workspace!(workspace)
    end
end

@testset "wait_workspace_idle! timeout returns while work runs" begin
    dir = _event_driven_dir(2)
    project = _event_driven_project(basename(dir); process_delay=1.0)
    workspace = MeasurementBrowser.open_workspace(
        project, test_source(project, dir); background_processing=true)
    try
        started = time()
        MeasurementBrowser.Workspace.wait_workspace_idle!(workspace; timeout=0.05)
        @test time() - started < 5
        # A second wait is not short-circuited by the first call's expired timer.
        MeasurementBrowser.Workspace.wait_workspace_idle!(workspace; timeout=30)
        @test !MeasurementBrowser.Workspace.engine_work_running(workspace)
    finally
        MeasurementBrowser.close_workspace!(workspace)
    end
end

@testset "close_workspace! never strands blocked waiters" begin
    dir = _event_driven_dir(4)
    project = _event_driven_project(basename(dir); process_delay=0.5)
    workspace = MeasurementBrowser.open_workspace(
        project, test_source(project, dir); background_processing=false)
    MeasurementBrowser.Workspace.wait_workspace_idle!(workspace; timeout=30)
    records = collect(values(workspace.index.items))
    waiters = [Threads.@spawn begin
        try
            (:ok, MeasurementBrowser.Workspace.materialize_items(workspace, [record]))
        catch error
            (:failed, error)
        end
    end for record in records]
    sleep(0.1)
    MeasurementBrowser.close_workspace!(workspace)
    for waiter in waiters
        outcome, _ = fetch(waiter)::Tuple
        @test outcome in (:ok, :failed)
    end
end

@testset "full API produces a plot without polling" begin
    dir = _event_driven_dir(2)
    project = _event_driven_project(basename(dir))
    workspace = MeasurementBrowser.open_workspace(
        project, test_source(project, dir); background_processing=false)
    try
        MeasurementBrowser.Workspace.wait_workspace_idle!(workspace; timeout=30)
        records = collect(values(workspace.index.items))
        MeasurementBrowser.select_items!(workspace, records)
        items = MeasurementBrowser.Workspace.materialize_items(workspace, records)
        plot_kind = only(MeasurementBrowser.registered_plot_kinds(project, :table))
        figure = MeasurementBrowser.setup_plot(workspace, plot_kind, items)
        @test figure isa Figure
        @test MeasurementBrowser.plot_data!(workspace, plot_kind, items, figure) === nothing
    finally
        MeasurementBrowser.close_workspace!(workspace)
    end
end
