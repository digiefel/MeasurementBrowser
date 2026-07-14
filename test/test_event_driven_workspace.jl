using Test
using CSV
using DataFrames: DataFrame, nrow
using GLMakie: Axis, Figure
using DataBrowser
using DataBrowser: display_label

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
        collection=(_data, _metadata) -> ["run"],
        label=(_data, metadata) -> metadata[:filename],
        process=function (data, _metadata)
            process_delay > 0 && sleep(process_delay)
            processed = DataFrame(data)
            processed.sum = data.x .+ data.y
            return processed
        end,
        analyze=(data, _metadata) -> Dict{Symbol,Any}(:rows => nrow(data)),
    )
    register_collection_analysis!(project, :table;
        analyze=(data, _metadata) -> Dict{Symbol,Any}(:members => length(data)),
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
    workspace = DataBrowser.open_workspace(
        project, test_source(project, dir); background_processing=true)
    try
        DataBrowserCore.Workspace.wait_workspace_idle!(workspace; timeout=30)
        @test !DataBrowserCore.Workspace.engine_work_running(workspace)
        @test workspace.scan.state in (:done, :unchanged)
        @test length(workspace.index.items) == 3
        @test isempty(workspace.index.analysis_errors)
        # Background processing ran item and collection analysis purely on completion events.
        @test all(md -> md[:rows] == 2, values(workspace.index.item_metadata))
        @test length(workspace.index.item_metadata) == 3
        collection_record = _registered_collection_record(workspace.index.collections, "run")
        @test collection_record.analysis[:members] == 3
    finally
        DataBrowser.close_workspace!(workspace)
    end
end

@testset "blocking materialize resolves from worker publishes" begin
    dir = _event_driven_dir(2)
    project = _event_driven_project(basename(dir))
    workspace = DataBrowser.open_workspace(
        project, test_source(project, dir); background_processing=false)
    try
        DataBrowserCore.Workspace.wait_workspace_idle!(workspace; timeout=30)
        records = sort!(collect(values(workspace.index.items)); by=record -> record.id)
        items = DataBrowserCore.Workspace.materialize_items(workspace, records)
        @test length(items) == 2
        @test all(item -> hasproperty(DataBrowser.item_data(item), :sum), items)
    finally
        DataBrowser.close_workspace!(workspace)
    end
end

@testset "status counts refresh after blocking materialize" begin
    dir = _event_driven_dir(1)
    project = _event_driven_project(basename(dir))
    workspace = DataBrowser.open_workspace(
        project, test_source(project, dir); background_processing=false)
    try
        DataBrowserCore.Workspace.wait_workspace_idle!(workspace; timeout=30)
        DataBrowserCore.Workspace.refresh_status!(workspace)
        @test workspace.status.counts.cache.processed == 0

        records = collect(values(workspace.index.items))
        DataBrowserCore.Workspace.materialize_items(workspace, records)
        DataBrowserCore.Workspace.wait_workspace_idle!(workspace; timeout=30)
        DataBrowserCore.Workspace.refresh_status!(workspace)
        @test workspace.status.counts.cache.processed == 1
        @test workspace.status.counts.cache.analyzed == 1
    finally
        DataBrowser.close_workspace!(workspace)
    end
end

@testset "wait_workspace_idle! timeout returns while work runs" begin
    dir = _event_driven_dir(2)
    project = _event_driven_project(basename(dir); process_delay=1.0)
    workspace = DataBrowser.open_workspace(
        project, test_source(project, dir); background_processing=true)
    try
        started = time()
        DataBrowserCore.Workspace.wait_workspace_idle!(workspace; timeout=0.05)
        @test time() - started < 5
        # A second wait is not short-circuited by the first call's expired timer.
        DataBrowserCore.Workspace.wait_workspace_idle!(workspace; timeout=30)
        @test !DataBrowserCore.Workspace.engine_work_running(workspace)
    finally
        DataBrowser.close_workspace!(workspace)
    end
end

@testset "close_workspace! never strands blocked waiters" begin
    dir = _event_driven_dir(4)
    project = _event_driven_project(basename(dir); process_delay=0.5)
    workspace = DataBrowser.open_workspace(
        project, test_source(project, dir); background_processing=false)
    DataBrowserCore.Workspace.wait_workspace_idle!(workspace; timeout=30)
    records = collect(values(workspace.index.items))
    waiters = [Threads.@spawn begin
        try
            (:ok, DataBrowserCore.Workspace.materialize_items(workspace, [record]))
        catch error
            (:failed, error)
        end
    end for record in records]
    sleep(0.1)
    DataBrowser.close_workspace!(workspace)
    for waiter in waiters
        outcome, _ = fetch(waiter)::Tuple
        @test outcome in (:ok, :failed)
    end
end

@testset "labels never observe torn parameters during publishes" begin
    dir = mktempdir()
    for index in 1:120
        write_test_source(joinpath(dir, "cycle$index.csv"), index)
    end
    project = define_project("EventDriven_$(basename(dir))")
    register_item!(project, :cycle;
        detect=file -> endswith(file.filename, ".csv"),
        read=file -> CSV.read(file.filepath, DataFrame),
        entries=(data, metadata) -> [(
            data=data,
            metadata=Dict{Symbol,Any}(
                :cycle => parse(Int, match(r"\d+", metadata[:filename]).match)),
        )],
        collection=(_data, _metadata) -> ["run"],
        analyze=(data, _metadata) -> Dict{Symbol,Any}(:rows => nrow(data)),
        label=(_data, metadata) -> "cycle $(metadata[:cycle])",
    )
    workspace = DataBrowser.open_workspace(
        project, test_source(project, dir); background_processing=true)
    try
        # Mimic the GUI: render labels continuously while stats publishes rebuild the hierarchy.
        failures = Base.Threads.Atomic{Int}(0)
        reader = Threads.@spawn while DataBrowserCore.Workspace.engine_work_running(workspace)
            for record in collect(values(workspace.index.items))
                try
                    display_label(project, record)
                catch
                    Base.Threads.atomic_add!(failures, 1)
                end
            end
            yield()
        end
        DataBrowserCore.Workspace.wait_workspace_idle!(workspace; timeout=60)
        wait(reader)
        @test failures[] == 0
        @test length(workspace.index.items) == 120
        labels = Set(display_label(project, record)
                     for record in values(workspace.index.items))
        @test labels == Set("cycle $index" for index in 1:120)
    finally
        DataBrowser.close_workspace!(workspace)
    end
end

"""Poll one predicate until it holds or the deadline passes; returns the final value."""
function _settles(predicate::Function; timeout::Real=5.0)::Bool
    deadline = time() + timeout
    while !predicate() && time() < deadline
        sleep(0.02)
    end
    return predicate()
end

"""Whether any processing callback is currently running (not merely queued)."""
function _process_running(workspace)::Bool
    graph = workspace.work
    return lock(graph.lock) do
        any(
            node -> node.key.kind === DataBrowserCore.Workspace.ITEM_PROCESS &&
                node.state === :running,
            values(graph.nodes),
        )
    end
end

@testset "aborting profile prep settles scan state" begin
    dir = _event_driven_dir(3)
    project = _event_driven_project(basename(dir); process_delay=1.0)
    workspace = DataBrowser.open_workspace(
        project, test_source(project, dir);
        background_processing=true, profile_internal=true)
    try
        # A running (uncancellable) processing callback holds the workspace busy, so the profile
        # prep is still pending — not already restarted — when the stop below aborts it.
        @test _settles(() -> _process_running(workspace))
        DataBrowserCore.Workspace.start_internal_profile!(workspace)
        @test workspace.profiler.state === :preparing
        DataBrowserCore.Workspace.stop_internal_profile!(workspace)
        @test workspace.profiler.state === :idle
        @test _settles(() -> !DataBrowserCore.Workspace.source_scan_running(workspace))
        @test workspace.scan.state != :canceling
    finally
        DataBrowser.close_workspace!(workspace)
    end
end

@testset "scan cancel settles promptly" begin
    dir = _event_driven_dir(3)
    project = _event_driven_project(basename(dir))
    workspace = DataBrowser.open_workspace(
        project, test_source(project, dir); background_processing=false)
    try
        DataBrowserCore.Workspace.cancel_scan!(workspace)
        @test _settles(() -> workspace.scan.state != :canceling)
        @test workspace.scan.state in (:canceled, :done, :unchanged)
    finally
        DataBrowser.close_workspace!(workspace)
    end
end

@testset "late cancel is a no-op on a settled scan" begin
    dir = _event_driven_dir(2)
    project = _event_driven_project(basename(dir))
    workspace = DataBrowser.open_workspace(
        project, test_source(project, dir); background_processing=false)
    try
        DataBrowserCore.Workspace.wait_workspace_idle!(workspace; timeout=30)
        settled = workspace.scan.state
        @test settled in (:done, :unchanged)
        DataBrowserCore.Workspace.cancel_scan!(workspace)
        @test workspace.scan.state === settled
        @test !DataBrowserCore.Workspace.source_scan_running(workspace)
    finally
        DataBrowser.close_workspace!(workspace)
    end
end

@testset "workspace_status survives concurrent analysis_errors mutations" begin
    dir = _event_driven_dir(1)
    project = _event_driven_project(basename(dir))
    workspace = DataBrowser.open_workspace(
        project, test_source(project, dir); background_processing=true)
    try
        DataBrowserCore.Workspace.wait_workspace_idle!(workspace; timeout=30)
        lock(workspace.publish_lock) do
            workspace.index.analysis_errors["a"] = "one"
            workspace.index.analysis_errors["b"] = "two"
        end
        status = DataBrowserCore.Workspace.workspace_status(workspace)
        @test status.errors == ["a" => "one", "b" => "two"]

        failures = Base.Threads.Atomic{Int}(0)
        done = Base.Threads.Atomic{Bool}(false)
        writer = Threads.@spawn while !done[]
            i = rand(1:100)
            lock(workspace.publish_lock) do
                workspace.index.analysis_errors["k$i"] = "v$i"
            end
            lock(workspace.work.lock) do
                for key in collect(keys(workspace.index.analysis_errors))[1:min(3, end)]
                    delete!(workspace.index.analysis_errors, key)
                end
            end
            yield()
        end
        reader = Threads.@spawn while !done[]
            try
                DataBrowserCore.Workspace.workspace_status(workspace)
            catch
                Base.Threads.atomic_add!(failures, 1)
            end
            yield()
        end
        sleep(0.5)
        done[] = true
        fetch(writer)
        fetch(reader)
        @test failures[] == 0
    finally
        DataBrowser.close_workspace!(workspace)
    end
end

@testset "full API produces a plot without polling" begin
    dir = _event_driven_dir(2)
    project = _event_driven_project(basename(dir))
    workspace = DataBrowser.open_workspace(
        project, test_source(project, dir); background_processing=false)
    try
        DataBrowserCore.Workspace.wait_workspace_idle!(workspace; timeout=30)
        records = collect(values(workspace.index.items))
        DataBrowser.select_items!(workspace, records)
        items = DataBrowserCore.Workspace.materialize_items(workspace, records)
        plot_kind = only(DataBrowser.registered_plot_kinds(project, :table))
        figure = DataBrowser.setup_plot(workspace, plot_kind, items)
        @test figure isa Figure
        @test DataBrowser.plot_data!(workspace, plot_kind, items, figure) === nothing
    finally
        DataBrowser.close_workspace!(workspace)
    end
end
