using MeasurementBrowser
using DataFrames: DataFrame, nrow
using Test

const MBP = MeasurementBrowser

"""
Build a fatigue-style project: one item per cycle CSV, item analyze counts rows, the collection
process folds a cumulative total back onto each member (data and metadata), and collection analyze
attaches a members count to the collection node only.
"""
function _pipeline_project(name::AbstractString; fail_process::Bool=false)
    project = MBP.define_project("Pipeline_$(name)")
    MBP.register_item!(project, :cycle;
        detect=file -> endswith(file.filename, ".csv"),
        read=file -> DataFrame(count=[parse(Int, strip(read(file.filepath, String)))]),
        entries=function (file, data)
            cycle = parse(Int, match(r"\d+", file.filename).match)
            return [MBP.DataItem(
                kind=:cycle,
                collection=["fatigue"],
                label="cycle $cycle",
                metadata=Dict{Symbol,Any}(:cycle => cycle),
                data=data)]
        end,
        process=item -> MBP.DataItem(item, DataFrame(count=item.data.count)),
        analyze=item -> Dict{Symbol,Any}(:events => item.data.count[1]),
    )
    MBP.register_collection_analysis!(project, :cycle;
        process=function (items)
            fail_process && error("collection process failed")
            ordered = sort(items; by=it -> MBP.metadata(it)[:cycle])
            running = 0
            outputs = MBP.AbstractDataItem[]
            for it in ordered
                running += it.data.count[1]
                data = DataFrame(count=it.data.count, cumulative=[running])
                push!(outputs, MBP.DataItem(
                    kind=MBP.kind(it),
                    collection=MBP.collection(it),
                    label=MBP.item_label(it),
                    metadata=merge(MBP.metadata(it), Dict{Symbol,Any}(:cumulative => running)),
                    data=data,
                    id=MBP.id(it)))
            end
            return outputs
        end,
        analyze=items -> Dict{Symbol,Any}(:members => length(items)),
    )
    return project
end

function _fatigue_dir(counts::Vector{Int})::String
    dir = mktempdir()
    for (cycle, count) in enumerate(counts)
        write(joinpath(dir, "cycle$cycle.csv"), string(count))
    end
    return dir
end

@testset "cumulative-history pipeline delivers both layers" begin
    dir = _fatigue_dir([2, 3, 5])
    project = _pipeline_project(basename(dir))
    workspace = MBP.open_workspace(
        project, test_source(project, dir); background_processing=true)
    try
        wait_workspace_idle!(workspace)
        @test workspace.scan.state in (:done, :unchanged)
        @test isempty(workspace.index.analysis_errors)

        records = sort(
            workspace.index.hierarchy.all_items; by=record -> record.metadata[:cycle])
        loaded = MBP.Workspace.materialize_items(workspace, records)
        loaded = sort(loaded; by=item -> item.metadata[:cycle])

        # Item analyze put :events on each item; collection process pushed cumulative totals down
        # and rewrote the data column — delivered items carry both.
        @test [item.metadata[:events] for item in loaded] == [2, 3, 5]
        @test [item.metadata[:cumulative] for item in loaded] == [2, 5, 10]
        @test [item.data.cumulative[1] for item in loaded] == [2, 5, 10]

        # Collection analyze landed on the collection node only, not on items.
        node = workspace.index.hierarchy.index[("fatigue",)]
        @test node.analysis[:members] == 3
        @test all(item -> !haskey(item.metadata, :members), loaded)
    finally
        MBP.close_workspace!(workspace)
    end
end

@testset "re-running item analysis replaces its layer" begin
    dir = _fatigue_dir([4])
    # A project with no collection recipe: the item's delivered layer is exactly its analyze output,
    # so re-running analyze must replace it wholesale (stale keys never accumulate).
    project = MBP.define_project("Relayer_$(basename(dir))")
    MBP.register_item!(project, :cycle;
        detect=file -> endswith(file.filename, ".csv"),
        read=file -> DataFrame(count=[parse(Int, strip(read(file.filepath, String)))]),
        entries=function (file, data)
            cycle = parse(Int, match(r"\d+", file.filename).match)
            return [MBP.DataItem(kind=:cycle, collection=["fatigue"],
                metadata=Dict{Symbol,Any}(:cycle => cycle), data=data)]
        end,
        analyze=item -> Dict{Symbol,Any}(:events => item.data.count[1]),
    )
    workspace = MBP.open_workspace(
        project, test_source(project, dir); background_processing=true)
    try
        wait_workspace_idle!(workspace)
        record = only(workspace.index.hierarchy.all_items)
        @test workspace.index.item_metadata[record.id][:events] == 4

        # Re-register analyze with a different key, then re-run the item's analyze stage: its output
        # replaces the previous computed layer wholesale, so stale keys never accumulate.
        MBP.register_item!(project, :cycle;
            detect=file -> endswith(file.filename, ".csv"),
            read=file -> DataFrame(count=[parse(Int, strip(read(file.filepath, String)))]),
            entries=function (file, data)
                cycle = parse(Int, match(r"\d+", file.filename).match)
                return [MBP.DataItem(kind=:cycle, collection=["fatigue"],
                    metadata=Dict{Symbol,Any}(:cycle => cycle), data=data)]
            end,
            analyze=item -> Dict{Symbol,Any}(:doubled => item.data.count[1] * 2),
        )
        recomputed = MBP.Workspace.run_item_analysis(workspace, record)
        lock(workspace.publish_lock) do
            workspace.index.item_metadata[record.id] =
                Dict{Symbol,Any}(k => v for (k, v) in recomputed)
        end
        @test workspace.index.item_metadata[record.id][:doubled] == 8
        @test !haskey(workspace.index.item_metadata[record.id], :events)
    finally
        MBP.close_workspace!(workspace)
    end
end

@testset "collection-process failure delivers processed payloads with the error surfaced" begin
    dir = _fatigue_dir([1, 2])
    project = _pipeline_project(basename(dir); fail_process=true)
    workspace = MBP.open_workspace(
        project, test_source(project, dir); background_processing=true)
    try
        wait_workspace_idle!(workspace)
        records = sort(
            workspace.index.hierarchy.all_items; by=record -> record.metadata[:cycle])
        loaded = MBP.Workspace.materialize_items(workspace, records)
        # The fold failed, so members deliver their :processed payloads (no cumulative column).
        @test length(loaded) == 2
        @test all(item -> !hasproperty(item.data, :cumulative), loaded)
        # The collection-process error is surfaced on the collection key.
        key = MBP.ItemIndex.collection_path_key(["fatigue"])
        @test haskey(workspace.index.analysis_errors, key)
    finally
        MBP.close_workspace!(workspace)
    end
end
