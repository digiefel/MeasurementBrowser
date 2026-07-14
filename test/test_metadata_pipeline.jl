using DataBrowser
using DataFrames: DataFrame, nrow
using Test

const MBP = DataBrowser

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
        entries=function (data, source_metadata)
            cycle = parse(Int, match(r"\d+", source_metadata[:filename]).match)
            return [(data=data, metadata=Dict{Symbol,Any}(:cycle => cycle))]
        end,
        collection=(data, metadata) -> ["fatigue"],
        label=(data, metadata) -> "cycle $(metadata[:cycle])",
        process=(data, metadata) -> DataFrame(count=data.count),
        analyze=(data, metadata) -> Dict{Symbol,Any}(:events => data.count[1]),
    )
    MBP.register_collection_analysis!(project, :cycle;
        process=function (data, metadata)
            fail_process && error("collection process failed")
            order = sortperm(metadata; by=values -> values[:cycle])
            running = 0
            outputs = copy(data)
            for position in order
                running += data[position].count[1]
                outputs[position] = DataFrame(
                    count=data[position].count,
                    cumulative=[running],
                )
            end
            return outputs
        end,
        analyze=(data, metadata) -> Dict{Symbol,Any}(:members => length(data)),
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
            collect(values(workspace.index.items)); by=record -> record.metadata[:cycle])
        loaded = DataBrowserCore.Workspace.materialize_items(workspace, records)
        loaded_by_id = Dict(record.id => item for (record, item) in zip(records, loaded))
        delivered = [DataBrowserCore.Workspace.delivered_metadata(
            workspace, record, workspace.index.collections) for record in records]

        @test [values[:events] for values in delivered] == [2, 3, 5]
        @test [MBP.item_data(loaded_by_id[record.id]).cumulative[1] for record in records] ==
            [2, 5, 10]

        # Collection analyze landed on the collection node only, not on items.
        collection_record = _registered_collection_record(
            workspace.index.collections, "fatigue")
        @test collection_record.analysis[:members] == 3
        @test all(values -> !haskey(values, :members), delivered)

        conflict = "metadata :events expected Int64, got Bool; value dropped"
        @test_logs (:warn, r"Workspace metadata conflict") begin
            DataBrowserCore.Workspace.publish_metadata_conflicts!(
                workspace, "manual-conflict", :test, [conflict])
        end
        @test workspace.index.analysis_errors["manual-conflict"] == conflict
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
        entries=function (data, source_metadata)
            cycle = parse(Int, match(r"\d+", source_metadata[:filename]).match)
            return [(data=data, metadata=Dict{Symbol,Any}(:cycle => cycle))]
        end,
        collection=(data, metadata) -> ["fatigue"],
        analyze=(data, metadata) -> Dict{Symbol,Any}(:events => data.count[1]),
    )
    workspace = MBP.open_workspace(
        project, test_source(project, dir); background_processing=true)
    try
        wait_workspace_idle!(workspace)
        record = only(collect(values(workspace.index.items)))
        @test workspace.index.item_metadata[record.id][:events] == 4

        # Re-register analyze with a different key, then re-run the item's analyze stage: its output
        # replaces the previous computed layer wholesale, so stale keys never accumulate.
        MBP.register_item!(project, :cycle;
            detect=file -> endswith(file.filename, ".csv"),
            read=file -> DataFrame(count=[parse(Int, strip(read(file.filepath, String)))]),
            entries=function (data, source_metadata)
                cycle = parse(Int, match(r"\d+", source_metadata[:filename]).match)
                return [(data=data, metadata=Dict{Symbol,Any}(:cycle => cycle))]
            end,
            collection=(data, metadata) -> ["fatigue"],
            analyze=(data, metadata) -> Dict{Symbol,Any}(:doubled => data.count[1] * 2),
        )
        recomputed = DataBrowserCore.Workspace.run_item_analysis(
            workspace, workspace.index.collections, record)
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

@testset "memory cache collection-process delivery after idle" begin
    dir = _fatigue_dir([2, 3, 5])
    project = _pipeline_project(basename(dir))
    workspace = MBP.open_workspace(
        project, test_source(project, dir); background_processing=true, cache=false)
    try
        wait_workspace_idle!(workspace)
        @test workspace.cache.db isa DataBrowserCache.MemoryCacheDB
        @test isempty(workspace.index.analysis_errors)

        records = sort(
            collect(values(workspace.index.items)); by=record -> record.metadata[:cycle])
        loaded = DataBrowserCore.Workspace.materialize_items(workspace, records)
        @test [MBP.item_data(item).cumulative[1] for item in loaded] == [2, 5, 10]
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
            collect(values(workspace.index.items)); by=record -> record.metadata[:cycle])
        loaded = DataBrowserCore.Workspace.materialize_items(workspace, records)
        # The fold failed, so members deliver their :processed payloads (no cumulative column).
        @test length(loaded) == 2
        @test all(item -> !hasproperty(MBP.item_data(item), :cumulative), loaded)
        # The collection-process error is surfaced on the collection key.
        key = _registered_collection_key(workspace.index.collections, "fatigue")
        @test haskey(workspace.index.analysis_errors, key)
    finally
        MBP.close_workspace!(workspace)
    end
end
