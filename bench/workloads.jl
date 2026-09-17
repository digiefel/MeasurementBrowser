using DataBrowserAPI
import DataBrowserAPI as API
import DataBrowserCore.Workspace as Workspace
using Statistics: median

const BENCH_TIMEOUT = 120.0
const PAYLOAD_ITEMS = 400
const PAYLOAD_ROWS = 10_000
const INDEX_ITEMS = 10_000
const SELECTION_SIZE = 10
const SELECTION_REPEATS = 5

"""Prepared user value: the benchmark times engine work, not generating or analyzing this data."""
struct BenchItem <: API.AbstractDataItem
    identity::String
    data::Any
    metadata::Dict{Symbol,Any}
end
API.id(item::BenchItem) = item.identity
API.label(item::BenchItem) = item.identity
API.item_data(item::BenchItem) = item.data
API.metadata(item::BenchItem) = item.metadata
API.reconstruct(::Type{BenchItem}, id::AbstractString, data, metadata::Dict) = BenchItem(id, data, metadata)

struct BenchSourceItem <: API.AbstractDataSourceItem
    item::BenchItem
end
API.id(item::BenchSourceItem) = API.id(item.item)
API.label(item::BenchSourceItem) = API.id(item)
API.fingerprint(::BenchSourceItem) = 1

struct PreparedSource <: API.AbstractDataSource
    items::Vector{BenchSourceItem}
end
API.source_id(::PreparedSource) = "prepared-benchmark"
API.source_label(::PreparedSource) = "Prepared benchmark data"
API.source_items(source::PreparedSource; kwargs...) = source.items

"""Callback counts are correctness checks, not performance measurements."""
Base.@kwdef struct BenchProject <: API.AbstractProject
    reads::Threads.Atomic{Int} = Threads.Atomic{Int}(0)
    processes::Threads.Atomic{Int} = Threads.Atomic{Int}(0)
    analyses::Threads.Atomic{Int} = Threads.Atomic{Int}(0)
end
function API.read(project::BenchProject, ::PreparedSource, item::BenchSourceItem)
    Threads.atomic_add!(project.reads, 1)
    return item.item
end
function API.process(project::BenchProject, item::BenchItem)
    Threads.atomic_add!(project.processes, 1)
    return item
end
function API.analyze(project::BenchProject, ::BenchItem)
    Threads.atomic_add!(project.analyses, 1)
    return Dict{Symbol,Any}(:analyzed => true)
end
callback_counts(project::BenchProject) = (project.reads[], project.processes[], project.analyses[])

"""Fixed four-column Float64 tables, prepared before timing and shared by source items."""
function prepared_source(n, rows)
    data = (time=Float64.(1:rows), voltage=sin.(Float64.(1:rows)),
        current=cos.(Float64.(1:rows)), charge=Float64.(0:rows-1) ./ max(rows-1, 1))
    items = [BenchSourceItem(BenchItem("item_$(lpad(i, 6, '0'))", data,
        Dict{Symbol,Any}(:number => i))) for i in 1:n]
    return PreparedSource(items)
end

"""Wait for engine completion and reject timeouts/errors. Close separately waits for disk writes."""
function settle!(ws)
    Workspace.wait_workspace_idle!(ws; timeout=BENCH_TIMEOUT)
    status = Workspace.workspace_status(ws)
    status.busy && error("Benchmark workspace timed out")
    status.level === :error && error("Benchmark workspace failed: $status")
    isempty(status.errors) || error("Benchmark stage failures: $(status.errors)")
    return nothing
end

"""Materialize exact public IDs; check data outside the timed call."""
function materialize(ws, ids)
    Workspace.select_items!(ws, ids)
    return Workspace.materialize_items(ws, ids)
end
function validate_items(items, ids, rows)
    API.id.(items) == ids || error("Materialization lost or reordered items")
    for item in items
        data = API.item_data(item)
        length(data.time) == rows || error("Incorrect payload row count")
        (data.time[1], data.time[end]) == (1.0, Float64(rows)) || error("Incorrect payload values")
    end
    return nothing
end
function timed_materialize(ws, ids, rows)
    sample = @timed materialize(ws, ids)
    validate_items(sample.value, ids, rows)
    return sample.time
end

"""Give each scenario a fresh data cache without changing the package depot used for loading code."""
function with_benchmark_cache(f)
    mktempdir() do depot
        pushfirst!(DEPOT_PATH, depot)
        try
            return f()
        finally
            popfirst!(DEPOT_PATH)
        end
    end
end

"""Per-item indexing overhead, including saving the index at close; callbacks only return prepared data."""
function measure_index(source)
    with_benchmark_cache() do
        ws = nothing
        project = BenchProject()
        try
            elapsed = @elapsed begin
                ws = Workspace.open_workspace(project, source)
                settle!(ws)
                sort(Workspace.query_items(ws)) == sort(API.id.(source.items)) || error("Incorrect indexed items")
                Workspace.close_workspace!(ws)
            end
            callback_counts(project) == (length(source.items), 0, 0) || error("Unexpected indexing callbacks")
            return elapsed * 1e6 / length(source.items)
        finally
            ws === nothing || Workspace.close_workspace!(ws)
        end
    end
end

"""Cache throughput and materialization latency, with identity processing and constant analysis."""
function measure_cache(source)
    with_benchmark_cache() do
        project = BenchProject()
        ids = API.id.(source.items)
        selected, remaining = ids[1:SELECTION_SIZE], ids[SELECTION_SIZE+1:end]
        ws = nothing
        bulk = nothing
        try
            # Persist the selected items before the write probe, so concurrent reads are cache hits.
            ws = Workspace.open_workspace(project, source)
            settle!(ws)
            timed_materialize(ws, selected, PAYLOAD_ROWS)
            settle!(ws)
            Workspace.close_workspace!(ws)
            ws = Workspace.open_workspace(project, source)
            settle!(ws)
            started = time_ns()
            bulk = Threads.@spawn begin
                result = Workspace.materialize_items(ws, remaining)
                settle!(ws)
                result
            end
            # Submit the same reads alongside every bulk write; either operation may finish first.
            concurrent_samples = [timed_materialize(ws, selected, PAYLOAD_ROWS) for _ in 1:SELECTION_REPEATS]
            written_items = fetch(bulk)
            bulk = nothing
            project.processes[] == length(ids) || error("Not every item was processed exactly once")
            project.analyses[] == length(ids) || error("Not every item was analyzed exactly once")
            Workspace.close_workspace!(ws)
            ws = nothing
            write_s = (time_ns() - started) / 1e9
            validate_items(written_items, remaining, PAYLOAD_ROWS)
            written_items = nothing

            saved_counts = callback_counts(project)
            reopen_s = @elapsed begin
                ws = Workspace.open_workspace(project, source)
                settle!(ws)
            end
            sort(Workspace.query_items(ws)) == sort(ids) || error("Restore lost items")
            read_s = timed_materialize(ws, ids, PAYLOAD_ROWS)
            cached_samples = [timed_materialize(ws, selected, PAYLOAD_ROWS) for _ in 1:SELECTION_REPEATS]
            settle!(ws)
            callback_counts(project) == saved_counts || error("Cache hits reran user callbacks")
            bytes_per_item = PAYLOAD_ROWS * 4 * sizeof(Float64)
            return (write_mib_s=length(remaining) * bytes_per_item / 1024^2 / write_s,
                read_mib_s=length(ids) * bytes_per_item / 1024^2 / read_s,
                cached_materialize_ms=median(cached_samples) * 1e3,
                concurrent_materialize_ms=median(concurrent_samples) * 1e3,
                reopen_ms=reopen_s * 1e3)
        finally
            ws === nothing || Workspace.close_workspace!(ws)
            bulk === nothing || wait(bulk; throw=false)
        end
    end
end

"""Warm each operation once, then measure identical workloads with fresh caches."""
function run_benchmark(repeats::Int)
    payload_source = prepared_source(PAYLOAD_ITEMS, PAYLOAD_ROWS)
    index_source = prepared_source(INDEX_ITEMS, 1)
    measure_cache(payload_source)
    measure_index(index_source)
    samples = NamedTuple[]
    for _ in 1:repeats
        GC.gc()
        cache = measure_cache(payload_source)
        overhead = measure_index(index_source)
        push!(samples, merge(cache, (index_us_per_item=overhead,)))
    end
    return samples
end
