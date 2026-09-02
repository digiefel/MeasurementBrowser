# Scaling sweep for the hot operations that must NOT grow with the item count.
#
#   julia --project=bench bench/scaling.jl [n1,n2,...]
#
# Default sizes (or MB_BENCH_SCALING_SIZES): 250,500,1000. The committed harness uses those.
# Each size is N aliases of bench/templates/kind1.csv.
#
#   status_refresh    refresh_status! every GUI frame
#   items_panel       per-frame gather+sort of the items panel
#   metadata_publish  source-metadata reconciliation
#   scan_build_per_item  cold scan + publication wall time / N
#
# Fits time ~ N^exponent (log-log). Exponent ~0 is flat, ~1 is linear, ~2 is quadratic.
# The sweep should include a size at or above the cache buffer row ceiling (~1000 items).
# Writes nothing on its own; bench/run.jl records exponent, R², and ms at each size.

using DataBrowserRecipes: Project, define_project, register_item!
using DataBrowserCore.Workspace:
    close_workspace!,
    open_workspace,
    reconcile_source_metadata_cache!,
    refresh_status!,
    wait_workspace_idle!,
    workspace_busy
using DataBrowserGUI.Browser: BrowserState, _items_of_selected_collections
using DataFrames: DataFrame
using Printf: @printf
using Statistics: mean, median

if !@isdefined(BenchSource)
    include(joinpath(@__DIR__, "custom_data_source.jl"))
end

const TEMPLATE = joinpath(@__DIR__, "templates", "kind1.csv")
const SCALING_OPS = (
    "scan_build_per_item",
    "status_refresh",
    "workspace_busy",
    "items_panel",
    "metadata_publish",
)

function scaling_sizes()
    if abspath(PROGRAM_FILE) == @__FILE__
        isempty(ARGS) || return parse.(Int, split(ARGS[1], ","))
    end
    return parse.(Int, split(get(ENV, "MB_BENCH_SCALING_SIZES", "250,500,1000"), ","))
end

"""Warm once, then median wall time in seconds. Repeats a fast `f` until the batch is ~1 ms."""
function _median_seconds(f; n::Int=5)
    sink = Ref{Any}(nothing)
    sink[] = f()
    t0 = time_ns()
    sink[] = f()
    dt = time_ns() - t0
    batch = dt < 1_000_000 ? Int(cld(1_000_000, max(dt, UInt64(1)))) : 1
    samples = Float64[]
    sizehint!(samples, n)
    for _ in 1:n
        t0 = time_ns()
        for _ in 1:batch
            sink[] = f()
        end
        push!(samples, (time_ns() - t0) / 1e9 / batch)
    end
    return median(samples)
end

"""One collection, one trivial item per file — the smallest project that still exercises the scan."""
function scaling_project(name::AbstractString)::Project
    project = define_project(name)
    register_item!(
        project,
        :row;
        detect=file -> endswith(file.filename, ".csv"),
        read=file -> DataFrame(v=[1]),
        collection=(_data, _metadata) -> ["batch"],
        label=(_data, metadata) -> metadata[:filename],
    )
    return project
end

"""Open a settled `n`-item workspace, run `probe(ws, build_seconds)`, and clean up."""
function with_workspace(probe::Function, n::Int)
    name = "scaling_" * string(n) * "_" * string(time_ns())
    cache_dir = joinpath(first(DEPOT_PATH), "databrowser", name)
    source = BenchSource(
        abspath(TEMPLATE) * "_" * name,
        alias_file(TEMPLATE, n, i -> "item_$i.csv"),
    )
    try
        workspace = nothing
        build_seconds = @elapsed begin
            workspace = open_workspace(scaling_project(name), source)
            wait_workspace_idle!(workspace; timeout=600)
        end
        result = probe(workspace, build_seconds)
        close_workspace!(workspace)
        return result
    finally
        rm(cache_dir; force=true, recursive=true)
    end
end

"""Median seconds of each hot operation on an `n`-item workspace."""
function measure(n::Int)::Dict{String,Float64}
    return with_workspace(n) do ws, build_seconds
        ws.selection.collection_ids = [only(values(ws.index.collections.records)).id]
        state = BrowserState(workspace=ws)
        Dict(
            "scan_build_per_item" => build_seconds / n,
            "status_refresh" => _median_seconds(() -> refresh_status!(ws)),
            "workspace_busy" => _median_seconds(() -> workspace_busy(ws)),
            "items_panel" => _median_seconds(() -> _items_of_selected_collections(state)),
            "metadata_publish" => _median_seconds(
                () -> reconcile_source_metadata_cache!(ws; refresh_hierarchy=true),
            ),
        )
    end
end

"""Least-squares fit of log(y) = intercept + exponent·log(x); `exponent` is the scaling order."""
function power_law(xs::Vector{<:Real}, ys::Vector{<:Real})
    all(iszero, ys) && return (exponent=0.0, r2=1.0)
    any(<=(0), ys) && error("power_law needs positive times, got $ys")
    lx = log.(xs)
    ly = log.(ys)
    x̄ = mean(lx)
    ȳ = mean(ly)
    exponent = sum((lx .- x̄) .* (ly .- ȳ)) / sum((lx .- x̄) .^ 2)
    intercept = ȳ - exponent * x̄
    ss_res = sum((ly .- (intercept .+ exponent .* lx)) .^ 2)
    ss_tot = sum((ly .- ȳ) .^ 2)
    return (exponent=exponent, r2=(ss_tot == 0 ? 1.0 : 1 - ss_res / ss_tot))
end

"""Run the sweep. Returns sizes, per-op milliseconds, and log-log fits."""
function run_scaling(sizes::Vector{Int}=scaling_sizes())
    println("scaling sweep over item counts: ", join(sizes, ", "))
    measure(min(first(sizes), 16))
    ms_by_op = Dict(op => Float64[] for op in SCALING_OPS)
    for n in sizes
        times = measure(n)
        for (op, seconds) in times
            push!(ms_by_op[op], 1e3 * seconds)
        end
        @printf(
            "  n=%-6d  build=%.4f ms/item  refresh=%.4f ms  busy=%.4f ms  items=%.4f ms  metadata=%.4f ms\n",
            n,
            1e3 * times["scan_build_per_item"],
            1e3 * times["status_refresh"],
            1e3 * times["workspace_busy"],
            1e3 * times["items_panel"],
            1e3 * times["metadata_publish"],
        )
    end
    fits = Dict{String,NamedTuple{(:exponent, :r2),Tuple{Float64,Float64}}}()
    println("\noperation          exponent   r2")
    for op in SCALING_OPS
        fit = power_law(Float64.(sizes), ms_by_op[op])
        fits[op] = fit
        @printf("%-18s %6.2f   %.3f\n", op, fit.exponent, fit.r2)
    end
    return (sizes=sizes, ops=collect(String, SCALING_OPS), fits=fits, ms_by_op=ms_by_op)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_scaling()
end
