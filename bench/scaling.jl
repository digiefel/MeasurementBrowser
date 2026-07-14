# Scaling sweep for the hot operations that must NOT grow with the item count.
#
#   julia --project=bench bench/scaling.jl [n1,n2,...]
#
# For each item count it builds one workspace, then uses BenchmarkTools to time real operations —
# nothing internal is fabricated, only public/engine functions are called:
#
#   status_refresh    refresh_status! every GUI frame (workspace_busy + optional status rebuild)
#   items_panel       the per-frame gather+sort the items panel runs (Browser._items_of_selected_collections)
#   metadata_publish  an explicit source-metadata reconciliation
#   scan_build_per_item  cold scan + publication wall time normalized by item count
#
# It fits time ~ N^exponent (a linear regression in log-log space) and reports the exponent and R²
# per operation. Exponent ~0 is flat, ~1 is linear, and ~2 is quadratic.
#
# The sweep must cross the cache buffer row ceiling (~1000 items): below it metadata reads can be
# served from the in-memory write buffer and look artificially flat. The normalized scan-build row
# catches cumulative publication costs: exponent ~0 means stable throughput, while exponent ~1
# means total build time is quadratic.
#
# Results land in bench/results/<timestamp>-scaling/scaling.csv.

using DataBrowserAPI: Project, define_project, register_item!
using DataBrowserCore.Workspace:
    close_workspace!,
    open_workspace,
    reconcile_source_metadata_cache!,
    refresh_status!,
    wait_workspace_idle!,
    workspace_busy
using DataBrowserGUI.Browser: BrowserState, _items_of_selected_collections
using BenchmarkTools
using DataFrames: DataFrame
using Dates: format, now
using Printf: @printf, @sprintf
using Statistics: mean

const SIZES = isempty(ARGS) ? [500, 1000, 2000, 4000] : parse.(Int, split(ARGS[1], ","))

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
    dir = mktempdir()
    for index in 1:n
        write(joinpath(dir, "item_$index.csv"), "v\n1\n")
    end
    name = "scaling_" * basename(dir)  # unique; the cache is keyed by project name
    cache_dir = joinpath(first(DEPOT_PATH), "databrowser", name)
    try
        workspace = nothing
        build_seconds = @elapsed begin
            workspace = open_workspace(scaling_project(name), dir)
            wait_workspace_idle!(workspace; timeout=600)
        end
        result = probe(workspace, build_seconds)
        close_workspace!(workspace)
        return result
    finally
        rm(dir; force=true, recursive=true)
        rm(cache_dir; force=true, recursive=true)
    end
end

"""Minimum seconds (via BenchmarkTools) of each hot operation on an `n`-item workspace."""
function measure(n::Int)::Dict{String,Float64}
    return with_workspace(n) do ws, build_seconds
        ws.selection.collection_ids = [only(values(ws.index.collections.records)).id]
        state = BrowserState(workspace=ws)
        Dict(
            "scan_build_per_item" => build_seconds / n,
            "status_refresh" => @belapsed(refresh_status!($ws)),
            "workspace_busy" => @belapsed(workspace_busy($ws)),
            "items_panel" => @belapsed(_items_of_selected_collections($state)),
            "metadata_publish" => @belapsed(
                reconcile_source_metadata_cache!($ws; refresh_hierarchy=true)),
        )
    end
end

"""Least-squares fit of log(y) = intercept + exponent·log(x); `exponent` is the scaling order."""
function power_law(xs::Vector{<:Real}, ys::Vector{<:Real})
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

function main()
    println("scaling sweep over item counts: ", join(SIZES, ", "))
    # Compile the workspace and panel paths before the sweep so the first size measures runtime
    # throughput rather than one-time Julia compilation.
    measure(min(first(SIZES), 16))
    ms_by_op = Dict(op => Float64[] for op in (
        "scan_build_per_item", "status_refresh", "workspace_busy", "items_panel",
        "metadata_publish"))
    for n in SIZES
        times = measure(n)
        for (op, seconds) in times
            push!(ms_by_op[op], 1e3 * seconds)
        end
        @printf("  n=%-6d  build=%.4f ms/item  refresh=%.4f ms  busy=%.4f ms  items=%.4f ms  metadata=%.4f ms\n",
            n, 1e3 * times["scan_build_per_item"], 1e3 * times["status_refresh"],
            1e3 * times["workspace_busy"],
            1e3 * times["items_panel"], 1e3 * times["metadata_publish"])
    end

    outdir = joinpath(@__DIR__, "results", format(now(), "yyyymmdd-HHMMSS") * "-scaling")
    mkpath(outdir)
    csv = joinpath(outdir, "scaling.csv")
    open(csv, "w") do io
        println(io, "operation,exponent,r2," * join(("ms_n$n" for n in SIZES), ","))
        println("\noperation          exponent   r2")
        for op in ("scan_build_per_item", "status_refresh", "workspace_busy", "items_panel",
                "metadata_publish")
            fit = power_law(Float64.(SIZES), ms_by_op[op])
            @printf("%-18s %6.2f   %.3f\n", op, fit.exponent, fit.r2)
            println(io, join(
                [op, @sprintf("%.4f", fit.exponent), @sprintf("%.4f", fit.r2),
                    (@sprintf("%.4f", ms) for ms in ms_by_op[op])...],
                ","))
        end
    end
    println("\nwrote ", csv)
end

main()
