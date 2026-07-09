# Scaling sweep for the hot operations that must NOT grow with the item count.
#
#   julia --project=bench bench/scaling.jl [n1,n2,...]
#
# For each item count it builds one workspace, then uses BenchmarkTools to time three real
# operations — nothing internal is fabricated, only public/engine functions are called:
#
#   status_refresh    refresh_status! every GUI frame (workspace_busy + optional status rebuild)
#   items_panel       the per-frame gather+sort the items panel runs (Browser._items_of_selected_collections)
#   metadata_publish  the per-publish source-metadata diff a scan runs once per source item
#
# It fits time ~ N^exponent (a linear regression in log-log space) and reports the exponent and R²
# per operation. Exponent ~0 is flat, ~1 is linear, ~2 is quadratic. The metadata call is O(N) per
# publish today and a scan makes N of them, so its per-call growth means an O(N^2) scan.
#
# The sweep must cross the cache buffer row ceiling (~1000 items): below it the metadata reads are
# served from the in-memory write buffer and look flat; above it they hit disk and the per-call cost
# grows sharply. Keep the largest size well past the ceiling, or the metadata cost is understated.
# The largest size also costs the most to *build* (the scan it profiles is the O(N^2) one), so a
# wide sweep is minutes, not seconds — pass a smaller list while iterating.
#
# Results land in bench/results/<timestamp>-scaling/scaling.csv.

using DataBrowserAPI
using DataBrowserCore.ItemIndex: DataItem
using DataBrowserCore.Workspace
using DataBrowserGUI.Browser
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
        entries=(file, data) -> [
            DataItem(
                kind=:row,
                collection=["batch"],
                label=file.filename,
                metadata=Dict{Symbol,Any}(),
                data=data,
            ),
        ],
    )
    return project
end

"""Open a settled `n`-item workspace, run `probe(ws)`, and clean up its data + cache."""
function with_workspace(probe::Function, n::Int)
    dir = mktempdir()
    for index in 1:n
        write(joinpath(dir, "item_$index.csv"), "v\n1\n")
    end
    name = "scaling_" * basename(dir)  # unique; the cache is keyed by project name
    cache_dir = joinpath(first(DEPOT_PATH), "measurementbrowser", name)
    try
        workspace = open_workspace(scaling_project(name), dir)
        wait_workspace_idle!(workspace; timeout=600)
        result = probe(workspace)
        close_workspace!(workspace)
        return result
    finally
        rm(dir; force=true, recursive=true)
        rm(cache_dir; force=true, recursive=true)
    end
end

"""Minimum seconds (via BenchmarkTools) of each hot operation on an `n`-item workspace."""
function measure(n::Int)::Dict{String,Float64}
    return with_workspace(n) do ws
        ws.selection.collection_paths = ["batch"]  # select the one collection, as the GUI would
        state = BrowserState(workspace=ws)
        Dict(
            "status_refresh" => @belapsed(refresh_status!($ws)),
            "workspace_busy" => @belapsed(workspace_busy($ws)),
            "items_panel" => @belapsed(_items_of_selected_collections($state)),
            "metadata_publish" => @belapsed(reconcile_source_metadata_cache!($ws)),
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
    ms_by_op = Dict(op => Float64[] for op in (
        "status_refresh", "workspace_busy", "items_panel", "metadata_publish"))
    for n in SIZES
        times = measure(n)
        for (op, seconds) in times
            push!(ms_by_op[op], 1e3 * seconds)
        end
        @printf("  n=%-6d  refresh=%.4f ms  busy=%.4f ms  items=%.4f ms  metadata=%.4f ms\n",
            n, 1e3 * times["status_refresh"], 1e3 * times["workspace_busy"],
            1e3 * times["items_panel"], 1e3 * times["metadata_publish"])
    end

    outdir = joinpath(@__DIR__, "results", format(now(), "yyyymmdd-HHMMSS") * "-scaling")
    mkpath(outdir)
    csv = joinpath(outdir, "scaling.csv")
    open(csv, "w") do io
        println(io, "operation,exponent,r2," * join(("ms_n$n" for n in SIZES), ","))
        println("\noperation          exponent   r2")
        for op in ("status_refresh", "workspace_busy", "items_panel", "metadata_publish")
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
