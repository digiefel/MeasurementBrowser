# Measure monolithic fresh publish substeps (HEAD scan_source path).
# Run in jmux after interactive_profile is loaded and Operations.jl marks are revised in:
#   include(".../measure_fresh_publish_substeps.jl")
#   DataBrowserFreshPublishMeasure.main()

module DataBrowserFreshPublishMeasure

using Dates: format, now
using DataBrowser
using DataBrowserCore.Workspace: close_workspace!, open_workspace, workspace_status
using DataBrowserProfiling: DebugTimings, with_debug_timings, write_debug_timings

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const CACHE_PATH = joinpath(
    first(DEPOT_PATH), "databrowser", "DataBrowserProfilingRuO2", "cache.duckdb")
const DATA_ROOT = get(
    ENV,
    "RUO2_DATA_ROOT",
    "/Users/davide/Library/CloudStorage/OneDrive-LundUniversity/projects/Borg/202501_RuO2test/electricaldata",
)

function duckdb_counts(path::AbstractString)
    ispath(path) || return (; exists=false, bytes=0, raw="missing")
    bytes = filesize(path)
    sql = """
    SELECT 'source_items' AS t, COUNT(*) AS n FROM source_items
    UNION ALL SELECT 'items', COUNT(*) FROM items
    UNION ALL SELECT 'collections', COUNT(*) FROM collections
    UNION ALL SELECT 'meta', COUNT(*) FROM meta;
    """
    raw = try
        strip(read(`duckdb $path -c $sql`, String))
    catch err
        "duckdb_cli_error: $(sprint(showerror, err))"
    end
    return (; exists=true, bytes, raw)
end

function status_snapshot(workspace)
    st = workspace_status(workspace)
    c = st.counts
    return (;
        busy=st.busy,
        label=st.label,
        detail=st.detail,
        sources_found=c.sources_found,
        sources_pending=c.sources_pending,
        cached_sources=c.cache.cached_sources,
        interpreted=c.cache.interpreted_items,
        scan_state=workspace.scan.state,
        discovered=workspace.scan.discovered[],
        cache_state=workspace.cache_state,
        cache_operation=workspace.cache.operation,
        index_items=length(workspace.index.items),
        index_collections=length(workspace.index.collections.records),
        index_source_isnothing=workspace.index.source === nothing,
        work_nodes=length(workspace.work.nodes),
    )
end

function main(; profile_project=nothing, discover_timeout::Float64=180.0, after_publish::Float64=3.0)
    outdir = joinpath(
        REPO_ROOT,
        "bench",
        "results",
        string(format(now(), "yyyymmdd-HHMMSS-sss"), "-ruo2-fresh-publish-substeps"),
    )
    mkpath(outdir)
    report = joinpath(outdir, "freshness_and_substeps.txt")
    io = open(report, "w")
    println(io, "outdir=$outdir")
    println(io, "cache_path=$CACHE_PATH")
    println(io, "threads=$(Threads.nthreads())")
    println(io, "commit=$(strip(read(`git -C $REPO_ROOT rev-parse --short HEAD`, String)))")
    flush(io)

    before_bytes = ispath(CACHE_PATH) ? filesize(CACHE_PATH) : 0
    println(io, "\n=== BEFORE explicit wipe ===")
    println(io, (; exists=ispath(CACHE_PATH), bytes=before_bytes))
    flush(io)

    for candidate in (CACHE_PATH, CACHE_PATH * ".wal")
        ispath(candidate) && rm(candidate; force=true)
    end
    println(io, "\n=== AFTER explicit rm of cache.duckdb(+wal) ===")
    println(io, (; exists=ispath(CACHE_PATH), bytes=0))
    flush(io)

    if profile_project === nothing
        isdefined(Main, :DataBrowserInteractiveProfile) || error(
            "Load DataBrowserInteractiveProfile first, or pass profile_project=")
        profile_project = Main.DataBrowserInteractiveProfile.project()
    end

    timings = DebugTimings()
    workspace = nothing
    with_debug_timings(timings) do
        open_seconds = @elapsed workspace = open_workspace(
            profile_project,
            DATA_ROOT;
            metadata_file="device_info.txt",
            rebuild=true,
            cache=true,
            background_processing=false,
        )
        t0 = status_snapshot(workspace)
        println(io, "\n=== t≈0 immediately after open_workspace(rebuild=true) ===")
        println(io, "open_workspace_seconds=$open_seconds")
        println(io, t0)
        flush(io)

        sleep(0.05)
        println(io, "\n=== t≈0.05s ===")
        println(io, status_snapshot(workspace))
        flush(io)

        # Wait until monolithic publish finished: scan leaves :discovering, or we see work queued
        # for ~all discovered sources after discovery completes (~6002).
        deadline = time() + discover_timeout
        last_print = 0.0
        while time() < deadline
            st = status_snapshot(workspace)
            if time() - last_print > 5.0
                println(io, "\n=== progress $(round(time() - (deadline - discover_timeout); digits=1))s ===")
                println(io, st)
                flush(io)
                last_print = time()
            end
            if st.scan_state === :idle
                println(io, "\n=== scan idle ===")
                println(io, st)
                flush(io)
                break
            end
            # Monolithic path: discovery walks to ~6002 then one big publish; after publish,
            # work_nodes jumps and scan eventually goes idle.
            if st.discovered >= 5900 && st.work_nodes >= 1000
                println(io, "\n=== post-publish signal (discovered>=5900, work_nodes>=1000) ===")
                println(io, st)
                flush(io)
                break
            end
            sleep(0.1)
        end

        sleep(after_publish)
        println(io, "\n=== before close_workspace ===")
        println(io, status_snapshot(workspace))
        flush(io)
        close_workspace!(workspace)
        workspace = nothing
    end

    # write_debug_timings treats the path as an output directory
    write_debug_timings(outdir, timings)
    println(io, "\n=== debug_timings.txt ===")
    write(io, read(joinpath(outdir, "debug_timings.txt"), String))
    println(io)

    after = duckdb_counts(CACHE_PATH)
    println(io, "\n=== AFTER close (DuckDB table counts) ===")
    println(io, after)
    close(io)
    println("Wrote $outdir")
    println(read(report, String))
    return outdir
end

end # module
