# Dump first-second index/tree evidence on a fresh-cache open.
# Usage (after interactive_profile loaded, or standalone with project()):
#   include(".../dump_fresh_first_second.jl")
#   DataBrowserFreshFirstSecond.dump!(; gui=true)

module DataBrowserFreshFirstSecond

using Dates: Dates, format, now, unix2datetime
using Printf
using DataBrowser
using DataBrowserCore.Workspace: close_workspace!, open_workspace, workspace_status

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const CACHE_PATH = joinpath(
    first(DEPOT_PATH), "databrowser", "DataBrowserProfilingRuO2", "cache.duckdb")
const DATA_ROOT = get(
    ENV,
    "RUO2_DATA_ROOT",
    "/Users/davide/Library/CloudStorage/OneDrive-LundUniversity/projects/Borg/202501_RuO2test/electricaldata",
)

function _cache_stat()
    path = CACHE_PATH
    wal = path * ".wal"
    return (;
        path,
        exists=ispath(path),
        bytes=ispath(path) ? filesize(path) : 0,
        mtime=ispath(path) ? string(unix2datetime(mtime(path))) : "missing",
        wal_exists=ispath(wal),
        wal_bytes=ispath(wal) ? filesize(wal) : 0,
    )
end

function _sample_labels(workspace; limit::Int=20)
    records = collect(values(workspace.index.collections.records))
    sort!(records; by=r -> r.key)
    labels = String[r.label for r in records]
    roots = String[
        r.label for r in records if r.parent_key === nothing
    ]
    return (;
        n_collections=length(records),
        n_items=length(workspace.index.items),
        sample_labels=labels[1:min(limit, length(labels))],
        root_labels=roots[1:min(limit, length(roots))],
    )
end

function _snap(workspace, t0::Float64)
    st = workspace_status(workspace)
    labels = _sample_labels(workspace)
    return (;
        t=round(time() - t0; digits=3),
        busy=st.busy,
        status_label=st.label,
        status_detail=st.detail,
        scan_state=workspace.scan.state,
        discovered=workspace.scan.discovered[],
        cache_state=workspace.cache_state,
        cache_operation=workspace.cache.operation,
        rebuild_flag_meaning="cache.operation === :rebuild => this open used rebuild=true",
        index_source_nothing=workspace.index.source === nothing,
        labels...,
        cache=_cache_stat(),
        work_nodes=length(workspace.work.nodes),
    )
end

"""Open fresh workspace (optional GUI) and dump status/index for the first second."""
function dump!(; gui::Bool=false, sample_hz::Float64=20.0, window_s::Float64=1.0)
    isdefined(Main, :DataBrowserInteractiveProfile) || error(
        "Load DataBrowserInteractiveProfile first")
    profile_project = Main.DataBrowserInteractiveProfile.project()
    outdir = joinpath(
        REPO_ROOT, "bench", "results",
        string(format(now(), "yyyymmdd-HHMMSS-sss"), "-fresh-first-second", gui ? "-gui" : "-headless"),
    )
    mkpath(outdir)
    report = joinpath(outdir, "first_second_dump.txt")

    open(report, "w") do io
        println(io, "outdir=$outdir")
        println(io, "gui=$gui")
        println(io, "threads=$(Threads.nthreads())")
        println(io, "commit=$(strip(read(`git -C $REPO_ROOT rev-parse --short HEAD`, String)))")
        println(io, "cache_before=$(_cache_stat())")
        flush(io)

        open_seconds = @elapsed workspace = open_workspace(
            profile_project,
            DATA_ROOT;
            metadata_file="device_info.txt",
            rebuild=true,
            cache=true,
            background_processing=false,
        )
        t0 = time()
        println(io, "\n=== immediately after open_workspace(rebuild=true) open_seconds=$open_seconds ===")
        println(io, _snap(workspace, t0))
        flush(io)

        browser = nothing
        first_frame_seconds = NaN
        if gui
            GLMakie = Main.GLMakie
            GLMakie.activate!()
            browser_started = time()
            browser = DataBrowser.open_browser(workspace; wait=false)
            deadline = time() + 120.0
            while isnan(browser.state.performance.first_frame_at) && time() < deadline
                sleep(0.005)
            end
            ff = browser.state.performance.first_frame_at
            first_frame_seconds = isnan(ff) ? NaN : ff - browser_started
            println(io, "\n=== first GUI frame at t=$(round(time() - t0; digits=3))s first_frame_seconds=$first_frame_seconds ===")
            println(io, _snap(workspace, t0))
            flush(io)
        end

        dt = 1.0 / sample_hz
        while (time() - t0) < window_s
            sleep(dt)
            snap = _snap(workspace, t0)
            println(io, "\n=== sample t=$(snap.t)s ===")
            println(io, snap)
            flush(io)
            if snap.n_collections > 0
                println(io, "\n*** COLLECTIONS PRESENT DURING FIRST WINDOW ***")
                println(io, "root_labels=$(snap.root_labels)")
                println(io, "sample_labels=$(snap.sample_labels)")
                flush(io)
            end
        end

        println(io, "\n=== end of $(window_s)s window ===")
        println(io, _snap(workspace, t0))
        flush(io)

        try
            browser isa DataBrowser.BrowserSession && DataBrowser.close_browser!(browser)
        finally
            close_workspace!(workspace)
        end
        println(io, "\n=== after close cache=$(_cache_stat()) ===")
    end
    println("Wrote $report")
    print(read(report, String))
    return report
end

end # module
