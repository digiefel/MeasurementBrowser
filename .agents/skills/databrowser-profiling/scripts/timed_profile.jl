#!/usr/bin/env julia
# Fresh-process boundary for the dedicated DataBrowserProfilingRuO2 cache. See SKILL.md.

const MODE = isempty(ARGS) ? "headless" : first(ARGS)
MODE in ("headless", "gui") || error("Unknown mode '$MODE'; use headless or gui")
option_seconds(name, default) = let m = findfirst(a -> startswith(a, "--$name="), ARGS)
    m === nothing ? default : parse(Float64, split(ARGS[m], '='; limit=2)[2])
end
option_string(name, default) = let m = findfirst(a -> startswith(a, "--$name="), ARGS)
    m === nothing ? default : split(ARGS[m], '='; limit=2)[2]
end
const BUDGET = option_seconds("budget", 60.0)
const CACHE_MODE = option_string("cache", "fresh")
CACHE_MODE in ("fresh", "resume") || error("--cache must be fresh or resume")
const FRESH_COMPILE = "--fresh-compile" in ARGS
const BACKGROUND_PROCESSING = "--background-processing" in ARGS
const PROFILER = let m = findfirst(a -> startswith(a, "--profile="), ARGS)
    m === nothing ? "none" : split(ARGS[m], '='; limit=2)[2]
end
# "cpu" is not offered: the Mach-based CPU sampler can wedge the process on macOS (see SKILL.md).
PROFILER in ("none", "wall", "allocs") ||
    error("Unknown --profile '$PROFILER'; use wall or allocs")

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const DEFINITIONS = normpath(joinpath(@__DIR__, "..", "project", "definitions.jl"))
const DATA_ROOT = get(
    ENV,
    "RUO2_DATA_ROOT",
    "/Users/davide/Library/CloudStorage/OneDrive-LundUniversity/projects/Borg/202501_RuO2test/electricaldata",
)
const RUN_LOCK = joinpath(tempdir(), "databrowser-profiling.pid")

using Dates: format, now
using FileWatching: trymkpidlock
using FileIO
using FlameGraphs
using PProf
using Printf: @printf, @sprintf
using Profile

include("profile_artifacts.jl")
using .ProfileArtifacts: save_profile

const RUN_LOCK_HANDLE = let run_lock = trymkpidlock(RUN_LOCK; stale_age=1)
    run_lock === false && error("Another profiling run is active: $RUN_LOCK")
    run_lock
end
atexit(() -> close(RUN_LOCK_HANDLE))

Base.cumulative_compile_timing(true)
compile_ns() = first(Base.cumulative_compile_time_ns())

const OUTDIR = abspath(joinpath(REPO_ROOT, "bench", "results",
    format(now(), "yyyymmdd-HHMMSS") * "-ruo2-v2-timed-" * MODE))
mkpath(OUTDIR)

const TIMELINE_EVENTS = Tuple{String,Float64}[]  # (event, seconds since process t0)
const T0 = time()
mark!(event::String) = push!(TIMELINE_EVENTS, (event, time() - T0))

# --- 1. fresh precompile of the code under test --------------------------------------------------
mark!("start")
using Pkg
if FRESH_COMPILE
    compiled = joinpath(first(DEPOT_PATH), "compiled", "v$(VERSION.major).$(VERSION.minor)")
    for entry in filter(startswith("DataBrowser"), readdir(compiled))
        rm(joinpath(compiled, entry); recursive=true, force=true)
    end
end
precompile_seconds = @elapsed Pkg.precompile()
mark!("precompile done")

# --- 2. package load ------------------------------------------------------------------------
compile_before_load = compile_ns()
load_seconds = @elapsed @eval using DataBrowser
load_compile_seconds = (compile_ns() - compile_before_load) / 1e9
mark!("using DataBrowser done")

pathof(DataBrowser) !== nothing && startswith(pathof(DataBrowser), REPO_ROOT) ||
    error("Loaded DataBrowser from $(pathof(DataBrowser)) instead of $REPO_ROOT")

# --- 3. project definitions (include + register! calls) --------------------------------------
include_seconds = @elapsed Base.include(Main, DEFINITIONS)
isdefined(Main, :PROJECT) || error("$DEFINITIONS did not define PROJECT")
project = getfield(Main, :PROJECT)
mark!("project definitions loaded")

using DataBrowserCore.Workspace: close_workspace!, open_workspace, workspace_status

# --- 4. open_workspace -----------------------------------------------------------------------
compile_before_open = compile_ns()
open_seconds = @elapsed workspace = open_workspace(
    project, DATA_ROOT; metadata_file="device_info.txt", cache=true,
    rebuild=CACHE_MODE == "fresh",
    background_processing=BACKGROUND_PROCESSING,
)
open_compile_seconds = (compile_ns() - compile_before_open) / 1e9
cache_path = workspace.cache.identity.cache_path
mark!("open_workspace returned")

browser_seconds = NaN
first_frame_seconds = NaN
browser = nothing
t_browser = time()
if MODE == "gui"
    browser = Main.DataBrowser.open_browser(workspace; wait=false)
    mark!("open_browser returned")
    browser_seconds = time() - t_browser
    # Wait until the first non-blank frame is submitted (startup surface or full UI).
    deadline = time() + 120.0
    while isnan(browser.state.performance.first_frame_at) && time() < deadline
        sleep(0.01)
    end
    first_frame_at = browser.state.performance.first_frame_at
    first_frame_seconds = isnan(first_frame_at) ? NaN : (first_frame_at - t_browser)
    mark!("first frame")
end

# --- 5. sample the scan for the budget window ------------------------------------------------
struct Sample
    t::Float64
    sources_found::Int
    sources_pending::Int
    cached_sources::Int
    interpreted::Int
    processed::Int
    analyzed::Int
    collection_processed::Int
    collection_analyzed::Int
    busy::Bool
end

samples = Sample[]
function sample_scan_window!(samples::Vector{Sample}, workspace, scan_t0::Float64)::Nothing
    while time() - scan_t0 < BUDGET
        status = workspace_status(workspace)
        c = status.counts
        push!(samples, Sample(time() - scan_t0, c.sources_found, c.sources_pending,
            c.cache.cached_sources, c.cache.interpreted_items, c.cache.processed,
            c.cache.analyzed, c.cache.collection_processed, c.cache.collection_analyzed,
            status.busy))
        sleep(0.1)  # ~10 Hz
    end
    return nothing
end

compile_before_scan = compile_ns()
scan_t0 = time()
profile_data = nothing
if PROFILER == "wall"
    Profile.clear()
    Profile.@profile_walltime sample_scan_window!(samples, workspace, scan_t0)
    profile_data = Profile.retrieve()
elseif PROFILER == "allocs"
    Profile.Allocs.clear()
    Profile.Allocs.@profile sample_rate = 0.01 sample_scan_window!(samples, workspace, scan_t0)
else
    sample_scan_window!(samples, workspace, scan_t0)
end
scan_window = time() - scan_t0
scan_compile_seconds = (compile_ns() - compile_before_scan) / 1e9
mark!("measurement period ended")

# --- 6. pipeline callback seconds vs wall capacity (app overhead) ----------------------------
phase_totals = Dict{Symbol,Float64}(
    :detect => 0.0, :read => 0.0, :entries => 0.0, :process => 0.0, :analyze => 0.0)
lock(project.profile_lock) do
    for entry in values(project.scan_profile)
        phase_totals[:detect] += entry.detect_seconds
        phase_totals[:read] += entry.read_seconds
        phase_totals[:entries] += entry.entries_seconds
        phase_totals[:process] += entry.process_seconds
        phase_totals[:analyze] += entry.analyze_seconds
    end
end
callback_seconds = sum(values(phase_totals))

# --- 7. shutdown ------------------------------------------------------------------------------
# Close the GUI first so it does not fight close_workspace! (its exit path also closes the workspace).
close_browser_seconds = NaN
if browser isa Main.DataBrowser.BrowserSession
    close_browser_seconds = @elapsed Main.DataBrowser.close_browser!(browser)
    mark!("browser closed")
end
close_seconds = @elapsed close_workspace!(workspace)
mark!("closed")

# --- report -----------------------------------------------------------------------------------
if profile_data !== nothing
    data, lidict = profile_data
    raw = joinpath(OUTDIR, "$(PROFILER)_profile.jls")
    jlprof = joinpath(OUTDIR, "$(PROFILER)_profile.jlprof")
    protobuf = joinpath(OUTDIR, "$(PROFILER)_profile.pb.gz")
    save_profile(raw, data, lidict)
    save(File{format"JLPROF"}(jlprof), data, lidict)
    PProf.pprof(
        data,
        lidict;
        web=false,
        out=protobuf,
        sampling_delay=ccall(:jl_profile_delay_nsec, UInt64, ()),
        full_signatures=true,
    )
elseif PROFILER == "allocs"
    open(joinpath(OUTDIR, "allocation_profile.txt"), "w") do io
        Profile.Allocs.print(io; format=:flat, sortedby=:count, mincount=2)
    end
end

open(joinpath(OUTDIR, "timeline.csv"), "w") do io
    println(io, "t_seconds,sources_found,sources_pending,cached_sources,interpreted,processed,analyzed,collection_processed,collection_analyzed,busy")
    for s in samples
        @printf(io, "%.3f,%d,%d,%d,%d,%d,%d,%d,%d,%s\n", s.t, s.sources_found,
            s.sources_pending, s.cached_sources, s.interpreted, s.processed, s.analyzed,
            s.collection_processed, s.collection_analyzed, s.busy)
    end
end

last_s = samples[end]
# Stage counts vs time (~10 Hz samples).
try
    using CairoMakie
    fig = CairoMakie.Figure(size = (900, 500))
    ax = CairoMakie.Axis(
        fig[1, 1]; xlabel = "t (s)", ylabel = "count", title = "scan timeline")
    ts = [s.t for s in samples]
    CairoMakie.lines!(ax, ts, [s.sources_found for s in samples]; label = "sources_found")
    CairoMakie.lines!(ax, ts, [s.cached_sources for s in samples]; label = "cached_sources")
    CairoMakie.lines!(ax, ts, [s.interpreted for s in samples]; label = "interpreted")
    CairoMakie.lines!(ax, ts, [s.processed for s in samples]; label = "processed")
    CairoMakie.lines!(ax, ts, [s.analyzed for s in samples]; label = "analyzed")
    CairoMakie.axislegend(ax; position = :lt)
    CairoMakie.save(joinpath(OUTDIR, "timeline.png"), fig)
catch err
    @warn "timeline.png not written" exception = err
end

commit = readchomp(`git -C $REPO_ROOT rev-parse --short HEAD`)
fmt(x) = isnan(x) ? "not reached" : @sprintf("%.1f s", x)
first_s = first(samples)
rate(last, first) = (last - first) / max(scan_window, eps())

println()
println("=== DataBrowser timed profile ($MODE, $commit, $(Threads.nthreads()) threads, " *
    "budget $(BUDGET) s, cache=$CACHE_MODE, background=$BACKGROUND_PROCESSING) ===")
println("Cache: $cache_path")
println()
println("Startup")
@printf("  precompile             %8.1f s%s\n", precompile_seconds,
    FRESH_COMPILE ? "" : "   (compiled package caches kept)")
@printf("  using DataBrowser      %8.1f s   (%.1f s aggregate compile)\n", load_seconds, load_compile_seconds)
@printf("  project script load    %8.1f s   (helpers + define_project/register_* calls)\n", include_seconds)
@printf("  open_workspace         %8.1f s   (%.1f s aggregate compile)\n", open_seconds, open_compile_seconds)
if MODE == "gui"
    @printf("  open_browser           %8.1f s\n", browser_seconds)
    @printf("  time to first frame    %8s\n", fmt(first_frame_seconds))
end
println()
println("Measurement period: $(round(scan_window; digits=1)) s " *
    (MODE == "gui" ? "after first frame" : "after open_workspace returned"))
@printf("  runtime compilation during the measurement, summed across Julia threads  %8.1f s\n",
    scan_compile_seconds)
println()
println("Progress and throughput during the measurement")
@printf("  cached sources         %8.1f sources/s      (%d → %d)\n",
    rate(last_s.cached_sources, first_s.cached_sources),
    first_s.cached_sources, last_s.cached_sources)
@printf("  interpreted items      %8.1f items/s        (%d → %d)\n",
    rate(last_s.interpreted, first_s.interpreted), first_s.interpreted, last_s.interpreted)
@printf("  processed items        %8.1f items/s        (%d → %d)\n",
    rate(last_s.processed, first_s.processed), first_s.processed, last_s.processed)
@printf("  analyzed items         %8.1f items/s        (%d → %d)\n",
    rate(last_s.analyzed, first_s.analyzed), first_s.analyzed, last_s.analyzed)
@printf("  processed collections  %8.1f collections/s  (%d → %d)\n",
    rate(last_s.collection_processed, first_s.collection_processed),
    first_s.collection_processed, last_s.collection_processed)
@printf("  analyzed collections   %8.1f collections/s  (%d → %d)\n",
    rate(last_s.collection_analyzed, first_s.collection_analyzed),
    first_s.collection_analyzed, last_s.collection_analyzed)
println()
if callback_seconds == 0
    println("Project callbacks: none ran during this profile")
else
    println("Summed elapsed time inside callbacks")
    for phase in (:detect, :read, :entries, :process, :analyze)
        @printf("  %-20s %8.1f s\n", "$(phase) callbacks", phase_totals[phase])
    end
    @printf("  total                   %8.1f s\n", callback_seconds)
end
println()
MODE == "gui" && !isnan(close_browser_seconds) &&
    @printf("Shutdown: close_browser! %.2f s, close_workspace! %.2f s\n", close_browser_seconds, close_seconds)
MODE != "gui" && @printf("Shutdown: close_workspace! %.2f s\n", close_seconds)
println()
println("Timeline events")
for (event, t) in TIMELINE_EVENTS
    @printf("  %7.1f s  %s\n", t, event)
end
println()
println("Artifacts: ", OUTDIR)
