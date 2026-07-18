#!/usr/bin/env julia
# Primary hot-loop harness: disposable cache (never touches databrowser/RuO2/cache.duckdb). See SKILL.md.

const MODE = isempty(ARGS) ? "headless" : first(ARGS)
MODE in ("headless", "gui") || error("Unknown mode '$MODE'; use headless or gui")
option_seconds(name, default) = let m = findfirst(a -> startswith(a, "--$name="), ARGS)
    m === nothing ? default : parse(Float64, split(ARGS[m], '='; limit=2)[2])
end
const BUDGET = option_seconds("budget", 30.0)
const WARMUP = option_seconds("warmup", 20.0)
const FRESH = "--fresh" in ARGS
const BACKGROUND_PROCESSING = !("--no-background-processing" in ARGS)
const PROFILER = let m = findfirst(a -> startswith(a, "--profile="), ARGS)
    m === nothing ? "none" : split(ARGS[m], '='; limit=2)[2]
end
PROFILER in ("none", "cpu", "wall", "allocs") ||
    error("Unknown --profile '$PROFILER'; use cpu, wall, or allocs")

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const DEFINITIONS = normpath(joinpath(@__DIR__, "..", "project", "definitions.jl"))
const DATA_ROOT = get(
    ENV,
    "RUO2_DATA_ROOT",
    "/Users/davide/Library/CloudStorage/OneDrive-LundUniversity/projects/Borg/202501_RuO2test/electricaldata",
)
const RUN_LOCK = joinpath(tempdir(), "databrowser-profiling.lock")
try
    mkdir(RUN_LOCK)
catch
    error("Another profiling run may be active; inspect processes before removing $RUN_LOCK")
end
atexit(() -> rm(RUN_LOCK; recursive=true, force=true))

using Dates: format, now
using Printf: @printf, @sprintf
using Profile

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
if FRESH
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
project = Base.invokelatest(() -> getfield(Main, :PROJECT))
mark!("project definitions loaded")

using DataBrowserCore.Workspace: close_workspace!, open_workspace, workspace_status

# --- 4. warmup scan: compile every callback and app path, then discard the result -------------
warmup_seconds = 0.0
warmup_jit_seconds = 0.0
if WARMUP > 0
    warmup_root = mktempdir()
    pushfirst!(DEPOT_PATH, warmup_root)
    jit_before_warmup = compile_ns()
    warmup_seconds = @elapsed begin
        warm = open_workspace(
            project, DATA_ROOT; metadata_file="device_info.txt", cache=true,
            background_processing=BACKGROUND_PROCESSING,
        )
        deadline = time() + WARMUP
        while time() < deadline && workspace_status(warm).busy
            sleep(0.25)
        end
        close_workspace!(warm)
    end
    warmup_jit_seconds = (compile_ns() - jit_before_warmup) / 1e9
    first(DEPOT_PATH) == warmup_root && popfirst!(DEPOT_PATH)
    rm(warmup_root; force=true, recursive=true)
    # The measured window must only aggregate its own callback timings.
    lock(project.profile_lock) do
        empty!(project.scan_profile)
    end
end
mark!("warmup done")

# --- 5. open_workspace with an isolated cache -----------------------------------------------
cache_root = mktempdir()
pushfirst!(DEPOT_PATH, cache_root)

compile_before_open = compile_ns()
open_seconds = @elapsed workspace = open_workspace(
    project, DATA_ROOT; metadata_file="device_info.txt", cache=true,
    background_processing=BACKGROUND_PROCESSING,
)
open_compile_seconds = (compile_ns() - compile_before_open) / 1e9
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
        status.busy || break
        sleep(0.1)  # ~10 Hz
    end
    return nothing
end

compile_before_scan = compile_ns()
scan_t0 = time()
if PROFILER == "cpu"
    Profile.clear()
    Profile.@profile sample_scan_window!(samples, workspace, scan_t0)
elseif PROFILER == "wall"
    Profile.clear()
    Profile.@profile_walltime sample_scan_window!(samples, workspace, scan_t0)
elseif PROFILER == "allocs"
    Profile.Allocs.clear()
    Profile.Allocs.@profile sample_rate = 0.01 sample_scan_window!(samples, workspace, scan_t0)
else
    sample_scan_window!(samples, workspace, scan_t0)
end
scan_window = time() - scan_t0
scan_compile_seconds = (compile_ns() - compile_before_scan) / 1e9
mark!("scan window closed")

# --- 6. pipeline callback seconds vs wall capacity (app overhead) ----------------------------
phase_totals = Dict{Symbol,Float64}(
    :detect => 0.0, :read => 0.0, :entries => 0.0, :process => 0.0, :analyze => 0.0)
worker_threads = Set{Int}()
lock(project.profile_lock) do
    for entry in values(project.scan_profile)
        phase_totals[:detect] += entry.detect_seconds
        phase_totals[:read] += entry.read_seconds
        phase_totals[:entries] += entry.entries_seconds
        phase_totals[:process] += entry.process_seconds
        phase_totals[:analyze] += entry.analyze_seconds
        union!(worker_threads, entry.thread_ids)
    end
end
callback_seconds = sum(values(phase_totals))
capacity_seconds = scan_window * max(length(worker_threads), 1)

# --- 7. shutdown ------------------------------------------------------------------------------
# Close the GUI first so it does not fight close_workspace! (its exit path also closes the workspace).
close_browser_seconds = NaN
if browser isa Main.DataBrowser.BrowserSession
    close_browser_seconds = @elapsed Main.DataBrowser.close_browser!(browser)
    mark!("browser closed")
end
close_seconds = @elapsed close_workspace!(workspace)
first(DEPOT_PATH) == cache_root && popfirst!(DEPOT_PATH)
rm(cache_root; force=true, recursive=true)
mark!("closed")

# --- report -----------------------------------------------------------------------------------
if PROFILER == "cpu"
    try
        @eval using ProfileSVG
        Base.invokelatest(ProfileSVG.save, joinpath(OUTDIR, "cpu_flame.svg"))
    catch err
        @warn "cpu_flame.svg not written" exception = err
        open(joinpath(OUTDIR, "cpu_profile.txt"), "w") do io
            Profile.print(io; format=:flat, sortedby=:count, mincount=2)
        end
    end
elseif PROFILER == "wall"
    open(joinpath(OUTDIR, "wall_profile.txt"), "w") do io
        Profile.print(io; format=:tree, mincount=2)
    end
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

function milestone(samples::Vector{Sample}, done::Function)::Float64
    i = findfirst(done, samples)
    return i === nothing ? NaN : samples[i].t
end

last_s = samples[end]
# Full discovery = first sample that already knows the final source count for this run.
discovery_done = milestone(samples, s -> s.sources_found == last_s.sources_found && last_s.sources_found > 0)
interpret_done = milestone(samples,
    s -> last_s.sources_found > 0 && s.cached_sources >= last_s.sources_found)
process_done = milestone(samples, s -> s.interpreted > 0 && s.processed >= s.interpreted)
analyze_done = milestone(samples, s -> s.interpreted > 0 && s.analyzed >= s.interpreted)
finished = !last_s.busy

# Stage counts vs time (~10 Hz samples).
try
    using GLMakie
    fig = Figure(size = (900, 500))
    ax = Axis(fig[1, 1]; xlabel = "t (s)", ylabel = "count", title = "scan timeline")
    ts = [s.t for s in samples]
    lines!(ax, ts, [s.sources_found for s in samples]; label = "sources_found")
    lines!(ax, ts, [s.cached_sources for s in samples]; label = "cached_sources")
    lines!(ax, ts, [s.interpreted for s in samples]; label = "interpreted")
    lines!(ax, ts, [s.processed for s in samples]; label = "processed")
    lines!(ax, ts, [s.analyzed for s in samples]; label = "analyzed")
    axislegend(ax; position = :lt)
    save(joinpath(OUTDIR, "timeline.png"), fig)
catch err
    @warn "timeline.png not written" exception = err
end

commit = readchomp(`git -C $REPO_ROOT rev-parse --short HEAD`)
fmt(x) = isnan(x) ? "not reached" : @sprintf("%.1f s", x)
rate(n) = n / max(scan_window, eps())

println()
println("=== DataBrowser timed profile ($MODE, $commit, $(Threads.nthreads()) threads, " *
    "budget $(BUDGET) s, fresh=$FRESH, background=$BACKGROUND_PROCESSING) ===")
println()
println("Startup")
@printf("  precompile             %8.1f s%s\n", precompile_seconds,
    FRESH ? "" : "   (caches kept — not a fresh measurement)")
@printf("  using DataBrowser      %8.1f s   (%.1f s compile)\n", load_seconds, load_compile_seconds)
@printf("  project include        %8.1f s\n", include_seconds)
@printf("  open_workspace         %8.1f s   (%.1f s compile)\n", open_seconds, open_compile_seconds)
if MODE == "gui"
    @printf("  open_browser           %8.1f s\n", browser_seconds)
    @printf("  time to first frame    %8s\n", fmt(first_frame_seconds))
end
println()
println("Scan window ($(round(scan_window; digits=1)) s sampled, $(finished ? "pipeline finished" : "budget hit"))")
@printf("  runtime JIT in window  %8.1f s\n", scan_compile_seconds)
println("  discovery done         ", fmt(discovery_done), "   ($(last_s.sources_found) sources)")
println("  interpretation done    ", fmt(interpret_done), "   ($(last_s.interpreted) items from $(last_s.cached_sources) sources)")
println("  processing done        ", fmt(process_done), "   ($(last_s.processed) processed)")
println("  analysis done          ", fmt(analyze_done), "   ($(last_s.analyzed) analyzed)")
println()
println("Throughput (items / scan-window second)")
@printf("  interpretation         %8.1f\n", rate(last_s.interpreted))
@printf("  processing             %8.1f\n", rate(last_s.processed))
println()
println("Pipeline callbacks vs capacity ($(length(worker_threads)) worker threads seen)")
for phase in (:detect, :read, :entries, :process, :analyze)
    @printf("  %-9s              %8.1f s\n", phase, phase_totals[phase])
end
@printf("  callbacks total        %8.1f s  of %.1f s worker capacity (%.0f%% callbacks, %.0f%% unattributed)\n",
    callback_seconds, capacity_seconds, 100 * callback_seconds / capacity_seconds,
    100 * (1 - callback_seconds / capacity_seconds))
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
