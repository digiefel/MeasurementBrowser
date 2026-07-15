#!/usr/bin/env julia

const USAGE = """
Profile DataBrowser with the real RuO2 v2 project using Julia's standard profilers.

Usage:
  julia --project=. --threads=auto ruo2_v2_profile.jl headless [options]
  julia -i --project=. --threads=auto ruo2_v2_profile.jl gui [options]
  julia --project=. ruo2_v2_profile.jl load-only

Options:
  --profile=cpu|allocs|none   Headless profiler (default: cpu)
  --background-processing    Process and analyze the complete workspace
  --sample-rate=RATE          Allocation sampling rate (default: 0.01)
  --timeout=SECONDS           Idle timeout (default: 1800)
  --output=PATH               Artifact directory
  --help                      Show this message
"""

if isempty(ARGS) || any(arg -> arg == "--help", ARGS)
    println(USAGE)
    exit(isempty(ARGS) ? 1 : 0)
end

const COMMAND = first(ARGS)
const OPTIONS = ARGS[2:end]
COMMAND in ("headless", "gui", "load-only") || error("Unknown command '$COMMAND'\n\n$USAGE")

using Dates: format, now
using Printf: @printf
using Profile

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const V2_DIR = get(
    ENV,
    "RUO2_V2_DIR",
    "/Users/davide/Library/CloudStorage/OneDrive-LundUniversity/projects/Borg/202501_RuO2test/analysis/v2",
)
const DATA_ROOT = get(ENV, "RUO2_DATA_ROOT", normpath(joinpath(V2_DIR, "..", "..", "electricaldata")))

option(name::String) = name in OPTIONS
function option_value(prefix::String, default::String)::String
    match = findfirst(arg -> startswith(arg, prefix * "="), OPTIONS)
    return match === nothing ? default : split(OPTIONS[match], '='; limit=2)[2]
end

function load_v2_project()
    browser_file = joinpath(V2_DIR, "browser.jl")
    isfile(browser_file) || error("RuO2 v2 project not found: $browser_file")
    source = read(browser_file, String)
    marker = "# finally, open the browser"
    occursin(marker, source) || error("RuO2 browser launcher marker changed in $browser_file")
    definitions = first(split(source, marker; limit=2))
    mktempdir() do loader_dir
        for relative_dir in ("data", "analysis", "plots")
            symlink(joinpath(V2_DIR, relative_dir), joinpath(loader_dir, relative_dir))
        end
        loader = joinpath(loader_dir, "browser.jl")
        write(loader, definitions)
        Base.include(Main, loader)
    end
    isdefined(Main, :PROJECT) || error("RuO2 browser did not define PROJECT")
    return Base.invokelatest(() -> getfield(Main, :PROJECT))
end

project = load_v2_project()
databrowser_source = pathof(Main.DataBrowser)
println("DataBrowser source: ", databrowser_source)
println("Active project: ", Base.active_project())
println("RuO2 project: ", V2_DIR)
println("Data root: ", DATA_ROOT)
startswith(String(databrowser_source), REPO_ROOT) || error(
    "Loaded DataBrowser from $databrowser_source instead of target checkout $REPO_ROOT",
)

COMMAND == "load-only" && exit(0)

using DataBrowserCore.Workspace:
    close_workspace!,
    engine_work_running,
    open_workspace,
    rebuild_cache!,
    wait_workspace_idle!

const TIMEOUT = parse(Float64, option_value("--timeout", "1800"))
const SAMPLE_RATE = parse(Float64, option_value("--sample-rate", "0.01"))
const TIMESTAMP = format(now(), "yyyymmdd-HHMMSS")
const OUTDIR = abspath(option_value(
    "--output",
    joinpath(REPO_ROOT, "bench", "results", TIMESTAMP * "-ruo2-v2-" * COMMAND),
))
mkpath(OUTDIR)

function wait_idle!(workspace)::Nothing
    wait_workspace_idle!(workspace; timeout=TIMEOUT)
    engine_work_running(workspace) && error("Workspace did not become idle within $TIMEOUT seconds")
    return nothing
end

function measured_rebuild!(workspace, profiler::String)
    GC.gc()
    bytes_before = Base.gc_bytes()
    elapsed = if profiler == "cpu"
        Profile.clear()
        @elapsed Profile.@profile begin
            rebuild_cache!(workspace)
            wait_idle!(workspace)
        end
    elseif profiler == "allocs"
        Profile.Allocs.clear()
        @elapsed Profile.Allocs.@profile sample_rate=SAMPLE_RATE begin
            rebuild_cache!(workspace)
            wait_idle!(workspace)
        end
    elseif profiler == "none"
        @elapsed begin
            rebuild_cache!(workspace)
            wait_idle!(workspace)
        end
    else
        error("Unknown profiler '$profiler'; use cpu, allocs, or none")
    end
    return elapsed, Base.gc_bytes() - bytes_before
end

function write_standard_profile(profiler::String)::Nothing
    if profiler == "cpu"
        open(joinpath(OUTDIR, "cpu_profile.txt"), "w") do io
            Profile.print(io; format=:flat, sortedby=:count, mincount=2)
        end
    elseif profiler == "allocs"
        open(joinpath(OUTDIR, "allocation_profile.txt"), "w") do io
            Profile.Allocs.print(io; format=:flat, sortedby=:count, mincount=2)
        end
    end
    return nothing
end

function run_headless()::Nothing
    profiler = option_value("--profile", "cpu")
    background = option("--background-processing")
    cache_root = mktempdir()
    pushfirst!(DEPOT_PATH, cache_root)
    workspace = nothing
    try
        workspace = open_workspace(
            project,
            DATA_ROOT;
            metadata_file="device_info.txt",
            cache=true,
            background_processing=background,
        )
        println("Warm-up build ...")
        wait_idle!(workspace)

        elapsed, allocated = measured_rebuild!(workspace, profiler)
        items = length(workspace.index.items)
        sources = length(workspace.index.items_by_source)
        throughput = items / elapsed
        close_seconds = @elapsed close_workspace!(workspace)
        workspace = nothing
        write_standard_profile(profiler)

        open(joinpath(OUTDIR, "run.csv"), "w") do io
            println(io, "commit,threads,background_processing,profiler,sources,items,seconds,items_per_s,allocated_mib,close_seconds")
            commit = readchomp(`git -C $REPO_ROOT rev-parse HEAD`)
            @printf(io, "%s,%d,%s,%s,%d,%d,%.6f,%.3f,%.3f,%.6f\n",
                commit, Threads.nthreads(), background, profiler, sources, items, elapsed,
                throughput, allocated / 1024^2, close_seconds)
        end
        @printf("Measured rebuild: %.3f s, %d items, %.1f items/s, %.1f MiB allocated\n",
            elapsed, items, throughput, allocated / 1024^2)
        println("Artifacts: ", OUTDIR)
    finally
        workspace === nothing || close_workspace!(workspace)
        first(DEPOT_PATH) == cache_root && popfirst!(DEPOT_PATH)
        rm(cache_root; force=true, recursive=true)
    end
    return nothing
end

function profile_gui_cpu(seconds::Real=30)::String
    path = joinpath(OUTDIR, "gui_cpu_profile.txt")
    Profile.clear()
    Profile.@profile sleep(Float64(seconds))
    open(path, "w") do io
        Profile.print(io; format=:flat, sortedby=:count, mincount=2)
    end
    println("GUI CPU profile: ", path)
    return path
end

function profile_gui_allocs(seconds::Real=15; rate::Real=0.01)::String
    path = joinpath(OUTDIR, "gui_allocation_profile.txt")
    Profile.Allocs.clear()
    Profile.Allocs.@profile sample_rate=rate sleep(Float64(seconds))
    open(path, "w") do io
        Profile.Allocs.print(io; format=:flat, sortedby=:count, mincount=2)
    end
    println("GUI allocation profile: ", path)
    return path
end

function run_gui()::Nothing
    isinteractive() || error("GUI mode requires `julia -i` so the native window remains attached")
    workspace = open_workspace(project, DATA_ROOT; metadata_file="device_info.txt")
    browser = Main.DataBrowser.open_browser(workspace)
    global DATABROWSER_PROFILE_WORKSPACE = workspace
    global DATABROWSER_PROFILE_BROWSER = browser
    println("Native DataBrowser launched. Confirm the window is visible before profiling.")
    println("CPU: profile_gui_cpu(30)")
    println("Allocations: profile_gui_allocs(15; rate=0.01)")
    println("Artifacts: ", OUTDIR)
    return nothing
end

COMMAND == "headless" ? run_headless() : run_gui()
