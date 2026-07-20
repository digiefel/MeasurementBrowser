module DataBrowserInteractiveProfile

using Dates: format, now
using FileWatching: trymkpidlock
using FileIO
using FlameGraphs
using PProf
using Printf: @printf, @sprintf
using Profile
import ProfileView
using DataBrowser
using GLMakie
import CairoMakie
using DataBrowserCore.Workspace:
    close_workspace!, open_workspace, rebuild_cache!, workspace_status

include("profile_artifacts.jl")
using .ProfileArtifacts: profile_sample_count, save_profile

Base.cumulative_compile_timing(true)

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const DATA_ROOT = get(
    ENV,
    "RUO2_DATA_ROOT",
    "/Users/davide/Library/CloudStorage/OneDrive-LundUniversity/projects/Borg/202501_RuO2test/electricaldata",
)
const RUN_LOCK = joinpath(tempdir(), "databrowser-profiling.pid")
const LAST_PROFILE = Ref{Any}(nothing)
const LAST_PPROF = Ref{Union{Nothing,String}}(nothing)
const LIVE_SESSION = Ref{Any}(nothing)

module DataBrowserProfilingProject
include(joinpath(@__DIR__, "..", "project", "definitions.jl"))
end

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

mutable struct LiveSession
    profile_project::Any
    workspace::Any
    browser::Any
    mode::Symbol
    cache_mode::Symbol
    cache_path::String
    background_processing::Bool
    open_seconds::Float64
    first_frame_seconds::Float64
    run_lock::Any
end

project() = DataBrowserProfilingProject.PROJECT

function save_profile_outputs(
    outdir::AbstractString,
    profiler::Symbol,
    profile_data,
)
    data, lidict = profile_data
    raw = joinpath(outdir, "$(profiler)_profile.jls")
    jlprof = joinpath(outdir, "$(profiler)_profile.jlprof")
    protobuf = joinpath(outdir, "$(profiler)_profile.pb.gz")
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
    LAST_PROFILE[] = profile_data
    LAST_PPROF[] = protobuf
    return (; raw, jlprof, protobuf)
end

function acquire_lock()
    run_lock = trymkpidlock(RUN_LOCK; stale_age=1)
    run_lock === false && error("Another profiling run is active: $RUN_LOCK")
    return run_lock
end

function clear_callback_profile!(profile_project)::Nothing
    lock(profile_project.profile_lock) do
        empty!(profile_project.scan_profile)
    end
    return nothing
end

function take_sample(workspace, started_at::Float64)::Sample
    status = workspace_status(workspace)
    counts = status.counts
    return Sample(
        time() - started_at,
        counts.sources_found,
        counts.sources_pending,
        counts.cache.cached_sources,
        counts.cache.interpreted_items,
        counts.cache.processed,
        counts.cache.analyzed,
        counts.cache.collection_processed,
        counts.cache.collection_analyzed,
        status.busy,
    )
end

function sample_scan_window!(
    samples::Vector{Sample},
    workspace,
    started_at::Float64,
    budget::Float64,
    stop_when_idle::Bool=true,
)::Nothing
    while true
        sample = take_sample(workspace, started_at)
        push!(samples, sample)
        (sample.t >= budget || (stop_when_idle && !sample.busy)) && break
        sleep(0.1)
    end
    return nothing
end

function open_profile_session(
    profile_project;
    mode::Symbol,
    cache_mode::Symbol,
    background_processing::Bool,
)
    workspace = nothing
    browser = nothing
    try
        open_seconds = @elapsed workspace = open_workspace(
            profile_project,
            DATA_ROOT;
            metadata_file="device_info.txt",
            rebuild=cache_mode === :fresh,
            cache=true,
            background_processing,
        )
        cache_path = workspace.cache.identity.cache_path
        first_frame_seconds = NaN
        if mode === :gui
            GLMakie.activate!()
            browser_started_at = time()
            browser = DataBrowser.open_browser(workspace; wait=false)
            deadline = time() + 120.0
            while isnan(browser.state.performance.first_frame_at) && time() < deadline
                sleep(0.01)
            end
            first_frame_at = browser.state.performance.first_frame_at
            first_frame_seconds =
                isnan(first_frame_at) ? NaN : first_frame_at - browser_started_at
        end
        return (; workspace, browser, cache_path, open_seconds, first_frame_seconds)
    catch
        browser isa DataBrowser.BrowserSession && DataBrowser.close_browser!(browser)
        workspace === nothing || close_workspace!(workspace)
        rethrow()
    end
end

function run_profile_window(
    profile_project;
    mode::Symbol,
    budget::Float64,
    cache_mode::Symbol,
    background_processing::Bool,
)
    workspace = nothing
    browser = nothing
    try
        opened = open_profile_session(
            profile_project; mode, cache_mode, background_processing)
        workspace = opened.workspace
        browser = opened.browser
        samples = Sample[]
        measurement_started_at = time()
        sample_scan_window!(samples, workspace, measurement_started_at, budget, false)
        measurement_seconds = time() - measurement_started_at
        return (;
            workspace,
            browser,
            cache_path=opened.cache_path,
            open_seconds=opened.open_seconds,
            first_frame_seconds=opened.first_frame_seconds,
            scan_seconds=measurement_seconds,
            samples,
        )
    catch
        browser isa DataBrowser.BrowserSession && DataBrowser.close_browser!(browser)
        workspace === nothing || close_workspace!(workspace)
        rethrow()
    end
end

function write_timeline(
    outdir::AbstractString,
    samples::Vector{Sample};
    restore_gl::Bool=false,
)::Nothing
    open(joinpath(outdir, "timeline.csv"), "w") do io
        println(io, "t_seconds,sources_found,sources_pending,cached_sources,interpreted,processed,analyzed,collection_processed,collection_analyzed,busy")
        for sample in samples
            @printf(
                io,
                "%.3f,%d,%d,%d,%d,%d,%d,%d,%d,%s\n",
                sample.t,
                sample.sources_found,
                sample.sources_pending,
                sample.cached_sources,
                sample.interpreted,
                sample.processed,
                sample.analyzed,
                sample.collection_processed,
                sample.collection_analyzed,
                sample.busy,
            )
        end
    end
    CairoMakie.activate!()
    try
        figure = CairoMakie.Figure(size=(900, 500))
        axis = CairoMakie.Axis(
            figure[1, 1];
            xlabel="t (s)",
            ylabel="count",
            title="scan timeline",
        )
        times = [sample.t for sample in samples]
        CairoMakie.lines!(
            axis, times, [sample.sources_found for sample in samples]; label="sources_found")
        CairoMakie.lines!(
            axis, times, [sample.cached_sources for sample in samples]; label="cached_sources")
        CairoMakie.lines!(
            axis, times, [sample.interpreted for sample in samples]; label="interpreted")
        CairoMakie.lines!(
            axis, times, [sample.processed for sample in samples]; label="processed")
        CairoMakie.lines!(
            axis, times, [sample.analyzed for sample in samples]; label="analyzed")
        CairoMakie.axislegend(axis; position=:lt)
        CairoMakie.save(joinpath(outdir, "timeline.png"), figure)
    finally
        restore_gl && GLMakie.activate!()
    end
    return nothing
end

function callback_totals(profile_project)
    phase_totals = Dict{Symbol,Float64}(
        :detect => 0.0,
        :read => 0.0,
        :entries => 0.0,
        :process => 0.0,
        :analyze => 0.0,
    )
    lock(profile_project.profile_lock) do
        for entry in values(profile_project.scan_profile)
            phase_totals[:detect] += entry.detect_seconds
            phase_totals[:read] += entry.read_seconds
            phase_totals[:entries] += entry.entries_seconds
            phase_totals[:process] += entry.process_seconds
            phase_totals[:analyze] += entry.analyze_seconds
        end
    end
    return phase_totals
end

function print_stage_row(
    io::IO,
    label::AbstractString,
    first_sample::Sample,
    last_sample::Sample,
    field::Symbol,
    scan_seconds::Float64,
)::Nothing
    first_value = getfield(first_sample, field)
    last_value = getfield(last_sample, field)
    change = last_value - first_value
    @printf(
        io,
        "  %-22s %8d %8d %+8d %10.1f\n",
        label,
        first_value,
        last_value,
        change,
        change / max(scan_seconds, eps()),
    )
    return nothing
end

function print_run_summary(
    io::IO;
    mode::Symbol,
    profiler::Symbol,
    cache_mode::Symbol,
    cache_path::AbstractString,
    background_processing::Bool,
    budget::Float64,
    profile_seconds::Float64,
    open_seconds::Float64,
    first_frame_seconds::Float64,
    scan_seconds::Float64,
    compile_seconds::Float64,
    samples::Vector{Sample},
    sample_count::Int,
    phases::Dict{Symbol,Float64},
    outdir::AbstractString,
    live::Bool=false,
)::Nothing
    first_sample, last_sample = first(samples), last(samples)
    fmt(value) = isnan(value) ? "not reached" : @sprintf("%.2f s", value)
    @printf(
        io,
        "DataBrowser %s profile: mode=%s profiler=%s budget=%.1f s threads=%d background_processing=%s\n",
        live ? "live" : "interactive",
        mode,
        profiler,
        budget,
        Threads.nthreads(),
        background_processing,
    )
    @printf(io, "cache: %s (%s)\n", cache_mode, cache_path)
    if !live
        @printf(io, "  open_workspace: %.2f s\n", open_seconds)
        mode == :gui && println(io, "  time to first frame: ", fmt(first_frame_seconds))
    end
    @printf(io, "measurement period: %.2f s\n", scan_seconds)
    @printf(
        io,
        "runtime compilation during the measurement, summed across Julia threads: %.2f s\n",
        compile_seconds,
    )
    println(io, "progress during this profile")
    println(io, "  stage                     start      end    delta     rate/s")
    print_stage_row(io, "sources found", first_sample, last_sample, :sources_found, scan_seconds)
    print_stage_row(io, "sources pending", first_sample, last_sample, :sources_pending, scan_seconds)
    print_stage_row(io, "cached sources", first_sample, last_sample, :cached_sources, scan_seconds)
    print_stage_row(io, "interpreted", first_sample, last_sample, :interpreted, scan_seconds)
    print_stage_row(io, "processed", first_sample, last_sample, :processed, scan_seconds)
    print_stage_row(io, "analyzed", first_sample, last_sample, :analyzed, scan_seconds)
    print_stage_row(
        io, "collection processed", first_sample, last_sample, :collection_processed, scan_seconds)
    print_stage_row(
        io, "collection analyzed", first_sample, last_sample, :collection_analyzed, scan_seconds)
    @printf(io, "profile samples: %d\n", sample_count)
    callback_seconds = sum(values(phases))
    if callback_seconds == 0
        println(io, "project callbacks: none ran during this profile")
    else
        println(io, "  summed elapsed time inside callbacks")
        println(io, "  callback                 seconds")
        for phase in (:detect, :read, :entries, :process, :analyze)
            @printf(io, "  %-22s %8.2f\n", "$(phase) callbacks", phases[phase])
        end
        @printf(io, "  %-22s %8.2f\n", "total", callback_seconds)
    end
    println(io, "artifacts: ", outdir)
    return nothing
end

revise!() = isdefined(Main, :Revise) ? Main.Revise.revise(; throw=true) : nothing

function validate_run_options(
    mode::Symbol,
    profiler::Symbol,
    cache_mode::Symbol,
    budget::Real,
)::Nothing
    mode in (:gui, :headless) || error("mode must be :gui or :headless")
    # :cpu (Profile.@profile) is not offered: the Mach-based CPU sampler in the Julia runtime can
    # suspend a thread mid-lock on macOS and wedge the whole process, and it does so intermittently.
    profiler in (:none, :wall) || error("profiler must be :none or :wall")
    cache_mode in (:fresh, :resume) || error("cache_mode must be :fresh or :resume")
    budget > 0 || error("budget must be positive")
    Threads.nthreads() > 1 || error(
        "Profiling requires multiple Julia threads; restart bench jmux with JULIA_NUM_THREADS=auto",
    )
    return nothing
end

"""
    start_live!(; mode=:gui, cache_mode=:resume, background_processing=false) -> LiveSession

Open one workspace and browser for repeated profiling and Revise-backed editing. The session keeps
the profiling lock until `stop_live!()` closes its browser, workspace, and lock.
"""
function start_live!(;
    mode::Symbol=:gui,
    cache_mode::Symbol=:resume,
    background_processing::Bool=false,
)
    LIVE_SESSION[] === nothing || error("A live profiling session is already open")
    validate_run_options(mode, :none, cache_mode, 1)
    revise!()
    profile_project = project()
    clear_callback_profile!(profile_project)
    run_lock = acquire_lock()
    opened = nothing
    try
        opened = open_profile_session(
            profile_project; mode, cache_mode, background_processing)
        session = LiveSession(
            profile_project,
            opened.workspace,
            opened.browser,
            mode,
            cache_mode,
            opened.cache_path,
            background_processing,
            opened.open_seconds,
            opened.first_frame_seconds,
            run_lock,
        )
        @printf(
            "DataBrowser live session: mode=%s cache=%s threads=%d background_processing=%s\n",
            mode,
            cache_mode,
            Threads.nthreads(),
            background_processing,
        )
        @printf("cache: %s\n", opened.cache_path)
        @printf("open_workspace: %.2f s\n", opened.open_seconds)
        mode === :gui && @printf("time to first frame: %.2f s\n", opened.first_frame_seconds)
        LIVE_SESSION[] = session
        return session
    catch
        try
            try
                opened !== nothing && opened.browser isa DataBrowser.BrowserSession &&
                    DataBrowser.close_browser!(opened.browser)
            finally
                opened === nothing || close_workspace!(opened.workspace)
            end
        finally
            close(run_lock)
        end
        rethrow()
    end
end

function live_session()::LiveSession
    session = LIVE_SESSION[]
    session isa LiveSession || error("No live profiling session; call start_live!() first")
    return session
end

"""Profile a fixed-duration window of the open live session and keep it running afterward."""
function profile_live!(; budget::Real=60, profiler::Symbol=:wall)
    session = live_session()
    validate_run_options(session.mode, profiler, session.cache_mode, budget)
    revise!()
    clear_callback_profile!(session.profile_project)
    outdir = abspath(joinpath(
        REPO_ROOT,
        "bench",
        "results",
        format(now(), "yyyymmdd-HHMMSS-sss") *
            "-ruo2-v2-live-$(session.mode)-$(profiler)-$(session.cache_mode)",
    ))
    mkpath(outdir)
    samples = Sample[]
    profile_data = nothing
    compile_before = first(Base.cumulative_compile_time_ns())
    profile_started_at = time()
    scan_started_at = time()
    if profiler === :wall
        Profile.clear()
        Profile.@profile_walltime sample_scan_window!(
            samples,
            session.workspace,
            scan_started_at,
            Float64(budget),
            false,
        )
        profile_data = Profile.retrieve()
    else
        sample_scan_window!(
            samples,
            session.workspace,
            scan_started_at,
            Float64(budget),
            false,
        )
    end
    scan_seconds = last(samples).t
    profile_seconds = time() - profile_started_at
    compile_seconds = (first(Base.cumulative_compile_time_ns()) - compile_before) / 1e9
    totals = callback_totals(session.profile_project)
    write_timeline(outdir, samples; restore_gl=session.mode === :gui)
    profile_outputs =
        profile_data === nothing ? nothing : save_profile_outputs(outdir, profiler, profile_data)
    count = profile_data === nothing ? 0 : profile_sample_count(first(profile_data))
    summary_path = joinpath(outdir, "summary.txt")
    open(summary_path, "w") do io
        print_run_summary(
            io;
            mode=session.mode,
            profiler,
            cache_mode=session.cache_mode,
            cache_path=session.cache_path,
            background_processing=session.background_processing,
            budget=Float64(budget),
            profile_seconds,
            open_seconds=0.0,
            first_frame_seconds=NaN,
            scan_seconds,
            compile_seconds,
            samples,
            sample_count=count,
            phases=totals,
            outdir,
            live=true,
        )
    end
    print(read(summary_path, String))
    return (;
        outdir,
        profile_outputs,
        cache_mode=session.cache_mode,
        cache_path=session.cache_path,
        background_processing=session.background_processing,
        profile_seconds,
        scan_seconds,
        compile_seconds,
        samples=count,
    )
end

"""Start a clean cache rebuild in the open live session without reopening Julia or the browser."""
function rebuild_live_cache!()::Nothing
    session = live_session()
    revise!()
    rebuild_cache!(session.workspace)
    session.cache_mode = :fresh
    return nothing
end

"""Replace the ProfileView window with the latest captured samples."""
function refresh_profileview!(; windowname::AbstractString="DataBrowser latest")
    profile_data = LAST_PROFILE[]
    profile_data === nothing && error("No saved profile; call profile_live!() first")
    data, lidict = profile_data
    ProfileView.closeall()
    return ProfileView.view(data; lidict, windowname)
end

"""Start or refresh the pprof web UI from the latest exported profile."""
function refresh_pprof!(; webhost::AbstractString="localhost", webport::Integer=57599)
    protobuf = LAST_PPROF[]
    protobuf === nothing && error("No saved profile; call profile_live!() first")
    return PProf.refresh(; webhost, webport, file=protobuf)
end

"""Close the live browser, workspace, and profiling lock."""
function stop_live!()::Nothing
    session = LIVE_SESSION[]
    session isa LiveSession || return nothing
    LIVE_SESSION[] = nothing
    try
        session.browser isa DataBrowser.BrowserSession &&
            DataBrowser.close_browser!(session.browser)
    finally
        try
            close_workspace!(session.workspace)
        finally
            close(session.run_lock)
        end
    end
    return nothing
end

"""
    profile_scan!(; mode=:gui, budget=60, profiler=:wall, cache_mode=:fresh,
                    background_processing=false) -> NamedTuple

Run a RuO2 scan in the persistent Revise session, write numeric/timeline/raw profile artifacts,
close all run-owned tasks, and preserve the dedicated profiling cache. `cache_mode=:fresh` passes
`rebuild=true`; `cache_mode=:resume` continues from the cache's current state.
"""
function profile_scan!(;
    mode::Symbol=:gui,
    budget::Real=60,
    profiler::Symbol=:wall,
    cache_mode::Symbol=:fresh,
    background_processing::Bool=false,
)
    validate_run_options(mode, profiler, cache_mode, budget)
    revise!()
    profile_project = project()
    clear_callback_profile!(profile_project)
    outdir = abspath(joinpath(
        REPO_ROOT,
        "bench",
        "results",
        format(now(), "yyyymmdd-HHMMSS-sss") *
            "-ruo2-v2-interactive-$(mode)-$(profiler)-$(cache_mode)",
    ))
    mkpath(outdir)
    run_lock = acquire_lock()
    workspace = nothing
    browser = nothing
    samples = Sample[]
    profile_data = nothing
    open_seconds = NaN
    first_frame_seconds = NaN
    scan_seconds = NaN
    compile_seconds = NaN
    profile_seconds = NaN
    cache_path = ""
    debug_timings = DataBrowser.DebugTimings(start_ns=time_ns())
    try
        compile_before = first(Base.cumulative_compile_time_ns())
        profile_started_at = time()
        result = DataBrowser.with_debug_timings(debug_timings) do
            profiled = nothing
            if profiler == :wall
                Profile.clear()
                Profile.@profile_walltime profiled = run_profile_window(
                    profile_project;
                    mode,
                    budget=Float64(budget),
                    cache_mode,
                    background_processing,
                )
                profile_data = Profile.retrieve()
            else
                profiled = run_profile_window(
                    profile_project;
                    mode,
                    budget=Float64(budget),
                    cache_mode,
                    background_processing,
                )
            end
            return profiled
        end
        profile_seconds = time() - profile_started_at
        compile_seconds = (first(Base.cumulative_compile_time_ns()) - compile_before) / 1e9
        workspace = result.workspace
        browser = result.browser
        cache_path = result.cache_path
        open_seconds = result.open_seconds
        first_frame_seconds = result.first_frame_seconds
        scan_seconds = result.scan_seconds
        samples = result.samples
        totals = callback_totals(profile_project)
        workspace = result.workspace
        browser = result.browser
        browser isa DataBrowser.BrowserSession && DataBrowser.close_browser!(browser)
        browser = nothing
        close_workspace!(workspace)
        workspace = nothing
        write_timeline(outdir, samples)
        DataBrowser.write_debug_timings(outdir, debug_timings)
        profile_outputs =
            profile_data === nothing ? nothing : save_profile_outputs(outdir, profiler, profile_data)
        count = profile_data === nothing ? 0 : profile_sample_count(first(profile_data))
        summary_path = joinpath(outdir, "summary.txt")
        open(summary_path, "w") do io
            print_run_summary(
                io;
                mode,
                profiler,
                cache_mode,
                cache_path,
                background_processing,
                budget=Float64(budget),
                profile_seconds,
                open_seconds,
                first_frame_seconds,
                scan_seconds,
                compile_seconds,
                samples,
                sample_count=count,
                phases=totals,
                outdir,
            )
        end
        print(read(summary_path, String))
        return (
            outdir,
            profile_outputs,
            debug_timings=(
                text=joinpath(outdir, "debug_timings.txt"),
                csv=joinpath(outdir, "debug_timings.csv"),
            ),
            cache_mode,
            cache_path,
            background_processing,
            profile_seconds,
            open_seconds,
            first_frame_seconds,
            scan_seconds,
            compile_seconds,
            samples=count,
        )
    finally
        browser isa DataBrowser.BrowserSession && DataBrowser.close_browser!(browser)
        workspace === nothing || close_workspace!(workspace)
        close(run_lock)
    end
end

end
