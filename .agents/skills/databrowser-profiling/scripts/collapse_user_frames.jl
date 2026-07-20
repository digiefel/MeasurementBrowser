"""
Collapse user-callback frames in a saved profile at the engine boundary.

Usage:

    julia --project=bench .agents/skills/databrowser-profiling/scripts/collapse_user_frames.jl \\
        <results-dir or *_profile.jls>

Reads the run's `*_profile.jls`, replaces every stretch of sampled frames that runs inside a
registered project callback with one opaque frame named after the engine call site, and writes
next to the input:

- `collapsed_<name>.pb.gz` — the collapsed profile for pprof; engine frames keep full detail,
  user code is a single frame per callback kind;
- `collapsed_<name>_summary.txt` — sample classification (user callback / app active / task
  waiting) overall and per thread, plus the app-side flat hotspots.

The collapse rule is deterministic: walking each stack from its root, once a frame of a known
callback-invoking engine function (`project_engine.jl`) is passed, every deeper frame whose file
is not under `lib/DataBrowser*` is user code and is folded into the nearest callback frame.
Engine frames reached from inside user code (callbacks calling back into the API) stay visible.
"""
module CollapseUserFrames

using Profile
using PProf

include(joinpath(@__DIR__, "profile_artifacts.jl"))
using .ProfileArtifacts: load_profile

# Engine functions in project_engine.jl that directly invoke registered project callbacks.
const BOUNDARY_FUNCS = (
    "_detect_recipe", "data_items", "_registered_item", "_processed_item", "process",
    "_process_collection", "_analyze_collection", "_analyze_item",
)

# Leaf-side functions that mean the sampled task was parked or its OS thread was blocked idle.
# Contended-lock spinning (_jl_mutex_wait) is deliberately NOT here: that is app waste to surface.
const WAIT_FUNCS = (
    "poptask", "wait", "task_done_hook", "jl_swap_fiber", "jl_start_fiber_swap",
    "jl_set_fiber", "start_task", "yieldto",
    "__psynch_cvwait", "kevent", "mach_msg2_trap", "usleep", "select",
)

const SYNTHETIC_BASE = 0xfeed_0000_0000_0000

_matches_func(name::AbstractString, target::AbstractString) =
    name == target || name == "#" * target ||
    startswith(name, target * "#") || startswith(name, "#" * target * "#")

function _boundary_func(frames::Vector{Base.StackTraces.StackFrame})::Union{Nothing,String}
    for frame in frames
        file = String(frame.file)
        endswith(file, "project_engine.jl") || continue
        name = String(frame.func)
        for target in BOUNDARY_FUNCS
            _matches_func(name, target) && return target
        end
    end
    return nothing
end

_is_engine(frames::Vector{Base.StackTraces.StackFrame})::Bool =
    any(frame -> occursin("lib/DataBrowser", String(frame.file)), frames)

_is_wait(frames::Vector{Base.StackTraces.StackFrame})::Bool =
    any(frame -> String(frame.func) in WAIT_FUNCS, frames)

struct Block
    ips::Vector{UInt64}      # leaf-first, as stored in the profile buffer
    meta::Vector{UInt64}     # the trailing meta words and null terminators, verbatim
    threadid::Int
end

"""Split a metadata-carrying profile buffer into per-sample blocks."""
function parse_blocks(data::Vector{UInt64})::Vector{Block}
    Profile.has_meta(data) || error(
        "profile buffer has no per-sample metadata; re-capture with a current Julia")
    blocks = Block[]
    start = 1
    for i in eachindex(data)
        Profile.is_block_end(data, i) || continue
        ips = data[start:i-6]
        meta = data[i-5:i]
        push!(blocks, Block(ips, meta, Int(data[i-Profile.META_OFFSET_THREADID])))
        start = i + 1
    end
    return blocks
end

"""Collapse one block's stack; returns the new leaf-first ips and the callback label hit."""
function collapse_ips(
    ips::Vector{UInt64},
    lidict::Profile.LineInfoDict,
    synthetic::Dict{String,UInt64},
)::Tuple{Vector{UInt64},Union{Nothing,String}}
    rooted = reverse(ips)
    frames_of = ip -> get(lidict, ip, Base.StackTraces.StackFrame[])
    boundary_index = findfirst(ip -> _boundary_func(frames_of(ip)) !== nothing, rooted)
    boundary_index === nothing && return ips, nothing
    out = rooted[1:boundary_index]
    current = _boundary_func(frames_of(rooted[boundary_index]))::String
    user_label = nothing
    run_open = false
    for j in boundary_index+1:length(rooted)
        ip = rooted[j]
        frames = frames_of(ip)
        inner = _boundary_func(frames)
        inner === nothing || (current = inner)
        if _is_engine(frames)
            push!(out, ip)
            run_open = false
        else
            user_label = current
            if !run_open
                key = "user code in " * current * " callback"
                fake = get!(synthetic, key) do
                    SYNTHETIC_BASE + UInt64(length(synthetic) + 1)
                end
                push!(out, fake)
                run_open = true
            end
        end
    end
    return reverse(out), user_label
end

function collapse(path::AbstractString)
    input = if isdir(path)
        candidates = filter(f -> endswith(f, "_profile.jls"), readdir(path))
        isempty(candidates) && error("no *_profile.jls in $path")
        joinpath(path, first(sort(candidates)))
    else
        path
    end
    payload = load_profile(input)
    blocks = parse_blocks(Vector{UInt64}(payload.data))
    lidict = payload.lidict

    synthetic = Dict{String,UInt64}()
    newdata = UInt64[]
    class_counts = Dict("user callback" => 0, "app active" => 0, "task waiting" => 0)
    per_thread = Dict{Int,Dict{String,Int}}()
    per_label = Dict{String,Int}()
    app_flat = Dict{String,Int}()

    for block in blocks
        ips, label = collapse_ips(block.ips, lidict, synthetic)
        append!(newdata, ips)
        append!(newdata, block.meta)
        leafward = [get(lidict, ip, Base.StackTraces.StackFrame[])
                    for ip in block.ips[1:min(4, length(block.ips))]]
        class = if label !== nothing
            per_label[label] = get(per_label, label, 0) + 1
            "user callback"
        elseif any(_is_wait, leafward)
            "task waiting"
        else
            if !isempty(leafward) && !isempty(leafward[1])
                leaf = first(leafward[1])
                key = String(leaf.func) * "  (" * basename(String(leaf.file)) * ")"
                app_flat[key] = get(app_flat, key, 0) + 1
            end
            "app active"
        end
        class_counts[class] += 1
        thread = get!(per_thread, block.threadid) do
            Dict("user callback" => 0, "app active" => 0, "task waiting" => 0)
        end
        thread[class] += 1
    end

    newlidict = Profile.LineInfoDict()
    for (ip, frames) in lidict
        newlidict[ip] = frames
    end
    for (label, fake) in synthetic
        newlidict[fake] = [Base.StackTraces.StackFrame(
            Symbol(label), Symbol("[registered project callback]"), 0,
            nothing, false, false, fake)]
    end

    dir = dirname(abspath(input))
    stem = replace(basename(input), r"\.jls$" => "")
    out_pb = joinpath(dir, "collapsed_" * stem * ".pb.gz")
    out_txt = joinpath(dir, "collapsed_" * stem * "_summary.txt")

    PProf.pprof(
        newdata, newlidict;
        web=false, out=out_pb, sampling_delay=UInt64(payload.sampling_delay_ns),
        full_signatures=true,
    )

    total = length(blocks)
    open(out_txt, "w") do io
        println(io, "Collapsed profile summary — $(basename(input)), $total samples")
        println(io)
        println(io, "Sample classes (samples are tasks in a wall profile: a waiting task on a")
        println(io, "busy thread wastes nothing, so these counts locate work, they are NOT")
        println(io, "capacity fractions — waste share comes from the callback timers):")
        for (name, count) in sort(collect(class_counts); by=last, rev=true)
            println(io, rpad("  " * name, 18), count, "  (",
                round(100 * count / max(total, 1); digits=1), "%)")
        end
        println(io)
        println(io, "User-callback samples by boundary:")
        for (name, count) in sort(collect(per_label); by=last, rev=true)
            println(io, rpad("  " * name, 24), count)
        end
        println(io)
        println(io, "Per thread (0-based ids as recorded by the profiler):")
        for threadid in sort(collect(keys(per_thread)))
            counts = per_thread[threadid]
            println(io, "  thread ", threadid, ": user=", counts["user callback"],
                " app=", counts["app active"], " waiting=", counts["task waiting"])
        end
        println(io)
        println(io, "App-side flat hotspots (leaf frames of app-active samples):")
        for (name, count) in first(sort(collect(app_flat); by=last, rev=true),
                min(25, length(app_flat)))
            println(io, rpad("  " * name, 60), count)
        end
    end
    return (; protobuf=out_pb, summary=out_txt)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("usage: collapse_user_frames.jl <results-dir or *_profile.jls>")
    paths = CollapseUserFrames.collapse(ARGS[1])
    println("wrote ", paths.protobuf)
    println("wrote ", paths.summary)
end
