"""Workspace-owned metrics, traces, CPU samples, and crash forensics."""
module Profiling

using JSON
using Profile
using Statistics: median, quantile

const MAX_PROFILE_EVENTS = 500_000
const COUNTER_INTERVAL_SECONDS = 0.1
const SAMPLING_TOTAL_BUFFER_SIZE = 2_000_000
const SAMPLING_DELAY_SECONDS = 0.01
const CURRENT_SPAN = Base.ScopedValues.ScopedValue(UInt64(0))
const CURRENT_SESSION = Base.ScopedValues.ScopedValue{Any}(nothing)

"""Darwin `mach_task_basic_info`; used only for current resident memory sampling."""
struct MachTaskBasicInfo
    virtual_size::UInt64
    resident_size::UInt64
    resident_size_max::UInt64
    user_seconds::Int32
    user_microseconds::Int32
    system_seconds::Int32
    system_microseconds::Int32
    policy::Int32
    suspend_count::Int32
end

"""Fixed trace dimensions retained without measurement values or metadata payloads."""
Base.@kwdef struct ProfileAttributes
    stage::Symbol = :none
    kind::Symbol = :none
    source_id::String = ""
    item_id::String = ""
    items::Int64 = 0
    rows::Int64 = 0
    bytes::Int64 = 0
    batch_size::Int64 = 0
    cache_items::Int64 = 0
    wait_ns::Int64 = 0
    service_ns::Int64 = 0
end

"""One completed structured span."""
struct ProfileEvent
    id::UInt64
    parent_id::UInt64
    category::Symbol
    operation::Symbol
    started_ns::UInt64
    duration_ns::UInt64
    task_id::UInt64
    start_thread::Int16
    end_thread::Int16
    status::Symbol
    attributes::ProfileAttributes
end

"""Process-wide counters sampled independently of overlapping spans."""
struct ProfileCounter
    sampled_ns::UInt64
    rss_bytes::Int64
    gc_total_ns::Int64
    gc_pause_count::Int64
    gc_full_sweeps::Int64
    gc_max_pause_ns::Int64
    gc_allocated_bytes::Int64
    scan_done::Int64
    scan_total::Int64
    processing_done::Int64
    processing_total::Int64
    queue_depth::Int64
    writer_busy_ns::Int64
    writer_wait_ns::Int64
end

"""One source line in a completed Julia CPU sampling profile."""
struct SamplingProfileRow
    samples::Int
    self_samples::Int
    function_name::String
    file::String
    line::Int
end

"""Reduced Julia CPU samples retained after raw profiler buffers are cleared."""
struct SamplingProfile
    total_samples::Int
    delay_seconds::Float64
    truncated::Bool
    rows::Vector{SamplingProfileRow}
end

"""Aggregated operation latency displayed by the internal profiler UI."""
struct ProfileSummaryRow
    category::Symbol
    operation::Symbol
    count::Int
    total_ms::Float64
    median_ms::Float64
    p90_ms::Float64
    p99_ms::Float64
    max_ms::Float64
    wait_ms::Float64
    service_ms::Float64
    mean_batch::Float64
    median_batch::Float64
    p90_batch::Float64
    max_batch::Float64
end

"""Immutable result of one stopped buffered profile session."""
struct ProfileReport
    started_ns::UInt64
    stopped_ns::UInt64
    events::Vector{ProfileEvent}
    counters::Vector{ProfileCounter}
    summary::Vector{ProfileSummaryRow}
    cpu::Union{Nothing,SamplingProfile}
    dropped_events::Int
end

"""Token used by manually instrumented spans whose final attributes are known at completion."""
struct ProfileSpanToken
    session::Any
    id::UInt64
    parent_id::UInt64
    category::Symbol
    operation::Symbol
    started_ns::UInt64
    task_id::UInt64
    start_thread::Int16
    attributes::ProfileAttributes
end

"""Mutable state for one workspace's opt-in internal profiler."""
mutable struct ProfileSession
    enabled::Bool
    cpu_enabled::Bool
    output_path::Union{Nothing,String}
    crash_path::Union{Nothing,String}
    crash_io::Union{Nothing,IO}
    crash_lock::ReentrantLock
    state::Symbol
    created_ns::UInt64
    started_ns::UInt64
    next_id::Base.Threads.Atomic{UInt64}
    next_slot::Base.Threads.Atomic{Int}
    dropped_events::Base.Threads.Atomic{Int}
    events::Vector{ProfileEvent}
    counters::Vector{ProfileCounter}
    counter_lock::ReentrantLock
    counter_stop::Base.Threads.Atomic{Bool}
    counter_task::Union{Nothing,Task}
    sampling_active::Bool
    report::Union{Nothing,ProfileReport}
    error::String
    export_sequence::Int
end

"""Parse one strict `0`/`1` environment flag."""
function environment_flag(name::AbstractString, default::Bool=false)::Bool
    value = get(ENV, String(name), nothing)
    value === nothing && return default
    value == "0" && return false
    value == "1" && return true
    throw(ArgumentError("Environment variable $name must be '0' or '1', got '$value'"))
end

"""Return an optional path from the environment without accepting an empty path."""
function environment_path(name::AbstractString)::Union{Nothing,String}
    value = get(ENV, String(name), nothing)
    value === nothing && return nothing
    isempty(value) && throw(ArgumentError("Environment variable $name cannot be empty"))
    return String(value)
end

"""Write one JSON object and flush it for crash durability."""
function _write_crash_line!(session::ProfileSession, value)::Nothing
    io = session.crash_io
    io === nothing && return nothing
    lock(session.crash_lock) do
        JSON.print(io, value)
        write(io, '\n')
        flush(io)
    end
    return nothing
end

"""Create one disabled or opt-in workspace profile session."""
function ProfileSession(
    enabled::Bool,
    cpu_enabled::Bool,
    output_path::Union{Nothing,AbstractString},
    crash_path::Union{Nothing,AbstractString},
)::ProfileSession
    cpu_enabled && !enabled && throw(ArgumentError(
        "profile_cpu=true requires profile_internal=true",
    ))
    output_path !== nothing && !enabled && throw(ArgumentError(
        "profile_output requires profile_internal=true",
    ))
    crash_string = crash_path === nothing ? nothing : String(crash_path)
    crash_io = if crash_string === nothing
        nothing
    else
        mkpath(dirname(crash_string))
        open(crash_string, "w")
    end
    now = time_ns()
    session = ProfileSession(
        enabled,
        cpu_enabled,
        output_path === nothing ? nothing : String(output_path),
        crash_string,
        crash_io,
        ReentrantLock(),
        enabled ? :idle : :disabled,
        now,
        UInt64(0),
        Base.Threads.Atomic{UInt64}(0),
        Base.Threads.Atomic{Int}(0),
        Base.Threads.Atomic{Int}(0),
        ProfileEvent[],
        ProfileCounter[],
        ReentrantLock(),
        Base.Threads.Atomic{Bool}(false),
        nothing,
        false,
        nothing,
        "",
        0,
    )
    _write_crash_line!(session, Dict(
        "type" => "session_start",
        "t_ns" => 0,
        "pid" => getpid(),
    ))
    return session
end

"""Return whether buffered spans or durable crash spans are currently requested."""
@inline function should_trace(session::ProfileSession)::Bool
    return session.state === :recording || session.crash_io !== nothing
end

@inline should_trace(::Nothing)::Bool = false

"""Return the profile session propagated through the current task scope."""
current_session()::Union{Nothing,ProfileSession} = CURRENT_SESSION[]

"""Return the number of complete buffered events recorded so far."""
function event_count(session::ProfileSession)::Int
    report = session.report
    report isa ProfileReport && session.state !== :recording && return length(report.events)
    return min(session.next_slot[], length(session.events))
end

"""Return a stable copy of the most recent completed buffered events."""
function recent_events(
    session::ProfileSession;
    limit::Integer=200,
)::Vector{ProfileEvent}
    limit >= 0 || throw(ArgumentError("Profile event limit must be nonnegative"))
    report = session.report
    source = report isa ProfileReport && session.state !== :recording ?
        report.events : session.events
    count = min(event_count(session), length(source))
    first_index = max(1, count - Int(limit) + 1)
    count == 0 && return ProfileEvent[]
    return copy(@view source[first_index:count])
end

"""Return the newest process counter from an active or completed capture."""
function latest_counter(session::ProfileSession)::Union{Nothing,ProfileCounter}
    report = session.report
    report isa ProfileReport && return isempty(report.counters) ? nothing : last(report.counters)
    return lock(session.counter_lock) do
        isempty(session.counters) ? nothing : last(session.counters)
    end
end

"""Start one span and emit its crash-durable begin record when configured."""
function start_span!(
    session::ProfileSession,
    category::Symbol,
    operation::Symbol,
    attributes::ProfileAttributes=ProfileAttributes(),
)::ProfileSpanToken
    id = Base.Threads.atomic_add!(session.next_id, UInt64(1)) + UInt64(1)
    parent_id = CURRENT_SPAN[]
    started_ns = time_ns()
    token = ProfileSpanToken(
        session,
        id,
        parent_id,
        category,
        operation,
        started_ns,
        UInt64(objectid(current_task())),
        Int16(Base.Threads.threadid()),
        attributes,
    )
    _write_crash_line!(session, Dict(
        "type" => "span_start",
        "id" => id,
        "parent_id" => parent_id,
        "category" => String(category),
        "operation" => String(operation),
        "t_ns" => Int64(started_ns - session.created_ns),
        "thread" => Int(token.start_thread),
    ))
    return token
end

"""Finish one span, preserving the supplied final status and attributes."""
function finish_span!(
    token::ProfileSpanToken;
    status::Symbol=:ok,
    attributes::ProfileAttributes=token.attributes,
)::Nothing
    session = token.session::ProfileSession
    stopped_ns = time_ns()
    duration_ns = stopped_ns - token.started_ns
    if session.state === :recording
        slot = Base.Threads.atomic_add!(session.next_slot, 1) + 1
        if slot <= length(session.events)
            session.events[slot] = ProfileEvent(
                token.id,
                token.parent_id,
                token.category,
                token.operation,
                token.started_ns,
                duration_ns,
                token.task_id,
                token.start_thread,
                Int16(Base.Threads.threadid()),
                status,
                attributes,
            )
        else
            Base.Threads.atomic_add!(session.dropped_events, 1)
        end
    end
    _write_crash_line!(session, Dict(
        "type" => "span_end",
        "id" => token.id,
        "status" => String(status),
        "t_ns" => Int64(stopped_ns - session.created_ns),
        "duration_ns" => Int64(duration_ns),
        "thread" => Base.Threads.threadid(),
    ))
    return nothing
end

"""Run `work` inside one nested structured span and preserve its exception."""
function with_profile_span(
    work::Function,
    session::ProfileSession,
    category::Symbol,
    operation::Symbol,
    attributes::ProfileAttributes,
)::Any
    token = start_span!(session, category, operation, attributes)
    status = :ok
    try
        return Base.ScopedValues.with(
            CURRENT_SPAN => token.id,
            CURRENT_SESSION => session,
        ) do
            work()
        end
    catch
        status = :error
        rethrow()
    finally
        finish_span!(token; status)
    end
end

"""Trace `expr` only when buffered or crash-durable tracing is active."""
macro profile_span(session, category, operation, attributes, expr)
    quote
        local _profile_session = $(esc(session))
        if Profiling.should_trace(_profile_session)
            Profiling.with_profile_span(
                _profile_session,
                $(esc(category)),
                $(esc(operation)),
                $(esc(attributes)),
            ) do
                $(esc(expr))
            end
        else
            $(esc(expr))
        end
    end
end

"""Start Julia's bounded all-thread CPU sampling profiler."""
function start_sampling!()::Nothing
    Profile.is_running() && error("Julia's CPU profiler is already running")
    buffer_per_thread = max(50_000,
        cld(SAMPLING_TOTAL_BUFFER_SIZE, Base.Threads.nthreads()))
    Profile.init(n=buffer_per_thread, delay=SAMPLING_DELAY_SECONDS)
    Profile.clear()
    Profile.start_timer()
    return nothing
end

"""Stop CPU sampling and reduce raw stacks to flat source-line counts."""
function stop_sampling!()::SamplingProfile
    running = Profile.is_running()
    truncated = Profile.is_buffer_full()
    (running || truncated) || error("Julia's CPU profiler is not running")
    running && Profile.stop_timer()
    data, line_info = Profile.retrieve(include_meta=false, limitwarn=false)
    flat_data, flat_info = Profile.flatten(data, line_info)
    counts = Dict{Tuple{String,String,Int},Int}()
    self_counts = Dict{Tuple{String,String,Int},Int}()
    frames = Tuple{String,String,Int}[]
    total_samples = 0
    for instruction_pointer in flat_data
        if instruction_pointer == 0
            isempty(frames) && continue
            total_samples += 1
            for frame in unique(frames)
                counts[frame] = get(counts, frame, 0) + 1
            end
            leaf = frames[1]
            self_counts[leaf] = get(self_counts, leaf, 0) + 1
            empty!(frames)
            continue
        end
        frame = flat_info[instruction_pointer]
        frame.from_c && continue
        push!(frames, (String(frame.func), String(frame.file), frame.line))
    end
    sample_rows = SamplingProfileRow[
        SamplingProfileRow(samples, get(self_counts, frame, 0), frame[1], frame[2], frame[3])
        for (frame, samples) in counts
    ]
    sort!(sample_rows; by=row -> (row.self_samples, row.samples), rev=true)
    result = SamplingProfile(total_samples, SAMPLING_DELAY_SECONDS, truncated, sample_rows)
    empty!(data)
    empty!(flat_data)
    empty!(line_info)
    empty!(flat_info)
    Profile.clear()
    GC.gc(true)
    return result
end

"""Cancel CPU sampling and release raw profiler buffers."""
function cancel_sampling!()::Nothing
    Profile.is_running() && Profile.stop_timer()
    Profile.clear()
    GC.gc(true)
    return nothing
end

"""Read current process resident memory without substituting peak RSS."""
function process_rss_bytes()::Int64
    if Sys.isapple()
        info = Ref(MachTaskBasicInfo(0, 0, 0, 0, 0, 0, 0, 0, 0))
        count = Ref{UInt32}(UInt32(div(sizeof(MachTaskBasicInfo), sizeof(Int32))))
        task = ccall(:mach_task_self, UInt32, ())
        result = ccall(
            :task_info,
            Int32,
            (UInt32, Int32, Ref{MachTaskBasicInfo}, Ref{UInt32}),
            task,
            20,
            info,
            count,
        )
        result == 0 || error("mach task_info failed while sampling RSS with code $result")
        return Int64(info[].resident_size)
    end
    if Sys.islinux()
        fields = split(read("/proc/self/statm", String))
        length(fields) >= 2 || error("/proc/self/statm does not contain a resident-page field")
        return parse(Int64, fields[2]) * Int64(Sys.page_size())
    end
    error("Current RSS sampling is unsupported on $(Sys.KERNEL)")
end

"""Capture one process-wide counter sample plus one workspace snapshot."""
function _capture_counter(session::ProfileSession, snapshot::Function)::Nothing
    session.state === :recording || return nothing
    values = snapshot()
    gc = Base.gc_num()
    sample = ProfileCounter(
        time_ns(),
        process_rss_bytes(),
        Int64(gc.total_time),
        Int64(gc.pause),
        Int64(gc.full_sweep),
        Int64(gc.max_pause),
        Int64(gc.total_allocd),
        Int64(values.scan_done),
        Int64(values.scan_total),
        Int64(values.processing_done),
        Int64(values.processing_total),
        Int64(values.queue_depth),
        Int64(values.writer_busy_ns),
        Int64(values.writer_wait_ns),
    )
    lock(session.counter_lock) do
        push!(session.counters, sample)
    end
    return nothing
end

"""Fail before any application work is canceled if a profile cannot start."""
function validate_start(session::ProfileSession)::Nothing
    session.enabled || error("Internal profiling is disabled for this workspace")
    session.state in (:idle, :complete) || error(
        "Cannot start internal profiling while session state is $(session.state)",
    )
    session.cpu_enabled && Profile.is_running() && error(
        "Cannot start internal profiling because Julia's CPU profiler is already running",
    )
    return nothing
end

"""Start a fresh bounded buffered profiling session."""
function start!(session::ProfileSession, snapshot::Function)::Nothing
    validate_start(session)
    session.events = Vector{ProfileEvent}(undef, MAX_PROFILE_EVENTS)
    empty!(session.counters)
    session.next_slot[] = 0
    session.dropped_events[] = 0
    session.report = nothing
    session.error = ""
    session.started_ns = time_ns()
    session.counter_stop[] = false
    session.state = :recording
    try
        if session.cpu_enabled
            start_sampling!()
            session.sampling_active = true
        end
        _capture_counter(session, snapshot)
        session.counter_task = Base.Threads.@spawn begin
            while !session.counter_stop[]
                sleep(COUNTER_INTERVAL_SECONDS)
                _capture_counter(session, snapshot)
            end
        end
    catch error
        session.sampling_active && cancel_sampling!()
        session.sampling_active = false
        session.error = "Profile start failed: " * sprint(showerror, error)
        session.state = :error
        rethrow()
    end
    return nothing
end

"""Aggregate completed events by category and operation."""
function _summarize(events::Vector{ProfileEvent})::Vector{ProfileSummaryRow}
    groups = Dict{Tuple{Symbol,Symbol},Vector{ProfileEvent}}()
    for event in events
        push!(get!(() -> ProfileEvent[], groups, (event.category, event.operation)), event)
    end
    rows = ProfileSummaryRow[]
    for ((category, operation), group) in groups
        durations = Float64[event.duration_ns / 1e6 for event in group]
        waits = sum(event.attributes.wait_ns for event in group) / 1e6
        service = sum(event.attributes.service_ns for event in group) / 1e6
        batches = Int64[event.attributes.batch_size for event in group
                        if event.attributes.batch_size > 0]
        push!(rows, ProfileSummaryRow(
            category,
            operation,
            length(group),
            sum(durations),
            median(durations),
            quantile(durations, 0.9),
            quantile(durations, 0.99),
            maximum(durations),
            waits,
            service,
            isempty(batches) ? 0.0 : sum(batches) / length(batches),
            isempty(batches) ? 0.0 : median(batches),
            isempty(batches) ? 0.0 : quantile(batches, 0.9),
            isempty(batches) ? 0.0 : maximum(batches),
        ))
    end
    sort!(rows; by=row -> row.total_ms, rev=true)
    return rows
end

"""Stop a recording and return its immutable report."""
function stop!(session::ProfileSession)::ProfileReport
    session.state === :recording || error(
        "Cannot stop internal profiling while session state is $(session.state)",
    )
    session.counter_stop[] = true
    task = session.counter_task
    try
        task === nothing || wait(task)
    catch error
        session.sampling_active && cancel_sampling!()
        session.sampling_active = false
        session.error = "Counter sampling failed: " * sprint(showerror, error)
        session.state = :error
        rethrow()
    end
    session.counter_task = nothing
    cpu = if session.sampling_active
        result = stop_sampling!()
        session.sampling_active = false
        result
    else
        nothing
    end
    count = event_count(session)
    events = copy(@view session.events[1:count])
    counters = lock(session.counter_lock) do
        copy(session.counters)
    end
    report = ProfileReport(
        session.started_ns,
        time_ns(),
        events,
        counters,
        _summarize(events),
        cpu,
        session.dropped_events[],
    )
    session.report = report
    session.events = ProfileEvent[]
    session.counters = ProfileCounter[]
    session.state = :complete
    if session.output_path !== nothing
        try
            export!(session, session.output_path; numbered=true)
        catch error
            session.error = "Profile export failed: " * sprint(showerror, error)
            session.state = :error
            rethrow()
        end
    end
    return report
end

"""Reset a stopped report while preserving workspace profiling configuration."""
function reset!(session::ProfileSession)::Nothing
    session.state === :recording && error("Stop internal profiling before resetting it")
    session.report = nothing
    session.events = ProfileEvent[]
    session.counters = ProfileCounter[]
    session.next_slot[] = 0
    session.dropped_events[] = 0
    session.error = ""
    session.state = session.enabled ? :idle : :disabled
    return nothing
end

"""Return one suffix-safe output path for repeated automatic exports."""
function _output_path(session::ProfileSession, path::String, numbered::Bool)::String
    numbered || return path
    session.export_sequence += 1
    session.export_sequence == 1 && return path
    stem, extension = splitext(path)
    return string(stem, "-", session.export_sequence, isempty(extension) ? ".json" : extension)
end

"""Convert fixed trace dimensions to Perfetto event arguments."""
function _event_arguments(event::ProfileEvent)::Dict{String,Any}
    attrs = event.attributes
    result = Dict{String,Any}(
        "span_id" => event.id,
        "parent_id" => event.parent_id,
        "task_id" => event.task_id,
        "end_thread" => event.end_thread,
        "status" => String(event.status),
    )
    attrs.stage === :none || (result["stage"] = String(attrs.stage))
    attrs.kind === :none || (result["kind"] = String(attrs.kind))
    isempty(attrs.source_id) || (result["source_id"] = attrs.source_id)
    isempty(attrs.item_id) || (result["item_id"] = attrs.item_id)
    for (name, value) in (
        "items" => attrs.items,
        "rows" => attrs.rows,
        "bytes" => attrs.bytes,
        "batch_size" => attrs.batch_size,
        "cache_items" => attrs.cache_items,
        "wait_ns" => attrs.wait_ns,
        "service_ns" => attrs.service_ns,
    )
        value == 0 || (result[name] = value)
    end
    return result
end

"""Build a Chrome/Perfetto trace object from one report."""
function _perfetto_trace(report::ProfileReport)::Dict{String,Any}
    pid = getpid()
    trace_events = Any[]
    threads = sort!(unique(Int(event.start_thread) for event in report.events))
    for thread in threads
        push!(trace_events, Dict(
            "ph" => "M", "pid" => pid, "tid" => thread,
            "name" => "thread_name", "args" => Dict("name" => "Julia thread $thread"),
        ))
    end
    for event in report.events
        push!(trace_events, Dict(
            "ph" => "X",
            "pid" => pid,
            "tid" => Int(event.start_thread),
            "cat" => String(event.category),
            "name" => String(event.operation),
            "ts" => (event.started_ns - report.started_ns) / 1e3,
            "dur" => event.duration_ns / 1e3,
            "args" => _event_arguments(event),
        ))
    end
    for counter in report.counters
        push!(trace_events, Dict(
            "ph" => "C", "pid" => pid, "tid" => 0,
            "cat" => "process", "name" => "workspace",
            "ts" => (counter.sampled_ns - report.started_ns) / 1e3,
            "args" => Dict(
                "rss_bytes" => counter.rss_bytes,
                "gc_total_ns" => counter.gc_total_ns,
                "gc_pause_count" => counter.gc_pause_count,
                "gc_full_sweeps" => counter.gc_full_sweeps,
                "gc_max_pause_ns" => counter.gc_max_pause_ns,
                "gc_allocated_bytes" => counter.gc_allocated_bytes,
                "scan_done" => counter.scan_done,
                "scan_total" => counter.scan_total,
                "processing_done" => counter.processing_done,
                "processing_total" => counter.processing_total,
                "queue_depth" => counter.queue_depth,
                "writer_busy_ns" => counter.writer_busy_ns,
                "writer_wait_ns" => counter.writer_wait_ns,
            ),
        ))
    end
    return Dict(
        "traceEvents" => trace_events,
        "displayTimeUnit" => "ms",
        "metadata" => Dict(
            "dropped_events" => report.dropped_events,
            "cpu_samples" => report.cpu === nothing ? 0 : report.cpu.total_samples,
        ),
    )
end

"""Export the stopped report as validated Chrome/Perfetto JSON."""
function export!(
    session::ProfileSession,
    path::AbstractString;
    numbered::Bool=false,
)::String
    report = session.report
    report isa ProfileReport || error("Stop internal profiling before exporting it")
    destination = _output_path(session, String(path), numbered)
    mkpath(dirname(destination))
    temporary = tempname(dirname(destination))
    try
        open(temporary, "w") do io
            JSON.print(io, _perfetto_trace(report))
        end
        parsed = JSON.parsefile(temporary)
        haskey(parsed, "traceEvents") || error(
            "Generated profiler output is missing traceEvents",
        )
        mv(temporary, destination; force=true)
    catch
        rm(temporary; force=true)
        rethrow()
    end
    return destination
end

"""Close CPU sampling and crash output owned by one workspace session."""
function close!(session::ProfileSession)::Nothing
    if session.state === :recording
        stop!(session)
    end
    _write_crash_line!(session, Dict(
        "type" => "session_end",
        "t_ns" => Int64(time_ns() - session.created_ns),
    ))
    io = session.crash_io
    session.crash_io = nothing
    io === nothing || close(io)
    return nothing
end

end # module Profiling
