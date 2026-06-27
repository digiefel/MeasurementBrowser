"""Result delivered to one selected caller waiting for shared processing work."""
struct ProcessingResult
    item::Union{Nothing,AbstractDataItem}
    failure::Union{Nothing,CapturedException}
end

const PROCESSING_BACKGROUND_BATCH_ITEMS = 8

"""Start the bounded worker set that owns processing for one workspace."""
function start_processing_workers!(workspace::Workspace)::Nothing
    worker_count = max(1, Base.Threads.nthreads() ÷ 2)
    for _ in 1:worker_count
        worker = Base.Threads.@spawn processing_worker!(workspace)
        push!(workspace.processing.workers, worker)
    end
    return nothing
end

"""Stop accepting processing work and wait for all workers to release workspace state."""
function stop_processing_workers!(workspace::Workspace)::Nothing
    queue = workspace.processing
    lock(queue.lock) do
        queue.closed = true
        notify(queue.condition; all=true)
    end
    foreach(wait, queue.workers)
    return nothing
end

"""Cancel every waiting processing job while allowing currently executing callbacks to finish."""
function cancel_waiting_processing!(workspace::Workspace)::Nothing
    queue = workspace.processing
    canceled = ProcessingJob[]
    lock(queue.lock) do
        for (id, job) in collect(queue.jobs)
            job.state === :waiting || continue
            delete!(queue.jobs, id)
            push!(canceled, job)
        end
        queue.completed += length(canceled)
        notify(queue.condition; all=true)
    end
    failure = ProcessingResult(
        nothing,
        CapturedException(JobCancelled(), Any[]),
    )
    for job in canceled, waiter in job.waiters
        put!(waiter, failure)
    end
    return nothing
end

"""Reset idle processing bookkeeping before a clean profiling rebuild."""
function reset_processing_queue!(workspace::Workspace)::Nothing
    queue = workspace.processing
    lock(queue.lock) do
        isempty(queue.jobs) || error("Cannot reset processing while jobs are active")
        empty!(queue.selected)
        empty!(queue.background)
        queue.background_index = 1
        empty!(queue.source_locks)
        empty!(queue.completed_items)
        queue.total = 0
        queue.completed = 0
    end
    return nothing
end

"""
Hand one processed value to the cache buffer's single write path.

Selected work (`write_payload`) is a priority deposit that bypasses backpressure so interaction never
stalls; background work only stages statistics and yields to the ceiling like any other producer.
"""
function enqueue_processed_write!(
    workspace::Workspace,
    record::ItemRecord,
    item::Union{Nothing,AbstractDataItem},
    stats::Union{Nothing,MetadataDict},
    write_payload::Bool,
)::Nothing
    stage_processed!(
        workspace.buffer, record, item, stats, write_payload; priority=write_payload)
    return nothing
end

"""Return the next valid waiting job, preferring selected items."""
function take_processing_job!(queue::ProcessingQueue)::Union{Nothing,ProcessingJob}
    while true
        # Closing abandons queued work at once: a worker finishes only its current item, never the rest
        # of the queue. Without this, close_workspace! drains the whole backlog before workers exit, so
        # the app cannot be shut down during a build.
        queue.closed && return nothing
        while !isempty(queue.selected)
            id = pop!(queue.selected)
            job = get(queue.jobs, id, nothing)
            if job !== nothing && job.state == :waiting
                job.state = :running
                return job
            end
        end
        while queue.background_index <= length(queue.background)
            id = queue.background[queue.background_index]
            queue.background_index += 1
            job = get(queue.jobs, id, nothing)
            if job !== nothing && job.state == :waiting
                job.state = :running
                return job
            end
        end
        queue.closed && return nothing
        wait(queue.condition)
    end
end

"""Return the next valid jobs, batching adjacent background work from one source item."""
function take_processing_jobs!(queue::ProcessingQueue)::Vector{ProcessingJob}
    job = take_processing_job!(queue)
    job === nothing && return ProcessingJob[]
    job.priority > 0 && return ProcessingJob[job]
    !isempty(job.waiters) && return ProcessingJob[job]

    jobs = ProcessingJob[job]
    source_item_id_value = job.record.source_item_id
    while length(jobs) < PROCESSING_BACKGROUND_BATCH_ITEMS &&
          queue.background_index <= length(queue.background)
        id = queue.background[queue.background_index]
        next = get(queue.jobs, id, nothing)
        if next === nothing || next.state != :waiting
            queue.background_index += 1
            continue
        end
        next.record.source_item_id == source_item_id_value || break
        next.priority == 0 || break
        next.state = :running
        queue.background_index += 1
        push!(jobs, next)
    end
    return jobs
end

"""
Queue one item once, optionally promoting it and attaching a selected waiter.

The queue is unbounded and this never blocks: it holds only data-less records, so producing work
cannot stall the upstream scan. Processing reads each item's interpreted data back from the cache.
"""
function enqueue_processing!(
    workspace::Workspace,
    record::ItemRecord;
    selected::Bool=false,
)::Union{Nothing,Channel{Any}}
    queue = workspace.processing
    waiter = selected ? Channel{Any}(1) : nothing
    scan_id = workspace.scan.id
    skipped = lock(queue.lock) do
        queue.closed && error("Cannot queue processing after the workspace has closed")
        if !selected && get(queue.completed_items, record.id, 0) == scan_id
            true
        else
            job = get(queue.jobs, record.id, nothing)
            if job === nothing
                job = ProcessingJob(
                    record, :waiting, selected ? 1 : 0, Channel{Any}[], time_ns(), scan_id)
                queue.jobs[record.id] = job
                queue.total += 1
                push!(selected ? queue.selected : queue.background, record.id)
            elseif selected && job.state == :waiting && job.priority == 0
                job.priority = 1
                push!(queue.selected, record.id)
            end
            waiter === nothing || push!(job.waiters, waiter)
            notify(queue.condition)
            false
        end
    end
    skipped && return nothing
    return waiter
end

"""Return the lock that deduplicates source fallback for one source item."""
function source_fallback_lock(queue::ProcessingQueue, source_item_id_value::String)::ReentrantLock
    return lock(queue.lock) do
        get!(() -> ReentrantLock(), queue.source_locks, source_item_id_value)
    end
end

"""Interpret one source item through the normal path and return the requested logical item."""
function source_fallback(workspace::Workspace, record::ItemRecord)::AbstractDataItem
    fallback_lock = source_fallback_lock(workspace.processing, record.source_item_id)
    return lock(fallback_lock) do
        cached = only(buffer_read_item_data(
            workspace.buffer, [record]; stage=:interpreted))
        cached === nothing || return cached

        discovered = source_items(workspace.source)
        position = findfirst(
            item -> source_item_id(item) == record.source_item_id,
            discovered,
        )
        position === nothing && error(
            "Cannot load item '$(record.id)': source item '$(record.source_item_id)' " *
            "is no longer present in source '$(source_id(workspace.source))'",
        )
        interpretation = interpret_source_item(
            workspace.project,
            workspace.source,
            discovered[position],
            workspace.profiler,
        )
        stage_interpreted!(
            workspace.buffer,
            interpretation.records,
            Any[item for item in interpretation.interpreted_items],
        )
        requested = findfirst(item -> item.id == record.id, interpretation.records)
        requested === nothing && error(
            "Source item '$(record.source_item_id)' no longer produces item '$(record.id)'",
        )
        return interpretation.interpreted_items[requested]
    end
end

"""Load interpreted data from DuckDB, or the shared source fallback."""
function interpreted_item(workspace::Workspace, job::ProcessingJob)::AbstractDataItem
    cached = only(buffer_read_item_data(
        workspace.buffer, [job.record]; stage=:interpreted))
    return cached === nothing ? source_fallback(workspace, job.record) : cached
end

"""Load interpreted data for same-source background jobs with one cache read."""
function interpreted_items(
    workspace::Workspace,
    jobs::Vector{ProcessingJob},
)::Vector{AbstractDataItem}
    records = ItemRecord[job.record for job in jobs]
    cached = buffer_read_item_data(workspace.buffer, records; stage=:interpreted)
    interpreted = Vector{AbstractDataItem}(undef, length(jobs))
    for (index, job) in pairs(jobs)
        item = cached[index]
        interpreted[index] = item === nothing ? source_fallback(workspace, job.record) :
            item::AbstractDataItem
    end
    return interpreted
end

"""Process one item, persist eligible data and stats, and return the view item."""
function process_item(
    workspace::Workspace,
    job::ProcessingJob,
)::AbstractDataItem
    needs_payload = job.priority > 0 || !isempty(job.waiters)
    if needs_payload
        cached = only(buffer_read_item_data(
            workspace.buffer, [job.record]; stage=:processed))
        cached === nothing || return cached
    end
    return process_item(workspace, job, interpreted_item(workspace, job))
end

"""
Process one item when its interpreted data has already been loaded.

The processed-cache short-circuit lives in the one-arg method (and only single, non-batched jobs ever
need a payload — [`take_processing_jobs!`](@ref) never batches selected or waited-on work), so by the
time we reach here the payload must be computed.
"""
function process_item(
    workspace::Workspace,
    job::ProcessingJob,
    interpreted::AbstractDataItem,
)::AbstractDataItem
    process_started = time_ns()
    processed = Profiling.@profile_span workspace.profiler :project :process Profiling.ProfileAttributes(
        kind=job.record.kind,
        source_id=job.record.source_item_id,
        item_id=job.record.id,
    ) begin
        process(workspace.project, workspace.source, interpreted)
    end
    record_scan_phase!(workspace.project, job.record.source_item_id, job.record.kind,
        :process, (time_ns() - process_started) / 1e9, Base.Threads.threadid())
    view_item = DataItem(job.record, item_data(processed))
    item_stats = copy(job.record.stats)
    stats_failure = try
        stats_started = time_ns()
        computed = Profiling.@profile_span workspace.profiler :project :item_stats Profiling.ProfileAttributes(
            kind=job.record.kind,
            source_id=job.record.source_item_id,
            item_id=job.record.id,
        ) begin
            metadata_dict(compute_item_stats(
                workspace.project, workspace.source, view_item))
        end
        record_scan_phase!(workspace.project, job.record.source_item_id, job.record.kind,
            :stats, (time_ns() - stats_started) / 1e9, Base.Threads.threadid())
        merge!(item_stats, computed)
        nothing
    catch error
        CapturedException(error, catch_backtrace())
    end
    needs_payload = job.priority > 0 || !isempty(job.waiters)
    cache_payload = cacheable(processed) && needs_payload
    enqueue_processed_write!(
        workspace,
        job.record,
        cache_payload ? view_item : nothing,
        stats_failure === nothing && !isempty(item_stats) ? item_stats : nothing,
        needs_payload,
    )
    if stats_failure === nothing
        put!(workspace.processing.events, (
            kind=:stats,
            item_id=job.record.id,
            stats=item_stats,
        ))
    else
        failure = stats_failure::CapturedException
        message = "stats: " * sprint(showerror, failure.ex)
        stage_failure!(workspace.buffer, job.record, message)
        put!(workspace.processing.events, (
            kind=:failure,
            item_id=job.record.id,
            source_item_id=job.record.source_item_id,
            message,
        ))
    end
    return DataItem(ItemRecord(job.record; stats=item_stats), item_data(processed))
end

"""Run queued processing jobs and deliver results without retaining completed item objects."""
function finish_processing_job!(
    queue::ProcessingQueue,
    job::ProcessingJob,
    result::ProcessingResult,
)::Nothing
    if result.failure !== nothing
        put!(queue.events, (
            kind=:failure,
            item_id=job.record.id,
            source_item_id=job.record.source_item_id,
            message="process: " * sprint(showerror, result.failure.ex),
        ))
    end
    waiters = lock(queue.lock) do
        delete!(queue.jobs, job.record.id)
        queue.completed_items[job.record.id] = job.scan_id
        queue.completed += 1
        notify(queue.condition; all=true)
        copy(job.waiters)
    end
    for waiter in waiters
        put!(waiter, result)
    end
    return nothing
end

function processing_worker!(workspace::Workspace)::Nothing
    queue = workspace.processing
    while true
        jobs = lock(queue.lock) do
            take_processing_jobs!(queue)
        end
        isempty(jobs) && return nothing
        interpreted = try
            length(jobs) == 1 ? nothing : interpreted_items(workspace, jobs)
        catch error
            failure = ProcessingResult(nothing, CapturedException(error, catch_backtrace()))
            for job in jobs
                stage_failure!(
                    workspace.buffer,
                    job.record,
                    "process: " * sprint(showerror, error),
                )
                finish_processing_job!(queue, job, failure)
            end
            continue
        end
        for (index, job) in pairs(jobs)
            result = try
                item = Profiling.@profile_span workspace.profiler :processing :item Profiling.ProfileAttributes(
                    kind=job.record.kind,
                    source_id=job.record.source_item_id,
                    item_id=job.record.id,
                    wait_ns=Int64(time_ns() - job.queued_ns),
                ) begin
                    if interpreted === nothing
                        process_item(workspace, job)
                    else
                        process_item(workspace, job, interpreted[index])
                    end
                end
                ProcessingResult(item, nothing)
            catch error
                stage_failure!(
                    workspace.buffer,
                    job.record,
                    "process: " * sprint(showerror, error),
                )
                ProcessingResult(nothing, CapturedException(error, catch_backtrace()))
            end
            finish_processing_job!(queue, job, result)
        end
    end
end

"""Return processed items, sharing and promoting any work not already committed."""
function request_processed_items(
    workspace::Workspace,
    records::Vector{ItemRecord},
)::Vector{AbstractDataItem}
    cached = buffer_read_item_data(workspace.buffer, records; stage=:processed)
    loaded = Vector{AbstractDataItem}(undef, length(records))
    waiters = Tuple{Int,Channel{Any}}[]
    for (position, record) in pairs(records)
        item = cached[position]
        if item === nothing
            waiter = enqueue_processing!(workspace, record; selected=true)
            push!(waiters, (position, waiter::Channel{Any}))
        else
            item isa AbstractDataItem || error(
                "Processed cache returned $(typeof(item)) for item '$(record.id)'",
            )
            stats = get(workspace.index.item_stats, record.id, record.stats)
            loaded[position] = DataItem(
                ItemRecord(record; stats=stats),
                item_data(item),
            )
        end
    end
    for (position, waiter) in waiters
        result = take!(waiter)::ProcessingResult
        result.failure === nothing || throw(result.failure)
        loaded[position] = result.item::AbstractDataItem
    end
    return loaded
end
