"""Result delivered to one selected caller waiting for shared processing work."""
struct ProcessingResult
    item::Union{Nothing,AbstractDataItem}
    failure::Union{Nothing,CapturedException}
end

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

"""Commit one processed value through the processing queue's shared cache writer."""
function commit_processed_item!(
    workspace::Workspace,
    record::ItemRecord,
    item::Union{Nothing,AbstractDataItem},
    stats::Union{Nothing,MetadataDict},
)::Union{Nothing,CapturedException}
    queue = workspace.processing
    request = ProcessedWriteRequest(
        record, item, stats, Channel{ProcessedWriteResult}(1))
    leader = lock(queue.lock) do
        push!(queue.pending_writes, request)
        if queue.write_active
            false
        else
            queue.write_active = true
            true
        end
    end
    if leader
        sleep(0.001)
        while true
            batch = lock(queue.lock) do
                ready = copy(queue.pending_writes)
                empty!(queue.pending_writes)
                ready
            end
            payload_failure = try
                write_started = time_ns()
                write_cached_item_data!(
                    workspace.cache.db,
                    ItemRecord[pending.record for pending in batch],
                    Any[pending.item for pending in batch];
                    stage=:processed,
                )
                record_cache_phase!(
                    workspace.monitor.processed_write_ns,
                    workspace.monitor.processed_writes, write_started)
                nothing
            catch error
                CapturedException(error, catch_backtrace())
            end
            stats_failure = if payload_failure === nothing
                batch_stats = Dict{String,MetadataDict}(
                    pending.record.id => (pending.stats::MetadataDict)
                    for pending in batch if pending.stats !== nothing
                )
                try
                    if !isempty(batch_stats)
                        write_started = time_ns()
                        persist_stats!(
                            workspace.cache.db,
                            batch_stats,
                            Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}(),
                        )
                        record_cache_phase!(
                            workspace.monitor.stats_write_ns,
                            workspace.monitor.stats_writes, write_started)
                    end
                    nothing
                catch error
                    CapturedException(error, catch_backtrace())
                end
            else
                nothing
            end
            result = ProcessedWriteResult(payload_failure, stats_failure)
            for pending in batch
                put!(pending.completion, result)
            end
            done = lock(queue.lock) do
                if isempty(queue.pending_writes)
                    queue.write_active = false
                    true
                else
                    false
                end
            end
            done && break
        end
    end
    result = take!(request.completion)
    result.payload_failure === nothing || throw(result.payload_failure)
    return result.stats_failure
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
    lock(queue.lock) do
        queue.closed && error("Cannot queue processing after the workspace has closed")
        job = get(queue.jobs, record.id, nothing)
        if job === nothing
            job = ProcessingJob(record, :waiting, selected ? 1 : 0, Channel{Any}[])
            queue.jobs[record.id] = job
            queue.total += 1
            push!(selected ? queue.selected : queue.background, record.id)
        elseif selected && job.state == :waiting && job.priority == 0
            job.priority = 1
            push!(queue.selected, record.id)
        end
        waiter === nothing || push!(job.waiters, waiter)
        notify(queue.condition)
    end
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
        cached = only(read_cached_item_data(
            workspace.cache.db, [record]; stage=:interpreted))
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
            workspace.project, workspace.source, discovered[position])
        stored = Any[
            cacheable(item) ? item : nothing
            for item in interpretation.interpreted_items
        ]
        reconcile_source_items!(
            workspace.cache.db,
            [interpretation.records],
            [stored];
            stage=:interpreted,
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
    cached = only(read_cached_item_data(
        workspace.cache.db, [job.record]; stage=:interpreted))
    return cached === nothing ? source_fallback(workspace, job.record) : cached
end

"""Process one item, persist eligible data and stats, and return the view item."""
function process_item(workspace::Workspace, job::ProcessingJob)::AbstractDataItem
    cached = only(read_cached_item_data(
        workspace.cache.db, [job.record]; stage=:processed))
    cached === nothing || return cached

    interpreted = interpreted_item(workspace, job)
    process_started = time_ns()
    processed = process(workspace.project, workspace.source, interpreted)
    record_scan_phase!(workspace.project, job.record.source_item_id, job.record.kind,
        :process, (time_ns() - process_started) / 1e9, Base.Threads.threadid())
    view_item = DataItem(job.record, item_data(processed))
    item_stats = copy(job.record.stats)
    stats_failure = try
        stats_started = time_ns()
        computed = metadata_dict(compute_item_stats(
            workspace.project, workspace.source, view_item))
        record_scan_phase!(workspace.project, job.record.source_item_id, job.record.kind,
            :stats, (time_ns() - stats_started) / 1e9, Base.Threads.threadid())
        merge!(item_stats, computed)
        nothing
    catch error
        CapturedException(error, catch_backtrace())
    end
    persisted_failure = commit_processed_item!(
        workspace,
        job.record,
        cacheable(processed) ? view_item : nothing,
        stats_failure === nothing && !isempty(item_stats) ? item_stats : nothing,
    )
    stats_failure === nothing && (stats_failure = persisted_failure)
    if stats_failure === nothing
        put!(workspace.processing.events, (
            kind=:stats,
            item_id=job.record.id,
            stats=item_stats,
        ))
    else
        failure = stats_failure::CapturedException
        message = "stats: " * sprint(showerror, failure.ex)
        persist_item_failure!(workspace.cache.db, job.record, message)
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
function processing_worker!(workspace::Workspace)::Nothing
    queue = workspace.processing
    while true
        job = lock(queue.lock) do
            take_processing_job!(queue)
        end
        job === nothing && return nothing
        result = try
            ProcessingResult(process_item(workspace, job), nothing)
        catch error
            persist_item_failure!(
                workspace.cache.db,
                job.record,
                "process: " * sprint(showerror, error),
            )
            ProcessingResult(nothing, CapturedException(error, catch_backtrace()))
        end
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
            queue.completed += 1
            notify(queue.condition; all=true)
            copy(job.waiters)
        end
        for waiter in waiters
            put!(waiter, result)
        end
    end
end

"""Return processed items, sharing and promoting any work not already committed."""
function request_processed_items(
    workspace::Workspace,
    records::Vector{ItemRecord},
)::Vector{AbstractDataItem}
    cached = read_cached_item_data(workspace.cache.db, records; stage=:processed)
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
