"""Result delivered to one caller waiting for processed item data."""
struct ProcessingResult
    item::Union{Nothing,AbstractDataItem}
    failure::Union{Nothing,CapturedException}
end

"""Return the revision of work derived from one indexed item."""
function item_revision(record::ItemRecord)::UInt64
    return UInt64(hash((
        record.source_item_id,
        record.source_item_fingerprint,
        record.id,
        record.item_fingerprint,
    )))
end

"""Return the revision of collection statistics derived from current member statistics."""
function collection_revision(workspace::Workspace, path::Tuple{Vararg{String}})::UInt64
    node = workspace.index.hierarchy.index[path]
    members = sort!(Tuple{String,UInt64}[
        (record.id, UInt64(hash(get(workspace.index.item_stats, record.id, record.stats))))
        for record in node.items
    ])
    return UInt64(hash((path, members)))
end

"""Start the bounded worker set shared by every work kind."""
function start_work_workers!(workspace::Workspace)::Nothing
    for _ in 1:max(1, Base.Threads.nthreads() ÷ 2)
        push!(workspace.work.workers, Base.Threads.@spawn work_worker!(workspace))
    end
    return nothing
end

"""Stop accepting work and wait for all running callbacks."""
function stop_work_workers!(workspace::Workspace)::Nothing
    graph = workspace.work
    lock(graph.lock) do
        graph.closed = true
        notify(graph.condition; all=true)
    end
    foreach(wait, graph.workers)
    return nothing
end

"""Cancel queued work while allowing callbacks already running to finish."""
function cancel_waiting_work!(workspace::Workspace)::Nothing
    graph = workspace.work
    canceled = WorkNode[]
    lock(graph.lock) do
        for node in values(graph.nodes)
            node.state === :queued || continue
            node.state = :missing
            append!(canceled, [node])
        end
        empty!(graph.queue)
        notify(graph.condition; all=true)
    end
    result = ProcessingResult(nothing, CapturedException(JobCancelled(), Any[]))
    for node in canceled, waiter in node.waiters
        put!(waiter, result)
    end
    return nothing
end

"""Reset idle work counters before a clean profiling rebuild."""
function reset_work_graph!(workspace::Workspace)::Nothing
    graph = workspace.work
    lock(graph.lock) do
        any(node -> node.state in (:queued, :running), values(graph.nodes)) &&
            error("Cannot reset the work graph while work is active")
        empty!(graph.nodes)
        empty!(graph.queue)
        empty!(graph.source_items)
        empty!(graph.source_locks)
        graph.total = 0
        graph.completed = 0
        graph.source_batch = 0
        graph.source_batch_open = false
    end
    return nothing
end

"""Queue one work revision once, promoting it when its requested priority increases."""
function enqueue_work!(
    workspace::Workspace,
    key::WorkKey,
    revision::UInt64;
    priority::Int,
    dependencies::Vector{WorkKey}=WorkKey[],
    waiter::Union{Nothing,Channel{Any}}=nothing,
)::WorkNode
    graph = workspace.work
    return lock(graph.lock) do
        graph.closed && error("Cannot queue work after the workspace has closed")
        node = get(graph.nodes, key, nothing)
        if node === nothing || node.revision != revision
            node = WorkNode(
                key,
                revision,
                :queued,
                priority,
                dependencies,
                Channel{Any}[],
                time_ns(),
            )
            graph.nodes[key] = node
            push!(graph.queue, (key, revision))
            graph.total += 1
        elseif node.state === :missing
            node.state = :queued
            node.priority = priority
            node.dependencies = dependencies
            node.queued_ns = time_ns()
            push!(graph.queue, (key, revision))
            graph.total += 1
        elseif node.state === :queued && priority > node.priority
            node.priority = priority
        end
        waiter === nothing || push!(node.waiters, waiter)
        notify(graph.condition)
        node
    end
end

"""Take the highest-priority current work revision."""
function take_work!(graph::WorkDependencyGraph)::Union{Nothing,WorkNode}
    while true
        graph.closed && return nothing
        best_position = 0
        best_priority = typemin(Int)
        for (position, (key, revision)) in pairs(graph.queue)
            node = get(graph.nodes, key, nothing)
            node === nothing && continue
            node.revision == revision && node.state === :queued || continue
            if node.priority > best_priority
                best_position = position
                best_priority = node.priority
            end
        end
        if best_position > 0
            key, revision = splice!(graph.queue, best_position)
            node = get(graph.nodes, key, nothing)
            node === nothing && continue
            node.revision == revision && node.state === :queued || continue
            node.state = :running
            return node
        end
        empty!(graph.queue) || empty!(graph.queue)
        wait(graph.condition)
    end
end

"""Return the source fallback lock shared by items from one source item."""
function source_fallback_lock(
    graph::WorkDependencyGraph,
    source_item_id_value::String,
)::ReentrantLock
    return lock(graph.lock) do
        get!(() -> ReentrantLock(), graph.source_locks, source_item_id_value)
    end
end

"""Interpret one source item and return the requested logical item."""
function source_fallback(workspace::Workspace, record::ItemRecord)::AbstractDataItem
    fallback_lock = source_fallback_lock(workspace.work, record.source_item_id)
    return lock(fallback_lock) do
        cached = only(read_item_data(workspace.cache.db, [record]; stage=:interpreted))
        cached === nothing || return cached

        source_item = lock(workspace.work.lock) do
            get(workspace.work.source_items, record.source_item_id, nothing)
        end
        if source_item === nothing
            discovered = source_items(workspace.source)
            position = findfirst(
                item -> source_item_id(item) == record.source_item_id,
                discovered,
            )
            position === nothing && error(
                "Cannot load item '$(record.id)': source item '$(record.source_item_id)' " *
                "is no longer present in source '$(source_id(workspace.source))'",
            )
            source_item = discovered[position]
        end
        interpretation = interpret_source_item(
            workspace.project, workspace.source, source_item, workspace.profiler)
        store_interpreted!(
            workspace.cache.db,
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

"""Run processing without computing or publishing statistics."""
function run_processing(workspace::Workspace, record::ItemRecord)::NamedTuple
    interpreted = only(read_item_data(
        workspace.cache.db, [record]; stage=:interpreted))
    interpreted === nothing && (interpreted = source_fallback(workspace, record))
    process_started = time_ns()
    processed = Profiling.@profile_span workspace.profiler :project :process Profiling.ProfileAttributes(
        kind=record.kind,
        source_id=record.source_item_id,
        item_id=record.id,
    ) begin
        process(workspace.project, workspace.source, interpreted)
    end
    record_scan_phase!(workspace.project, record.source_item_id, record.kind,
        :process, (time_ns() - process_started) / 1e9, Base.Threads.threadid())
    return (item=DataItem(record, item_data(processed)), cacheable=cacheable(processed))
end

"""Run item statistics from the already-published processed result."""
function run_item_stats(workspace::Workspace, record::ItemRecord)::MetadataDict
    processed = only(read_item_data(workspace.cache.db, [record]; stage=:processed))
    processed === nothing && error(
        "Cannot compute statistics for item '$(record.id)': processed data is missing",
    )
    stats_started = time_ns()
    computed = Profiling.@profile_span workspace.profiler :project :item_stats Profiling.ProfileAttributes(
        kind=record.kind,
        source_id=record.source_item_id,
        item_id=record.id,
    ) begin
        metadata_dict(compute_item_stats(
            workspace.project,
            workspace.source,
            DataItem(record, item_data(processed)),
        ))
    end
    record_scan_phase!(workspace.project, record.source_item_id, record.kind,
        :stats, (time_ns() - stats_started) / 1e9, Base.Threads.threadid())
    return merge(copy(record.stats), computed)
end

"""Run one collection-stat callback against a stable snapshot of published item statistics."""
function run_collection_stats(workspace::Workspace, key::String)::MetadataDict
    path = collection_path_tuple(key)
    records, stats = lock(workspace.work.lock) do
        node = get(workspace.index.hierarchy.index, path, nothing)
        node === nothing && error("Cannot summarize missing collection '$key'")
        (copy(node.items), deepcopy(workspace.index.item_stats))
    end
    items = AbstractDataItem[
        DataItem(ItemRecord(record; stats=get(stats, record.id, record.stats)), nothing)
        for record in records
    ]
    return metadata_dict(collection_stats(
        workspace.project,
        workspace.source,
        collect(path),
        items,
    ))
end

"""Execute one work node and emit a typed completion event."""
function execute_work!(workspace::Workspace, node::WorkNode)::Nothing
    key = node.key
    result = try
        if key.kind === INTERPRET_SOURCE
            source_item = lock(workspace.work.lock) do
                get(workspace.work.source_items, key.entity, nothing)
            end
            source_item === nothing &&
                error("Cannot interpret removed source item '$(key.entity)'")
            interpret_source_item(
                workspace.project, workspace.source, source_item, workspace.profiler)
        elseif key.kind === PROCESS_ITEM
            record = lock(workspace.work.lock) do
                get(workspace.index.items, key.entity, nothing)
            end
            record === nothing && error("Cannot process removed item '$(key.entity)'")
            run_processing(workspace, record)
        elseif key.kind === ITEM_STATS
            record = lock(workspace.work.lock) do
                get(workspace.index.items, key.entity, nothing)
            end
            record === nothing && error("Cannot summarize removed item '$(key.entity)'")
            run_item_stats(workspace, record)
        else
            run_collection_stats(workspace, key.entity)
        end
    catch error
        CapturedException(error, catch_backtrace())
    end
    put!(workspace.work.events, (
        kind=:work_complete,
        key,
        revision=node.revision,
        result,
    ))
    return nothing
end

"""Run work until workspace shutdown."""
function work_worker!(workspace::Workspace)::Nothing
    graph = workspace.work
    while true
        node = lock(graph.lock) do
            take_work!(graph)
        end
        node === nothing && return nothing
        execute_work!(workspace, node)
    end
end

"""Queue one item-processing revision and optionally attach a selected waiter."""
function enqueue_processing!(
    workspace::Workspace,
    record::ItemRecord;
    selected::Bool=false,
)::Union{Nothing,Channel{Any}}
    waiter = selected ? Channel{Any}(1) : nothing
    enqueue_work!(
        workspace,
        WorkKey(PROCESS_ITEM, record.id),
        item_revision(record);
        priority=selected ? 4 : 2,
        dependencies=WorkKey[WorkKey(INTERPRET_SOURCE, record.source_item_id)],
        waiter,
    )
    return waiter
end

"""Return processed items, joining and promoting shared work when required."""
function request_processed_items(
    workspace::Workspace,
    records::Vector{ItemRecord},
)::Vector{AbstractDataItem}
    cached = read_item_data(workspace.cache.db, records; stage=:processed)
    loaded = Vector{AbstractDataItem}(undef, length(records))
    waiters = Tuple{Int,Channel{Any}}[]
    for (position, record) in pairs(records)
        item = cached[position]
        if item === nothing
            waiter = enqueue_processing!(workspace, record; selected=true)
            push!(waiters, (position, waiter::Channel{Any}))
        else
            loaded[position] = DataItem(
                ItemRecord(record; stats=get(
                    workspace.index.item_stats, record.id, record.stats)),
                item_data(item),
            )
        end
    end
    for (position, waiter) in waiters
        while !isready(waiter)
            poll_workspace!(workspace)
            yield()
        end
        result = take!(waiter)::ProcessingResult
        result.failure === nothing || throw(result.failure)
        loaded[position] = result.item::AbstractDataItem
    end
    return loaded
end
