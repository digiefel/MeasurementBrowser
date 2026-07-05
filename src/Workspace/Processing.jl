"""Result delivered to one caller waiting for processed item data."""
struct ProcessingResult
    item::Union{Nothing,AbstractDataItem}
    failure::Union{Nothing,CapturedException}
end

"""
Bump one work key's monotonic revision counter and return the new value.

Freshness propagates through the work graph, not through input hashes: an event bumps the key's
counter, which discards stale in-flight completions and dedups enqueues. A key with no node yet
starts at revision 1.
"""
function bump_revision!(graph::WorkDependencyGraph, key::WorkKey)::UInt64
    node = get(graph.nodes, key, nothing)
    return node === nothing ? UInt64(1) : node.revision + UInt64(1)
end

"""Return one work key's current revision, or 1 when no node exists yet."""
function current_revision(graph::WorkDependencyGraph, key::WorkKey)::UInt64
    node = get(graph.nodes, key, nothing)
    return node === nothing ? UInt64(1) : node.revision
end

"""Start one work-conserving worker pool shared by every work kind."""
function start_work_workers!(workspace::Workspace)::Nothing
    for _ in 1:max(2, Base.Threads.nthreads())
        push!(workspace.work.workers, Base.Threads.@spawn work_worker!(workspace))
    end
    return nothing
end

"""Stop accepting work, wait for all running callbacks, and fail any still-attached waiters."""
function stop_work_workers!(workspace::Workspace)::Nothing
    graph = workspace.work
    lock(graph.lock) do
        graph.closed = true
        notify(graph.condition; all=true)
    end
    foreach(wait, graph.workers)
    result = ProcessingResult(nothing, CapturedException(JobCancelled(), Any[]))
    lock(graph.lock) do
        for node in values(graph.nodes), waiter in node.waiters
            put!(waiter, result)
        end
        foreach(node -> empty!(node.waiters), values(graph.nodes))
    end
    return nothing
end

"""Cancel queued work while allowing callbacks already running to finish."""
function cancel_waiting_work!(workspace::Workspace)::Nothing
    graph = workspace.work
    canceled = WorkNode[]
    lock(graph.lock) do
        for node in values(graph.nodes)
            node.state in (:queued, :waiting) || continue
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
    lock(workspace.publish_lock) do
        notify(workspace.idle_condition; all=true)
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
        empty!(graph.dependents)
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

work_state_terminal(state::Symbol)::Bool = state in (:ready, :failed)

function dependency_ready_for(kind::WorkKind, dependency::WorkNode)::Bool
    kind in (COLLECTION_PROCESS, COLLECTION_ANALYZE) &&
        return work_state_terminal(dependency.state)
    return dependency.state === :ready
end

function dependencies_ready(graph::WorkDependencyGraph, node::WorkNode)::Bool
    return all(node.dependencies) do dependency_key
        dependency = get(graph.nodes, dependency_key, nothing)
        dependency !== nothing && dependency_ready_for(node.key.kind, dependency)
    end
end

function remove_dependency_edges!(graph::WorkDependencyGraph, node::WorkNode)::Nothing
    for dependency in node.dependencies
        dependents = get(graph.dependents, dependency, nothing)
        dependents === nothing && continue
        delete!(dependents, node.key)
        isempty(dependents) && delete!(graph.dependents, dependency)
    end
    return nothing
end

function set_dependencies!(
    graph::WorkDependencyGraph,
    node::WorkNode,
    dependencies::Vector{WorkKey},
)::Nothing
    remove_dependency_edges!(graph, node)
    node.dependencies = dependencies
    for dependency in dependencies
        push!(get!(() -> Set{WorkKey}(), graph.dependents, dependency), node.key)
    end
    return nothing
end

function replace_work_node!(graph::WorkDependencyGraph, node::WorkNode)::Nothing
    previous = get(graph.nodes, node.key, nothing)
    previous === nothing || remove_dependency_edges!(graph, previous)
    graph.nodes[node.key] = node
    set_dependencies!(graph, node, node.dependencies)
    return nothing
end

"""Append one queue entry to its priority bucket."""
function push_queue_entry!(
    graph::WorkDependencyGraph,
    priority::Int,
    entry::Tuple{WorkKey,UInt64},
)::Nothing
    push!(get!(() -> Tuple{WorkKey,UInt64}[], graph.queue, priority), entry)
    return nothing
end

function queue_ready_node!(graph::WorkDependencyGraph, node::WorkNode)::Nothing
    node.state = :queued
    node.queued_ns = time_ns()
    push_queue_entry!(graph, node.priority, (node.key, node.revision))
    # One queued node needs one worker. Waking the whole pool here made every bulk enqueue a
    # thundering herd serialized behind the caller's locks.
    notify(graph.condition; all=false)
    return nothing
end

function wake_ready_dependents!(graph::WorkDependencyGraph, dependency_key::WorkKey)::Nothing
    for dependent_key in collect(get(graph.dependents, dependency_key, Set{WorkKey}()))
        dependent = get(graph.nodes, dependent_key, nothing)
        dependent === nothing && continue
        dependent.state === :waiting || continue
        dependencies_ready(graph, dependent) || continue
        queue_ready_node!(graph, dependent)
    end
    return nothing
end

"""Bump one work key's revision and record it missing without queuing it (event-driven staleness)."""
function bump_work_missing!(
    workspace::Workspace,
    key::WorkKey;
    dependencies::Vector{WorkKey}=WorkKey[],
)::Nothing
    graph = workspace.work
    lock(graph.lock) do
        revision = bump_revision!(graph, key)
        replace_work_node!(graph, WorkNode(
            key,
            revision,
            :missing,
            0,
            dependencies,
            Channel{Any}[],
            0,
        ))
    end
    return nothing
end

"""Seed a completed work revision loaded from a fresh cache index."""
function seed_work_node!(
    workspace::Workspace,
    key::WorkKey,
    revision::UInt64,
    state::Symbol,
)::Nothing
    state in (:ready, :failed) ||
        throw(ArgumentError("cached work state must be :ready or :failed, got $state"))
    graph = workspace.work
    lock(graph.lock) do
        node = get(graph.nodes, key, nothing)
        if node === nothing || node.revision != revision || node.state ∉ (:queued, :running)
            replace_work_node!(graph, WorkNode(
                key,
                revision,
                state,
                0,
                WorkKey[],
                Channel{Any}[],
                0,
            ))
        end
        wake_ready_dependents!(graph, key)
    end
    return nothing
end

"""
Queue one work revision once, promoting it when its requested priority increases.

A closed graph drops the request: in-flight completions publishing during shutdown may legitimately
try to queue follow-up work, which shutdown discards. An attached waiter fails immediately.
"""
function enqueue_work!(
    workspace::Workspace,
    key::WorkKey,
    revision::UInt64;
    priority::Int,
    dependencies::Vector{WorkKey}=WorkKey[],
    waiter::Union{Nothing,Channel{Any}}=nothing,
)::Union{Nothing,WorkNode}
    graph = workspace.work
    return lock(graph.lock) do
        if graph.closed
            waiter === nothing ||
                put!(waiter, ProcessingResult(nothing, CapturedException(JobCancelled(), Any[])))
            return nothing
        end
        node = get(graph.nodes, key, nothing)
        if node !== nothing && node.revision == revision && node.state in (:ready, :failed)
            # Terminal at this revision: no future publish fires waiters attached now. Answer the
            # waiter immediately — an empty success means "already published, read the cache".
            if waiter !== nothing
                failure = node.state === :failed ?
                    CapturedException(ErrorException(get(
                        workspace.index.analysis_errors, key.entity, "work failed")), Any[]) :
                    nothing
                put!(waiter, ProcessingResult(nothing, failure))
            end
            return node
        end
        if node === nothing || node.revision != revision
            state = :waiting
            node = WorkNode(
                key,
                revision,
                state,
                priority,
                dependencies,
                Channel{Any}[],
                time_ns(),
            )
            replace_work_node!(graph, node)
            dependencies_ready(graph, node) && queue_ready_node!(graph, node)
            graph.total += 1
        elseif node.state === :missing
            node.priority = priority
            set_dependencies!(graph, node, dependencies)
            node.queued_ns = time_ns()
            if dependencies_ready(graph, node)
                queue_ready_node!(graph, node)
            else
                node.state = :waiting
            end
            graph.total += 1
        elseif node.state === :waiting
            node.priority = max(node.priority, priority)
            set_dependencies!(graph, node, dependencies)
            if dependencies_ready(graph, node)
                queue_ready_node!(graph, node)
            end
        elseif node.state === :queued && priority > node.priority
            # Promote by re-queuing at the higher priority; the stale lower-bucket entry is
            # skipped on take because the node is no longer :queued by the time it surfaces.
            node.priority = priority
            push_queue_entry!(graph, priority, (key, revision))
        end
        waiter === nothing || push!(node.waiters, waiter)
        node
    end
end

"""Pop the highest-priority queued current revision; caller must hold `graph.lock`."""
function pop_queued_node!(graph::WorkDependencyGraph)::Union{Nothing,WorkNode}
    for priority in sort!(collect(keys(graph.queue)); rev=true)
        bucket = graph.queue[priority]
        while !isempty(bucket)
            key, revision = popfirst!(bucket)
            node = get(graph.nodes, key, nothing)
            node === nothing && continue
            node.revision == revision && node.state === :queued || continue
            node.state = :running
            return node
        end
    end
    return nothing
end

"""Take the highest-priority current work revision, blocking until one is ready or the graph closes."""
function take_work!(graph::WorkDependencyGraph)::Union{Nothing,WorkNode}
    lock(graph.lock)
    try
        while true
            graph.closed && return nothing
            node = pop_queued_node!(graph)
            node === nothing || return node
            wait(graph.condition)
        end
    finally
        unlock(graph.lock)
    end
end

"""Whether this worker still owns the current running node."""
function work_node_current(workspace::Workspace, node::WorkNode)::Bool
    return lock(workspace.work.lock) do
        get(workspace.work.nodes, node.key, nothing) === node && node.state === :running
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
            source_item,
            interpretation.records,
            interpretation.interpreted_items,
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
    materialized_record = effective_record(workspace.index.hierarchy, record)
    input = DataItem(materialized_record, item_data(interpreted))
    processed = Profiling.@profile_span workspace.profiler :project :process Profiling.ProfileAttributes(
        kind=record.kind,
        source_id=record.source_item_id,
        item_id=record.id,
    ) begin
        process(workspace.project, workspace.source, input)
    end
    record_scan_phase!(workspace.project, record.source_item_id, record.kind,
        :process, (time_ns() - process_started) / 1e9, Base.Threads.threadid())
    return (
        item=DataItem(input, item_data(processed)),
        record=materialized_record,
        cacheable=cacheable(processed),
    )
end

"""Return the delivered metadata for one record: inherited ⊕ entries ⊕ computed layers."""
function delivered_metadata(workspace::Workspace, record::ItemRecord)::MetadataDict
    effective = effective_metadata(workspace.index.hierarchy, record)
    computed = get(workspace.index.item_metadata, record.id, nothing)
    computed === nothing || merge!(effective, metadata_dict(computed))
    return effective
end

"""Return records carrying their delivered metadata, for materializing a collection's members."""
function delivered_records(workspace::Workspace, records::Vector{ItemRecord})::Vector{ItemRecord}
    return ItemRecord[
        ItemRecord(record; metadata=delivered_metadata(workspace, record)) for record in records
    ]
end

"""Run item analysis from the already-published processed result; merge output over the item's layer."""
function run_item_analysis(workspace::Workspace, record::ItemRecord)::MetadataDict
    delivered_record = ItemRecord(record; metadata=delivered_metadata(workspace, record))
    processed = only(read_item_data(
        workspace.cache.db, [delivered_record]; stage=:processed))
    processed === nothing && error(
        "Cannot analyze item '$(record.id)': processed data is missing",
    )
    analyze_started = time_ns()
    computed = Profiling.@profile_span workspace.profiler :project :analyze Profiling.ProfileAttributes(
        kind=record.kind,
        source_id=record.source_item_id,
        item_id=record.id,
    ) begin
        metadata_dict(analyze_item(
            workspace.project,
            workspace.source,
            DataItem(
                delivered_record,
                item_data(processed),
            ),
        ))
    end
    record_scan_phase!(workspace.project, record.source_item_id, record.kind,
        :analyze, (time_ns() - analyze_started) / 1e9, Base.Threads.threadid())
    return computed
end

"""
Materialize and rewrite one collection's members through registered collection `process`.

Members are materialized from their processed payloads, folded, and each rewritten member (whose
`data` is not `===` the input) is re-cached at the `:collection_processed` stage. Returns the
rewritten members and the ids whose payload was rewritten, persisted worker-side before publish.
"""
function run_collection_process(workspace::Workspace, key::String)::NamedTuple
    path = collection_path_tuple(key)
    records = lock(workspace.work.lock) do
        node = get(workspace.index.hierarchy.index, path, nothing)
        node === nothing && error("Cannot process missing collection '$key'")
        copy(node.items)
    end
    delivered = delivered_records(workspace, records)
    payloads = read_item_data(workspace.cache.db, delivered; stage=:processed)
    inputs = AbstractDataItem[
        DataItem(delivered[index],
            payloads[index] === nothing ? nothing : item_data(payloads[index]))
        for index in eachindex(delivered)
    ]
    outputs = process_collection(
        workspace.project, workspace.source, collect(path), inputs)
    by_id = Dict(id(input) => input for input in inputs)
    metadata_by_id = Dict{String,MetadataDict}()
    rewritten_ids = String[]
    for output in outputs
        input = by_id[id(output)]
        metadata_by_id[id(output)] = metadata_dict(Projects.metadata(output))
        item_data(output) === item_data(input) && continue
        push!(rewritten_ids, id(output))
    end
    return (
        records=records,
        outputs=outputs,
        metadata_by_id=metadata_by_id,
        rewritten_ids=Set(rewritten_ids),
    )
end

"""Run one collection's analyze folds over the post-process members."""
function run_collection_analysis(workspace::Workspace, key::String)::MetadataDict
    path = collection_path_tuple(key)
    records = lock(workspace.work.lock) do
        node = get(workspace.index.hierarchy.index, path, nothing)
        node === nothing && error("Cannot summarize missing collection '$key'")
        copy(node.items)
    end
    delivered = delivered_records(workspace, records)
    # A member's kind decides its payload stage: a collection process rewrites only the payloads it
    # changed, so rewritten members read `:collection_processed` and everything else `:processed`.
    payloads = Vector{Any}(nothing, length(delivered))
    rewritten = [
        index for index in eachindex(records)
        if has_collection_process(workspace.project, records[index].kind)
    ]
    folded = read_item_data(workspace.cache.db, delivered[rewritten]; stage=:collection_processed)
    for (position, index) in pairs(rewritten)
        payloads[index] = folded[position]
    end
    remaining = [index for index in eachindex(records) if payloads[index] === nothing]
    base = read_item_data(workspace.cache.db, delivered[remaining]; stage=:processed)
    for (position, index) in pairs(remaining)
        payloads[index] = base[position]
    end
    items = AbstractDataItem[
        DataItem(delivered[index],
            payloads[index] === nothing ? nothing : item_data(payloads[index]))
        for index in eachindex(delivered)
    ]
    return metadata_dict(analyze_collection(
        workspace.project,
        workspace.source,
        collect(path),
        items,
    ))
end

"""Execute one work node and publish its completion immediately."""
function execute_work!(workspace::Workspace, node::WorkNode)::Nothing
    key = node.key
    result = try
        if key.kind === SOURCE_INTERPRET
            source_item = lock(workspace.work.lock) do
                get(workspace.work.source_items, key.entity, nothing)
            end
            source_item === nothing &&
                error("Cannot interpret removed source item '$(key.entity)'")
            interpretation = interpret_source_item(
                workspace.project, workspace.source, source_item, workspace.profiler)
            # The worker persists what it computed before publication: cache writes are part of
            # the work, so publication stays pure in-memory state under the publish lock.
            conflicts = String[]
            if work_node_current(workspace, node)
                conflicts = store_interpreted!(
                    workspace.cache.db,
                    source_item,
                    interpretation.records,
                    interpretation.interpreted_items,
                )
            end
            (records=interpretation.records, failures=interpretation.failures, conflicts)
        elseif key.kind === ITEM_PROCESS
            record = lock(workspace.work.lock) do
                get(workspace.index.items, key.entity, nothing)
            end
            record === nothing && error("Cannot process removed item '$(key.entity)'")
            processing = run_processing(workspace, record)
            # The payload store can block on write backpressure for whole flush cycles; running it
            # here throttles processing production without stalling any publication. It must land
            # before the node publishes as :ready, because late joiners then read the cache.
            if work_node_current(workspace, node)
                store_processed!(workspace.cache.db, processing.record, processing.item)
            end
            processing
        elseif key.kind === ITEM_ANALYZE
            record = lock(workspace.work.lock) do
                get(workspace.index.items, key.entity, nothing)
            end
            record === nothing && error("Cannot analyze removed item '$(key.entity)'")
            run_item_analysis(workspace, record)
        elseif key.kind === COLLECTION_PROCESS
            processing = run_collection_process(workspace, key.entity)
            # Persist rewritten payloads and each member's metadata worker-side before publication.
            if work_node_current(workspace, node)
                by_id = Dict(record.id => record for record in processing.records)
                for output in processing.outputs
                    id(output) in processing.rewritten_ids || continue
                    store_processed!(
                        workspace.cache.db, by_id[id(output)], output;
                        stage=:collection_processed)
                end
                store_collection_process_result!(workspace.cache.db, key.entity)
            end
            processing
        else
            run_collection_analysis(workspace, key.entity)
        end
    catch error
        CapturedException(error, catch_backtrace())
    end
    publish_work_completion!(workspace, key, node.revision, result)
    return nothing
end

"""Run pooled work of any kind until workspace shutdown."""
function work_worker!(workspace::Workspace)::Nothing
    graph = workspace.work
    while true
        node = take_work!(graph)
        node === nothing && return nothing
        execute_work!(workspace, node)
    end
end

"""Whether any current node of one kind is queued or running."""
function work_kind_running(workspace::Workspace, kind::WorkKind)::Bool
    return lock(workspace.work.lock) do
        any(
            node -> node.key.kind === kind && node.state in (:queued, :running),
            values(workspace.work.nodes),
        )
    end
end

"""Return whether one work key's current revision has already published successfully."""
function work_ready(workspace::Workspace, key::WorkKey)::Bool
    return lock(workspace.work.lock) do
        node = get(workspace.work.nodes, key, nothing)
        node !== nothing && node.state === :ready
    end
end

"""Return completed, total, and active work counters for status/profiling."""
function work_counts(workspace::Workspace)::Tuple{Int,Int,Int}
    return lock(workspace.work.lock) do
        active = count(
            node -> node.state in (:queued, :running),
            values(workspace.work.nodes),
        )
        (workspace.work.completed, workspace.work.total, active)
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
        WorkKey(ITEM_PROCESS, record.id),
        current_revision(workspace.work, WorkKey(ITEM_PROCESS, record.id));
        priority=selected ? 4 : 2,
        waiter,
    )
    return waiter
end

"""
The delivery gate for one record: the work node whose readiness makes the delivered payload current.

Kinds with a registered collection `process` are gated on their collection's COLLECTION_PROCESS and
read the `:collection_processed` payload (falling back to `:processed`); others gate on ITEM_PROCESS.
"""
function delivery_gate(workspace::Workspace, record::ItemRecord)::Tuple{WorkKey,Symbol}
    if has_collection_process(workspace.project, record.kind)
        return WorkKey(COLLECTION_PROCESS, collection_path_key(record.collection)), :collection
    end
    return WorkKey(ITEM_PROCESS, record.id), :item
end

"""Read one record's delivered payload from the cache, honoring its gate's payload stage."""
function _delivered_payload(workspace::Workspace, record::ItemRecord, mode::Symbol)::Any
    delivered = ItemRecord(record; metadata=delivered_metadata(workspace, record))
    if mode === :collection
        payload = only(read_item_data(
            workspace.cache.db, [delivered]; stage=:collection_processed))
        payload === nothing || return payload
    end
    return only(read_item_data(workspace.cache.db, [delivered]; stage=:processed))
end

"""Promote a collection-process gate and its member dependencies to selected priority, with a waiter."""
function enqueue_collection_delivery!(
    workspace::Workspace,
    gate::WorkKey,
)::Union{Nothing,Channel{Any}}
    path = collection_path_tuple(gate.entity)
    node = get(workspace.index.hierarchy.index, path, nothing)
    node === nothing && return nothing
    lock(workspace.work.lock) do
        for record in node.items
            enqueue_work!(
                workspace, WorkKey(ITEM_PROCESS, record.id),
                current_revision(workspace.work, WorkKey(ITEM_PROCESS, record.id));
                priority=4)
        end
    end
    waiter = Channel{Any}(1)
    enqueue_work!(
        workspace, gate, current_revision(workspace.work, gate);
        priority=3, waiter,
        dependencies=WorkKey[WorkKey(ITEM_ANALYZE, record.id) for record in node.items],
    )
    return waiter
end

"""Return processed items, joining and promoting shared work when required."""
function request_processed_items(
    workspace::Workspace,
    records::Vector{ItemRecord},
)::Vector{AbstractDataItem}
    loaded_item(record, item) = DataItem(
        ItemRecord(record; metadata=delivered_metadata(workspace, record)),
        item_data(item),
    )
    # The work graph decides whether a stored result is current: only keys it has published as
    # :ready read the cache; everything else joins or starts the work.
    loaded = Vector{AbstractDataItem}(undef, length(records))
    waiters = Tuple{Int,Symbol,Channel{Any}}[]
    for (position, record) in pairs(records)
        gate, mode = delivery_gate(workspace, record)
        if work_ready(workspace, gate)
            payload = _delivered_payload(workspace, record, mode)
            payload === nothing && error(
                "Delivered data for item '$(record.id)' is missing from the cache")
            loaded[position] = loaded_item(record, payload)
        elseif mode === :collection
            waiter = enqueue_collection_delivery!(workspace, gate)
            push!(waiters, (position, :collection, waiter::Channel{Any}))
        else
            waiter = enqueue_processing!(workspace, record; selected=true)
            push!(waiters, (position, :item, waiter::Channel{Any}))
        end
    end
    for (position, mode, waiter) in waiters
        record = records[position]
        result = take!(waiter)::ProcessingResult
        # A failed collection process is terminal, not fatal: members still deliver their
        # `:processed` payloads with the fold error surfaced. Item-stage failures do surface.
        if result.failure !== nothing
            (mode === :collection && !is_job_cancelled(result.failure.ex)) ||
                throw(result.failure)
        elseif mode === :item && result.item !== nothing
            loaded[position] = result.item::AbstractDataItem
            continue
        end
        payload = _delivered_payload(workspace, record, mode)
        payload === nothing && error(
            "Delivered data for item '$(record.id)' is missing from the cache")
        loaded[position] = loaded_item(record, payload)
    end
    return loaded
end
