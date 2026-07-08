"""Result delivered to one caller waiting for processed item data."""
struct ProcessingResult
    item::Union{Nothing,AbstractDataItem}
    failure::Union{Nothing,CapturedException}
end

"""Map one work key to its cache `result_states` kind."""
function _work_key_cache_kind(key::WorkKey)::CacheResultKind
    return key.kind === ITEM_PROCESS ? PROCESSING_RESULT :
        key.kind === ITEM_ANALYZE ? ITEM_ANALYSIS_RESULT :
        key.kind === COLLECTION_PROCESS ? COLLECTION_PROCESS_RESULT : COLLECTION_ANALYSIS_RESULT
end

"""
Return `:ready`, `:failed`, or `:absent` for one finished work key.

Finished work lives in the cache ledger; live jobs are detected via `haskey(work.nodes, key)`.
"""
function cache_work_status(workspace::Workspace, key::WorkKey)::Symbol
    cachedb = workspace.cache.db
    entity = String(key.entity)
    memory = !(cachedb isa CacheDB)
    if key.kind === SOURCE_INTERPRET
        if memory
            return lock(cachedb.lock) do
                haskey(cachedb.failures, (entity, entity)) && return :failed
                entity in cachedb.source_items ? :ready : :absent
            end
        end
        failures = read(cachedb.failures)
        haskey(failures, (entity, entity)) && return :failed
        return haskey(read(cachedb.source_items), entity) ? :ready : :absent
    end
    kind = _work_key_cache_kind(key)
    states = memory ? lock(() -> copy(cachedb.result_states), cachedb.lock) :
        read(cachedb.result_states)
    state = get(states, (Int8(kind), entity), nothing)
    state === nothing && return :absent
    return CacheResultStatus(state.status) === RESULT_READY ? :ready : :failed
end

"""Drop one work key's finished-state ledger row so invalidation must rerun it."""
function clear_work_result_state!(workspace::Workspace, key::WorkKey)::Nothing
    cachedb = workspace.cache.db
    entity = String(key.entity)
    if key.kind === SOURCE_INTERPRET
        clear_cached_source_state!(cachedb, entity)
        return nothing
    end
    clear_cached_result_state!(cachedb, _work_key_cache_kind(key), entity)
    return nothing
end

"""Bump one work key's monotonic revision and return the new value."""
function bump_revision!(graph::WorkDependencyGraph, key::WorkKey)::UInt16
    node = get(graph.nodes, key, nothing)
    return node === nothing ? UInt16(1) : node.revision + UInt16(1)
end

"""Return one work key's current revision, or 1 when no live node exists."""
function current_revision(graph::WorkDependencyGraph, key::WorkKey)::UInt16
    node = get(graph.nodes, key, nothing)
    return node === nothing ? UInt16(1) : node.revision
end

"""Register reverse edges and count live unmet dependencies for one new node."""
function seed_node_dependencies!(
    graph::WorkDependencyGraph,
    node::WorkNode,
    dependencies::Vector{WorkKey},
)::Nothing
    node.pending = 0
    for dependency in dependencies
        dependency_node = get(graph.nodes, dependency, nothing)
        dependency_node === nothing && continue
        push!(dependency_node.dependents, node.key)
        node.pending += 1
    end
    return nothing
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
        for node in collect(values(graph.nodes))
            node.state in (:queued, :waiting) || continue
            node.state === :queued && (graph.active -= 1)
            delete!(graph.nodes, node.key)
            push!(canceled, node)
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
        empty!(graph.queue)
        empty!(graph.source_items)
        empty!(graph.source_locks)
        graph.total = 0
        graph.completed = 0
        graph.active = 0
        graph.source_batch = 0
        graph.source_batch_open = false
    end
    return nothing
end

dependencies_ready(node::WorkNode)::Bool = node.pending == UInt64(0)

"""Append one queue entry to its priority bucket."""
function push_queue_entry!(
    graph::WorkDependencyGraph,
    priority::Int,
    entry::Tuple{WorkKey,UInt16},
)::Nothing
    push!(get!(() -> Tuple{WorkKey,UInt16}[], graph.queue, priority), entry)
    return nothing
end

function queue_ready_node!(graph::WorkDependencyGraph, node::WorkNode)::Nothing
    node.state === :queued && return nothing
    node.state = :queued
    graph.active += 1
    node.queued_ns = time_ns()
    push_queue_entry!(graph, node.priority, (node.key, node.revision))
    notify(graph.condition; all=false)
    return nothing
end

"""Decrement dependents' pending counters and queue any that become ready."""
function wake_ready_dependents!(graph::WorkDependencyGraph, node::WorkNode)::Nothing
    for dependent_key in node.dependents
        dependent = get(graph.nodes, dependent_key, nothing)
        dependent === nothing && continue
        dependent.pending -= 1
        dependent.pending == 0 &&
            dependent.state === :waiting &&
            queue_ready_node!(graph, dependent)
    end
    return nothing
end

"""Remove one finished live node and wake its dependents."""
function finish_work_node!(workspace::Workspace, node::WorkNode)::Vector{Channel{Any}}
    return lock(workspace.work.lock) do
        current = get(workspace.work.nodes, node.key, nothing)
        current === node || return Channel{Any}[]
        node.state in (:queued, :running) && (workspace.work.active -= 1)
        wake_ready_dependents!(workspace.work, node)
        delete!(workspace.work.nodes, node.key)
        workspace.work.completed += 1
        waiters = copy(node.waiters)
        empty!(node.waiters)
        waiters
    end
end

"""
Queue one work revision once, promoting it when its requested priority increases.

Re-enqueue replaces the node with a fresh revision and re-seeds `pending`/`dependents` from
`dependencies`; it never patches an existing node in place.
"""
function enqueue_work!(
    workspace::Workspace,
    key::WorkKey,
    revision::UInt16;
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
        if node !== nothing && node.revision == revision && node.state === :queued &&
                priority > node.priority
            node.priority = priority
            push_queue_entry!(graph, priority, (key, revision))
            waiter === nothing || push!(node.waiters, waiter)
            return node
        end
        if node === nothing || node.revision != revision
            if node !== nothing
                # Superseding a still-live revision reuses its work slot: the old revision never
                # completes, so counting the replacement as new work would leave total > completed
                # forever and make the progress denominator drift up as invalidations pile in.
                node.state in (:queued, :running) && (graph.active -= 1)
                previous_waiters = node.waiters
            else
                previous_waiters = Channel{Any}[]
                graph.total += 1
            end
            node = WorkNode(
                key,
                revision,
                :waiting,
                priority,
                Set{WorkKey}(),
                UInt64(0),
                previous_waiters,
                time_ns(),
            )
            graph.nodes[key] = node
            seed_node_dependencies!(graph, node, dependencies)
            dependencies_ready(node) && queue_ready_node!(graph, node)
        elseif node.state === :waiting
            node.priority = max(node.priority, priority)
            seed_node_dependencies!(graph, node, dependencies)
            dependencies_ready(node) && queue_ready_node!(graph, node)
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
        metadata_by_id[id(output)] = metadata_dict(DataBrowserAPI.metadata(output))
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
    publish_work_completion!(workspace, node, result)
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

"""Return completed, total, and active work counters for status/profiling."""
function work_counts(workspace::Workspace)::Tuple{Int,Int,Int}
    return lock(workspace.work.lock) do
        (workspace.work.completed, workspace.work.total, workspace.work.active)
    end
end

"""Queue one item's background processing at its current revision."""
function enqueue_processing!(workspace::Workspace, record::ItemRecord)::Nothing
    key = WorkKey(ITEM_PROCESS, record.id)
    enqueue_work!(workspace, key, current_revision(workspace.work, key); priority=2)
    return nothing
end

"""
The delivery gate for one record: the work key whose readiness makes the delivered payload current.

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

"""
Block until one work key is current: wait on a live node, read a ready cache row, or enqueue work.
"""
function ensure_uptodate!(
    workspace::Workspace,
    key::WorkKey;
    dependencies::Vector{WorkKey}=WorkKey[],
    priority::Int=2,
)::ProcessingResult
    waiter = Channel{Any}(1)
    should_wait = lock(workspace.work.lock) do
        node = get(workspace.work.nodes, key, nothing)
        node === nothing && return false
        push!(node.waiters, waiter)
        return true
    end
    if should_wait
        result = take!(waiter)::ProcessingResult
        return result
    end
    status = cache_work_status(workspace, key)
    status === :ready && return ProcessingResult(nothing, nothing)
    status === :failed && return ProcessingResult(
        nothing,
        CapturedException(
            ErrorException(get(workspace.index.analysis_errors, key.entity, "work failed")),
            Any[],
        ),
    )
    enqueue_work!(
        workspace,
        key,
        current_revision(workspace.work, key);
        priority,
        dependencies,
        waiter,
    )
    return take!(waiter)::ProcessingResult
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
    loaded = Vector{AbstractDataItem}(undef, length(records))
    for (position, record) in pairs(records)
        gate, mode = delivery_gate(workspace, record)
        if mode === :collection
            path = collection_path_tuple(gate.entity)
            hierarchy_node = get(workspace.index.hierarchy.index, path, nothing)
            hierarchy_node === nothing && error(
                "Cannot deliver item '$(record.id)': collection '$(gate.entity)' is missing")
            lock(workspace.work.lock) do
                for member in hierarchy_node.items
                    process_key = WorkKey(ITEM_PROCESS, member.id)
                    get(workspace.work.nodes, process_key, nothing) === nothing && continue
                    enqueue_work!(
                        workspace,
                        process_key,
                        current_revision(workspace.work, process_key);
                        priority=4,
                    )
                end
            end
            result = ensure_uptodate!(
                workspace, gate;
                priority=3,
                dependencies=WorkKey[
                    WorkKey(ITEM_ANALYZE, member.id) for member in hierarchy_node.items],
            )
            if result.failure !== nothing
                is_job_cancelled(result.failure.ex) && throw(result.failure)
            end
        else
            result = ensure_uptodate!(workspace, gate; priority=4)
            if result.failure !== nothing
                throw(result.failure)
            elseif result.item !== nothing
                loaded[position] = result.item::AbstractDataItem
                continue
            end
        end
        payload = _delivered_payload(workspace, record, mode)
        payload === nothing && error(
            "Delivered data for item '$(record.id)' is missing from the cache")
        loaded[position] = loaded_item(record, payload)
    end
    return loaded
end
