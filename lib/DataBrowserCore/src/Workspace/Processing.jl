"""Result delivered to one caller waiting for processed item data."""
struct ProcessingResult
    item::Union{Nothing,AbstractDataItem}
    failure::Union{Nothing,CapturedException}
end

"""Answer one waiter with workspace-shutdown cancellation."""
function _workspace_cancelled_result(workspace::Workspace)::ProcessingResult
    return ProcessingResult(
        nothing,
        CapturedException(
            OperationCanceledException(get_token(workspace.cancel_source)),
            Any[],
        ),
    )
end

"""Answer one waiter dropped from the queue without shutting the workspace down."""
function _dropped_work_result()::ProcessingResult
    source = CancellationTokenSource()
    cancel(source)
    return ProcessingResult(
        nothing,
        CapturedException(OperationCanceledException(get_token(source)), Any[]),
    )
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
    memory = !(cachedb isa CacheDB)
    if key.kind === SOURCE_INTERPRET
        entity = key.entity::String
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
    if key.kind in (COLLECTION_PROCESS, COLLECTION_ANALYZE)
        collection_key = key.entity::Int64
        states = memory ? lock(() -> copy(cachedb.collection_result_states), cachedb.lock) :
            read(cachedb.collection_result_states)
        state = get(states, (Int8(kind), collection_key), nothing)
        state === nothing && return :absent
        return CacheResultStatus(state.status) === RESULT_READY ? :ready : :failed
    end
    states = memory ? lock(() -> copy(cachedb.result_states), cachedb.lock) :
        read(cachedb.result_states)
    entity = key.entity::String
    state = get(states, (Int8(kind), entity), nothing)
    state === nothing && return :absent
    return CacheResultStatus(state.status) === RESULT_READY ? :ready : :failed
end

"""Drop one work key's finished-state ledger row so invalidation must rerun it."""
function clear_work_result_state!(workspace::Workspace, key::WorkKey)::Nothing
    cachedb = workspace.cache.db
    if key.kind === SOURCE_INTERPRET
        entity = key.entity::String
        clear_cached_source_state!(cachedb, entity)
        return nothing
    end
    cache_entity = key.kind in (COLLECTION_PROCESS, COLLECTION_ANALYZE) ?
        key.entity::Int64 : key.entity::String
    clear_cached_result_state!(cachedb, _work_key_cache_kind(key), cache_entity)
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
    result = _workspace_cancelled_result(workspace)
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
    result = _dropped_work_result()
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

"""Start one work-conserving worker pool shared by every work kind."""
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
                put!(waiter, _workspace_cancelled_result(workspace))
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

"""Remove one finished live node and wake its dependents."""
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
            discovered = source_items(
                workspace.source;
                cancel_token=get_token(workspace.cancel_source),
            )
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
            workspace.project,
            workspace.source,
            source_item,
        )
        resolved_records = ItemRecord[
            get(workspace.index.items, candidate.id, candidate)
            for candidate in interpretation.records
        ]
        store_interpreted!(
            workspace.cache.db,
            source_item,
            resolved_records,
            interpretation.interpreted_items,
        )
        requested = findfirst(item -> item.id == record.id, interpretation.records)
        requested === nothing && error(
            "Source item '$(record.source_item_id)' no longer produces item '$(record.id)'",
        )
        return interpretation.interpreted_items[requested]
    end
end

"""Materialize a registered item with its user-facing registration string path."""
function registered_data_item(
    collections::CollectionIndex,
    record::ItemRecord,
    data,
)::RegisteredDataItem
    names = registration_names(collections, record.collection_key)
    names === nothing && error(
        "Cannot materialize registered item '$(record.id)' from a typed collection path",
    )
    return RegisteredDataItem(
        record,
        data,
        names,
    )
end

"""Run processing without computing or publishing statistics."""
function run_processing(
    workspace::Workspace,
    collections::CollectionIndex,
    record::ItemRecord,
)::NamedTuple
    interpreted = only(read_item_data(
        workspace.cache.db, [record]; stage=:interpreted))
    interpreted === nothing && (interpreted = source_fallback(workspace, record))
    process_started = time_ns()
    materialized_record = effective_record(collections, record)
    input = interpreted isa RegisteredDataItem ?
        registered_data_item(
            collections, materialized_record, item_data(interpreted)) : interpreted
    processed = @timed_dbg process(workspace.project, workspace.source, input)
    record_scan_phase!(workspace.project, record.source_item_id, record.kind,
        :process, (time_ns() - process_started) / 1e9, Base.Threads.threadid())
    return (
        item=processed,
        record=materialized_record,
        cacheable=cacheable(processed),
    )
end

"""Return the delivered metadata for one record: inherited ⊕ entries ⊕ computed layers."""
function delivered_metadata(
    workspace::Workspace,
    record::ItemRecord,
    collections::CollectionIndex,
)::MetadataDict
    effective = effective_metadata(collections, record)
    computed = get(workspace.index.item_metadata, record.id, nothing)
    computed === nothing || merge!(effective, metadata_dict(computed))
    return effective
end

"""Return records carrying their delivered metadata, for materializing a collection's members."""
function delivered_records(
    workspace::Workspace,
    collections::CollectionIndex,
    records::Vector{ItemRecord},
)::Vector{ItemRecord}
    return ItemRecord[
        ItemRecord(record; metadata=delivered_metadata(workspace, record, collections))
        for record in records
    ]
end

"""Run item analysis from the already-published processed result; merge output over the item's layer."""
function run_item_analysis(
    workspace::Workspace,
    collections::CollectionIndex,
    record::ItemRecord,
)::MetadataDict
    delivered_record = ItemRecord(
        record; metadata=delivered_metadata(workspace, record, collections))
    processed = only(read_item_data(
        workspace.cache.db, [delivered_record]; stage=:processed))
    processed === nothing && error(
        "Cannot analyze item '$(record.id)': processed data is missing",
    )
    analyze_started = time_ns()
    computed = @timed_dbg "analyze" begin
        input = processed isa RegisteredDataItem ?
            registered_data_item(
                collections, delivered_record, item_data(processed)) : processed
        metadata_dict(_analyze_item(
            workspace.project,
            workspace.source,
            input,
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
function run_collection_process(workspace::Workspace, collection_key::Int64)::NamedTuple
    index = workspace.index
    collections = index.collections
    records = lock(workspace.work.lock) do
        haskey(collections.records, collection_key) ||
            error("Cannot process missing collection '$collection_key'")
        ItemRecord[
            index.items[id]
            for id in collection_item_ids(collections, collection_key)
            if haskey(index.items, id)
        ]
    end
    delivered = delivered_records(workspace, collections, records)
    payloads = read_item_data(workspace.cache.db, delivered; stage=:processed)
    any(isnothing, payloads) && error(
        "Cannot process collection '$collection_key': one or more processed members are missing",
    )
    inputs = AbstractDataItem[
        payload isa RegisteredDataItem ?
            registered_data_item(
                collections, delivered[index], item_data(payload)) :
            payload::AbstractDataItem
        for (index, payload) in pairs(payloads)
    ]
    outputs = _process_collection(
        workspace.project,
        workspace.source,
        inputs,
    )
    by_id = Dict(id(input) => input for input in inputs)
    rewritten_ids = String[]
    for output in outputs
        input = by_id[id(output)]
        item_data(output) === item_data(input) && continue
        push!(rewritten_ids, id(output))
    end
    return (
        records=records,
        outputs=outputs,
        rewritten_ids=Set(rewritten_ids),
    )
end

"""Run one collection's analyze folds over the post-process members."""
function run_collection_analysis(workspace::Workspace, collection_key::Int64)::MetadataDict
    index = workspace.index
    collections = index.collections
    records = lock(workspace.work.lock) do
        haskey(collections.records, collection_key) ||
            error("Cannot summarize missing collection '$collection_key'")
        ItemRecord[
            index.items[id]
            for id in collection_item_ids(collections, collection_key)
            if haskey(index.items, id)
        ]
    end
    delivered = delivered_records(workspace, collections, records)
    payloads = Vector{Any}(nothing, length(delivered))
    rewritten = [
        index for index in eachindex(records)
        if _has_collection_process(workspace.project, records[index].kind)
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
    any(isnothing, payloads) && error(
        "Cannot analyze collection '$collection_key': one or more processed members are missing",
    )
    items = AbstractDataItem[
        payload isa RegisteredDataItem ?
            registered_data_item(
                collections, delivered[index], item_data(payload)) :
            payload::AbstractDataItem
        for (index, payload) in pairs(payloads)
    ]
    return metadata_dict(_analyze_collection(
        workspace.project,
        workspace.source,
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
                workspace.project, workspace.source, source_item)
            (
                source_item=source_item,
                interpretation=interpretation,
            )
        elseif key.kind === ITEM_PROCESS
            index = workspace.index
            record = get(index.items, key.entity, nothing)
            record === nothing && error("Cannot process removed item '$(key.entity)'")
            processing = run_processing(workspace, index.collections, record)
            if work_node_current(workspace, node)
                store_processed!(workspace.cache.db, processing.record, processing.item)
            end
            processing
        elseif key.kind === ITEM_ANALYZE
            index = workspace.index
            record = get(index.items, key.entity, nothing)
            record === nothing && error("Cannot analyze removed item '$(key.entity)'")
            run_item_analysis(workspace, index.collections, record)
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
                store_collection_process_result!(
                    workspace.cache.db,
                    key.entity::Int64,
                )
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
    if _has_collection_process(workspace.project, record.kind)
        return WorkKey(
            COLLECTION_PROCESS,
            record.collection_key::Int64,
        ), :collection
    end
    return WorkKey(ITEM_PROCESS, record.id), :item
end

"""Read one record's delivered payload from the cache, honoring its gate's payload stage."""
function _delivered_payload(
    workspace::Workspace,
    collections::CollectionIndex,
    record::ItemRecord,
    mode::Symbol,
)::Any
    delivered = ItemRecord(
        record; metadata=delivered_metadata(workspace, record, collections))
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
    index = workspace.index
    collections = index.collections
    loaded_item(record, item) = item isa RegisteredDataItem ?
        registered_data_item(
            collections,
            ItemRecord(record; metadata=delivered_metadata(workspace, record, collections)),
            item_data(item),
        ) : item
    loaded = Vector{AbstractDataItem}(undef, length(records))
    for (position, record) in pairs(records)
        gate, mode = delivery_gate(workspace, record)
        if mode === :collection
            collection_key = gate.entity::Int64
            haskey(collections.records, collection_key) || error(
                "Cannot deliver item '$(record.id)': collection '$(gate.entity)' is missing")
            members = ItemRecord[
                index.items[id]
                for id in collection_item_ids(collections, collection_key)
                if haskey(index.items, id)
            ]
            lock(workspace.work.lock) do
                for member in members
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
                    WorkKey(ITEM_ANALYZE, member.id) for member in members],
            )
            if result.failure !== nothing
                result.failure.ex isa OperationCanceledException && throw(result.failure)
            end
        else
            result = ensure_uptodate!(workspace, gate; priority=4)
            if result.failure !== nothing
                throw(result.failure)
            elseif result.item !== nothing
                loaded[position] = loaded_item(record, result.item::AbstractDataItem)
                continue
            end
        end
        payload = _delivered_payload(workspace, collections, record, mode)
        payload === nothing && error(
            "Delivered data for item '$(record.id)' is missing from the cache")
        loaded[position] = loaded_item(record, payload::AbstractDataItem)
    end
    return loaded
end
