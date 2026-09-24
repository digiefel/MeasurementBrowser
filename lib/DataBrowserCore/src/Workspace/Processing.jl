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

"""
Return `:ready`, `:failed`, or `:absent` for one finished work key.

Finished work lives in the cache ledger; live jobs are detected via `haskey(work.nodes, key)`.
"""
function cache_work_status(workspace::Workspace, key::WorkKey)::Symbol
    cachedb = workspace.cache.db
    kind = key.kind
    entity = key.entity
    state = cached_result_state(cachedb, kind, entity)
    if kind in (SOURCE_READ, ITEM_PROCESS)
        state !== nothing && CacheResultStatus(state.status) === RESULT_FAILED && return :failed
        state === nothing && return :absent
        return has_payload(cachedb, entity; stage=kind) ? :ready : :absent
    end
    state === nothing && return :absent
    return CacheResultStatus(state.status) === RESULT_READY ? :ready : :failed
end

"""Drop one work key's finished-state ledger row so invalidation must rerun it."""
function clear_work_result_state!(workspace::Workspace, key::WorkKey)::Nothing
    cachedb = workspace.cache.db
    clear_cached_result_state!(cachedb, key.kind, key.entity)
    return nothing
end

"""
Pull a source interpretation's read dependency, reusing a resident read result when available.

Interpretation has higher priority than reading so completed reads are consumed promptly instead
of accumulating inputs behind the remaining source reads. The dependent pins its input until it
finishes; the cache independently retains the value for later replay.
"""
function enqueue_source_interpretation!(workspace::Workspace, source_key::Int64;
        priority::Int=3, waiter::Union{Nothing,Channel{Any}}=nothing, supersede::Bool=false)
    return lock(workspace.work.lock) do
        revision(key) = supersede ? bump_revision!(workspace.work, key) :
            current_revision(workspace.work, key)
        read_key = WorkKey(SOURCE_READ, source_key)
        interpret_key = WorkKey(SOURCE_INTERPRET, source_key)
        cached = cache_work_status(workspace, read_key) === :ready ?
            read_payload(workspace.cache.db, source_key) : nothing
        dependencies = WorkKey[]
        if haskey(workspace.work.nodes, read_key) || cached === nothing
            enqueue_work!(workspace, read_key, revision(read_key); priority)
            push!(dependencies, read_key)
        end
        node = enqueue_work!(workspace, interpret_key,
            revision(interpret_key);
            priority=priority + 1, dependencies, waiter)
        node === nothing || isempty(dependencies) && (node.input = cached)
        return node
    end
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
                nothing,
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
    source_item_key_value::Int64,
)::ReentrantLock
    return lock(graph.lock) do
        get!(() -> ReentrantLock(), graph.source_locks, source_item_key_value)
    end
end

"""
Obtain the read input during source replay, with the caller holding the source fallback lock.

Replay runs inline because an item worker must not block waiting for another job in the same pool.
It uses the same cache boundary as scheduled source work and never reruns a saved read failure.
"""
function read_source_input!(workspace::Workspace, source_item::AbstractDataSourceItem,
        source_key::Int64)
    state = cached_result_state(workspace.cache.db, SOURCE_READ, source_key)
    if state !== nothing && CacheResultStatus(state.status) === RESULT_FAILED
        error(something(state.message, "Source read failed"))
    end
    cached = state === nothing ? nothing : read_payload(workspace.cache.db, source_key)
    cached === nothing || return something(cached)
    loaded = read(workspace.project, workspace.source, source_item)
    store_source_read!(workspace.cache.db, source_key, loaded)
    return loaded
end

"""Interpret one source item and return the requested logical item."""
function source_fallback(workspace::Workspace, record::ItemRecord)::AbstractDataItem
    fallback_lock = source_fallback_lock(workspace.work, record.source_item_key)
    return lock(fallback_lock) do
        source_ref = source_item_id(workspace, record.source_item_key)
        source_item = lock(workspace.work.lock) do
            get(workspace.work.source_items, record.source_item_key, nothing)
        end
        if source_item === nothing
            discovered = source_items(
                workspace.source;
                cancel_token=get_token(workspace.cancel_source),
            )
            position = findfirst(
                item -> id(item) == source_ref,
                discovered,
            )
            position === nothing && error(
                "Cannot load item '$(record.id)': source item '$source_ref' " *
                "is no longer present in source '$(source_id(workspace.source))'",
            )
            source_item = discovered[position]
        end
        loaded = read_source_input!(workspace, source_item, record.source_item_key)
        interpretation = interpret_source_item(
            workspace.project, workspace.source, source_item, loaded;
            source_item_key=record.source_item_key)
        store_interpreted_data!(
            workspace.cache.db, interpretation.records,
            item_data.(interpretation.interpreted_items))
        requested = findfirst(item -> item.id == record.id, interpretation.records)
        requested === nothing && error(
            "Source item '$source_ref' no longer produces item '$(record.id)'",
        )
        return interpretation.interpreted_items[requested]
    end
end

"""
Prepare one interpreted item for `process`.

A cached payload is reconstructed with interpretation-time metadata. When reconstruction declines
or the resident payload is missing, source fallback reruns `entries` with the cached read input.
It calls `read` only if that input is missing. The caller still runs `process` exactly once.
"""
function interpreted_item(
    workspace::Workspace,
    collections::CollectionIndex,
    record::ItemRecord,
    cached::Union{Nothing,Some},
)::AbstractDataItem
    path = indexed_collection_path(workspace, collections, record)
    effective = effective_record(collections, record)
    if cached !== nothing
        rebuilt = reconstruct(
            record.type, record.id, something(cached), effective.metadata)
        rebuilt === nothing || return attach_record(rebuilt, effective, path)
    end
    return attach_record(source_fallback(workspace, record), effective, path)
end

"""
The collection path the index holds for one record, ancestor to self.

Each level is rebuilt from its stored type, identity, and own metadata. Unlike items, a collection
has no rerun-from-source fallback, so `reconstruct` is required rather than optional.
"""
indexed_collection_path(
    ::Workspace,
    collections::CollectionIndex,
    record::ItemRecord,
)::Vector{AbstractCollection} = collection_value_path(collections, record.collection_key)

"""The concrete collection value for one indexed collection key."""
function collection_value(
    ::Workspace,
    collections::CollectionIndex,
    collection_key::Int64,
)::AbstractCollection
    path = collection_value_path(collections, collection_key)
    isempty(path) && error("Collection '$collection_key' has no indexed path")
    return last(path)
end

"""Rebuild an item through `entries` and `process`, obtaining the read input from cache or source."""
function reprocess_item(
    workspace::Workspace,
    collections::CollectionIndex,
    record::ItemRecord,
)::AbstractDataItem
    upstream_record = get(workspace.index.items, record.id, record)
    cached = only(read_payload(
        workspace.cache.db, [upstream_record]; stage=SOURCE_INTERPRET))
    input = interpreted_item(workspace, collections, upstream_record, cached)
    processed = process(workspace.project, input)
    processed isa AbstractDataItem || error(
        "process(::$(typeof(workspace.project)), ::$(typeof(input))) must return an " *
        "AbstractDataItem; got $(typeof(processed))",
    )
    return processed
end

"""
Turn one cached payload into an item ready for the next project stage.

Every payload is rebuilt through `reconstruct` on the record's own concrete type, including a
payload that is itself an `AbstractDataItem`. When that method declines, the engine reruns
`entries` and `process`, reading the source only if the read input is no longer cached. A rebuilt item then adopts its record and the index's collection path, so
identity comes from the engine and never from stale cached state.
"""
function materialized_item(
    workspace::Workspace,
    collections::CollectionIndex,
    record::ItemRecord,
    payload,
)::AbstractDataItem
    path = indexed_collection_path(workspace, collections, record)
    rebuilt = reconstruct(
        record.type, record.id, payload, effective_metadata(collections, record))
    rebuilt === nothing && return reprocess_item(workspace, collections, record)
    return attach_record(rebuilt, effective_record(collections, record), path)
end

"""Run processing without computing or publishing statistics."""
function run_processing(
    workspace::Workspace,
    collections::CollectionIndex,
    record::ItemRecord,
)::NamedTuple
    cached = only(read_payload(
        workspace.cache.db, [record]; stage=SOURCE_INTERPRET))
    materialized_record = effective_record(collections, record)
    input = interpreted_item(workspace, collections, record, cached)
    processed = @timed_dbg process(workspace.project, input)
    processed isa AbstractDataItem || error(
        "process(::$(typeof(workspace.project)), ::$(typeof(input))) must return an " *
        "AbstractDataItem; got $(typeof(processed))",
    )
    processed_metadata = merge(
        copy(materialized_record.metadata), metadata_dict(metadata(processed)))
    return (item=processed, record=ItemRecord(materialized_record; metadata=processed_metadata))
end

"""
Return an owned snapshot of inherited, entries and computed metadata under the publication lock.
Callers can reconstruct user values after releasing the lock while other items publish analysis.
"""
function delivered_metadata(
    workspace::Workspace,
    record::ItemRecord,
    collections::CollectionIndex,
)::MetadataDict
    return lock(workspace.publish_lock) do
        effective = effective_metadata(collections, record)
        computed = get(workspace.index.item_metadata, record.id, nothing)
        computed === nothing || merge!(effective, metadata_dict(computed))
        effective
    end
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
    processed = only(read_payload(
        workspace.cache.db, [delivered_record]; stage=ITEM_PROCESS))
    processed === nothing && error(
        "Cannot analyze item '$(record.id)': processed data is missing",
    )
    return @timed_dbg "analyze" begin
        metadata_dict(analyze(
            workspace.project,
            materialized_item(
                workspace, collections, delivered_record, something(processed))))
    end
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
    payloads = read_payload(workspace.cache.db, delivered; stage=ITEM_PROCESS)
    any(isnothing, payloads) && error(
        "Cannot process collection '$collection_key': one or more processed members are missing",
    )
    inputs = AbstractDataItem[
        materialized_item(workspace, collections, delivered[index], something(payload))
        for (index, payload) in pairs(payloads)
    ]
    outputs = process(
        workspace.project, collection_value(workspace, collections, collection_key), inputs)
    length(outputs) == length(inputs) || error(
        "Collection process for '$collection_key' returned $(length(outputs)) items for " *
        "$(length(inputs)) members; it must return one output per input, in order",
    )
    # Paired by position, not by id: a typed item carries no id of its own, and the contract
    # already requires one output per input.
    rewritten = [
        position for position in eachindex(inputs)
        if item_data(outputs[position]) !== item_data(inputs[position])
    ]
    return (
        records=records,
        outputs=outputs,
        rewritten=rewritten,
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
    # A member has a collection-processed payload only when a fold rewrote it; everyone else
    # analyzes from their own processed payload.
    payloads = read_payload(
        workspace.cache.db, delivered; stage=COLLECTION_PROCESS)
    remaining = [index for index in eachindex(records) if payloads[index] === nothing]
    base = read_payload(workspace.cache.db, delivered[remaining]; stage=ITEM_PROCESS)
    for (position, index) in pairs(remaining)
        payloads[index] = base[position]
    end
    any(isnothing, payloads) && error(
        "Cannot analyze collection '$collection_key': one or more processed members are missing",
    )
    items = AbstractDataItem[
        materialized_item(workspace, collections, delivered[index], something(payload))
        for (index, payload) in pairs(payloads)
    ]
    return metadata_dict(analyze(
        workspace.project, collection_value(workspace, collections, collection_key), items))
end

"""Execute one work node and publish its completion immediately."""
function execute_work!(workspace::Workspace, node::WorkNode)::Nothing
    key = node.key
    result = try
        if key.kind in (SOURCE_READ, SOURCE_INTERPRET)
            source_item = lock(workspace.work.lock) do
                get(workspace.work.source_items, key.entity, nothing)
            end
            source_item === nothing && error("Cannot load removed source item '$(key.entity)'")
            if key.kind === SOURCE_READ
                (loaded=read(workspace.project, workspace.source, source_item),)
            else
                loaded = something(node.input)
                node.input = nothing
                interpretation = interpret_source_item(
                    workspace.project, workspace.source, source_item, loaded;
                    source_item_key=key.entity::Int64)
                (source_item=source_item, interpretation=interpretation)
            end
        elseif key.kind === ITEM_PROCESS
            index = workspace.index
            record = get(index.items, key.entity, nothing)
            record === nothing && error("Cannot process removed item '$(key.entity)'")
            processing = run_processing(workspace, index.collections, record)
            if work_node_current(workspace, node)
                store_processed!(
                    workspace.cache.db, processing.record, item_data(processing.item))
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
                for position in processing.rewritten
                    store_processed!(
                        workspace.cache.db,
                        processing.records[position],
                        item_data(processing.outputs[position]);
                        stage=COLLECTION_PROCESS)
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

Always ITEM_PROCESS. A collection `process` that rewrites this member supersedes the delivered
payload when its fold lands, through the normal publish path — delivery does not wait for it, so a
project with collection stages still shows item-processed data while the fold runs.
"""
delivery_gate(::Workspace, record::ItemRecord)::WorkKey = WorkKey(ITEM_PROCESS, record.id)

"""
Read one record's delivered payload from the cache: the deepest stage present.

A collection `process` that rewrote this member left a `:collection_processed` payload, and that
supersedes the member's own. Absent one, the member's `:processed` payload is what it delivers —
including while a fold is still running.
"""
function _delivered_payload(
    workspace::Workspace,
    collections::CollectionIndex,
    record::ItemRecord,
)::Any
    delivered = ItemRecord(
        record; metadata=delivered_metadata(workspace, record, collections))
    folded = only(read_payload(
        workspace.cache.db, [delivered]; stage=COLLECTION_PROCESS))
    folded === nothing || return folded
    return only(read_payload(workspace.cache.db, [delivered]; stage=ITEM_PROCESS))
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
    if key.kind === SOURCE_INTERPRET
        upstream = cached_result_state(workspace.cache.db, SOURCE_READ, key.entity)
        if upstream !== nothing && CacheResultStatus(upstream.status) === RESULT_FAILED
            return ProcessingResult(nothing, CapturedException(
                ErrorException(something(upstream.message, "Source read failed")), Any[]))
        end
        enqueue_source_interpretation!(workspace, key.entity; priority, waiter)
    else
        enqueue_work!(workspace, key, current_revision(workspace.work, key);
            priority, dependencies, waiter)
    end
    return take!(waiter)::ProcessingResult
end

"""Return processed items, joining and promoting shared work when required."""
function request_processed_items(
    workspace::Workspace,
    records::Vector{ItemRecord},
)::Vector{AbstractDataItem}
    index = workspace.index
    collections = index.collections
    loaded_item(record, item) = materialized_item(
        workspace,
        collections,
        ItemRecord(record; metadata=delivered_metadata(workspace, record, collections)),
        item,
    )
    loaded = Vector{AbstractDataItem}(undef, length(records))
    for (position, record) in pairs(records)
        result = ensure_uptodate!(workspace, delivery_gate(workspace, record); priority=4)
        result.failure === nothing || throw(result.failure)
        cached = _delivered_payload(workspace, collections, record)
        cached === nothing && error(
            "Delivered data for item '$(record.id)' is missing from the cache")
        loaded[position] = loaded_item(record, something(cached))
    end
    return loaded
end
