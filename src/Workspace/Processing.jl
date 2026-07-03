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

"""Start fixed worker pools for each work kind."""
function start_work_workers!(workspace::Workspace)::Nothing
    worker_count = max(1, Base.Threads.nthreads() ÷ 4)
    for kind in (INTERPRET_SOURCE, PROCESS_ITEM, ITEM_STATS, COLLECTION_STATS)
        for _ in 1:worker_count
            push!(workspace.work.workers, Base.Threads.@spawn work_worker!(workspace, kind))
        end
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
    lock(workspace.poll_lock) do
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
    kind === COLLECTION_STATS && return work_state_terminal(dependency.state)
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

function queue_ready_node!(graph::WorkDependencyGraph, node::WorkNode)::Nothing
    node.state = :queued
    node.queued_ns = time_ns()
    push!(graph.queue, (node.key, node.revision))
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

"""Record a current missing work revision without queuing it."""
function mark_work_missing!(
    workspace::Workspace,
    key::WorkKey,
    revision::UInt64;
    dependencies::Vector{WorkKey}=WorkKey[],
)::Nothing
    graph = workspace.work
    lock(graph.lock) do
        node = get(graph.nodes, key, nothing)
        if node === nothing || node.revision != revision || node.state !== :running
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
        notify(graph.condition; all=true)
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
        notify(graph.condition; all=true)
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
            dependencies_ready(graph, node) && (node.state = :queued)
            replace_work_node!(graph, node)
            node.state === :queued && push!(graph.queue, (key, revision))
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
            node.priority = priority
        end
        waiter === nothing || push!(node.waiters, waiter)
        notify(graph.condition; all=true)
        node
    end
end

"""Take the highest-priority current work revision."""
function take_work!(
    graph::WorkDependencyGraph,
    kind::WorkKind,
)::Union{Nothing,WorkNode}
    while true
        graph.closed && return nothing
        best_position = 0
        best_priority = typemin(Int)
        for (position, (key, revision)) in pairs(graph.queue)
            node = get(graph.nodes, key, nothing)
            node === nothing && continue
            node.key.kind === kind && node.revision == revision && node.state === :queued ||
                continue
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
        wait(graph.condition)
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

"""Execute one work node and publish its completion immediately."""
function execute_work!(workspace::Workspace, node::WorkNode)::Nothing
    key = node.key
    result = try
        if key.kind === INTERPRET_SOURCE
            source_item = lock(workspace.work.lock) do
                get(workspace.work.source_items, key.entity, nothing)
            end
            source_item === nothing &&
                error("Cannot interpret removed source item '$(key.entity)'")
            interpretation = interpret_source_item(
                workspace.project, workspace.source, source_item, workspace.profiler)
            if work_node_current(workspace, node)
                store_interpreted_data!(
                    workspace.cache.db,
                    interpretation.records,
                    interpretation.interpreted_items,
                )
            end
            (records=interpretation.records, failures=interpretation.failures)
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
    publish_work_completion!(workspace, key, node.revision, result)
    return nothing
end

"""Run one fixed work-kind queue until workspace shutdown."""
function work_worker!(workspace::Workspace, kind::WorkKind)::Nothing
    graph = workspace.work
    while true
        node = lock(graph.lock) do
            take_work!(graph, kind)
        end
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
        WorkKey(PROCESS_ITEM, record.id),
        item_revision(record);
        priority=selected ? 4 : 2,
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
        result = take!(waiter)::ProcessingResult
        result.failure === nothing || throw(result.failure)
        loaded[position] = result.item::AbstractDataItem
    end
    return loaded
end
