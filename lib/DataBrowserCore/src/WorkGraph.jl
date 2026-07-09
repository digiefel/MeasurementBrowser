"""Live background job graph and priority scheduling queue."""
module WorkGraph

import DataBrowserAPI: AbstractDataSourceItem

@enum WorkKind begin
    SOURCE_INTERPRET
    ITEM_PROCESS
    ITEM_ANALYZE
    COLLECTION_PROCESS
    COLLECTION_ANALYZE
end

"""Stable identity of one result-producing operation."""
struct WorkKey
    kind::WorkKind
    entity::String
end

"""One revision of a live background job (`:waiting`, `:queued`, or `:running` only)."""
mutable struct WorkNode
    key::WorkKey
    revision::UInt16
    state::Symbol
    priority::Int
    dependents::Set{WorkKey}
    pending::UInt64
    waiters::Vector{Channel{Any}}
    queued_ns::UInt64
end

"""Workspace-owned live job graph and scheduling queue (FIFO buckets keyed by priority)."""
mutable struct WorkDependencyGraph
    lock::ReentrantLock
    condition::Base.Threads.Condition
    nodes::Dict{WorkKey,WorkNode}
    queue::Dict{Int,Vector{Tuple{WorkKey,UInt16}}}
    source_items::Dict{String,AbstractDataSourceItem}
    source_locks::Dict{String,ReentrantLock}
    workers::Vector{Task}
    total::Int
    completed::Int
    active::Int
    source_batch::Int
    source_batch_open::Bool
    closed::Bool
end

function WorkDependencyGraph()::WorkDependencyGraph
    work_lock = ReentrantLock()
    return WorkDependencyGraph(
        work_lock,
        Base.Threads.Condition(work_lock),
        Dict{WorkKey,WorkNode}(),
        Dict{Int,Vector{Tuple{WorkKey,UInt16}}}(),
        Dict{String,AbstractDataSourceItem}(),
        Dict{String,ReentrantLock}(),
        Task[],
        0,
        0,
        0,
        0,
        false,
        false,
    )
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

end
