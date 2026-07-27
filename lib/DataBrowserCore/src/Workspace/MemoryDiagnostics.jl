"""A cheap point-in-time memory attribution snapshot for benchmarks and internal profiling."""
struct WorkspaceMemorySnapshot
    sampled_ns::UInt64
    gc_live_bytes::Int64
    gc_allocated_bytes::Int64
    index_items::Int64
    index_collections::Int64
    item_metadata::Int64
    analysis_errors::Int64
    processing_jobs::Int64
    pending_writes::Int64
    pending_write_rows::Int64
    selected_queue::Int64
    background_waiting::Int64
end

"""Return an attributed workspace memory sample without walking large object graphs."""
function workspace_memory_snapshot(workspace::Workspace)::WorkspaceMemorySnapshot
    gc = Base.gc_num()
    gc_live_bytes = Int64(Base.gc_live_bytes())
    collections = workspace.index.collections
    _completed, _total, active = work_counts(workspace)
    processing = (jobs=Int64(active), selected_queue=Int64(0), background_waiting=Int64(0))
    # Writes live in the cache buffers, not in the work graph.
    staged = cache_pending_counts(workspace.cache.db)
    return WorkspaceMemorySnapshot(
        time_ns(),
        gc_live_bytes,
        Int64(gc.total_allocd),
        Int64(length(workspace.index.items)),
        Int64(length(collections.records)),
        Int64(length(workspace.index.item_metadata)),
        Int64(length(workspace.index.analysis_errors)),
        processing.jobs,
        Int64(staged.items),
        Int64(staged.rows),
        processing.selected_queue,
        processing.background_waiting,
    )
end
