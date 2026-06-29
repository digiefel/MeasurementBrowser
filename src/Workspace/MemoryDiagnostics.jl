"""A cheap point-in-time memory attribution snapshot for benchmarks and internal profiling."""
struct WorkspaceMemorySnapshot
    sampled_ns::UInt64
    rss_bytes::Int64
    gc_live_bytes::Int64
    gc_allocated_bytes::Int64
    rss_minus_gc_live_bytes::Int64
    index_items::Int64
    index_collections::Int64
    item_stats::Int64
    analysis_errors::Int64
    processing_jobs::Int64
    pending_writes::Int64
    pending_write_rows::Int64
    selected_queue::Int64
    background_waiting::Int64
end

"""Return an attributed workspace memory sample without walking large object graphs."""
function workspace_memory_snapshot(workspace::Workspace)::WorkspaceMemorySnapshot
    rss_bytes = Profiling.process_rss_bytes()
    gc = Base.gc_num()
    gc_live_bytes = Int64(Base.gc_live_bytes())
    hierarchy = workspace.index.hierarchy
    processing = lock(workspace.processing.lock) do
        (
            jobs=Int64(length(workspace.processing.jobs)),
            selected_queue=Int64(length(workspace.processing.selected)),
            background_waiting=Int64(max(
                length(workspace.processing.background) - workspace.processing.background_index + 1,
                0,
            )),
        )
    end
    # Writes now live in the cache buffer, not the processing queue.
    staged = buffer_pending_counts(workspace.cache.db.buffer)
    return WorkspaceMemorySnapshot(
        time_ns(),
        rss_bytes,
        gc_live_bytes,
        Int64(gc.total_allocd),
        max(Int64(0), rss_bytes - gc_live_bytes),
        Int64(length(workspace.index.items)),
        Int64(length(hierarchy.index)),
        Int64(length(workspace.index.item_stats)),
        Int64(length(workspace.index.analysis_errors)),
        processing.jobs,
        Int64(staged.items),
        Int64(staged.rows),
        processing.selected_queue,
        processing.background_waiting,
    )
end
