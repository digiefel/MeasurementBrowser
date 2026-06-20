"""
Materialize the loaded items for selected records.

Each record is served from the in-memory LRU when fresh, else from the DuckDB item-data cache (one
snapshot read, no origin access), else by reloading from the origin. The interactive path performs no
writes: source scanning populates the durable cache while each loaded item is already in memory, so a
plot never pays a database write. Within a session the LRU keeps reselected items hot.
"""
function materialize_items(
    workspace::Workspace,
    records::Vector{ItemRecord},
)::Vector{Any}
    cache = workspace.loaded_items
    loaded = Vector{Any}(undef, length(records))
    misses = ItemRecord[]
    miss_positions = Int[]
    for (position, record) in pairs(records)
        fingerprint = (record.source_item_fingerprint, record.item_fingerprint)
        entry = lookup_item(cache, record.id, nothing)
        if entry !== nothing && first(entry) == fingerprint
            loaded[position] = last(entry)
        else
            push!(misses, record)
            push!(miss_positions, position)
        end
    end
    isempty(misses) && return loaded

    # DuckDB may already hold the loaded item data, sparing a source read.
    cached_data = read_cached_item_data(workspace.cache.db, misses)
    for (index, record) in pairs(misses)
        raw = cached_data[index]
        if raw === nothing
            # FIXME: The low-level contract reloads by source-item id and id, so sources may need
            # to re-find context they had during discovery. Pass a richer handle if this gets painful.
            raw = @timeit_debug TIMER "read/load_origin" load_data_item(
                workspace.project, workspace.source, record.source_item_id, record.id)
        end
        item = _effective_loaded_item(record, raw)
        fingerprint = (record.source_item_fingerprint, record.item_fingerprint)
        cache_item!(cache, record.id, fingerprint, item)
        loaded[miss_positions[index]] = item
    end
    return loaded
end

"""
Return loaded item data for callers that inspect raw tabular values.
"""
function read_item_data(
    workspace::Workspace,
    records::Vector{ItemRecord},
)::Vector{Any}
    return Any[item_data(item) for item in materialize_items(workspace, records)]
end
