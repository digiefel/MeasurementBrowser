"""
Materialize the loaded items for selected records.

Each record is served from the in-memory LRU when fresh, else from the DuckDB payload cache (no
origin access), else by reloading from the origin — whose payload is then written back so later opens
and evictions need not re-read it. Only items the project marks [`cacheable`](@ref) are persisted.
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

    # The DuckDB payload cache may hold the raw loaded item, sparing an origin read.
    payloads = read_item_payloads(workspace.cache.db, misses)
    fresh_items = ItemRecord[]
    fresh_payloads = Any[]
    for (index, record) in pairs(misses)
        raw = payloads[index]
        if raw === nothing
            # FIXME: The low-level contract reloads by source-item id and id, so sources may need
            # to re-find context they had during discovery. Pass a richer handle if this gets painful.
            raw = load_data_item(
                workspace.project, workspace.source, record.source_item_id, record.id)
            if cacheable(raw)
                push!(fresh_items, record)
                push!(fresh_payloads, raw)
            end
        end
        item = _effective_loaded_item(record, raw)
        fingerprint = (record.source_item_fingerprint, record.item_fingerprint)
        cache_item!(cache, record.id, fingerprint, item)
        loaded[miss_positions[index]] = item
    end
    isempty(fresh_items) ||
        write_item_payloads!(workspace.cache.db, fresh_items, fresh_payloads)
    return loaded
end

"""
Return loaded item payloads for callers that inspect raw tabular data.
"""
function read_item_data(
    workspace::Workspace,
    records::Vector{ItemRecord},
)::Vector{Any}
    return Any[item_data(item) for item in materialize_items(workspace, records)]
end
