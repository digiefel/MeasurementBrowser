"""
Materialize the loaded items for selected records.
"""
function materialize_items(
    workspace::Workspace,
    records::Vector{ItemRecord},
)::Vector{Any}
    loaded = Any[]
    sizehint!(loaded, length(records))
    for record in records
        key = item_record_key(record)
        fingerprint = (record.source_item_fingerprint, record.item_fingerprint)
        cached = get(workspace.loaded_items, key, nothing)
        if cached !== nothing && first(cached) == fingerprint
            push!(loaded, last(cached))
            continue
        end
        # FIXME: The low-level contract reloads by source-item id and item id, so sources may need
        # to re-find context they had during discovery. Pass a richer handle if this becomes painful.
        item = load_data_item(workspace.source, record.source_item_id, record.item_id)
        workspace.loaded_items[key] = (fingerprint, item)
        push!(loaded, item)
    end
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
