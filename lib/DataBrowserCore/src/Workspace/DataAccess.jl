"""
Materialize processed items for selected records through the workspace work graph.

Committed processed data comes from DuckDB. Missing results promote the existing background job or
create one whose interpreted-data dependency may use DuckDB or the normal source fallback.
"""
function materialize_items(
    workspace::Workspace,
    records::Vector{ItemRecord},
)::Vector{AbstractDataItem}
    return request_processed_items(workspace, records)
end

"""Materialize processed items by stable workspace ID, preserving the requested order."""
function materialize_items(
    workspace::Workspace,
    ids::AbstractVector{<:AbstractString},
)::Vector{AbstractDataItem}
    records = ItemRecord[]
    for id_value in ids
        item_id = String(id_value)
        haskey(workspace.index.items, item_id) || error(
            "Cannot materialize item id '$item_id': no indexed item with that id exists in this workspace",
        )
        push!(records, workspace.index.items[item_id])
    end
    return materialize_items(workspace, records)
end

"""Materialize the workspace's current item selection."""
materialize_items(workspace::Workspace)::Vector{AbstractDataItem} =
    materialize_items(workspace, workspace.selection.item_ids)

"""
Return loaded item data for callers that inspect raw tabular values.
"""
function read_item_data(
    workspace::Workspace,
    records::Vector{ItemRecord},
)::Vector{Any}
    return Any[item_data(item) for item in materialize_items(workspace, records)]
end

"""Return processed item payloads by stable workspace ID."""
read_item_data(workspace::Workspace, ids::AbstractVector{<:AbstractString})::Vector{Any} =
    Any[item_data(item) for item in materialize_items(workspace, ids)]

"""Return processed payloads for the workspace's current item selection."""
read_item_data(workspace::Workspace)::Vector{Any} =
    Any[item_data(item) for item in materialize_items(workspace)]
