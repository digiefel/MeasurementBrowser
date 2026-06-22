"""
Materialize processed items for selected records through the shared processing queue.

Committed processed data comes from DuckDB. Missing results promote the existing background job or
create one whose interpreted-data dependency may use DuckDB or the normal source fallback.
"""
function materialize_items(
    workspace::Workspace,
    records::Vector{ItemRecord},
)::Vector{AbstractDataItem}
    return request_processed_items(workspace, records)
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
