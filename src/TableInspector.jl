module TableInspector

using DataFrames: DataFrame, names, nrow

export InspectorTable, merge_item_tables

"""
Merged table built from one or more materialized items' `.data` DataFrames.

`columns` is the union of all item columns in stable order.
`row_item[r]` is the 1-based index into `item_labels` for row `r`.
`getcell(row, col)` returns display text for any cell.

When only one item is selected, `item_labels` has one entry and all `row_item` values are 1;
the provenance chrome is suppressed at render time.
"""
struct InspectorTable
    columns::Vector{String}
    rows::Int
    row_item::Vector{Int}
    item_labels::Vector{String}
    getcell::Function   # (row::Int, col::Int) -> String
end

"""Return display text for one table cell value."""
function _cell_text(value::Any)::String
    text = sprint(show, value)
    return length(text) > 90 ? first(text, 87) * "..." : text
end

"""
Build an `InspectorTable` from a list of `(item, data)` pairs, where `data` is a DataFrame.

Multiple items are merged by column union (missing columns render blank) with per-row provenance.
Non-DataFrame data is skipped and noted in `warnings` (returned separately).

Returns `(table::InspectorTable, warnings::Vector{String})`.
"""
function merge_item_tables(
    pairs::Vector{Tuple{Any,Any}},    # (item, data)
)::Tuple{InspectorTable,Vector{String}}
    warnings = String[]
    col_set = Set{String}()
    columns = String[]
    dfs = DataFrame[]
    labels = String[]

    for (item, data) in pairs
        if !(data isa DataFrame)
            push!(warnings, "Item '$(item)' has non-tabular data ($(typeof(data))); skipped.")
            continue
        end
        cols = String.(names(data))
        for c in cols
            if c ∉ col_set
                push!(col_set, c)
                push!(columns, c)
            end
        end
        push!(dfs, data)
        push!(labels, string(item))
    end

    if isempty(dfs)
        empty_table = InspectorTable(columns, 0, Int[], labels, (_, _) -> "")
        return empty_table, warnings
    end

    row_item = Int[]
    row_offsets = Int[]
    for (i, df) in enumerate(dfs)
        for r in 1:nrow(df)
            push!(row_item, i)
            push!(row_offsets, r)
        end
    end

    total_rows = length(row_item)
    col_indices = [
        Dict(c => j for (j, c) in enumerate(String.(names(df))))
        for df in dfs
    ]

    function getcell(row::Int, col::Int)::String
        item_i = row_item[row]
        row_in_item = row_offsets[row]
        df = dfs[item_i]
        col_name = columns[col]
        ci = get(col_indices[item_i], col_name, nothing)
        ci === nothing && return ""
        return _cell_text(df[row_in_item, ci])
    end

    table = InspectorTable(columns, total_rows, row_item, labels, getcell)
    return table, warnings
end

"""
Convenience: build an InspectorTable from items whose `item_data` is accessed via a callback.

`get_data(item) -> Any` extracts the loaded data.
`get_label(item) -> String` extracts the display label.
"""
function merge_item_tables(
    items::Vector,
    get_data::Function,
    get_label::Function,
)::Tuple{InspectorTable,Vector{String}}
    pairs = [(item, get_data(item)) for item in items]
    labeled_pairs = Tuple{Any,Any}[(item, data) for (item, data) in pairs]
    warnings = String[]
    col_set = Set{String}()
    columns = String[]
    dfs = DataFrame[]
    labels = String[]

    for (item, data) in labeled_pairs
        if !(data isa DataFrame)
            push!(warnings, "Item '$(get_label(item))' has non-tabular data ($(typeof(data))); skipped.")
            continue
        end
        cols = String.(names(data))
        for c in cols
            if c ∉ col_set
                push!(col_set, c)
                push!(columns, c)
            end
        end
        push!(dfs, data)
        push!(labels, get_label(item))
    end

    if isempty(dfs)
        return InspectorTable(columns, 0, Int[], labels, (_, _) -> ""), warnings
    end

    row_item = Int[]
    row_offsets = Int[]
    for (i, df) in enumerate(dfs)
        for r in 1:nrow(df)
            push!(row_item, i)
            push!(row_offsets, r)
        end
    end

    total_rows = length(row_item)
    col_indices = [
        Dict(c => j for (j, c) in enumerate(String.(names(df))))
        for df in dfs
    ]

    function getcell(row::Int, col::Int)::String
        item_i = row_item[row]
        row_in_item = row_offsets[row]
        df = dfs[item_i]
        col_name = columns[col]
        ci = get(col_indices[item_i], col_name, nothing)
        ci === nothing && return ""
        return _cell_text(df[row_in_item, ci])
    end

    return InspectorTable(columns, total_rows, row_item, labels, getcell), warnings
end

end
