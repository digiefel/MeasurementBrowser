using Tables

using DataBrowserAPI: item_data, label
using DataBrowserAPI.ItemIndex: ItemRecord
import DataBrowserCore.Workspace

export ItemTable, merge_item_tables, materialize_item_table

"""
Table view over one or more items' Tables-compatible payloads.

`columns` is the union of all item columns in stable order.
`row_item[r]` is the 1-based index into `item_labels` for row `r`.
`getcell(row, col)` returns display text for any cell; `getvalue(row, col)` returns the raw
typed value (`missing` when the row's item lacks that column). Consumers that compute — plotting,
fitting, export — must read `getvalue`, never parse display text.

When only one item is selected, `item_labels` has one entry and all `row_item` values are 1;
the provenance chrome is suppressed at render time.
"""
struct ItemTable
    columns::Vector{String}
    rows::Int
    row_item::Vector{Int}
    item_labels::Vector{String}
    getcell::Function   # (row::Int, col::Int) -> String
    getvalue::Function  # (row::Int, col::Int) -> Any
end

"""Return display text for one table cell value."""
function _cell_text(value::Any)::String
    text = sprint(show, value)
    return length(text) > 90 ? first(text, 87) * "..." : text
end

function _table_rowcount(table)::Int
    rowcount = Tables.rowcount(table)
    rowcount === nothing && return count(_ -> true, Tables.rows(table))
    return rowcount
end

function _column_name_strings(table)::Vector{String}
    return [string(name) for name in Tables.columnnames(table)]
end

function _append_table!(
    col_set::Set{String},
    columns::Vector{String},
    tables::Vector,
    labels::Vector{String},
    label::AbstractString,
    table,
)::Nothing
    for c in _column_name_strings(table)
        if c ∉ col_set
            push!(col_set, c)
            push!(columns, c)
        end
    end
    push!(tables, table)
    push!(labels, String(label))
    return nothing
end

function _item_table_from_tables(
    columns::Vector{String},
    tables::Vector,
    labels::Vector{String},
)::ItemTable
    isempty(tables) &&
        return ItemTable(columns, 0, Int[], labels, (_, _) -> "", (_, _) -> missing)

    row_item = Int[]
    row_offsets = Int[]
    for (i, table) in enumerate(tables)
        for r in 1:_table_rowcount(table)
            push!(row_item, i)
            push!(row_offsets, r)
        end
    end

    total_rows = length(row_item)
    col_indices = [
        Dict(c => j for (j, c) in enumerate(_column_name_strings(table)))
        for table in tables
    ]

    function getcell(row::Int, col::Int)::String
        item_i = row_item[row]
        row_in_item = row_offsets[row]
        table = tables[item_i]
        col_name = columns[col]
        ci = get(col_indices[item_i], col_name, nothing)
        ci === nothing && return ""
        column = Tables.getcolumn(table, ci)
        return _cell_text(column[row_in_item])
    end

    function getvalue(row::Int, col::Int)::Any
        item_i = row_item[row]
        ci = get(col_indices[item_i], columns[col], nothing)
        ci === nothing && return missing
        return Tables.getcolumn(tables[item_i], ci)[row_offsets[row]]
    end

    return ItemTable(columns, total_rows, row_item, labels, getcell, getvalue)
end

"""
Build an `ItemTable` from a list of `(label, table)` pairs.

Pairs whose data does not satisfy `Tables.istable` are skipped with a warning, so one non-tabular
item never hides its tabular siblings. Multiple items are merged by column union (missing columns
render blank) with per-row provenance.

Returns `(table::ItemTable, warnings::Vector{String})`.
"""
function merge_item_tables(pairs)::Tuple{ItemTable,Vector{String}}
    col_set = Set{String}()
    columns = String[]
    tables = Any[]
    labels = String[]
    warnings = String[]
    for (label, table) in pairs
        if !Tables.istable(table)
            push!(warnings, "Item '$(label)' has non-tabular data; skipped.")
            continue
        end
        _append_table!(col_set, columns, tables, labels, string(label), Tables.columns(table))
    end
    return _item_table_from_tables(columns, tables, labels), warnings
end

"""
    materialize_item_table(workspace, records) -> (table, warnings)

Load the processed items for `records` and combine their Tables-compatible payloads into an
`ItemTable`. Record labels identify each item's rows. Payload extraction failures and non-tabular
payloads are reported in `warnings` and skipped; materialization and table construction errors
propagate to the caller. This function does not select items or retain view state.
"""
function materialize_item_table(
    workspace::Workspace.Workspace,
    records::Vector{ItemRecord},
)::Tuple{ItemTable,Vector{String}}
    items = Workspace.materialize_items(workspace, records)
    pairs = Tuple{String,Any}[]
    warnings = String[]
    for (record, item) in zip(records, items)
        item_label = label(record)
        data = try
            item_data(item)
        catch err
            @error "Failed to load item data for table" label=item_label exception=(err, catch_backtrace())
            push!(warnings, "Item '$item_label': failed to load data; skipped.")
            continue
        end
        push!(pairs, (item_label, data))
    end
    table, merge_warnings = merge_item_tables(pairs)
    append!(warnings, merge_warnings)
    return table, warnings
end
