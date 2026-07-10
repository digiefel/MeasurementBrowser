using Tables

export InspectorTable, merge_item_tables

"""
Merged table built from one or more materialized items' tabular `.data` payloads.

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

function _inspector_table_from_tables(
    columns::Vector{String},
    tables::Vector,
    labels::Vector{String},
)::InspectorTable
    isempty(tables) && return InspectorTable(columns, 0, Int[], labels, (_, _) -> "")

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

    return InspectorTable(columns, total_rows, row_item, labels, getcell)
end

"""
Build an `InspectorTable` from a list of `(label, table)` pairs.

Multiple items are merged by column union (missing columns render blank) with per-row provenance.
"""
function merge_item_tables(pairs::Vector{Tuple{Any,Any}})::InspectorTable
    col_set = Set{String}()
    columns = String[]
    tables = Any[]
    labels = String[]
    for (label, table) in pairs
        _append_table!(col_set, columns, tables, labels, string(label), table)
    end
    return _inspector_table_from_tables(columns, tables, labels)
end

"""
Build an `InspectorTable` from items whose data and labels are extracted via callbacks.
"""
function merge_item_tables(
    items::Vector,
    get_data::Function,
    get_label::Function,
)::InspectorTable
    pairs = Tuple{Any,Any}[(get_label(item), get_data(item)) for item in items]
    return merge_item_tables(pairs)
end
