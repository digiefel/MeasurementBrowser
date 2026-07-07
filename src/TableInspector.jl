module TableInspector

using CSV
using DataFrames: DataFrame, names, nrow, hcat

export TablePreview, inspect_table
export InspectorTable, merge_item_tables

const DEFAULT_PREVIEW_ROWS = 500

"""
Preview of one arbitrary delimited text table.

This is intentionally not item data. It keeps the file's own columns and only records how
the table was detected so the browser can show the user what it found.
"""
struct TablePreview
    path::String
    delimiter::Char
    header_row::Union{Nothing,Int}
    data_start_row::Int
    columns::Vector{String}
    row_count::Int
    preview_rows::Int
    table::DataFrame
    warnings::Vector{String}
end

"""Inspect a delimited text file and return a bounded preview plus detected table metadata."""
function inspect_table(
    path::AbstractString;
    max_rows::Integer=DEFAULT_PREVIEW_ROWS,
)::TablePreview
    max_rows > 0 || throw(ArgumentError("max_rows must be positive"))
    filepath = normpath(String(path))
    isfile(filepath) || throw(ArgumentError("Not a file: $filepath"))

    lines = readlines(filepath)
    layout = detect_table_layout(lines)
    table = read_preview_table(filepath, layout; max_rows=Int(max_rows))
    row_count = count_data_lines(lines, layout.data_start_row)
    warnings = row_count > max_rows ?
        ["Showing first $(max_rows) rows of approximately $(row_count)."] :
        String[]

    return TablePreview(
        filepath,
        layout.delimiter,
        layout.header_row,
        layout.data_start_row,
        String.(names(table)),
        row_count,
        nrow(table),
        table,
        warnings,
    )
end

"""Detect delimiter, optional header row, and first data row from file text."""
function detect_table_layout(lines::Vector{String})::NamedTuple
    delimiter = choose_delimiter(lines)
    for (row, line) in enumerate(lines)
        fields = table_fields(line, delimiter)
        length(fields) >= 2 || continue
        numeric = count(field -> tryparse(Float64, field) !== nothing, fields)
        has_header = numeric < length(fields)
        return (
            delimiter=delimiter,
            header_row=has_header ? row : nothing,
            data_start_row=has_header ? row + 1 : row,
        )
    end
    throw(ArgumentError("Could not find a table with at least two columns"))
end

"""Choose the delimiter that creates the most multi-column rows near the top of the file."""
function choose_delimiter(lines::Vector{String})::Char
    candidates = (',', '\t', ';')
    scores = Dict(
        delimiter => sum(max(length(table_fields(line, delimiter)) - 1, 0)
                         for line in Iterators.take(lines, 80))
        for delimiter in candidates
    )
    delimiter = first(candidates)
    best_score = -1
    for candidate in candidates
        score = scores[candidate]
        if score > best_score
            delimiter = candidate
            best_score = score
        end
    end
    best_score > 0 || throw(ArgumentError("Could not detect a delimited table"))
    return delimiter
end

"""Split a possible table row into stripped fields, ignoring blank and comment lines."""
function table_fields(line::AbstractString, delimiter::Char)::Vector{String}
    text = strip(line)
    (isempty(text) || startswith(text, "#")) && return String[]
    return strip.(String.(split(text, delimiter; keepempty=true)))
end

"""Read the detected table with CSV.jl while preserving source column names."""
function read_preview_table(
    filepath::AbstractString,
    layout::NamedTuple;
    max_rows::Int,
)::DataFrame
    options = (
        delim=layout.delimiter,
        skipto=layout.data_start_row,
        limit=max_rows,
        normalizenames=false,
        silencewarnings=true,
    )
    if layout.header_row === nothing
        return CSV.read(filepath, DataFrame; options..., header=false)
    end
    return CSV.read(filepath, DataFrame; options..., header=layout.header_row)
end

"""Count non-empty source lines after the detected table start."""
function count_data_lines(lines::Vector{String}, data_start_row::Int)::Int
    data_start_row <= length(lines) || return 0
    return count(line -> !isempty(strip(line)), @view(lines[data_start_row:end]))
end

# ---------------------------------------------------------------------------
# Item-data merge (pure, no GUI)
# ---------------------------------------------------------------------------

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
    # Collect column union in stable order (first-seen wins for ordering)
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

    # Build per-row item index and a combined row-offset list
    row_item = Int[]
    row_offsets = Int[]   # (item_index, row_in_item)
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
    # Reuse the pair-based overload but inject labels properly
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
