using CSV
using DataFrames: DataFrame, names, nrow

const DEFAULT_PREVIEW_ROWS = 500

"""
Preview of one arbitrary delimited text table.

This is intentionally not item data. It keeps the file's own columns and only records how
the table was detected so the browser can show the user what it found.
"""
struct TabularFileSource
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
)::TabularFileSource
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

    return TabularFileSource(
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
