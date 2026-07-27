"""
Preview of one arbitrary delimited text table.

This is intentionally not item data. It keeps the file's own columns and only records how the table
was detected so the browser can show the user what it found. Reading a file *as data* is a project
concern, served by the recipes in `DataBrowserRecipes`; a source only has to discover files and show
enough of one for a person to recognize it, so this deliberately owns no format parser.

`table` holds one vector per entry in `columns`, positionally aligned.
"""
struct TabularFileSource
    path::String
    delimiter::Char
    header_row::Union{Nothing,Int}
    data_start_row::Int
    columns::Vector{String}
    row_count::Int
    preview_rows::Int
    table::Vector{AbstractVector}
    warnings::Vector{String}
end

"""Inspect a delimited text file and return table metadata plus loaded rows."""
function inspect_table(
    path::AbstractString;
    max_rows::Integer=typemax(Int),
)::TabularFileSource
    max_rows > 0 || throw(ArgumentError("max_rows must be positive"))
    filepath = normpath(String(path))
    isfile(filepath) || throw(ArgumentError("Not a file: $filepath"))

    lines = readlines(filepath)
    layout = detect_table_layout(lines)
    column_names, table = read_preview_table(lines, layout; max_rows=Int(max_rows))
    row_count = count_data_lines(lines, layout.data_start_row)
    warnings = row_count > max_rows ?
        ["Showing first $(max_rows) rows of approximately $(row_count)."] :
        String[]

    return TabularFileSource(
        filepath,
        layout.delimiter,
        layout.header_row,
        layout.data_start_row,
        column_names,
        row_count,
        isempty(table) ? 0 : length(first(table)),
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

"""
Read the detected table into column vectors, preserving the file's own column names.

Splitting is the same `table_fields` the layout detection already used, so a preview needs no
parser beyond it. Ragged rows are padded rather than rejected: the point is to show the user what
is in the file, including that it is ragged.
"""
function read_preview_table(
    lines::Vector{String},
    layout::NamedTuple;
    max_rows::Int,
)::Tuple{Vector{String},Vector{AbstractVector}}
    header = layout.header_row === nothing ? String[] :
        table_fields(lines[layout.header_row], layout.delimiter)
    rows = Vector{String}[]
    for line in @view(lines[min(layout.data_start_row, length(lines) + 1):end])
        length(rows) >= max_rows && break
        fields = table_fields(line, layout.delimiter)
        isempty(fields) && continue
        push!(rows, fields)
    end
    width = maximum(length, rows; init=length(header))
    names = String[
        index <= length(header) && !isempty(header[index]) ? header[index] : "column_$index"
        for index in 1:width
    ]
    columns = AbstractVector[
        _typed_column(String[index <= length(row) ? row[index] : "" for row in rows])
        for index in 1:width
    ]
    return names, columns
end

"""
Narrow one column of raw fields to the most specific type its every value supports.

Empty fields become `missing`, which is what makes a numeric column with gaps stay numeric.
"""
function _typed_column(fields::Vector{String})::AbstractVector
    present = [field for field in fields if !isempty(field)]
    for type in (Int64, Float64)
        all(field -> tryparse(type, field) !== nothing, present) || continue
        any(isempty, fields) || return type[parse(type, field) for field in fields]
        return Union{Missing,type}[
            isempty(field) ? missing : parse(type, field) for field in fields]
    end
    any(isempty, fields) || return fields
    return Union{Missing,String}[isempty(field) ? missing : field for field in fields]
end

"""Count non-empty source lines after the detected table start."""
function count_data_lines(lines::Vector{String}, data_start_row::Int)::Int
    data_start_row <= length(lines) || return 0
    return count(line -> !isempty(strip(line)), @view(lines[data_start_row:end]))
end
