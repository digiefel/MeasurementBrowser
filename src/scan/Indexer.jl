function _read_csv_header_summary(path::AbstractString; max_lines::Int=50)
    summary = Dict{String,String}()
    open(path, "r") do io
        line_count = 0
        for line in eachline(io)
            line_count += 1
            stripped = strip(line)
            if startswith(stripped, "#")
                m = match(r"^#\s*([^:]+):\s*(.*)$", stripped)
                m !== nothing && (summary[strip(m.captures[1])] = strip(m.captures[2]))
            elseif startswith(line, "Setup title,")
                parts = split(line, ','; limit=2)
                length(parts) > 1 && (summary["Setup title"] = strip(parts[2], '"'))
            elseif startswith(line, "Test date,")
                parts = split(line, ','; limit=2)
                length(parts) > 1 && (summary["Test date"] = strip(parts[2]))
            elseif startswith(line, "Test time,")
                parts = split(line, ','; limit=2)
                length(parts) > 1 && (summary["Test time"] = strip(parts[2]))
            end
            line_count >= max_lines && break
        end
    end
    return summary
end

function file_fingerprint(path::AbstractString)
    normalized = normpath(abspath(expanduser(String(path))))
    stat_info = stat(normalized)
    return FileFingerprint(normalized, Int64(stat_info.size), Int64(stat_info.mtime * 1_000_000_000))
end

function index_source_file(path::AbstractString)
    normalized = normpath(abspath(expanduser(String(path))))
    filename = basename(normalized)
    return SourceFile(
        normalized,
        normalized,
        filename,
        parse_timestamp(filename),
        _read_csv_header_summary(normalized),
        file_fingerprint(normalized),
        [],
    )
end

"""Return whether a directory entry is a visible CSV source file."""
function is_source_filename(name::AbstractString)::Bool
    return endswith(lowercase(name), ".csv") && !startswith(name, ".")
end

function walk_source_files(
    root_path::AbstractString;
    on_file::Function,
)
    for (root, _, names) in walkdir(root_path)
        for name in names
            _check_cancel()
            is_source_filename(name) || continue
            on_file(index_source_file(joinpath(root, name)))
        end
    end
    return nothing
end

function collect_source_files(
    root_path::AbstractString;
    on_file::Union{Nothing,Function}=nothing,
)
    source_files = SourceFile[]
    walk_source_files(root_path; on_file=file -> begin
        push!(source_files, file)
        on_file !== nothing && on_file(file, length(source_files))
    end)
    return source_files
end
