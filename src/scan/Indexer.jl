function _read_csv_header_summary(path::AbstractString; max_lines::Int=50)
    summary = Dict{Symbol,Any}()
    open(path, "r") do io
        line_count = 0
        for line in eachline(io)
            line_count += 1
            if startswith(line, "Setup title,")
                parts = split(line, ',')
                length(parts) > 1 && (summary[:setup_title] = strip(parts[2], '"'))
            elseif startswith(line, "Test date,")
                parts = split(line, ',')
                length(parts) > 1 && (summary[:test_date] = strip(parts[2]))
            elseif startswith(line, "Test time,")
                parts = split(line, ',')
                length(parts) > 1 && (summary[:test_time] = strip(parts[2]))
            elseif startswith(line, "Device ID,")
                parts = split(line, ',')
                length(parts) > 1 && (summary[:device_id] = strip(parts[2]))
            end
            line_count >= max_lines && break
        end
    end
    return summary
end

function index_csv_file(path::AbstractString)
    filename = basename(String(path))
    return IndexedCsvFile(
        file_id(path),
        String(path),
        filename,
        parse_timestamp(filename),
        _read_csv_header_summary(path),
    )
end

function walk_indexed_csv_files(
    root_path::AbstractString;
    on_file::Function,
    should_cancel::Union{Nothing,Function}=nothing,
)
    for (root, _, names) in walkdir(root_path)
        for name in names
            should_cancel !== nothing && should_cancel() && throw(ScanCancelled())
            endswith(lowercase(name), ".csv") || continue
            on_file(index_csv_file(joinpath(root, name)))
        end
    end
    return nothing
end

function collect_indexed_csv_files(
    root_path::AbstractString;
    should_cancel::Union{Nothing,Function}=nothing,
)
    indexed_files = IndexedCsvFile[]
    walk_indexed_csv_files(root_path; should_cancel=should_cancel, on_file=file -> push!(indexed_files, file))
    return indexed_files
end
