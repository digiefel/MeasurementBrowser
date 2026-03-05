"""
DeviceParser.jl - Parse device hierarchy from measurement filenames
"""

using Dates

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
const MAX_HEADER_LINES = 50

# ---------------------------------------------------------------------------
# Project dispatch types
# ---------------------------------------------------------------------------

abstract type AbstractProject end
struct RuO2Project <: AbstractProject end
struct TASEProject <: AbstractProject end

const KNOWN_PROJECTS = AbstractProject[]
const _default_project = Ref{Union{AbstractProject,Nothing}}(nothing)

# Interface (implemented by each project via multiple dispatch):
#   accepts_file(::P, filename) → Bool
#   parse_device_info(::P, filename) → DeviceInfo
#   detect_kind(::P, filename) → Symbol
#   kind_label(::P, kind) → String
#   display_label(::P, meas) → String
#   expand_measurement(::P, meas) → Vector{MeasurementInfo}
#   figure_for_file(::P, path, kind; kwargs...) → Union{Figure,Nothing}
#   figure_for_files(::P, paths, combined_kind; kwargs...) → Union{Figure,Nothing}
#   combined_plot_types(::P) → Vector{Tuple}
#   compatible_kinds(::P, combined_kind) → Vector{Symbol}

# ---------------------------------------------------------------------------
# DeviceInfo
# ---------------------------------------------------------------------------

struct DeviceInfo
    location::Vector{String}              # variable-length hierarchy
    parameters::Dict{Symbol,Any}          # device-level metadata
end

DeviceInfo(location::Vector{String}) = DeviceInfo(location, Dict{Symbol,Any}())

# ---------------------------------------------------------------------------
# Measurement related structs
# ---------------------------------------------------------------------------
struct MeasurementInfo
    filename::String
    filepath::String
    clean_title::String
    measurement_kind::Symbol
    timestamp::Union{DateTime,Nothing}
    device_info::DeviceInfo
    parameters::Dict{Symbol,Any}
    wakeup_pulse_count::Union{Int,Nothing}
end

struct HierarchyNode
    name::String
    kind::Symbol
    children::Vector{HierarchyNode}
    measurements::Vector{MeasurementInfo}
end

struct MeasurementHierarchy
    root::HierarchyNode
    all_measurements::Vector{MeasurementInfo}
    root_path::String
    index::Dict{Tuple{Vararg{String}},HierarchyNode}
    has_device_metadata::Bool
    project::AbstractProject
    skipped_count::Int          # CSV files rejected by accepts_file in last scan
end

HierarchyNode(name::String, kind::Symbol) = HierarchyNode(name, kind, HierarchyNode[], MeasurementInfo[])

# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

function MeasurementInfo(filepath::AbstractString, project::AbstractProject)
    filename = basename(filepath)
    device_info = parse_device_info(project, filename)
    measurement_kind = detect_kind(project, filename)
    timestamp = parse_timestamp(filename)
    parameters = parse_parameters(filename)
    file_info = extract_file_info(filepath)
    exp_label = kind_label(project, measurement_kind)

    device_label = join(device_info.location, "_")

    date_str = ""
    d = get(file_info, "test_date", "")
    if !isempty(d)
        date_str = try
            parts = split(d)
            if length(parts) >= 3
                month = parts[2]
                day = try
                    parse(Int, parts[3])
                catch
                    nothing
                end
                if day !== nothing
                    "$(month)$(day)"
                else
                    d
                end
            else
                d
            end
        catch
            d
        end
    end

    parts = filter(!isempty, (exp_label == "Unknown" ? "" : exp_label, device_label, date_str))
    clean_title = isempty(parts) ? strip(replace(filename, r"\.csv$" => "")) : join(parts, " ")

    # Cache wakeup pulse count once (avoid per-frame I/O later)
    wakeup_count = nothing
    if measurement_kind == :wakeup
        try
            lines = readlines(filepath)
            data_start = 1
            for (i, line) in enumerate(lines)
                if occursin("Time,MeasResult1_value,MeasResult2_value", line)
                    data_start = i + 1
                    break
                end
            end
            cnt = 0
            for line in lines[data_start:end]
                if !isempty(strip(line)) && occursin(',', line)
                    cnt += 1
                end
            end
            wakeup_count = cnt
        catch
            # ignore, leave as nothing
        end
    end

    return MeasurementInfo(filename, filepath, clean_title, measurement_kind, timestamp, device_info, parameters, wakeup_count)
end

# Backwards-compat: use the registered default project
MeasurementInfo(filepath::AbstractString) = MeasurementInfo(filepath, _default_project[])

function extract_file_info(path::AbstractString)
    file_stat = stat(path)
    size_bytes = file_stat.size
    setup_title = test_date = test_time = device_id = ""
    line_count = 0
    open(path, "r") do io
        for line in eachline(io)
            line_count += 1
            if startswith(line, "Setup title,")
                parts = split(line, ',')
                if length(parts) > 1
                    setup_title = strip(parts[2], '"')
                end
            elseif startswith(line, "Test date,")
                parts = split(line, ',')
                length(parts) > 1 && (test_date = parts[2])
            elseif startswith(line, "Test time,")
                parts = split(line, ',')
                length(parts) > 1 && (test_time = parts[2])
            elseif startswith(line, "Device ID,")
                parts = split(line, ',')
                device_id = length(parts) > 1 ? parts[2] : ""
            end
            line_count >= MAX_HEADER_LINES && break
        end
    end
    return Dict(
        "setup_title" => setup_title,
        "test_date" => test_date,
        "test_time" => test_time,
        "device_id" => device_id,
        "size_bytes" => size_bytes,
    )
end

function parse_timestamp(filename::String)
    if (m = match(r"; (\d{4}-\d{2}-\d{2}) (\d{2})_(\d{2})_(\d{2})\]", filename)) !== nothing
        date_str, hour, minute, second = m.captures
        try
            return DateTime("$date_str $hour:$minute:$second", "yyyy-mm-dd HH:MM:SS")
        catch
            return nothing
        end
    elseif (m = match(r"_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})_", filename)) !== nothing
        year, month, day, hour, minute, second = m.captures
        try
            return DateTime("$year-$month-$day $hour:$minute:$second", "yyyy-mm-dd HH:MM:SS")
        catch
            return nothing
        end
    end
    return nothing
end

function parse_parameters(filename::String)
    params = Dict{Symbol,Any}()
    if (m = match(r"(\d+(?:\.\d+)?)V", filename)) !== nothing
        params[:voltage_V] = parse(Float64, m.captures[1])
    end
    if (m = match(r"(\d+(?:\.\d+)?)(khz|hz)", lowercase(filename))) !== nothing
        val = parse(Float64, m.captures[1])
        unit = m.captures[2]
        params[:frequency_Hz] = unit == "khz" ? val * 1e3 : val
    end
    if (m = match(r"\((\d+)\)", filename)) !== nothing
        params[:count] = parse(Int, m.captures[1])
    end
    if (m = match(r"(\d+)K", filename)) !== nothing
        params[:temperature_K] = parse(Int, m.captures[1])
    end
    return params
end

# ---------------------------------------------------------------------------
# Sorting helpers
# ---------------------------------------------------------------------------
function roman_value(s::AbstractString)
    ROMAN_MAP = Dict('I' => 1, 'V' => 5, 'X' => 10, 'L' => 50)
    isempty(s) && return nothing
    total = 0
    prev = 0
    for c in reverse(uppercase(s))
        v = get(ROMAN_MAP, c, 0)
        v == 0 && return nothing
        if v < prev
            total -= v
        else
            total += v
            prev = v
        end
    end
    return total
end

function natural_key(s::AbstractString)
    toks = eachmatch(r"\d+|\D+", String(s))
    parts = Any[]
    for t in toks
        seg = t.match
        if all(isdigit, seg)
            push!(parts, (1, parse(Int, seg)))
        else
            push!(parts, (0, lowercase(seg)))
        end
    end
    return Tuple(parts)
end

function Base.sort!(node::HierarchyNode)
    function _roman_sortable(children::Vector{HierarchyNode})
        isempty(children) && return false
        for ch in children
            roman_value(ch.name) === nothing && return false
        end
        return true
    end
    for ch in node.children
        sort!(ch)
    end
    if _roman_sortable(node.children)
        sort!(node.children, by=c -> roman_value(c.name))
    else
        sort!(node.children, by=c -> natural_key(c.name))
    end
    for ch in node.children
        sort!(ch.measurements, by=m -> m.timestamp === nothing ? DateTime(Dates.year(typemax(Date))) : m.timestamp)
    end
    return node
end

function Base.sort!(mh::MeasurementHierarchy)
    sort!(mh.root)
    return mh
end

# ---------------------------------------------------------------------------
# Hierarchy construction
# ---------------------------------------------------------------------------
function MeasurementHierarchy(measurements::Vector{MeasurementInfo}, root_path::String, has_dev_metadata::Bool, project::AbstractProject, skipped_count::Int=0)
    root = HierarchyNode("/", :root)
    index = Dict{Tuple{Vararg{String}},HierarchyNode}()
    function ensure_child(parent::HierarchyNode, name::String, kind::Symbol, path_tuple::Tuple{Vararg{String}})
        for ch in parent.children
            if ch.name == name
                return ch
            end
        end
        node = HierarchyNode(name, kind)
        push!(parent.children, node)
        index[path_tuple] = node
        return node
    end
    for m in measurements
        segs = m.device_info.location
        parent = root
        for (i, seg) in enumerate(segs)
            kind = i == length(segs) ? :leaf : :level
            path_tuple = Tuple(segs[1:i])
            parent = ensure_child(parent, seg, kind, path_tuple)
        end
        push!(parent.measurements, m)
    end
    mh = MeasurementHierarchy(root, measurements, root_path, index, has_dev_metadata, project, skipped_count)
    sort!(mh)
    return mh
end

# Backwards-compat convenience constructors
MeasurementHierarchy(measurements::Vector{MeasurementInfo}, root_path::String) =
    MeasurementHierarchy(measurements, root_path, false, _default_project[], 0)

children(node::HierarchyNode) = node.children
isleaf(node::HierarchyNode) = isempty(node.children)

# ---------------------------------------------------------------------------
# Scanning
# ---------------------------------------------------------------------------

struct ScanCancelled <: Exception end

function _check_cancel(should_cancel::Union{Nothing,Function})
    should_cancel !== nothing && should_cancel() && throw(ScanCancelled())
end

function _emit_progress(
    on_progress::Union{Nothing,Function};
    phase::Symbol,
    total_csv::Int,
    processed_csv::Int,
    loaded_measurements::Int,
    skipped_csv::Int,
    current_path::String="",
)
    on_progress === nothing && return
    on_progress((
        phase=phase,
        total_csv=total_csv,
        processed_csv=processed_csv,
        loaded_measurements=loaded_measurements,
        skipped_csv=skipped_csv,
        current_path=current_path,
    ))
end

function _count_csv(root_path::String; should_cancel::Union{Nothing,Function}=nothing, on_progress::Union{Nothing,Function}=nothing)
    total = 0
    for (root, _, files) in walkdir(root_path)
        for file in files
            _check_cancel(should_cancel)
            endswith(lowercase(file), ".csv") || continue
            total += 1
            _emit_progress(on_progress;
                phase=:counting,
                total_csv=total,
                processed_csv=total,
                loaded_measurements=0,
                skipped_csv=0,
                current_path=joinpath(root, file),
            )
        end
    end
    return total
end

"infer primitive types from string"
function _infer_meta_value(s::AbstractString)
    v = strip(s)
    isempty(v) && return nothing
    low = lowercase(v)
    if low in ("true", "t", "yes", "y", "1")
        return true
    end
    if low in ("false", "f", "no", "n", "0")
        return false
    end
    try
        return parse(Int, v)
    catch
    end
    try
        return parse(Float64, v)
    catch
    end
    # simple date/datetime patterns
    for fmt in (dateformat"yyyy-mm-dd", dateformat"yyyy/mm/dd")
        try
            return Date(v, fmt)
        catch
        end
    end
    for fmt in (dateformat"yyyy-mm-dd HH:MM:SS", dateformat"yyyy-mm-ddTHH:MM:SS")
        try
            return DateTime(v, fmt)
        catch
        end
    end
    return v
end

"""
    _load_device_info_txt(path) -> Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}

Reads device_info.txt where first column is `device_path` (slash-separated),
remaining columns are arbitrary device-level parameters.
"""
function _load_device_info_txt(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && return Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}()
    header = strip.(split(lines[1], ','))
    length(header) >= 2 || return Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}()
    devcol = header[1]
    param_cols = header[2:end]
    meta = Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}()
    for ln in lines[2:end]
        isempty(strip(ln)) && continue
        parts = split(ln, ',')
        length(parts) < 1 && continue
        raw_path = strip(parts[1])
        isempty(raw_path) && continue
        segs = Tuple(filter(!isempty, split(raw_path, '/')))
        params = Dict{Symbol,Any}()
        for (i, col) in enumerate(param_cols)
            idx = i + 1
            idx > length(parts) && continue
            cell = strip(parts[idx])
            val = _infer_meta_value(cell)
            val === nothing && continue
            params[Symbol(strip(col))] = val
        end
        meta[segs] = params
    end
    return meta
end

function _lookup_device_params(meta::Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}, loc::Vector{String})
    isempty(loc) && return nothing

    candidates = Tuple{Vararg{String}}[]
    push!(candidates, (last(loc),))           # geometry-only (e.g., L100W2)
    push!(candidates, (loc[1],))              # chip
    for k in reverse(1:length(loc)-1)
        push!(candidates, Tuple(loc[1:k]))    # chip/block...
    end
    push!(candidates, Tuple(loc))             # full path

    merged = Dict{Symbol,Any}()
    seen = Set{Tuple{Vararg{String}}}()
    for cand in candidates
        cand in seen && continue
        push!(seen, cand)
        if haskey(meta, cand)
            merge!(merged, meta[cand])
        end
    end

    return isempty(merged) ? nothing : merged
end

function _auto_detect_project(root_path::String)
    for (r, _, files) in walkdir(root_path)
        for file in files
            endswith(lowercase(file), ".csv") || continue
            for proj in KNOWN_PROJECTS
                accepts_file(proj, file) && return proj
            end
        end
    end
    return nothing
end

function _resolve_scan_project(root_path::String, project::Union{AbstractProject,Nothing})
    return if project !== nothing
        project
    else
        p = _auto_detect_project(root_path)
        p !== nothing ? p : _default_project[]
    end
end

function _load_scan_metadata(root_path::String)
    meta_path = joinpath(root_path, "device_info.txt")
    return isfile(meta_path) ? _load_device_info_txt(meta_path) : nothing
end

function _scan_accepted_csv!(measurements::Vector{MeasurementInfo}, proj::AbstractProject, filepath::String, meta)
    measurement_info = MeasurementInfo(filepath, proj)
    expanded = expand_measurement(proj, measurement_info)
    for m in expanded
        if meta !== nothing
            dev_params = _lookup_device_params(meta, m.device_info.location)
            if dev_params !== nothing
                merge!(m.device_info.parameters, dev_params)
            end
        end
        push!(measurements, m)
    end
end

function scan_directory(
    root_path::String;
    project::Union{AbstractProject,Nothing}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
    should_cancel::Union{Nothing,Function}=nothing,
    count_first::Bool=false,
)::MeasurementHierarchy
    proj = _resolve_scan_project(root_path, project)
    measurements = MeasurementInfo[]
    skipped_count = 0
    meta = _load_scan_metadata(root_path)
    processed_csv = 0
    total_csv = count_first ? _count_csv(root_path; should_cancel, on_progress) : 0
    _emit_progress(on_progress;
        phase=:scanning,
        total_csv=total_csv,
        processed_csv=processed_csv,
        loaded_measurements=length(measurements),
        skipped_csv=skipped_count,
    )

    for (root, dirs, files) in walkdir(root_path)
        for file in files
            _check_cancel(should_cancel)
            if endswith(lowercase(file), ".csv")
                filepath = joinpath(root, file)
                processed_csv += 1
                if !accepts_file(proj, file)
                    skipped_count += 1
                    _emit_progress(on_progress;
                        phase=:scanning,
                        total_csv=total_csv,
                        processed_csv=processed_csv,
                        loaded_measurements=length(measurements),
                        skipped_csv=skipped_count,
                        current_path=filepath,
                    )
                    continue
                end
                try
                    _scan_accepted_csv!(measurements, proj, filepath, meta)
                catch e
                    @warn "Could not parse measurement file $filepath" error = e
                end
                _emit_progress(on_progress;
                    phase=:scanning,
                    total_csv=total_csv,
                    processed_csv=processed_csv,
                    loaded_measurements=length(measurements),
                    skipped_csv=skipped_count,
                    current_path=filepath,
                )
                yield()
            end
        end
    end
    return MeasurementHierarchy(measurements, root_path, meta !== nothing, proj, skipped_count)
end

"""
Get statistics about a set of measurements.
"""
function get_measurements_stats(measurements::Vector{MeasurementInfo}, project::AbstractProject)
    stats = Dict{Symbol,Any}()
    stats[:total_measurements] = length(measurements)
    stats[:measurement_types] = unique([kind_label(project, m.measurement_kind) for m in measurements])
    timestamps = [m.timestamp for m in measurements if m.timestamp !== nothing]
    if !isempty(timestamps)
        stats[:first_measurement] = minimum(timestamps)
        stats[:last_measurement] = maximum(timestamps)
    end
    all_params = Dict{Symbol,Vector{Any}}()
    for measurement in measurements
        for (key, value) in measurement.parameters
            if !haskey(all_params, key)
                all_params[key] = Any[]
            end
            push!(all_params[key], value)
        end
    end
    stats[:parameter_ranges] = Dict{Symbol,Any}()
    for (param, values) in all_params
        if eltype(values) <: Number && !isempty(values)
            stats[:parameter_ranges][param] = (minimum(values), maximum(values))
        end
    end
    return stats
end
