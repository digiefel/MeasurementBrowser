"""
DeviceParser.jl - Parse device hierarchy from measurement filenames
"""

using Dates

# ---------------------------------------------------------------------------
# Constants / Regex patterns
# ---------------------------------------------------------------------------
const MAX_HEADER_LINES = 50
const REGEX_DEVICE = r"RuO2test_([A-Z0-9]+)_([A-Z0-9]+)_([A-Z0-9]+(?:W[0-9]+)?)"
const REGEX_DEVICE_NEW = r"^RuO2test_([^_]+)_([^_]+)_([^_]+)_([^_]+)_"

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
end

HierarchyNode(name::String, kind::Symbol) = HierarchyNode(name, kind, HierarchyNode[], MeasurementInfo[])

# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

function MeasurementInfo(filepath::AbstractString)
    filename = basename(filepath)
    device_info = parse_device_info(filename)
    measurement_kind = detect_measurement_kind(filename)
    timestamp = parse_timestamp(filename)
    parameters = parse_parameters(filename)
    file_info = extract_file_info(filepath)
    exp_label = measurement_label(measurement_kind)
    
    device_label = ""
    if (m = match(REGEX_DEVICE_NEW, filename)) !== nothing
        device_label = join(m.captures, "_")
    elseif (m = match(REGEX_DEVICE, filename)) !== nothing
        device_label = join(m.captures, "_")
    end

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

function parse_device_info(filename::String)
    if (m = match(REGEX_DEVICE_NEW, filename)) !== nothing
        caps = String[]
        for c in m.captures
            c === nothing && continue
            push!(caps, String(c))
        end
        loc = vcat("RuO2test_" * caps[1], caps[2:end])
        return DeviceInfo(loc)
    elseif (m = match(REGEX_DEVICE, filename)) !== nothing
        caps = String[]
        for c in m.captures
            c === nothing && continue
            push!(caps, String(c))
        end
        loc = vcat("RuO2test_" * caps[1], caps[2:end])
        return DeviceInfo(loc)
    end
    error("Unrecognized device filename format: $filename")
end

function detect_measurement_kind(filename::String)::Symbol
    lower = lowercase(filename)
    if occursin("fe pund", lower) || occursin("fepund", lower)
        return :pund
    elseif occursin("i_v sweep", lower) || occursin("iv sweep", lower)
        return :iv
    elseif occursin("tlm_4p", lower) || occursin("tlm", lower)
        return :tlm4p
    elseif occursin("break", lower) || occursin("breakdown", lower)
        return :breakdown
    elseif occursin("wakeup", lower)
        return :wakeup
    else
        return :unknown
    end
end

function measurement_label(kind::Symbol)::String
    kind === :pund && return "FE PUND"
    kind === :iv && return "I-V Sweep"
    kind === :tlm4p && return "TLM 4-Point"
    kind === :breakdown && return "Breakdown"
    kind === :wakeup && return "Wakeup"
    return "Unknown"
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

function display_label(meas::MeasurementInfo)
    label = measurement_label(meas.measurement_kind)
    
    temp_str = ""
    if haskey(meas.parameters, :temperature_K)
        temp_str = " $(meas.parameters[:temperature_K])K"
    end

    if meas.measurement_kind == :wakeup
        if meas.wakeup_pulse_count !== nothing && meas.wakeup_pulse_count > 0
            return "$(meas.timestamp) $(label) $(meas.wakeup_pulse_count)×$(temp_str)"
        end
    elseif meas.measurement_kind == :pund
        try
            amplitude_match = match(r"(\d+(?:\.\d+)?)V", meas.filename)
            if amplitude_match !== nothing
                voltage = parse(Float64, amplitude_match.captures[1])
                voltage_str = voltage == floor(voltage) ? "$(Int(voltage))V" : "$(voltage)V"
                return "$(meas.timestamp) $(label) $(voltage_str)$(temp_str)"
            end
        catch
            # ignore, fall through to default
        end
    end
    return "$(meas.timestamp) $(label)$(temp_str)"
end



# ---------------------------------------------------------------------------
# Multi-device (Breakdown) expansion
# ---------------------------------------------------------------------------
function expand_multi_device(meas::MeasurementInfo)::Vector{MeasurementInfo}
    meas.measurement_kind == :breakdown || return [meas]
    dev = last(meas.device_info.location)
        if (m = match(r"^([A-Z][0-9]+)([A-Z][0-9]+)$", dev)) === nothing
        return [meas]
    end
    parts = m.captures
    loc = copy(meas.device_info.location)
    return [MeasurementInfo(
        meas.filename,
        meas.filepath,
        replace(meas.clean_title, dev => p),
        meas.measurement_kind,
        meas.timestamp,
        DeviceInfo(vcat(loc[1:end-1], [p]), deepcopy(meas.device_info.parameters)),
        deepcopy(meas.parameters),
        meas.wakeup_pulse_count,
    ) for p in parts]
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
function MeasurementHierarchy(measurements::Vector{MeasurementInfo}, root_path::String, has_dev_metadata::Bool)
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
    mh = MeasurementHierarchy(root, measurements, root_path, index, has_dev_metadata)
    sort!(mh)
    return mh
end
# Backwards-compatible convenience constructor (assumes no metadata)
MeasurementHierarchy(measurements::Vector{MeasurementInfo}, root_path::String) =
    MeasurementHierarchy(measurements, root_path, false)

children(node::HierarchyNode) = node.children
isleaf(node::HierarchyNode) = isempty(node.children)

# ---------------------------------------------------------------------------
# Scanning
# ---------------------------------------------------------------------------
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
remaining columns are arbitrary device-level parameters. Keys may optionally
include the `RuO2test_` chip prefix (e.g. `RuO2test_A9/...`) to disambiguate
chip IDs like `A9` from other devices; both prefixed and unprefixed keys are
honored during lookup.

Example device_info.txt format for TLM measurements (including oxide and electrode thickness):
```
device_path,length_um,width_um,area_um2,t_HZO_nm,t_RuO2_nm,notes
TLML800W2,800,2,1600,10,30,TLM structure L800 W2
TLML800W4,800,4,3200,10,30,TLM structure L800 W4
TLML400W2,400,2,800,10,30,TLM structure L400 W2
TLML400W4,400,4,1600,10,30,TLM structure L400 W4
```

For TLM combined analysis, the length_um and width_um parameters are required
to calculate width-normalized resistance and extract sheet resistance.
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

# ---------------------------------------------------------------------------
# Scanning
# ---------------------------------------------------------------------------
function _lookup_device_params(meta::Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}, loc::Vector{String})
    isempty(loc) && return nothing

    candidates = Tuple{Vararg{String}}[]
    push!(candidates, (last(loc),))           # geometry-only (e.g., L100W2)
    push!(candidates, (loc[1],))              # chip (includes prefix if present)
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

function scan_directory(root_path::String)::MeasurementHierarchy
    measurements = MeasurementInfo[]

    # load device_info.txt at root level if present
    meta_path = joinpath(root_path, "device_info.txt")
    meta = isfile(meta_path) ? _load_device_info_txt(meta_path) : nothing

    for (root, dirs, files) in walkdir(root_path)
        for file in files
            if endswith(lowercase(file), ".csv")
                filepath = joinpath(root, file)
                try
                    measurement_info = MeasurementInfo(filepath)
                    # expand (may duplicate)
                    for m in expand_multi_device(measurement_info)
                        if meta !== nothing
                            dev_params = _lookup_device_params(meta, m.device_info.location)
                            if dev_params !== nothing
                                merge!(m.device_info.parameters, dev_params)
                            end
                        end
                        push!(measurements, m)
                    end
                catch e
                    @warn "Could not parse measurement file $filepath" error = e
                end
            end
        end
    end
    return MeasurementHierarchy(measurements, root_path, meta !== nothing)
end

"""
Get statistics about a set of measurements.
"""
function get_measurements_stats(measurements::Vector{MeasurementInfo})
    stats = Dict{Symbol,Any}()
    stats[:total_measurements] = length(measurements)
    stats[:measurement_types] = unique([measurement_label(m.measurement_kind) for m in measurements])
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
