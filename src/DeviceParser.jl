"""
DeviceParser.jl - Parse device hierarchy from measurement filenames
"""

# TODO: invert ownership — `DeviceInfo` is currently nested inside every
# `MeasurementInfo`, but devices own measurements, not the other way around.
# `DeviceInfo` should live on the hierarchy leaf; `MeasurementInfo` should
# carry only the path/leaf reference. This will simplify per-device state
# (parameters, tags, notes, coords) and remove redundant `DeviceInfo` copies.

using Dates

# ---------------------------------------------------------------------------
# Project dispatch types
# ---------------------------------------------------------------------------

abstract type AbstractProject end
struct RuO2Project <: AbstractProject end
struct TASEProject <: AbstractProject end

const KNOWN_PROJECTS = AbstractProject[]
const _default_project = Ref{Union{AbstractProject,Nothing}}(nothing)

# Interface (implemented by each project via multiple dispatch):
#   parse_device_info(::P, source::SourceFile) → DeviceInfo
#   detect_kind(::P, filename) → Symbol
#   interpret_file(::P, source) → Vector{MeasurementItem}
#   kind_label(::P, kind) → String
#   display_label(::P, meas) → String
#   expand_measurement(::P, meas) → Vector{MeasurementInfo}
#   load_plot_for_file(::P, path, kind; kwargs...) → Any
#   analyze_plot_for_file(::P, kind, loaded; kwargs...) → Any
#   draw_plot_for_file(::P, kind, analyzed; kwargs...) → Union{Figure,Nothing}
#   load_plot_for_files(::P, paths, combined_kind; kwargs...) → Any
#   analyze_plot_for_files(::P, combined_kind, loaded; kwargs...) → Any
#   draw_plot_for_files(::P, combined_kind, analyzed; kwargs...) → Union{Figure,Nothing}
#   available_analyses(::P, measurements) → Vector{NamedTuple}
#   run_analysis(::P, key, measurements; kwargs...) → Union{AnalysisResult,Nothing}
#   draw_analysis_view(result, view) → Union{Figure,Nothing}
#   combined_plot_types(::P) → Vector{Tuple}
#   compatible_kinds(::P, combined_kind) → Vector{Symbol}
#   annotate_measurements!(::P, measurements) → nothing
#   augment_measurements_stats!(::P, stats, measurements) → stats

load_plot_for_file(::AbstractProject, path::AbstractString, kind::Union{Symbol,Nothing}; kwargs...) = nothing
analyze_plot_for_file(::AbstractProject, kind::Union{Symbol,Nothing}, loaded; kwargs...) = loaded
draw_plot_for_file(::AbstractProject, kind::Union{Symbol,Nothing}, analyzed; kwargs...) = nothing
load_plot_for_files(::AbstractProject, paths, combined_kind; kwargs...) = nothing
analyze_plot_for_files(::AbstractProject, combined_kind, loaded; kwargs...) = loaded
draw_plot_for_files(::AbstractProject, combined_kind, analyzed; kwargs...) = nothing
available_analyses(::AbstractProject, measurements) = NamedTuple[]
run_analysis(::AbstractProject, key::Symbol, measurements; kwargs...) = nothing
draw_analysis_view(result::AnalysisResult, view::NamedTuple) = nothing
interpret_file(::AbstractProject, source::SourceFile; kwargs...) = MeasurementItem[]

"""
Attach project-specific derived metadata to measurements after the complete measurement list is known.
Use this for values that depend on device-level measurement history and should reach plots or cache metadata.
"""
annotate_measurements!(::AbstractProject, measurements::Vector) = nothing

"""
Add project-specific aggregate values to the stats shown for a selected device.
Use this for information-panel summaries computed from the selected device's measurements.
"""
augment_measurements_stats!(::AbstractProject, stats::Dict{Symbol,Any}, measurements::Vector) = stats

struct PlotCancelled <: Exception end

function _check_plot_cancel(should_cancel::Union{Nothing,Function})
    should_cancel !== nothing && should_cancel() && throw(PlotCancelled())
end

# ---------------------------------------------------------------------------
# DeviceInfo
# ---------------------------------------------------------------------------

struct DeviceInfo
    location::Vector{String}              # variable-length hierarchy
    parameters::Dict{Symbol,Any}          # device-level metadata
end

DeviceInfo(location::Vector{String}) = DeviceInfo(location, Dict{Symbol,Any}())
device_path_label(::AbstractProject, device_info::DeviceInfo) = join(device_info.location, "_")
device_path_key(location::AbstractVector{<:AbstractString}) = join(location, "/")
device_path_key(device_info::DeviceInfo) = device_path_key(device_info.location)

"""
Parse a stored device path key back into the tuple form used by the hierarchy index.
"""
function device_path_tuple(key::AbstractString)
    stripped = strip(String(key))
    isempty(stripped) && error("Device path key cannot be empty")
    segs = split(stripped, '/')
    any(isempty, segs) && error("Invalid device path key '$key'")
    return Tuple(String.(segs))
end

# ---------------------------------------------------------------------------
# Measurement related structs
# ---------------------------------------------------------------------------
struct MeasurementInfo
    id::String
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
    skipped_count::Int          # CSV files not interpreted into visible measurements in last scan
end

HierarchyNode(name::String, kind::Symbol) = HierarchyNode(name, kind, HierarchyNode[], MeasurementInfo[])

struct SourceScan
    root_path::String
    project::AbstractProject
    files::Vector{SourceFile}
    hierarchy::MeasurementHierarchy
end

# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

"""
Construct `MeasurementInfo` directly from one source file.
This is mainly for single-file code paths; full scans use `interpret_measurements` first.
"""
function MeasurementInfo(filepath::AbstractString, project::AbstractProject)
    indexed = index_source_file(filepath)
    filename = indexed.filename
    device_info = parse_device_info(project, indexed)
    measurement_kind = detect_kind(project, filename)
    timestamp = indexed.timestamp
    parameters = parse_parameters(filename)
    exp_label = kind_label(project, measurement_kind)

    device_label = device_path_label(project, device_info)

    date_str = ""
    d = get(indexed.header_summary, :test_date, "")
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

    return MeasurementInfo(indexed.id, filename, filepath, clean_title, measurement_kind, timestamp, device_info, parameters, wakeup_count)
end

"""
Extract a measurement timestamp from supported filename conventions.
Returns `nothing` when the filename does not carry a parseable timestamp.
"""
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

"""
Extract lightweight measurement parameters encoded in the filename.
These are generic hints used before project-specific parsing or device metadata is applied.
"""
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

"""
Return the chronological sort key for measurement lists, placing missing timestamps last.
"""
measurement_timestamp_key(m::MeasurementInfo) =
    m.timestamp === nothing ? DateTime(Dates.year(typemax(Date))) : m.timestamp

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
        sort!(ch.measurements, by=measurement_timestamp_key)
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
"""
Build the device tree and lookup index from a flat measurement list.
The constructor also runs the project annotation hook, so measurement parameters may be finalized here.
"""
function MeasurementHierarchy(measurements::Vector{MeasurementInfo}, root_path::String, has_dev_metadata::Bool, project::AbstractProject, skipped_count::Int=0)
    annotate_measurements!(project, measurements)
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

function MeasurementHierarchy(root_path::String, has_dev_metadata::Bool, project::AbstractProject, skipped_count::Int=0)
    return MeasurementHierarchy(
        HierarchyNode("/", :root),
        MeasurementInfo[],
        root_path,
        Dict{Tuple{Vararg{String}},HierarchyNode}(),
        has_dev_metadata,
        project,
        skipped_count,
    )
end

function _ensure_hierarchy_child!(
    parent::HierarchyNode,
    index::Dict{Tuple{Vararg{String}},HierarchyNode},
    name::String,
    kind::Symbol,
    path_tuple::Tuple{Vararg{String}},
)
    existing = get(index, path_tuple, nothing)
    if existing !== nothing
        return existing
    end
    node = HierarchyNode(name, kind)
    push!(parent.children, node)
    index[path_tuple] = node
    return node
end

function insert_measurement!(mh::MeasurementHierarchy, measurement::MeasurementInfo)
    parent = mh.root
    for (i, seg) in enumerate(measurement.device_info.location)
        kind = i == length(measurement.device_info.location) ? :leaf : :level
        path_tuple = Tuple(measurement.device_info.location[1:i])
        parent = _ensure_hierarchy_child!(parent, mh.index, seg, kind, path_tuple)
    end
    push!(parent.measurements, measurement)
    push!(mh.all_measurements, measurement)
    return measurement
end

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

Reads device_info.txt where the first column contains an exact unit path
(`device`, `device_path`, etc.; slash-separated for multi-unit scopes) and the
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

"""
Find device metadata entries that apply to a parsed device location.
Exact path fragments can match at any level, with more specific matches overriding broader ones.
"""
function _lookup_device_params(meta::Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}, loc::Vector{String})
    isempty(loc) && return nothing

    merged = Dict{Symbol,Any}()
    # Match exact slash-separated unit sequences anywhere in the parsed path.
    # More specific matches win because longer sequences are merged later.
    for width in 1:length(loc)
        for start_idx in 1:(length(loc) - width + 1)
            candidate = Tuple(loc[start_idx:start_idx + width - 1])
            if haskey(meta, candidate)
                merge!(merged, meta[candidate])
            end
        end
    end

    return isempty(merged) ? nothing : merged
end

_resolve_scan_project(project::Union{AbstractProject,Nothing}) = project !== nothing ? project : _default_project[]

function _load_scan_metadata(root_path::String)
    meta_path = joinpath(root_path, "device_info.txt")
    return isfile(meta_path) ? _load_device_info_txt(meta_path) : nothing
end

"""
Create the short title shown for a measurement in panels and plots.
The title combines project label, device label, and header date when available, with the filename as fallback.
"""
function build_clean_title(
    project::AbstractProject,
    filename::String,
    measurement_kind::Symbol,
    device_info::DeviceInfo,
    header_summary::Dict{Symbol,Any},
)
    exp_label = kind_label(project, measurement_kind)
    device_label = device_path_label(project, device_info)
    date_str = ""
    d = get(header_summary, :test_date, "")
    if d isa AbstractString && !isempty(d)
        date_str = try
            parts = split(d)
            if length(parts) >= 3
                month = parts[2]
                day = try
                    parse(Int, parts[3])
                catch
                    nothing
                end
                day !== nothing ? "$(month)$(day)" : d
            else
                d
            end
        catch
            d
        end
    end
    parts = filter(!isempty, (exp_label == "Unknown" ? "" : exp_label, device_label, date_str))
    return isempty(parts) ? strip(replace(filename, r"\.csv$" => "")) : join(parts, " ")
end

function _measurement_info_from_item(item::MeasurementItem)
    wakeup_count = get(item.parameters, :wakeup_pulse_count, nothing)
    wakeup_count isa Int || (wakeup_count = nothing)
    return MeasurementInfo(
        item.id,
        basename(item.filepath),
        item.filepath,
        item.title,
        item.kind,
        item.timestamp,
        DeviceInfo(copy(item.device_path), deepcopy(item.device_parameters)),
        deepcopy(item.parameters),
        wakeup_count,
    )
end

"""
Interpret a source CSV into project measurement items.
The project parser handles filename/header semantics; this wrapper overlays device metadata loaded
from `device_info.txt` when present. Source scans and cache builds both use this path before
constructing `MeasurementInfo`.
"""
function interpret_measurements(
    project::AbstractProject,
    source::SourceFile,
    meta::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}};
    should_cancel::Union{Nothing,Function}=nothing,
)
    items = interpret_file(project, source; should_cancel=should_cancel)
    if meta !== nothing
        for item in items
            dev_params = _lookup_device_params(meta, item.device_path)
            dev_params !== nothing && merge!(item.device_parameters, dev_params)
        end
    end
    return items
end

function _is_scan_cancel_error(err)
    err isa ScanCancelled && return true
    err isa CompositeException || return false
    return any(_is_scan_cancel_error, err.exceptions)
end

"""
Interpret source CSV files concurrently while preserving the original file order in callbacks.
Progress counts files attempted, measurements loaded, and files that produced no visible measurements.
"""
function _interpret_source_files(
    project::AbstractProject,
    source_files::Vector{SourceFile},
    meta::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}};
    should_cancel::Union{Nothing,Function}=nothing,
    on_result::Function,
    on_progress::Union{Nothing,Function}=nothing,
)
    processed_csv = Base.Threads.Atomic{Int}(0)
    loaded_measurements = Base.Threads.Atomic{Int}(0)
    skipped_csv = Base.Threads.Atomic{Int}(0)
    progress_lock = ReentrantLock()
    worker_limit = Base.Semaphore(max(1, Base.Threads.nthreads()))

    @sync for (index, source) in pairs(source_files)
        Base.acquire(worker_limit)
        Base.Threads.@spawn begin
            try
                items = interpret_measurements(project, source, meta; should_cancel=should_cancel)
                measurements = [_measurement_info_from_item(item) for item in items]
                scanned = source_file_with_measurements(source, measurements)
                on_result(index, scanned)
                isempty(measurements) && Base.Threads.atomic_add!(skipped_csv, 1)
                isempty(measurements) || Base.Threads.atomic_add!(loaded_measurements, length(measurements))
                processed_now = Base.Threads.atomic_add!(processed_csv, 1) + 1

                if on_progress !== nothing
                    progress = (
                        total_csv=length(source_files),
                        processed_csv=processed_now,
                        loaded_measurements=loaded_measurements[],
                        skipped_csv=skipped_csv[],
                        current_path=source.filepath,
                    )
                    lock(progress_lock) do
                        on_progress(progress)
                    end
                end
            finally
                Base.release(worker_limit)
            end
        end
    end

    return (
        processed_csv=processed_csv[],
        loaded_measurements=loaded_measurements[],
        skipped_csv=skipped_csv[],
        total_csv=length(source_files),
    )
end

"""
Scan a project directory into the measurement hierarchy used by the browser.
The scan indexes CSV files, interprets them with the selected project, applies device metadata, and builds the tree.
"""
function scan_source(
    root_path::String;
    project::Union{AbstractProject,Nothing}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
    should_cancel::Union{Nothing,Function}=nothing,
    count_first::Bool=false,
)::SourceScan
    proj = _resolve_scan_project(project)
    root = source_path(root_path)
    meta = _load_scan_metadata(root)
    count_first && _count_csv(root; should_cancel, on_progress)
    source_files = collect_source_files(root; should_cancel=should_cancel)
    scanned_files = Vector{SourceFile}(undef, length(source_files))
    _emit_progress(on_progress;
        phase=:scanning,
        total_csv=length(source_files),
        processed_csv=0,
        loaded_measurements=0,
        skipped_csv=0,
    )

    summary = _interpret_source_files(
        proj,
        source_files,
        meta;
        should_cancel=should_cancel,
        on_result=(index, source) -> (scanned_files[index] = source),
        on_progress=(progress) -> _emit_progress(
            on_progress;
            phase=:scanning,
            total_csv=progress.total_csv,
            processed_csv=progress.processed_csv,
            loaded_measurements=progress.loaded_measurements,
            skipped_csv=progress.skipped_csv,
            current_path=progress.current_path,
        ),
    )

    measurements = MeasurementInfo[]
    sizehint!(measurements, summary.loaded_measurements)
    for source in scanned_files
        append!(measurements, source.measurements)
    end

    hierarchy = MeasurementHierarchy(measurements, root, meta !== nothing, proj, summary.skipped_csv)
    return SourceScan(root, proj, scanned_files, hierarchy)
end

"""
Collect generic device-summary fields for the information panel.
Projects can extend the returned dictionary with measurements-derived stats.
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
    return augment_measurements_stats!(project, stats, measurements)
end
