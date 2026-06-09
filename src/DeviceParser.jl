"""
DeviceParser.jl - Parse device hierarchy from measurement filenames
"""

# TODO: invert ownership — `DeviceInfo` is currently nested inside every
# `MeasurementInfo`, but devices own measurements, not the other way around.
# `DeviceInfo` should live on the hierarchy leaf; `MeasurementInfo` should
# carry only the path/leaf reference. This will simplify per-device state
# (parameters, tags, notes, coords) and remove redundant `DeviceInfo` copies.

using Dates
using GLMakie: Figure

# ---------------------------------------------------------------------------
# Project dispatch types
# ---------------------------------------------------------------------------

abstract type AbstractProject end
struct RuO2Project <: AbstractProject end
struct TASEProject <: AbstractProject end

abstract type PlotKind end

const KNOWN_PROJECTS = AbstractProject[]
const _default_project = Ref{Union{AbstractProject,Nothing}}(nothing)

struct MeasurementAnalysisFailure
    filepath::String
    measurement_id::String
    message::String
end

# Interface (implemented by each project via multiple dispatch):
#   parse_device_info(::P, file::SourceFile) → DeviceInfo
#   detect_kind(::P, filename) → Symbol
#   interpret_file(::P, file) → Vector{MeasurementInfo}
#   kind_label(::P, kind) → String
#   display_label(::P, meas) → String
#   load_source_data(::P, source_file; measurement=nothing) → DataFrame
#   process_measurement_data(::P, measurement, data) → DataFrame
#   setup_plot(::P, plot_kind, measurements) → Figure
#   plot_data!(::P, plot_kind, measurements, figure) → nothing
#   debug_plot(::P, measurements, loaded; kwargs...) → Union{Figure,Nothing}
#   available_analyses(::P, measurements) → Vector{NamedTuple}
#   run_analysis(::P, key, measurements; kwargs...) → Union{AnalysisResult,Nothing}
#   draw_analysis_view(result, view) → Union{Figure,Nothing}
#   compute_and_add_measurement_stats!(::P, measurements, files) → Vector{MeasurementAnalysisFailure}

available_analyses(::AbstractProject, measurements) = NamedTuple[]
run_analysis(::AbstractProject, key::Symbol, measurements; kwargs...) = nothing
draw_analysis_view(result::AnalysisResult, view::NamedTuple) = nothing
interpret_file(::AbstractProject, file::SourceFile) = MeasurementInfo[]

"""Compute per-measurement stats after the complete measurement list and source headers are known."""
compute_and_add_measurement_stats!(
    ::AbstractProject,
    measurements::Vector,
    files::Vector{SourceFile},
    ;
    on_progress::Union{Nothing,Function}=nothing,
)::Vector{MeasurementAnalysisFailure} = MeasurementAnalysisFailure[]

struct JobCancelled <: Exception end

const _CANCEL_CALLBACK_KEY = :MeasurementBrowser_cancel_requested

function _with_cancel(f::Function, cancel_requested::Union{Nothing,Function})
    cancel_requested === nothing && return f()
    return task_local_storage(f, _CANCEL_CALLBACK_KEY, cancel_requested)
end

function _check_cancel()
    storage = task_local_storage()
    cancel_requested = get(storage, _CANCEL_CALLBACK_KEY, nothing)
    cancel_requested !== nothing && cancel_requested() && throw(JobCancelled())
    return nothing
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
"""
MeasurementInfo is the core struct representing a single measurement in the browser. 
It may be derived from one or more source files, and multiple MeasurementInfo entries 
may be derived from the same file (e.g. one per device or cycle).

`clean_title` is the short title shown in the GUI.
`device_info` carries metadata about the device under test.

`parameters` are measurement/acquisition "settings" parsed from filename/header.
`stats` are computed from file contents or measurement history later on.
The same key may appear in both; the container is part of the meaning.
"""
struct MeasurementInfo
    unique_id::String
    filename::String
    filepath::String
    clean_title::String
    measurement_kind::Symbol
    timestamp::Union{DateTime,Nothing}
    device_info::DeviceInfo
    parameters::Dict{Symbol,Any}
    stats::Dict{Symbol,Any}
end

function MeasurementInfo(;
    filepath,
    measurement_kind,
    device_info,
    clean_title,
    filename=basename(String(filepath)),
    unique_id=filepath,
    timestamp=nothing,
    parameters=Dict{Symbol,Any}(),
    stats=Dict{Symbol,Any}(),
)
    filepath === nothing && error("MeasurementInfo requires filepath")
    filename === nothing && error("MeasurementInfo requires filename")
    unique_id === nothing && error("MeasurementInfo requires unique_id")
    clean_title === nothing && error("MeasurementInfo requires clean_title")
    measurement_kind === nothing && error("MeasurementInfo requires measurement_kind")
    device_info === nothing && error("MeasurementInfo requires device_info")
    return MeasurementInfo(
        String(unique_id),
        String(filename),
        String(filepath),
        String(clean_title),
        measurement_kind,
        timestamp,
        device_info,
        parameters,
        stats,
    )
end

function MeasurementInfo(
    measurement::MeasurementInfo;
    filepath=measurement.filepath,
    filename=measurement.filename,
    unique_id=measurement.unique_id,
    clean_title=measurement.clean_title,
    measurement_kind=measurement.measurement_kind,
    timestamp=measurement.timestamp,
    device_info=DeviceInfo(
        copy(measurement.device_info.location),
        deepcopy(measurement.device_info.parameters),
    ),
    parameters=deepcopy(measurement.parameters),
    stats=deepcopy(measurement.stats),
)
    return MeasurementInfo(;
        filepath,
        filename,
        unique_id,
        clean_title,
        measurement_kind,
        timestamp,
        device_info,
        parameters,
        stats,
    )
end

setup_plot(
    project::AbstractProject,
    kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
)::Figure =
    error("No setup_plot implementation for $(project_name(project)) plot kind '$kind'")

plot_data!(
    project::AbstractProject,
    kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
    figure::Figure,
)::Nothing =
    error("No plot_data! implementation for $(project_name(project)) plot kind '$kind'")

function _append_plot_kinds!(kinds::Vector{Type{<:PlotKind}}, parent::Type)
    for name in names(@__MODULE__; all=true)
        isdefined(@__MODULE__, name) || continue
        child = getfield(@__MODULE__, name)
        child isa DataType || continue
        supertype(child) === parent || continue
        if isabstracttype(child)
            _append_plot_kinds!(kinds, child)
        elseif child <: PlotKind
            push!(kinds, child)
        end
    end
    return kinds
end

function plot_kinds()::Vector{Type{<:PlotKind}}
    kinds = _append_plot_kinds!(Type{<:PlotKind}[], PlotKind)
    sort!(kinds; by=kind -> String(nameof(kind)))
    return kinds
end

function _plot_job_data(
    project::AbstractProject,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo};
    debug::Bool=false,
)
    return nothing
end

function _plot_job_figure(
    project::AbstractProject,
    plot_kind::Type{<:PlotKind},
    measurements::Vector{MeasurementInfo},
    data;
    debug::Bool=false,
    device_params::Vector{Dict{Symbol,Any}}=Dict{Symbol,Any}[],
)
    if debug
        return debug_plot(project, measurements, data; device_params, plot_kind)
    end
    fig = setup_plot(project, plot_kind, measurements)
    plot_data!(project, plot_kind, measurements, fig)
    return fig
end

debug_plot(
    ::AbstractProject,
    measurements::Vector{MeasurementInfo},
    loaded;
    kwargs...,
) = error("Debug plots are not implemented for this project")

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
    analysis_failures::Vector{MeasurementAnalysisFailure}
end

function SourceScan(
    root_path::String,
    project::AbstractProject,
    files::Vector{SourceFile},
    hierarchy::MeasurementHierarchy,
)
    return SourceScan(root_path, project, files, hierarchy, MeasurementAnalysisFailure[])
end

# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

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
"""
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

function _count_csv(root_path::String; on_progress::Union{Nothing,Function}=nothing)
    total = 0
    for (root, _, files) in walkdir(root_path)
        for file in files
            _check_cancel()
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

const DEVICE_INFO_FILENAME = "device_info.txt"

_device_info_path(root_path::AbstractString) = joinpath(root_path, DEVICE_INFO_FILENAME)
_has_device_metadata(root_path::AbstractString)::Bool = isfile(_device_info_path(root_path))

function _load_scan_metadata(root_path::String)
    meta_path = _device_info_path(root_path)
    return isfile(meta_path) ? _load_device_info_txt(meta_path) : nothing
end

"""
Create the short title shown for a measurement in panels and plots.
The title combines project label and device label, with the filename as fallback.
"""
function build_clean_title(
    project::AbstractProject,
    filename::String,
    measurement_kind::Symbol,
    device_info::DeviceInfo,
    header_summary::Dict{String,String},
)
    exp_label = kind_label(project, measurement_kind)
    device_label = device_path_label(project, device_info)
    parts = filter(!isempty, (exp_label == "Unknown" ? "" : exp_label, device_label))
    return isempty(parts) ? strip(replace(filename, r"\.csv$" => "")) : join(parts, " ")
end

"""
Interpret a file into project measurements.
The project parser handles filename/header semantics; this wrapper overlays device metadata loaded
from `device_info.txt` when present.
"""
function interpret_measurements(
    project::AbstractProject,
    file::SourceFile,
    meta::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}},
)
    measurements = interpret_file(project, file)
    if meta !== nothing
        for measurement in measurements
            dev_params = _lookup_device_params(meta, measurement.device_info.location)
            dev_params !== nothing && merge!(measurement.device_info.parameters, dev_params)
        end
    end
    return measurements
end

"""
Interpret a single file into all its expanded `MeasurementInfo` entries.
Equivalent to the scan pipeline for one file; useful in tests and single-file code paths.
"""
function measurements_for_file(
    project::AbstractProject,
    filepath::AbstractString;
    meta::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}}=nothing,
)
    file = index_source_file(filepath)
    measurements = interpret_measurements(project, file, meta)
    compute_and_add_measurement_stats!(project, measurements, [file])
    return measurements
end

"""
Interpret files concurrently.
`on_result` receives each interpreted `SourceFile` with its original file index so callers can
assemble the final scan in stable order. `on_measurements` receives interpreted measurements as
soon as each file finishes, allowing the UI to update while the complete scan is still running.
Progress counts files attempted, measurements loaded, and files that produced no visible
measurements.
"""
function _interpret_source_files(
    project::AbstractProject,
    files::Vector{SourceFile},
    meta::Union{Nothing,Dict{Tuple{Vararg{String}},Dict{Symbol,Any}}};
    on_result::Function,
    on_measurements::Union{Nothing,Function}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
)
    processed_csv = Base.Threads.Atomic{Int}(0)
    loaded_measurements = Base.Threads.Atomic{Int}(0)
    skipped_csv = Base.Threads.Atomic{Int}(0)
    progress_lock = ReentrantLock()
    worker_limit = Base.Semaphore(max(1, Base.Threads.nthreads()))
    cancel_requested = get(task_local_storage(), _CANCEL_CALLBACK_KEY, nothing)

    @sync for (index, file) in pairs(files)
        Base.acquire(worker_limit)
        Base.Threads.@spawn begin
            try
                measurements = _with_cancel(cancel_requested) do
                    interpret_measurements(project, file, meta)
                end
                scanned = SourceFile(file, measurements)
                on_result(index, scanned)
                if on_measurements !== nothing && !isempty(measurements)
                    lock(progress_lock) do
                        on_measurements(measurements)
                    end
                end
                isempty(measurements) && Base.Threads.atomic_add!(skipped_csv, 1)
                isempty(measurements) || Base.Threads.atomic_add!(loaded_measurements, length(measurements))
                processed_now = Base.Threads.atomic_add!(processed_csv, 1) + 1

                if on_progress !== nothing
                    progress = (
                        total_csv=length(files),
                        processed_csv=processed_now,
                        loaded_measurements=loaded_measurements[],
                        skipped_csv=skipped_csv[],
                        current_path=file.filepath,
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
        total_csv=length(files),
    )
end

"""
Scan a project directory into the measurement hierarchy used by the browser.
The scan indexes CSV files, interprets them with the selected project, applies device metadata,
streams interpreted measurements when requested, computes project stats, and builds the final tree.
"""
function scan_source(
    root_path::String;
    project::Union{AbstractProject,Nothing}=nothing,
    on_progress::Union{Nothing,Function}=nothing,
    on_measurements::Union{Nothing,Function}=nothing,
    count_first::Bool=false,
)::SourceScan
    proj = _resolve_scan_project(project)
    root = normpath(abspath(expanduser(String(root_path))))
    meta = _load_scan_metadata(root)
    count_first && _count_csv(root; on_progress)
    _emit_progress(on_progress;
        phase=:discovering,
        total_csv=0,
        processed_csv=0,
        loaded_measurements=0,
        skipped_csv=0,
    )
    files = collect_source_files(
        root;
        on_file=(file, count) -> _emit_progress(
            on_progress;
            phase=:discovering,
            total_csv=0,
            processed_csv=count,
            loaded_measurements=0,
            skipped_csv=0,
            current_path=file.filepath,
        ),
    )
    scanned_files = Vector{SourceFile}(undef, length(files))
    _emit_progress(on_progress;
        phase=:scanning,
        total_csv=length(files),
        processed_csv=0,
        loaded_measurements=0,
        skipped_csv=0,
    )

    summary = _interpret_source_files(
        proj,
        files,
        meta;
        on_result=(index, file) -> (scanned_files[index] = file),
        on_measurements=on_measurements,
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
    for file in scanned_files
        append!(measurements, file.measurements)
    end

    analysis_progress_seen = Ref(false)
    _emit_progress(on_progress;
        phase=:analyzing,
        total_csv=0,
        processed_csv=0,
        loaded_measurements=length(measurements),
        skipped_csv=summary.skipped_csv,
        current_path=root,
    )
    analysis_failures = compute_and_add_measurement_stats!(
        proj,
        measurements,
        scanned_files;
        on_progress=(progress) -> begin
            analysis_progress_seen[] = true
            _emit_progress(on_progress;
                phase=:analyzing,
                total_csv=progress.total,
                processed_csv=progress.processed,
                loaded_measurements=length(measurements),
                skipped_csv=summary.skipped_csv,
                current_path=progress.current_path,
            )
        end,
    )
    analysis_progress_seen[] || _emit_progress(on_progress;
        phase=:analyzing,
        total_csv=1,
        processed_csv=1,
        loaded_measurements=length(measurements),
        skipped_csv=summary.skipped_csv,
        current_path=root,
    )
    hierarchy = MeasurementHierarchy(measurements, root, meta !== nothing, proj, summary.skipped_csv)
    return SourceScan(root, proj, scanned_files, hierarchy, analysis_failures)
end

"""Collect generic device-summary fields for the information panel."""
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
