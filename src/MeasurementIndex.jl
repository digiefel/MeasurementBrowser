module MeasurementIndex

using Dates

import ..Projects:
    AbstractProject,
    compute_and_add_measurement_stats!,
    device_path_label,
    interpret_file,
    kind_label,
    project_name,
    reset_scan_profile!

"""
Failure produced while interpreting or analyzing one physical source file.
"""
struct MeasurementAnalysisFailure
    filepath::String
    measurement_id::String
    message::String
end

"""
Cancellation raised by package-owned background work.
"""
struct JobCancelled <: Exception end

const CANCEL_CALLBACK_KEY = :MeasurementBrowser_cancel_requested

"""
Return whether an exception represents cancellation of all contained work.
"""
function is_job_cancelled(error::Exception)::Bool
    error isa JobCancelled && return true
    error isa CompositeException || return false
    return !isempty(error.exceptions) && all(is_job_cancelled, error.exceptions)
end

"""
Run work with a task-local cancellation callback.
"""
function with_cancel(
    work::Function,
    cancel_requested::Union{Nothing,Function},
)::Any
    cancel_requested === nothing && return work()
    return task_local_storage(work, CANCEL_CALLBACK_KEY, cancel_requested)
end

"""
Stop the current package job when its task-local cancellation callback requests it.
"""
function check_cancel()::Nothing
    cancel_requested = get(task_local_storage(), CANCEL_CALLBACK_KEY, nothing)
    cancel_requested !== nothing && cancel_requested() && throw(JobCancelled())
    return nothing
end

"""
Device identity and metadata shared by measurements from the same device path.
"""
struct DeviceInfo
    location::Vector{String}
    parameters::Dict{Symbol,Any}
end

DeviceInfo(location::Vector{String}) = DeviceInfo(location, Dict{Symbol,Any}())

device_path_label(::AbstractProject, device_info::DeviceInfo)::String =
    join(device_info.location, "_")

device_path_key(location::AbstractVector{<:AbstractString})::String = join(location, "/")
device_path_key(device_info::DeviceInfo)::String = device_path_key(device_info.location)

"""
Parse a stored slash-separated device key into the tuple used by the hierarchy index.
"""
function device_path_tuple(key::AbstractString)::Tuple{Vararg{String}}
    stripped = strip(String(key))
    isempty(stripped) && error("Device path key cannot be empty")
    segments = split(stripped, '/')
    any(isempty, segments) && error("Invalid device path key '$key'")
    return Tuple(String.(segments))
end

"""
One logical measurement discovered inside a physical source file.

`parameters` describe acquisition settings known while interpreting the file. `stats` contain values
computed after the required measurement context is available.
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

"""
Construct a logical measurement while normalizing its string fields.
"""
function MeasurementInfo(;
    filepath::AbstractString,
    measurement_kind::Symbol,
    device_info::DeviceInfo,
    clean_title::AbstractString,
    filename::AbstractString=basename(filepath),
    unique_id::AbstractString=filepath,
    timestamp::Union{DateTime,Nothing}=nothing,
    parameters::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    stats::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)::MeasurementInfo
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

"""
Copy a measurement while replacing selected fields.
"""
function MeasurementInfo(
    measurement::MeasurementInfo;
    filepath::AbstractString=measurement.filepath,
    filename::AbstractString=measurement.filename,
    unique_id::AbstractString=measurement.unique_id,
    clean_title::AbstractString=measurement.clean_title,
    measurement_kind::Symbol=measurement.measurement_kind,
    timestamp::Union{DateTime,Nothing}=measurement.timestamp,
    device_info::DeviceInfo=DeviceInfo(
        copy(measurement.device_info.location),
        deepcopy(measurement.device_info.parameters),
    ),
    parameters::Dict{Symbol,Any}=deepcopy(measurement.parameters),
    stats::Dict{Symbol,Any}=deepcopy(measurement.stats),
)::MeasurementInfo
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

"""
Extract a measurement timestamp from the supported filename conventions.
"""
function parse_timestamp(filename::AbstractString)::Union{DateTime,Nothing}
    if (result = match(
        r"; (\d{4}-\d{2}-\d{2}) (\d{2})_(\d{2})_(\d{2})\]",
        filename,
    )) !== nothing
        date, hour, minute, second = result.captures
        return try
            DateTime("$date $hour:$minute:$second", "yyyy-mm-dd HH:MM:SS")
        catch
            nothing
        end
    end
    result = match(r"_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})_", filename)
    result === nothing && return nothing
    year, month, day, hour, minute, second = result.captures
    return try
        DateTime("$year-$month-$day $hour:$minute:$second", "yyyy-mm-dd HH:MM:SS")
    catch
        nothing
    end
end

include("MeasurementIndex/SourceFiles.jl")

"""
One node in the device hierarchy.
"""
struct HierarchyNode
    name::String
    kind::Symbol
    children::Vector{HierarchyNode}
    measurements::Vector{MeasurementInfo}
end

HierarchyNode(name::String, kind::Symbol) =
    HierarchyNode(name, kind, HierarchyNode[], MeasurementInfo[])

"""
The complete device tree and its indexes for one source root.
"""
struct MeasurementHierarchy
    root::HierarchyNode
    all_measurements::Vector{MeasurementInfo}
    root_path::String
    index::Dict{Tuple{Vararg{String}},HierarchyNode}
    has_device_metadata::Bool
    project::AbstractProject
    skipped_count::Int
end

"""
The authoritative result of one completed filesystem scan.
"""
struct SourceScan
    root_path::String
    project::AbstractProject
    files::Vector{SourceFile}
    hierarchy::MeasurementHierarchy
    analysis_failures::Vector{MeasurementAnalysisFailure}
end

"""Construct a successful scan with no recorded analysis failures."""
function SourceScan(
    root_path::String,
    project::AbstractProject,
    files::Vector{SourceFile},
    hierarchy::MeasurementHierarchy,
)::SourceScan
    return SourceScan(root_path, project, files, hierarchy, MeasurementAnalysisFailure[])
end

"""
Convert a Roman-numeral path segment to its integer value.
"""
function roman_value(text::AbstractString)::Union{Nothing,Int}
    values = Dict('I' => 1, 'V' => 5, 'X' => 10, 'L' => 50)
    isempty(text) && return nothing
    total = 0
    previous = 0
    for character in reverse(uppercase(text))
        value = get(values, character, 0)
        value == 0 && return nothing
        total += value < previous ? -value : value
        previous = max(previous, value)
    end
    return total
end

"""
Return a tuple that sorts text segments and embedded integers naturally.
"""
function natural_key(text::AbstractString)::Tuple
    parts = Any[]
    for result in eachmatch(r"\d+|\D+", String(text))
        segment = result.match
        push!(parts, all(isdigit, segment) ?
            (1, parse(Int, segment)) :
            (0, lowercase(segment)))
    end
    return Tuple(parts)
end

measurement_timestamp_key(measurement::MeasurementInfo)::DateTime =
    measurement.timestamp === nothing ?
    DateTime(Dates.year(typemax(Date))) :
    measurement.timestamp

"""
Sort one hierarchy node recursively using natural names and measurement time.
"""
function Base.sort!(node::HierarchyNode)::HierarchyNode
    foreach(sort!, node.children)
    roman = !isempty(node.children) &&
        all(child -> roman_value(child.name) !== nothing, node.children)
    sort!(
        node.children;
        by=roman ?
            child -> roman_value(child.name) :
            child -> natural_key(child.name),
    )
    foreach(
        child -> sort!(child.measurements; by=measurement_timestamp_key),
        node.children,
    )
    return node
end

Base.sort!(hierarchy::MeasurementHierarchy)::MeasurementHierarchy = (
    sort!(hierarchy.root);
    hierarchy
)

"""
Create an empty hierarchy ready for progressive scan results.
"""
function MeasurementHierarchy(
    root_path::String,
    has_device_metadata::Bool,
    project::AbstractProject,
    skipped_count::Int=0,
)::MeasurementHierarchy
    return MeasurementHierarchy(
        HierarchyNode("/", :root),
        MeasurementInfo[],
        root_path,
        Dict{Tuple{Vararg{String}},HierarchyNode}(),
        has_device_metadata,
        project,
        skipped_count,
    )
end

"""
Insert one measurement into an existing hierarchy.
"""
function insert_measurement!(
    hierarchy::MeasurementHierarchy,
    measurement::MeasurementInfo,
)::MeasurementInfo
    parent = hierarchy.root
    for (depth, segment) in enumerate(measurement.device_info.location)
        path = Tuple(measurement.device_info.location[1:depth])
        child = get(hierarchy.index, path, nothing)
        if child === nothing
            kind = depth == length(measurement.device_info.location) ? :leaf : :level
            child = HierarchyNode(segment, kind)
            push!(parent.children, child)
            hierarchy.index[path] = child
        end
        parent = child
    end
    push!(parent.measurements, measurement)
    push!(hierarchy.all_measurements, measurement)
    return measurement
end

"""
Build a complete device tree from a flat measurement list.
"""
function MeasurementHierarchy(
    measurements::Vector{MeasurementInfo},
    root_path::String,
    has_device_metadata::Bool,
    project::AbstractProject,
    skipped_count::Int=0,
)::MeasurementHierarchy
    hierarchy = MeasurementHierarchy(
        root_path,
        has_device_metadata,
        project,
        skipped_count,
    )
    foreach(measurement -> insert_measurement!(hierarchy, measurement), measurements)
    return sort!(hierarchy)
end

children(node::HierarchyNode)::Vector{HierarchyNode} = node.children

compute_and_add_measurement_stats!(
    ::AbstractProject,
    ::Vector{MeasurementInfo},
    ::Vector{SourceFile};
    on_progress::Union{Nothing,Function}=nothing,
)::Vector{MeasurementAnalysisFailure} = MeasurementAnalysisFailure[]

include("MeasurementIndex/Scanning.jl")

end
