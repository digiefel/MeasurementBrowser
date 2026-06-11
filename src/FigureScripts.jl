"""
FigureScripts.jl - Deprecated helpers for publication script generation.

This file is kept for old generated scripts while the app moves toward simpler live workflows.
Do not build new plotting or analysis behavior here.
"""

"""
    MeasurementFilterClause(; source_file=nothing, measurement_kind=nothing,
                             device_path=String[],
                             parameter_conditions=Pair{Symbol,Any}[],
                             timestamp_range=nothing)

One exact filter clause over logical measurements. All populated conditions must match.
`timestamp_range` is an optional `(lo, hi)` `DateTime` interval (inclusive on both ends);
measurements with `timestamp === nothing` never match a clause that sets it.
"""
struct MeasurementFilterClause
    source_file::Union{Nothing,String}
    measurement_kind::Union{Nothing,Symbol}
    device_path::Vector{String}
    parameter_conditions::Vector{Pair{Symbol,Any}}
    timestamp_range::Union{Nothing,Tuple{DateTime,DateTime}}
end

function MeasurementFilterClause(;
    source_file::Union{Nothing,AbstractString}=nothing,
    measurement_kind::Union{Nothing,Symbol}=nothing,
    device_path::AbstractVector{<:AbstractString}=String[],
    parameter_conditions::AbstractVector{<:Pair}=Pair{Symbol,Any}[],
    timestamp_range::Union{Nothing,Tuple{DateTime,DateTime}}=nothing,
)
    source = source_file === nothing ? nothing : String(source_file)
    path = String.(device_path)
    normalized_conditions = _normalize_parameter_conditions(parameter_conditions)
    isempty(path) || any(isempty, path) && error("device_path cannot contain empty segments")
    if timestamp_range !== nothing && timestamp_range[1] > timestamp_range[2]
        error("timestamp_range lower bound must not exceed upper bound")
    end
    return MeasurementFilterClause(source, measurement_kind, path, normalized_conditions, timestamp_range)
end

"""
    MeasurementGroupFilter(clauses)

An OR of exact `MeasurementFilterClause` values.
"""
struct MeasurementGroupFilter
    clauses::Vector{MeasurementFilterClause}
    function MeasurementGroupFilter(clauses::Vector{MeasurementFilterClause})
        isempty(clauses) && throw(FigureScriptValidationError("Measurement group filters must contain at least one clause"))
        return new(clauses)
    end
end

function MeasurementGroupFilter(clauses::AbstractVector{<:MeasurementFilterClause})
    return MeasurementGroupFilter(collect(clauses))
end

"""
    NamedMeasurementGroup(name, filter; exclude_bad=nothing)

One fixed named group in a generated figure script.

`exclude_bad` is an optional `Annotations.Tags.TagState`. When set, the group's
matching pipeline first rejects any measurement carrying the "bad" tag (directly
or via an ancestor device path), then evaluates `filter` against the survivors.
This acts as an implicit first clause that keeps known-bad measurements out of
generated figures unless one of them was explicitly selected.
"""
struct NamedMeasurementGroup
    name::String
    filter::MeasurementGroupFilter
    exclude_bad::Union{Nothing,Annotations.Tags.TagState}
    function NamedMeasurementGroup(
        name::String,
        filter::MeasurementGroupFilter;
        exclude_bad::Union{Nothing,Annotations.Tags.TagState}=nothing,
    )
        normalized_name = strip(name)
        isempty(normalized_name) && throw(FigureScriptValidationError("Measurement group names cannot be empty"))
        return new(normalized_name, filter, exclude_bad)
    end
end

function NamedMeasurementGroup(
    name::AbstractString,
    filter::MeasurementGroupFilter;
    exclude_bad::Union{Nothing,Annotations.Tags.TagState}=nothing,
)
    return NamedMeasurementGroup(String(name), filter; exclude_bad=exclude_bad)
end

"""
    FigureMeasurement

One logical measurement inside figure-script `data`.
"""
struct FigureMeasurement
    measurement::MeasurementInfo
    parameters::Dict{Symbol,Any}
    loaded::Any
    analyzed::Any
end

"""
    FigureScriptData

Typed container backing generated figure scripts.
"""
struct FigureScriptData
    source_files::Vector{String}
    measurements::Vector{FigureMeasurement}
    groups::Dict{String,NamedMeasurementGroup}
    group_names::Vector{String}
    group_matches::Dict{String,Vector{Int}}
end

struct FigureScriptProfileSection
    key::Symbol
    calls::Int
    duration_ms::Float64
    alloc_bytes::Int
end

struct FigureScriptInferenceProfile
    group_name::String
    selected_count::Int
    measurement_count::Int
    total_ms::Float64
    sections::Vector{FigureScriptProfileSection}
    counters::Dict{Symbol,Int}
end

struct FigureScriptProfiledError <: Exception
    cause::Exception
    bt
    profile::FigureScriptInferenceProfile
end

Base.showerror(io::IO, err::FigureScriptProfiledError) = showerror(io, err.cause)

struct FigureScriptValidationError <: Exception
    message::String
end

Base.showerror(io::IO, err::FigureScriptValidationError) = print(io, err.message)

struct FigureScriptResolutionError <: Exception
    message::String
end

Base.showerror(io::IO, err::FigureScriptResolutionError) = print(io, err.message)

struct FigureScriptExistsError <: Exception
    path::String
end

Base.showerror(io::IO, err::FigureScriptExistsError) = print(io, "Figure script already exists: ", err.path)

struct FigureScriptIOError <: Exception
    path::String
    operation::Symbol
    cause::Exception
end

function Base.showerror(io::IO, err::FigureScriptIOError)
    print(io, "Could not $(String(err.operation)) $(err.path): ")
    showerror(io, err.cause)
end

function Base.haskey(data::FigureScriptData, name::AbstractString)
    return haskey(data.groups, String(name))
end

function Base.keys(data::FigureScriptData)
    return copy(data.group_names)
end

function Base.getindex(data::FigureScriptData, name::AbstractString)
    group_name = String(name)
    haskey(data.groups, group_name) || throw(FigureScriptResolutionError(
        "Figure-script group '$group_name' is not defined",
    ))
    indices = get(data.group_matches, group_name, nothing)
    indices === nothing && throw(FigureScriptResolutionError(
        "Figure-script group '$group_name' has no computed matches",
    ))
    return [data.measurements[index] for index in indices]
end

function _normalize_parameter_conditions(parameter_conditions::AbstractVector{<:Pair})
    normalized = Pair{Symbol,Any}[]
    seen = Set{Symbol}()
    for condition in parameter_conditions
        name = Symbol(condition.first)
        name in seen && error("Duplicate parameter condition for '$name'")
        push!(seen, name)
        push!(normalized, name => condition.second)
    end
    sort!(normalized; by=condition -> String(condition.first))
    return normalized
end

function _validate_named_measurement_groups(groups::Vector{NamedMeasurementGroup})
    isempty(groups) && throw(FigureScriptValidationError("At least one measurement group is required"))

    seen_names = Set{String}()
    for group in groups
        group.name in seen_names && throw(FigureScriptValidationError("Duplicate measurement group name '$(group.name)'"))
        push!(seen_names, group.name)
        isempty(group.filter.clauses) && throw(FigureScriptValidationError("Measurement group '$(group.name)' is empty"))
    end
    return nothing
end

function _figure_script_output_directory(output_directory::AbstractString)
    directory = strip(String(output_directory))
    isempty(directory) && throw(FigureScriptValidationError("Choose an output directory before generating a figure script"))
    isabspath(directory) || throw(FigureScriptValidationError("Figure script output directory must be an absolute path"))
    return normpath(abspath(directory))
end

function _figure_script_basename(script_name::AbstractString)
    base = strip(String(script_name))
    isempty(base) && throw(FigureScriptValidationError("Script name cannot be empty"))
    occursin('/', base) && throw(FigureScriptValidationError("Script name must not contain '/'"))
    occursin('\\', base) && throw(FigureScriptValidationError("Script name must not contain '\\'"))
    if endswith(lowercase(base), ".jl")
        base = base[1:end-3]
    end
    isempty(base) && throw(FigureScriptValidationError("Script name cannot be empty"))
    return base
end

function figure_script_path(output_directory::AbstractString, script_name::AbstractString)
    return joinpath(
        _figure_script_output_directory(output_directory),
        _figure_script_basename(script_name) * ".jl",
    )
end

function _script_project_binding(::TASEProject)
    return "MeasurementBrowser.TASE_PROJECT"
end

function _path_within_root(root_path::AbstractString, path::AbstractString)
    root_abs = normpath(abspath(String(root_path)))
    path_abs = normpath(abspath(String(path)))
    root_parts = splitpath(root_abs)
    path_parts = splitpath(path_abs)
    length(path_parts) >= length(root_parts) || return false
    return path_parts[1:length(root_parts)] == root_parts
end

function _validate_source_file(root_path::AbstractString, source_file::AbstractString)
    normalized = normpath(abspath(String(source_file)))
    _path_within_root(root_path, normalized) || throw(FigureScriptResolutionError(
        "Figure script source file is outside root path '$root_path': $normalized",
    ))
    isfile(normalized) || throw(FigureScriptResolutionError(
        "Figure script source file does not exist: $normalized",
    ))
    return normalized
end

function _normalize_source_files(root_path::AbstractString, source_files::Vector{String})
    isempty(source_files) && throw(FigureScriptValidationError("At least one source file is required"))
    seen = Set{String}()
    normalized = String[]
    for source_file in source_files
        file = _validate_source_file(root_path, source_file)
        file in seen && continue
        push!(seen, file)
        push!(normalized, file)
    end
    return normalized
end

function _is_filterable_parameter_value(value)
    return value isa Bool ||
        value isa Integer ||
        value isa AbstractFloat ||
        value isa AbstractString ||
        value isa Symbol ||
        value isa Dates.TimeType
end

function _measurement_parameters(
    device_parameters::Dict{Symbol,Any},
    measurement_parameters::Dict{Symbol,Any},
)
    return merge(device_parameters, measurement_parameters)
end

function _measurement_parameters(measurement::MeasurementInfo)
    return merge(
        _measurement_parameters(measurement.device_info.parameters, measurement.parameters),
        measurement.stats,
    )
end

function _measurement_parameter_conditions(parameters::Dict{Symbol,Any})
    conditions = Pair{Symbol,Any}[]
    for (name, value) in sort!(collect(parameters); by=entry -> String(first(entry)))
        _is_filterable_parameter_value(value) || continue
        push!(conditions, name => value)
    end
    return conditions
end

function _measurement_matches_clause(measurement::MeasurementInfo, parameters::Dict{Symbol,Any}, clause::MeasurementFilterClause)
    clause.source_file !== nothing && measurement.filepath != clause.source_file && return false
    clause.measurement_kind !== nothing && measurement.measurement_kind != clause.measurement_kind && return false

    if !isempty(clause.device_path)
        length(measurement.device_info.location) >= length(clause.device_path) || return false
        measurement.device_info.location[1:length(clause.device_path)] == clause.device_path || return false
    end

    for condition in clause.parameter_conditions
        haskey(parameters, condition.first) || return false
        isequal(parameters[condition.first], condition.second) || return false
    end

    if clause.timestamp_range !== nothing
        measurement.timestamp === nothing && return false
        lo, hi = clause.timestamp_range
        (lo <= measurement.timestamp <= hi) || return false
    end
    return true
end

function _measurement_matches_filter(measurement::MeasurementInfo, parameters::Dict{Symbol,Any}, filter::MeasurementGroupFilter)
    for clause in filter.clauses
        _measurement_matches_clause(measurement, parameters, clause) && return true
    end
    return false
end

function _ancestor_device_keys(dev_key::AbstractString)
    parts = split(dev_key, '/'; keepempty=false)
    return [join(parts[1:i], '/') for i in 1:(length(parts) - 1)]
end

function _measurement_is_bad(tag_state::Annotations.Tags.TagState, measurement::MeasurementInfo)
    dev_key = device_path_key(measurement.device_info)
    ancestors = _ancestor_device_keys(dev_key)
    return "bad" in Annotations.Tags.effective(tag_state, measurement.unique_id, [dev_key; ancestors])
end

function _measurement_matches_group(
    measurement::MeasurementInfo,
    parameters::Dict{Symbol,Any},
    group::NamedMeasurementGroup,
)
    group.exclude_bad !== nothing && _measurement_is_bad(group.exclude_bad, measurement) && return false
    return _measurement_matches_filter(measurement, parameters, group.filter)
end

function _matching_measurements(measurements::Vector{MeasurementInfo}, group::NamedMeasurementGroup)
    return _matching_measurements(measurements, _measurement_parameter_maps(measurements), group)
end

function _matching_measurements(
    measurements::Vector{MeasurementInfo},
    parameter_maps::Vector{Dict{Symbol,Any}},
    group::NamedMeasurementGroup,
)
    length(measurements) == length(parameter_maps) || error("Measurement and parameter map lengths must match")
    matched = MeasurementInfo[]
    for index in eachindex(measurements, parameter_maps)
        measurement = measurements[index]
        _measurement_matches_group(measurement, parameter_maps[index], group) && push!(matched, measurement)
    end
    return matched
end

function _load_measurements_from_source_files(
    root_path::AbstractString,
    project::AbstractProject,
    source_files::Vector{String},
)
    source_files = _normalize_source_files(root_path, source_files)
    meta = load_scan_metadata(String(root_path))
    measurements = MeasurementInfo[]
    for source_file in source_files
        indexed = index_source_file(source_file)
        check_cancel()
        items = interpret_measurements(project, indexed, meta)
        isempty(items) && throw(FigureScriptResolutionError(
            "Measurement source file '$source_file' did not produce any measurements for project $(project_name(project))",
        ))
        append!(measurements, [_measurement_info_from_item(item) for item in items])
    end
    return measurements
end

function _group_matches(groups::Vector{NamedMeasurementGroup}, measurements::Vector{MeasurementInfo})
    return _group_matches(groups, measurements, _measurement_parameter_maps(measurements))
end

function _group_matches(
    groups::Vector{NamedMeasurementGroup},
    measurements::Vector{MeasurementInfo},
    parameter_maps::Vector{Dict{Symbol,Any}},
)
    matches = Dict{String,Vector{MeasurementInfo}}()
    for group in groups
        matched = _matching_measurements(measurements, parameter_maps, group)
        isempty(matched) && throw(FigureScriptResolutionError(
            "Measurement group '$(group.name)' matches no measurements in the current dataset",
        ))
        matches[group.name] = matched
    end
    return matches
end

function _figure_script_plot_kind(project::AbstractProject, measurement::MeasurementInfo)::Type{<:PlotKind}
    throw(FigureScriptResolutionError(
        "Figure script plot data is not implemented for $(project_name(project))",
    ))
end

function _build_figure_measurements(
    workspace::Workspace.Workspace,
    measurements::Vector{MeasurementInfo},
)
    project = workspace.project
    records = FigureMeasurement[]
    for measurement in measurements
        parameters = _measurement_parameters(measurement)
        plot_kind = _figure_script_plot_kind(project, measurement)
        df = only(read_measurement_data(workspace, [measurement]))
        loaded = _ruo2_plot_data(measurement, df)
        analyzed = _analyze_ruo2_file_plot(
            project,
            plot_kind,
            loaded;
            device_params=parameters,
        )
        push!(records, FigureMeasurement(measurement, parameters, loaded, analyzed))
    end
    return records
end

"""
    prepare_figure_script_data(root_path, project, source_files, groups)

Load the source-file union for a figure script and expose the fixed named groups as typed data.
"""
function prepare_figure_script_data(
    root_path::AbstractString,
    project::AbstractProject,
    source_files::Vector{String},
    groups::Vector{NamedMeasurementGroup},
)
    _validate_named_measurement_groups(groups)
    normalized_source_files = _normalize_source_files(root_path, source_files)
    measurements = _load_measurements_from_source_files(root_path, project, normalized_source_files)
    _group_matches(groups, measurements)
    workspace = Workspace.Workspace(project, root_path)
    records = _build_figure_measurements(workspace, measurements)

    group_map = Dict(group.name => group for group in groups)
    match_indices = Dict{String,Vector{Int}}()
    for group in groups
        indices = Int[]
        for (index, record) in enumerate(records)
            _measurement_matches_group(record.measurement, record.parameters, group) && push!(indices, index)
        end
        isempty(indices) && throw(FigureScriptResolutionError(
            "Measurement group '$(group.name)' matches no logical measurements after loading",
        ))
        match_indices[group.name] = indices
    end

    return FigureScriptData(
        copy(normalized_source_files),
        records,
        group_map,
        [group.name for group in groups],
        match_indices,
    )
end

function _clause_condition_count(clause::MeasurementFilterClause)
    return (clause.source_file === nothing ? 0 : 1) +
        (clause.measurement_kind === nothing ? 0 : 1) +
        (isempty(clause.device_path) ? 0 : 1) +
        length(clause.parameter_conditions) +
        (clause.timestamp_range === nothing ? 0 : 1)
end

function _measurement_parameter_maps(measurements::Vector{MeasurementInfo})
    return [_measurement_parameters(measurement) for measurement in measurements]
end

mutable struct _FigureScriptProfileAccumulator
    group_name::String
    selected_count::Int
    measurement_count::Int
    started_ns::UInt64
    sections::Dict{Symbol,FigureScriptProfileSection}
    counters::Dict{Symbol,Int}
end

function _new_figure_script_profile_accumulator(
    group_name::AbstractString,
    selected_count::Int,
    measurement_count::Int,
)
    return _FigureScriptProfileAccumulator(
        String(group_name),
        selected_count,
        measurement_count,
        time_ns(),
        Dict{Symbol,FigureScriptProfileSection}(),
        Dict{Symbol,Int}(),
    )
end

struct _FigureScriptFactIndex
    measurements::Vector{MeasurementInfo}
    parameter_maps::Vector{Dict{Symbol,Any}}
end

function _build_figure_script_fact_index(
    measurements::Vector{MeasurementInfo},
    profiler::Union{Nothing,_FigureScriptProfileAccumulator}=nothing,
)
    parameter_maps = _profile_section!(profiler, :collect_parameter_maps) do
        _measurement_parameter_maps(measurements)
    end
    return _FigureScriptFactIndex(measurements, parameter_maps)
end

function _matching_measurements(index::_FigureScriptFactIndex, group::NamedMeasurementGroup)
    return _matching_measurements(index.measurements, index.parameter_maps, group)
end

function _profile_section!(f::Function, ::Nothing, key::Symbol)
    return f()
end

function _profile_section!(f::Function, acc::_FigureScriptProfileAccumulator, key::Symbol)
    value = nothing
    t0 = time_ns()
    alloc_bytes = @allocated value = f()
    duration_ms = (time_ns() - t0) / 1e6
    current = get(acc.sections, key, FigureScriptProfileSection(key, 0, 0.0, 0))
    acc.sections[key] = FigureScriptProfileSection(
        key,
        current.calls + 1,
        current.duration_ms + duration_ms,
        current.alloc_bytes + alloc_bytes,
    )
    return value
end

_profile_counter!(::Nothing, key::Symbol, delta::Int=1) = nothing

function _profile_counter!(acc::_FigureScriptProfileAccumulator, key::Symbol, delta::Int=1)
    acc.counters[key] = get(acc.counters, key, 0) + delta
    return nothing
end

function _finalize_figure_script_profile(acc::_FigureScriptProfileAccumulator)
    total_ms = (time_ns() - acc.started_ns) / 1e6
    sections = sort!(collect(values(acc.sections)); by=section -> String(section.key))
    return FigureScriptInferenceProfile(
        acc.group_name,
        acc.selected_count,
        acc.measurement_count,
        total_ms,
        sections,
        copy(acc.counters),
    )
end

function _selected_measurement_ids(selected_measurements::Vector{MeasurementInfo})
    selected_ids = Set(measurement.unique_id for measurement in selected_measurements)
    length(selected_ids) == length(selected_measurements) || throw(FigureScriptValidationError(
        "Figure-script group selection contains duplicate logical measurements",
    ))
    return selected_ids
end

function _selected_measurement_indices(selected_ids::Set{String}, all_measurements::Vector{MeasurementInfo})
    selected_indices = findall(measurement -> measurement.unique_id in selected_ids, all_measurements)
    length(selected_indices) == length(selected_ids) || throw(FigureScriptResolutionError(
        "Figure-script group selection references measurements that are not present in the current scan",
    ))
    return selected_indices
end

function _subset_mask(indices::Vector{Int}, measurement_count::Int)
    mask = falses(measurement_count)
    @inbounds for i in indices
        mask[i] = true
    end
    return mask
end

function _check_clause_mask!(
    mask::BitVector,
    measurements::Vector{MeasurementInfo},
    parameter_maps::Vector{Dict{Symbol,Any}},
    clause::MeasurementFilterClause,
    profiler::Union{Nothing,_FigureScriptProfileAccumulator}=nothing,
)
    _profile_counter!(profiler, :full_mask_checks)
    @inbounds for i in eachindex(measurements)
        mask[i] = _measurement_matches_clause(measurements[i], parameter_maps[i], clause)
    end
    return mask
end

function _exact_selector_clause(measurement::MeasurementInfo, parameters::Dict{Symbol,Any})
    return MeasurementFilterClause(
        source_file=measurement.filepath,
        measurement_kind=measurement.measurement_kind,
        device_path=measurement.device_info.location,
        parameter_conditions=_measurement_parameter_conditions(parameters),
    )
end

function _longest_common_prefix(locations::Vector{Vector{String}})
    isempty(locations) && return String[]
    min_length = minimum(length, locations)
    common = String[]
    @inbounds for i in 1:min_length
        first_value = locations[1][i]
        all_equal = true
        for j in 2:length(locations)
            if locations[j][i] != first_value
                all_equal = false
                break
            end
        end
        all_equal || break
        push!(common, first_value)
    end
    return common
end

function _shared_parameter_conditions(
    subset::Vector{Int},
    parameter_maps::Vector{Dict{Symbol,Any}},
    key_filter::Function,
    skip_names::Set{Symbol},
)
    isempty(subset) && return Pair{Symbol,Any}[]
    first_map = parameter_maps[subset[1]]
    shared = Pair{Symbol,Any}[]
    for (name, value) in first_map
        name in skip_names && continue
        key_filter(name, subset[1]) || continue
        _is_filterable_parameter_value(value) || continue
        all_match = true
        for i in 2:length(subset)
            other = parameter_maps[subset[i]]
            if !key_filter(name, subset[i]) || !haskey(other, name) || !isequal(other[name], value)
                all_match = false
                break
            end
        end
        all_match && push!(shared, name => value)
    end
    sort!(shared; by=condition -> String(condition.first))
    return shared
end

function _describe_device_stage(
    measurements::Vector{MeasurementInfo},
    parameter_maps::Vector{Dict{Symbol,Any}},
    subset::Vector{Int},
)
    locations = [measurements[i].device_info.location for i in subset]
    common_prefix = _longest_common_prefix(locations)
    shared_device = _shared_parameter_conditions(
        subset, parameter_maps,
        (name, i) -> haskey(measurements[i].device_info.parameters, name),
        Set{Symbol}(),
    )
    return MeasurementFilterClause(
        device_path=common_prefix,
        parameter_conditions=shared_device,
    )
end

function _add_measurement_stage(
    base::MeasurementFilterClause,
    measurements::Vector{MeasurementInfo},
    parameter_maps::Vector{Dict{Symbol,Any}},
    subset::Vector{Int},
)
    isempty(subset) && return base
    first_kind = measurements[subset[1]].measurement_kind
    kind = all(measurements[i].measurement_kind == first_kind for i in subset) ? first_kind : nothing
    existing_names = Set{Symbol}(condition.first for condition in base.parameter_conditions)
    shared_measurement = _shared_parameter_conditions(
        subset, parameter_maps,
        (name, i) -> haskey(measurements[i].parameters, name),
        existing_names,
    )
    combined = vcat(base.parameter_conditions, shared_measurement)
    return MeasurementFilterClause(
        source_file=base.source_file,
        measurement_kind=kind,
        device_path=base.device_path,
        parameter_conditions=combined,
        timestamp_range=base.timestamp_range,
    )
end

function _add_timestamp_stage(
    base::MeasurementFilterClause,
    measurements::Vector{MeasurementInfo},
    subset::Vector{Int},
)
    timestamps = DateTime[]
    for i in subset
        ts = measurements[i].timestamp
        ts === nothing && return nothing
        push!(timestamps, ts)
    end
    isempty(timestamps) && return nothing
    return MeasurementFilterClause(
        source_file=base.source_file,
        measurement_kind=base.measurement_kind,
        device_path=base.device_path,
        parameter_conditions=base.parameter_conditions,
        timestamp_range=(minimum(timestamps), maximum(timestamps)),
    )
end

function _with_source_file(base::MeasurementFilterClause, source_file::AbstractString)
    return MeasurementFilterClause(
        source_file=source_file,
        measurement_kind=base.measurement_kind,
        device_path=base.device_path,
        parameter_conditions=base.parameter_conditions,
        timestamp_range=base.timestamp_range,
    )
end

function _split_by_next_device_level(
    measurements::Vector{MeasurementInfo},
    subset::Vector{Int},
    shared_prefix::Vector{String},
)
    next_pos = length(shared_prefix) + 1
    groups = Dict{Union{Nothing,String},Vector{Int}}()
    order = Union{Nothing,String}[]
    for i in subset
        location = measurements[i].device_info.location
        key = length(location) >= next_pos ? location[next_pos] : nothing
        if !haskey(groups, key)
            groups[key] = Int[]
            push!(order, key)
        end
        push!(groups[key], i)
    end
    return [groups[k] for k in order]
end

function _split_by_kind(measurements::Vector{MeasurementInfo}, subset::Vector{Int})
    groups = Dict{Symbol,Vector{Int}}()
    order = Symbol[]
    for i in subset
        kind = measurements[i].measurement_kind
        if !haskey(groups, kind)
            groups[kind] = Int[]
            push!(order, kind)
        end
        push!(groups[kind], i)
    end
    return [groups[k] for k in order]
end

function _split_by_source_file(measurements::Vector{MeasurementInfo}, subset::Vector{Int})
    groups = Dict{String,Vector{Int}}()
    order = String[]
    for i in subset
        file = measurements[i].filepath
        if !haskey(groups, file)
            groups[file] = Int[]
            push!(order, file)
        end
        push!(groups[file], i)
    end
    return [groups[k] for k in order]
end

function _mask_equals(left::BitVector, right::BitVector)
    length(left) == length(right) || return false
    @inbounds for i in eachindex(left)
        left[i] == right[i] || return false
    end
    return true
end

function _minimize_clause(
    clause::MeasurementFilterClause,
    measurements::Vector{MeasurementInfo},
    parameter_maps::Vector{Dict{Symbol,Any}},
    subset_mask::BitVector,
    scratch_mask::BitVector,
    profiler::Union{Nothing,_FigureScriptProfileAccumulator}=nothing,
)
    candidate = clause
    if candidate.source_file !== nothing
        trial = MeasurementFilterClause(
            source_file=nothing,
            measurement_kind=candidate.measurement_kind,
            device_path=candidate.device_path,
            parameter_conditions=candidate.parameter_conditions,
            timestamp_range=candidate.timestamp_range,
        )
        _check_clause_mask!(scratch_mask, measurements, parameter_maps, trial, profiler)
        _mask_equals(scratch_mask, subset_mask) && (candidate = trial)
    end
    if candidate.timestamp_range !== nothing
        trial = MeasurementFilterClause(
            source_file=candidate.source_file,
            measurement_kind=candidate.measurement_kind,
            device_path=candidate.device_path,
            parameter_conditions=candidate.parameter_conditions,
            timestamp_range=nothing,
        )
        _check_clause_mask!(scratch_mask, measurements, parameter_maps, trial, profiler)
        _mask_equals(scratch_mask, subset_mask) && (candidate = trial)
    end
    if !isempty(candidate.parameter_conditions)
        i = length(candidate.parameter_conditions)
        while i >= 1
            without = vcat(candidate.parameter_conditions[1:i-1], candidate.parameter_conditions[i+1:end])
            trial = MeasurementFilterClause(
                source_file=candidate.source_file,
                measurement_kind=candidate.measurement_kind,
                device_path=candidate.device_path,
                parameter_conditions=without,
                timestamp_range=candidate.timestamp_range,
            )
            _check_clause_mask!(scratch_mask, measurements, parameter_maps, trial, profiler)
            if _mask_equals(scratch_mask, subset_mask)
                candidate = trial
            end
            i -= 1
        end
    end
    if candidate.measurement_kind !== nothing
        trial = MeasurementFilterClause(
            source_file=candidate.source_file,
            measurement_kind=nothing,
            device_path=candidate.device_path,
            parameter_conditions=candidate.parameter_conditions,
            timestamp_range=candidate.timestamp_range,
        )
        _check_clause_mask!(scratch_mask, measurements, parameter_maps, trial, profiler)
        _mask_equals(scratch_mask, subset_mask) && (candidate = trial)
    end
    if !isempty(candidate.device_path)
        trial = MeasurementFilterClause(
            source_file=candidate.source_file,
            measurement_kind=candidate.measurement_kind,
            device_path=String[],
            parameter_conditions=candidate.parameter_conditions,
            timestamp_range=candidate.timestamp_range,
        )
        _check_clause_mask!(scratch_mask, measurements, parameter_maps, trial, profiler)
        _mask_equals(scratch_mask, subset_mask) && (candidate = trial)
    end
    return candidate
end

function _try_emit_stage!(
    clauses::Vector{MeasurementFilterClause},
    clause::MeasurementFilterClause,
    measurements::Vector{MeasurementInfo},
    parameter_maps::Vector{Dict{Symbol,Any}},
    subset_mask::BitVector,
    scratch_mask::BitVector,
    profiler::Union{Nothing,_FigureScriptProfileAccumulator},
)
    _check_clause_mask!(scratch_mask, measurements, parameter_maps, clause, profiler)
    if _mask_equals(scratch_mask, subset_mask)
        minimized = _minimize_clause(clause, measurements, parameter_maps, subset_mask, scratch_mask, profiler)
        push!(clauses, minimized)
        return true
    end
    return false
end

function _infer_clauses_for_subset!(
    clauses::Vector{MeasurementFilterClause},
    measurements::Vector{MeasurementInfo},
    parameter_maps::Vector{Dict{Symbol,Any}},
    subset::Vector{Int},
    scratch_mask::BitVector,
    profiler::Union{Nothing,_FigureScriptProfileAccumulator},
)
    isempty(subset) && return
    _profile_counter!(profiler, :recursive_subsets)

    subset_mask = _subset_mask(subset, length(measurements))

    clause_device = _profile_section!(profiler, :stage_device) do
        _describe_device_stage(measurements, parameter_maps, subset)
    end
    _try_emit_stage!(clauses, clause_device, measurements, parameter_maps, subset_mask, scratch_mask, profiler) && return

    clause_measurement = _profile_section!(profiler, :stage_measurement) do
        _add_measurement_stage(clause_device, measurements, parameter_maps, subset)
    end
    _try_emit_stage!(clauses, clause_measurement, measurements, parameter_maps, subset_mask, scratch_mask, profiler) && return

    clause_timestamp = _profile_section!(profiler, :stage_timestamp) do
        _add_timestamp_stage(clause_measurement, measurements, subset)
    end
    if clause_timestamp !== nothing
        _try_emit_stage!(clauses, clause_timestamp, measurements, parameter_maps, subset_mask, scratch_mask, profiler) && return
    end

    # Splits, in priority order
    kind_parts = _split_by_kind(measurements, subset)
    if length(kind_parts) >= 2
        for part in kind_parts
            _infer_clauses_for_subset!(clauses, measurements, parameter_maps, part, scratch_mask, profiler)
        end
        return
    end

    device_parts = _split_by_next_device_level(measurements, subset, clause_device.device_path)
    if length(device_parts) >= 2
        for part in device_parts
            _infer_clauses_for_subset!(clauses, measurements, parameter_maps, part, scratch_mask, profiler)
        end
        return
    end

    _profile_counter!(profiler, :stage_fallback_entries)
    source_parts = _split_by_source_file(measurements, subset)
    if length(source_parts) >= 2
        for part in source_parts
            _infer_clauses_for_subset!(clauses, measurements, parameter_maps, part, scratch_mask, profiler)
        end
        return
    end

    shared_source = measurements[subset[1]].filepath
    clause_with_source = _with_source_file(clause_measurement, shared_source)
    if _try_emit_stage!(clauses, clause_with_source, measurements, parameter_maps, subset_mask, scratch_mask, profiler)
        return
    end

    _profile_counter!(profiler, :exact_selector_emits, length(subset))
    for i in subset
        push!(clauses, _exact_selector_clause(measurements[i], parameter_maps[i]))
    end
end

function _prepare_group_inference(
    selected_measurements::Vector{MeasurementInfo},
    all_measurements::Vector{MeasurementInfo},
    profiler::Union{Nothing,_FigureScriptProfileAccumulator},
)
    isempty(selected_measurements) && throw(FigureScriptValidationError(
        "Select one or more measurements before creating a figure-script group",
    ))
    isempty(all_measurements) && throw(FigureScriptValidationError(
        "No scanned measurements are available for figure-script inference",
    ))

    selected_ids = _profile_section!(profiler, :selected_ids) do
        _selected_measurement_ids(selected_measurements)
    end
    selected_indices = _profile_section!(profiler, :selected_indices) do
        _selected_measurement_indices(selected_ids, all_measurements)
    end
    _profile_counter!(profiler, :selected_count, length(selected_indices))
    _profile_counter!(profiler, :measurement_count, length(all_measurements))

    parameter_maps = _profile_section!(profiler, :collect_parameter_maps) do
        _measurement_parameter_maps(all_measurements)
    end
    return selected_ids, selected_indices, parameter_maps
end

function _finalize_group_inference(
    name::AbstractString,
    selected_ids::Set{String},
    measurements::Vector{MeasurementInfo},
    parameter_maps::Vector{Dict{Symbol,Any}},
    filter::MeasurementGroupFilter,
    exclude_bad::Union{Nothing,Annotations.Tags.TagState},
    profiler::Union{Nothing,_FigureScriptProfileAccumulator},
)
    _profile_counter!(profiler, :final_clause_count, length(filter.clauses))
    group = NamedMeasurementGroup(name, filter; exclude_bad=exclude_bad)
    matched_ids = _profile_section!(profiler, :verify) do
        Set(measurement.unique_id for measurement in _matching_measurements(measurements, parameter_maps, group))
    end
    matched_ids == selected_ids || throw(FigureScriptResolutionError(
        "Inferred figure-script group '$(group.name)' does not reproduce the selected logical measurements exactly",
    ))
    return group
end

function _effective_exclude_bad(
    tag_state::Union{Nothing,Annotations.Tags.TagState},
    selected_measurements::Vector{MeasurementInfo},
)
    tag_state === nothing && return nothing
    # If any selected measurement is itself bad, leave bad measurements in the
    # match pool so the explicit selection is honoured.
    any(m -> _measurement_is_bad(tag_state, m), selected_measurements) && return nothing
    return tag_state
end

function _infer_measurement_group_staged(
    name::AbstractString,
    selected_measurements::Vector{MeasurementInfo},
    all_measurements::Vector{MeasurementInfo},
    profiler::Union{Nothing,_FigureScriptProfileAccumulator};
    tag_state::Union{Nothing,Annotations.Tags.TagState}=nothing,
)
    exclude_bad = _effective_exclude_bad(tag_state, selected_measurements)
    active_measurements = if exclude_bad === nothing
        all_measurements
    else
        filter(m -> !_measurement_is_bad(exclude_bad, m), all_measurements)
    end

    selected_ids, selected_indices, parameter_maps = _prepare_group_inference(
        selected_measurements, active_measurements, profiler,
    )
    scratch_mask = falses(length(active_measurements))
    clauses = MeasurementFilterClause[]
    _profile_section!(profiler, :infer_clauses) do
        _infer_clauses_for_subset!(
            clauses, active_measurements, parameter_maps, selected_indices, scratch_mask, profiler,
        )
    end
    isempty(clauses) && throw(FigureScriptResolutionError(
        "Could not infer any figure-script filter clauses from the current selection",
    ))
    filter_obj = MeasurementGroupFilter(clauses)
    # Verify against the full measurement list; exclude_bad on the group will
    # filter out bad measurements again so the result reproduces the selection.
    verify_parameter_maps = if exclude_bad === nothing
        parameter_maps
    else
        _measurement_parameter_maps(all_measurements)
    end
    return _finalize_group_inference(
        name, selected_ids, all_measurements, verify_parameter_maps, filter_obj, exclude_bad, profiler,
    )
end

function infer_measurement_group(
    name::AbstractString,
    selected_measurements::Vector{MeasurementInfo},
    all_measurements::Vector{MeasurementInfo};
    tag_state::Union{Nothing,Annotations.Tags.TagState}=nothing,
)
    try
        group, _ = _infer_measurement_group_profiled(
            name, selected_measurements, all_measurements; tag_state=tag_state,
        )
        return group
    catch err
        err isa FigureScriptProfiledError || rethrow()
        rethrow(err.cause)
    end
end

function _infer_measurement_group_profiled(
    name::AbstractString,
    selected_measurements::Vector{MeasurementInfo},
    all_measurements::Vector{MeasurementInfo};
    tag_state::Union{Nothing,Annotations.Tags.TagState}=nothing,
)
    profiler = _new_figure_script_profile_accumulator(name, length(selected_measurements), length(all_measurements))
    try
        group = _infer_measurement_group_staged(
            name,
            selected_measurements,
            all_measurements,
            profiler;
            tag_state=tag_state,
        )
        return group, _finalize_figure_script_profile(profiler)
    catch err
        throw(FigureScriptProfiledError(err, catch_backtrace(), _finalize_figure_script_profile(profiler)))
    end
end

function _group_source_files(groups::Vector{NamedMeasurementGroup}, measurements::Vector{MeasurementInfo})
    matched = _group_matches(groups, measurements)
    needed_files = Set{String}()
    source_files = String[]
    for measurement in measurements
        needed = false
        for group in groups
            measurement in matched[group.name] && (needed = true; break)
        end
        if needed && measurement.filepath ∉ needed_files
            push!(needed_files, measurement.filepath)
            push!(source_files, measurement.filepath)
        end
    end
    return source_files
end

function _render_clause(io::IO, clause::MeasurementFilterClause)
    println(io, "            MeasurementBrowser.MeasurementFilterClause(")
    clause.source_file !== nothing && println(io, "                source_file = ", repr(clause.source_file), ",")
    clause.measurement_kind !== nothing && println(io, "                measurement_kind = ", repr(clause.measurement_kind), ",")
    if !isempty(clause.device_path)
        println(io, "                device_path = ", repr(clause.device_path), ",")
    end
    if !isempty(clause.parameter_conditions)
        println(io, "                parameter_conditions = [")
        for condition in clause.parameter_conditions
            println(io, "                    ", repr(condition.first), " => ", repr(condition.second), ",")
        end
        println(io, "                ],")
    end
    if clause.timestamp_range !== nothing
        lo, hi = clause.timestamp_range
        println(io, "                timestamp_range = (", repr(lo), ", ", repr(hi), "),")
    end
    println(io, "            ),")
end

function _render_group(io::IO, group::NamedMeasurementGroup)
    println(io, "    MeasurementBrowser.NamedMeasurementGroup(")
    println(io, "        ", repr(group.name), ",")
    println(io, "        MeasurementBrowser.MeasurementGroupFilter([")
    for clause in group.filter.clauses
        _render_clause(io, clause)
    end
    if group.exclude_bad === nothing
        println(io, "        ]),")
    else
        println(io, "        ]);")
        println(io, "        exclude_bad = tag_state,")
    end
    println(io, "    ),")
end

function _render_figure_script(
    root_path::AbstractString,
    project::AbstractProject,
    source_files::Vector{String},
    groups::Vector{NamedMeasurementGroup},
    measurements::Vector{MeasurementInfo},
)
    _validate_named_measurement_groups(groups)
    normalized_source_files = _normalize_source_files(root_path, source_files)
    _group_matches(groups, measurements)
    root_abs = normpath(abspath(String(root_path)))
    project_binding = _script_project_binding(project)

    io = IOBuffer()
    println(io, "using MeasurementBrowser")
    needs_dates = any(group -> any(clause -> clause.timestamp_range !== nothing, group.filter.clauses), groups)
    needs_dates && println(io, "using Dates")
    needs_tag_state = any(group -> group.exclude_bad !== nothing, groups)
    needs_tag_state && println(io, "using Annotations")
    println(io)
    println(io, "# Read this file: ", normpath(joinpath(@__DIR__, "..", "docs", "figure_scripts.md")))
    println(io)
    println(io, "root_path = ", repr(root_abs))
    println(io, "project = ", project_binding)
    needs_tag_state && println(io, "tag_state = Annotations.Tags.load(root_path)")
    println(io, "source_files = joinpath.(root_path, [")
    for source_file in normalized_source_files
        rel = relpath(source_file, root_abs)
        println(io, "    ", repr(rel), ",")
    end
    println(io, "])")
    println(io)
    println(io, "groups = [")
    for group in groups
        _render_group(io, group)
    end
    println(io, "]")
    println(io)
    println(io, "data = MeasurementBrowser.prepare_figure_script_data(root_path, project, source_files, groups)")
    return String(take!(io))
end

"""
    write_figure_script(output_directory, root_path, project, script_name, groups, measurements; overwrite=false)

Render and write a figure-preparation script.
"""
function write_figure_script(
    output_directory::AbstractString,
    root_path::AbstractString,
    project::AbstractProject,
    script_name::AbstractString,
    groups::Vector{NamedMeasurementGroup},
    measurements::Vector{MeasurementInfo};
    overwrite::Bool=false,
)
    _validate_named_measurement_groups(groups)
    isempty(measurements) && throw(FigureScriptValidationError("Current scan contains no measurements for figure-script generation"))
    path = figure_script_path(output_directory, script_name)
    if isfile(path) && !overwrite
        throw(FigureScriptExistsError(path))
    end

    source_files = _group_source_files(groups, measurements)
    script_contents = _render_figure_script(root_path, project, source_files, groups, measurements)
    directory = dirname(path)
    try
        mkpath(directory)
        write(path, script_contents)
    catch err
        (err isa IOError || err isa SystemError) && throw(FigureScriptIOError(path, :write, err))
        rethrow()
    end
    return path
end
