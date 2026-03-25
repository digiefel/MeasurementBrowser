"""
FigureScripts.jl - Reusable helpers for publication script generation
"""

"""
    MeasurementFilterClause(; source_file=nothing, measurement_kind=nothing,
                             device_path=String[], device_path_mode=:exact,
                             parameter_conditions=Pair{Symbol,Any}[])

One exact filter clause over logical measurements. All populated conditions must match.
"""
struct MeasurementFilterClause
    source_file::Union{Nothing,String}
    measurement_kind::Union{Nothing,Symbol}
    device_path_mode::Symbol
    device_path::Vector{String}
    parameter_conditions::Vector{Pair{Symbol,Any}}
end

function MeasurementFilterClause(;
    source_file::Union{Nothing,AbstractString}=nothing,
    measurement_kind::Union{Nothing,Symbol}=nothing,
    device_path::AbstractVector{<:AbstractString}=String[],
    device_path_mode::Symbol=:exact,
    parameter_conditions::AbstractVector{<:Pair}=Pair{Symbol,Any}[],
)
    source = source_file === nothing ? nothing : String(source_file)
    path = String.(device_path)
    mode = isempty(path) ? :none : device_path_mode
    mode in (:none, :exact, :prefix) || error("Invalid device_path_mode '$mode'")
    mode == :none && !isempty(path) && error("device_path_mode=:none is only valid for an empty device_path")
    normalized_conditions = _normalize_parameter_conditions(parameter_conditions)
    isempty(path) || any(isempty, path) && error("device_path cannot contain empty segments")
    return MeasurementFilterClause(source, measurement_kind, mode, path, normalized_conditions)
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
    NamedMeasurementGroup(name, filter)

One fixed named group in a generated figure script.
"""
struct NamedMeasurementGroup
    name::String
    filter::MeasurementGroupFilter
    function NamedMeasurementGroup(name::String, filter::MeasurementGroupFilter)
        normalized_name = strip(name)
        isempty(normalized_name) && throw(FigureScriptValidationError("Measurement group names cannot be empty"))
        return new(normalized_name, filter)
    end
end

function NamedMeasurementGroup(name::AbstractString, filter::MeasurementGroupFilter)
    return NamedMeasurementGroup(String(name), filter)
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

function _script_project_binding(::RuO2Project)
    return "MeasurementBrowser.RUO2_PROJECT"
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

function _measurement_parameters(measurement::MeasurementInfo)
    return merge(deepcopy(measurement.device_info.parameters), deepcopy(measurement.parameters))
end

function _measurement_parameter_conditions(parameters::Dict{Symbol,Any})
    conditions = Pair{Symbol,Any}[]
    for (name, value) in sort!(collect(parameters); by=entry -> String(first(entry)))
        _is_filterable_parameter_value(value) || continue
        push!(conditions, name => value)
    end
    return conditions
end

function _measurement_parameter_conditions(measurement::MeasurementInfo)
    return _measurement_parameter_conditions(_measurement_parameters(measurement))
end

function _measurement_matches_clause(measurement::MeasurementInfo, parameters::Dict{Symbol,Any}, clause::MeasurementFilterClause)
    clause.source_file !== nothing && measurement.filepath != clause.source_file && return false
    clause.measurement_kind !== nothing && measurement.measurement_kind != clause.measurement_kind && return false

    if clause.device_path_mode == :exact
        measurement.device_info.location == clause.device_path || return false
    elseif clause.device_path_mode == :prefix
        length(measurement.device_info.location) >= length(clause.device_path) || return false
        measurement.device_info.location[1:length(clause.device_path)] == clause.device_path || return false
    end

    for condition in clause.parameter_conditions
        haskey(parameters, condition.first) || return false
        isequal(parameters[condition.first], condition.second) || return false
    end
    return true
end

function _measurement_matches_filter(measurement::MeasurementInfo, parameters::Dict{Symbol,Any}, filter::MeasurementGroupFilter)
    for clause in filter.clauses
        _measurement_matches_clause(measurement, parameters, clause) && return true
    end
    return false
end

function _matching_measurements(measurements::Vector{MeasurementInfo}, group::NamedMeasurementGroup)
    matched = MeasurementInfo[]
    for measurement in measurements
        params = _measurement_parameters(measurement)
        _measurement_matches_filter(measurement, params, group.filter) && push!(matched, measurement)
    end
    return matched
end

function _load_measurements_from_source_files(
    root_path::AbstractString,
    project::AbstractProject,
    source_files::Vector{String};
    should_cancel::Union{Nothing,Function}=nothing,
)
    source_files = _normalize_source_files(root_path, source_files)
    meta = _load_scan_metadata(String(root_path))
    measurements = MeasurementInfo[]
    for source_file in source_files
        indexed = index_csv_file(source_file)
        items = interpret_measurements(project, indexed, meta; should_cancel=should_cancel)
        isempty(items) && throw(FigureScriptResolutionError(
            "Measurement source file '$source_file' did not produce any measurements for project $(project_name(project))",
        ))
        append!(measurements, [_measurement_info_from_item(item) for item in items])
    end
    return measurements
end

function _group_matches(groups::Vector{NamedMeasurementGroup}, measurements::Vector{MeasurementInfo})
    matches = Dict{String,Vector{MeasurementInfo}}()
    for group in groups
        matched = _matching_measurements(measurements, group)
        isempty(matched) && throw(FigureScriptResolutionError(
            "Measurement group '$(group.name)' matches no measurements in the current dataset",
        ))
        matches[group.name] = matched
    end
    return matches
end

function _build_figure_measurements(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo};
    should_cancel::Union{Nothing,Function}=nothing,
)
    records = FigureMeasurement[]
    for measurement in measurements
        parameters = _measurement_parameters(measurement)
        loaded = load_plot_for_file(
            project,
            measurement.filepath,
            measurement.measurement_kind;
            device_params=parameters,
            should_cancel=should_cancel,
        )
        analyzed = analyze_plot_for_file(
            project,
            measurement.measurement_kind,
            loaded;
            device_params=parameters,
            should_cancel=should_cancel,
        )
        push!(records, FigureMeasurement(measurement, parameters, loaded, analyzed))
    end
    return records
end

"""
    prepare_figure_script_data(root_path, project, source_files, groups; should_cancel=nothing)

Load the source-file union for a figure script and expose the fixed named groups as typed data.
"""
function prepare_figure_script_data(
    root_path::AbstractString,
    project::AbstractProject,
    source_files::Vector{String},
    groups::Vector{NamedMeasurementGroup};
    should_cancel::Union{Nothing,Function}=nothing,
)
    _validate_named_measurement_groups(groups)
    normalized_source_files = _normalize_source_files(root_path, source_files)
    measurements = _load_measurements_from_source_files(root_path, project, normalized_source_files; should_cancel=should_cancel)
    _group_matches(groups, measurements)
    records = _build_figure_measurements(project, measurements; should_cancel=should_cancel)

    group_map = Dict(group.name => group for group in groups)
    match_indices = Dict{String,Vector{Int}}()
    for group in groups
        indices = Int[]
        for (index, record) in enumerate(records)
            _measurement_matches_filter(record.measurement, record.parameters, group.filter) && push!(indices, index)
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

struct _CandidateClause
    clause::MeasurementFilterClause
    cover::BitVector
    condition_count::Int
    exact_source_selector::Int
end

function _clause_condition_count(clause::MeasurementFilterClause)
    return (clause.source_file === nothing ? 0 : 1) +
        (clause.measurement_kind === nothing ? 0 : 1) +
        (clause.device_path_mode == :none ? 0 : 1) +
        length(clause.parameter_conditions)
end

function _is_exact_source_selector_clause(clause::MeasurementFilterClause)
    clause.source_file === nothing && return 0
    clause.measurement_kind === nothing &&
        clause.device_path_mode == :none &&
        isempty(clause.parameter_conditions) && return 0
    return 1
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

function _clause_key(clause::MeasurementFilterClause)
    return (
        clause.source_file,
        clause.measurement_kind,
        clause.device_path_mode,
        Tuple(clause.device_path),
        Tuple((condition.first, condition.second) for condition in clause.parameter_conditions),
    )
end

function _local_cover(mask::BitVector, selected_indices::Vector{Int})
    return BitVector(mask[index] for index in selected_indices)
end

function _selected_measurement_ids(selected_measurements::Vector{MeasurementInfo})
    selected_ids = Set(measurement.id for measurement in selected_measurements)
    length(selected_ids) == length(selected_measurements) || throw(FigureScriptValidationError(
        "Figure-script group selection contains duplicate logical measurements",
    ))
    return selected_ids
end

function _selected_measurement_indices(selected_ids::Set{String}, all_measurements::Vector{MeasurementInfo})
    selected_indices = findall(measurement -> measurement.id in selected_ids, all_measurements)
    length(selected_indices) == length(selected_ids) || throw(FigureScriptResolutionError(
        "Figure-script group selection references measurements that are not present in the current scan",
    ))
    return selected_indices
end

function _selected_measurement_mask(selected_indices::Vector{Int}, measurement_count::Int)
    selected_mask = falses(measurement_count)
    selected_mask[selected_indices] .= true
    return selected_mask
end

function _count_common_true(left::BitVector, right::BitVector)
    length(left) == length(right) || error("BitVector lengths must match")
    total = 0
    @inbounds for index in eachindex(left, right)
        left[index] && right[index] && (total += 1)
    end
    return total
end

function _clause_relaxations(clause::MeasurementFilterClause)
    relaxations = MeasurementFilterClause[]
    seen = Set{Tuple}()

    function push_relaxation(relaxed::MeasurementFilterClause)
        key = _clause_key(relaxed)
        key in seen && return
        push!(seen, key)
        push!(relaxations, relaxed)
    end

    if clause.source_file !== nothing
        push_relaxation(MeasurementFilterClause(
            measurement_kind=clause.measurement_kind,
            device_path=clause.device_path,
            device_path_mode=clause.device_path_mode,
            parameter_conditions=clause.parameter_conditions,
        ))
    end

    if clause.measurement_kind !== nothing
        push_relaxation(MeasurementFilterClause(
            source_file=clause.source_file,
            device_path=clause.device_path,
            device_path_mode=clause.device_path_mode,
            parameter_conditions=clause.parameter_conditions,
        ))
    end

    if clause.device_path_mode == :exact
        if length(clause.device_path) > 1
            push_relaxation(MeasurementFilterClause(
                source_file=clause.source_file,
                measurement_kind=clause.measurement_kind,
                device_path=clause.device_path[1:end-1],
                device_path_mode=:prefix,
                parameter_conditions=clause.parameter_conditions,
            ))
        elseif length(clause.device_path) == 1
            push_relaxation(MeasurementFilterClause(
                source_file=clause.source_file,
                measurement_kind=clause.measurement_kind,
                parameter_conditions=clause.parameter_conditions,
            ))
        end
    elseif clause.device_path_mode == :prefix
        if length(clause.device_path) > 1
            push_relaxation(MeasurementFilterClause(
                source_file=clause.source_file,
                measurement_kind=clause.measurement_kind,
                device_path=clause.device_path[1:end-1],
                device_path_mode=:prefix,
                parameter_conditions=clause.parameter_conditions,
            ))
        elseif length(clause.device_path) == 1
            push_relaxation(MeasurementFilterClause(
                source_file=clause.source_file,
                measurement_kind=clause.measurement_kind,
                parameter_conditions=clause.parameter_conditions,
            ))
        end
    end

    for index in eachindex(clause.parameter_conditions)
        conditions = copy(clause.parameter_conditions)
        deleteat!(conditions, index)
        push_relaxation(MeasurementFilterClause(
            source_file=clause.source_file,
            measurement_kind=clause.measurement_kind,
            device_path=clause.device_path,
            device_path_mode=clause.device_path_mode,
            parameter_conditions=conditions,
        ))
    end

    return relaxations
end

"""
    _collect_exact_candidates(measurements, selected_indices)

Build a small set of exact clauses by starting from exact selectors and greedily
relaxing them in a fixed order while exactness is preserved.
"""
function _collect_exact_candidates(
    measurements::Vector{MeasurementInfo},
    selected_indices::Vector{Int},
    profiler::Union{Nothing,_FigureScriptProfileAccumulator}=nothing,
)
    parameter_maps = _profile_section!(profiler, :collect_parameter_maps) do
        _measurement_parameter_maps(measurements)
    end
    selected_mask = _selected_measurement_mask(selected_indices, length(measurements))
    remaining_cover_mask = trues(length(selected_indices))
    candidates = Dict{Any,_CandidateClause}()

    while any(remaining_cover_mask)
        _profile_counter!(profiler, :anchor_count)
        anchor_position = findfirst(remaining_cover_mask)
        anchor_position === nothing && break
        anchor_index = selected_indices[anchor_position]
        candidate = _candidate_from_clause(
            measurements,
            parameter_maps,
            selected_mask,
            selected_indices,
            _exact_selector_clause(measurements[anchor_index]),
            profiler,
        )
        candidate === nothing && throw(FigureScriptResolutionError(
            "Could not build an exact selector for the current figure-script selection",
        ))

        while true
            relaxed = false
            relaxed_clauses = _profile_section!(profiler, :candidate_relaxations) do
                _clause_relaxations(candidate.clause)
            end
            _profile_counter!(profiler, :relaxation_clause_count, length(relaxed_clauses))
            for relaxed_clause in relaxed_clauses
                _profile_counter!(profiler, :relaxation_attempt_count)
                relaxed_candidate = _candidate_from_clause(
                    measurements,
                    parameter_maps,
                    selected_mask,
                    selected_indices,
                    relaxed_clause,
                    profiler,
                )
                relaxed_candidate === nothing && continue
                candidate = relaxed_candidate
                _profile_counter!(profiler, :relaxation_accept_count)
                relaxed = true
                break
            end
            relaxed || break
        end

        key = _clause_key(candidate.clause)
        haskey(candidates, key) || (candidates[key] = candidate)
        remaining_cover_mask .&= .!candidate.cover
    end

    isempty(candidates) && throw(FigureScriptResolutionError(
        "Could not infer any exact figure-script filter clauses from the current selection",
    ))
    _profile_counter!(profiler, :candidate_count, length(candidates))
    return collect(values(candidates))
end

function _clause_match_mask(
    measurements::Vector{MeasurementInfo},
    parameter_maps::Vector{Dict{Symbol,Any}},
    clause::MeasurementFilterClause,
)
    length(measurements) == length(parameter_maps) || error("Measurement and parameter map lengths must match")
    mask = falses(length(measurements))
    for index in eachindex(measurements)
        _measurement_matches_clause(measurements[index], parameter_maps[index], clause) && (mask[index] = true)
    end
    return mask
end

function _candidate_from_clause(
    measurements::Vector{MeasurementInfo},
    parameter_maps::Vector{Dict{Symbol,Any}},
    selected_mask::BitVector,
    selected_indices::Vector{Int},
    clause::MeasurementFilterClause,
    profiler::Union{Nothing,_FigureScriptProfileAccumulator}=nothing,
)
    _profile_counter!(profiler, :candidate_evaluation_count)
    mask = _profile_section!(profiler, :candidate_match_mask) do
        _clause_match_mask(measurements, parameter_maps, clause)
    end
    any(mask) || return nothing
    any(mask .& .!selected_mask) && return nothing
    cover = _local_cover(mask, selected_indices)
    any(cover) || return nothing
    return _CandidateClause(
        clause,
        cover,
        _clause_condition_count(clause),
        _is_exact_source_selector_clause(clause),
    )
end

function _prune_redundant_choice(chosen::Vector{_CandidateClause}, selected_count::Int)
    keep = copy(chosen)
    changed = true
    while changed
        changed = false
        for index in eachindex(keep)
            covered = falses(selected_count)
            for other_index in eachindex(keep)
                other_index == index && continue
                covered .|= keep[other_index].cover
            end
            if all(covered)
                deleteat!(keep, index)
                changed = true
                break
            end
        end
    end
    return keep
end

function _choose_group_candidates_greedy(candidates::Vector{_CandidateClause}, selected_count::Int)
    uncovered = trues(selected_count)
    chosen = _CandidateClause[]

    while any(uncovered)
        best_candidate = nothing
        best_score = nothing
        for candidate in candidates
            newly_covered = _count_common_true(candidate.cover, uncovered)
            newly_covered == 0 && continue
            score = (-newly_covered, candidate.condition_count, candidate.exact_source_selector)
            if best_score === nothing || score < best_score
                best_candidate = candidate
                best_score = score
            end
        end

        best_candidate === nothing && throw(FigureScriptResolutionError(
            "Could not greedily choose an exact figure-script filter for the current selection",
        ))
        push!(chosen, best_candidate)
        uncovered .&= .!best_candidate.cover
    end

    chosen = _prune_redundant_choice(chosen, selected_count)
    return MeasurementGroupFilter([candidate.clause for candidate in chosen])
end

function _exact_selector_clause(measurement::MeasurementInfo)
    return MeasurementFilterClause(
        source_file=measurement.filepath,
        measurement_kind=measurement.measurement_kind,
        device_path=measurement.device_info.location,
        device_path_mode=:exact,
        parameter_conditions=_measurement_parameter_conditions(measurement),
    )
end

function infer_measurement_group(
    name::AbstractString,
    selected_measurements::Vector{MeasurementInfo},
    all_measurements::Vector{MeasurementInfo},
)
    try
        group, _ = _infer_measurement_group_profiled(name, selected_measurements, all_measurements)
        return group
    catch err
        err isa FigureScriptProfiledError || rethrow()
        rethrow(err.cause)
    end
end

function _infer_measurement_group_impl(
    name::AbstractString,
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
    candidates = _profile_section!(profiler, :collect_candidates) do
        _collect_exact_candidates(all_measurements, selected_indices, profiler)
    end
    filter = _profile_section!(profiler, :choose_filter) do
        _choose_group_candidates_greedy(candidates, length(selected_indices))
    end
    _profile_counter!(profiler, :final_clause_count, length(filter.clauses))
    group = NamedMeasurementGroup(name, filter)

    matched_ids = _profile_section!(profiler, :verify) do
        Set(measurement.id for measurement in _matching_measurements(all_measurements, group))
    end
    matched_ids == selected_ids || throw(FigureScriptResolutionError(
        "Inferred figure-script group '$(group.name)' does not reproduce the selected logical measurements exactly",
    ))
    return group
end

function _infer_measurement_group_profiled(
    name::AbstractString,
    selected_measurements::Vector{MeasurementInfo},
    all_measurements::Vector{MeasurementInfo},
)
    profiler = _new_figure_script_profile_accumulator(name, length(selected_measurements), length(all_measurements))
    try
        group = _infer_measurement_group_impl(name, selected_measurements, all_measurements, profiler)
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
    if clause.device_path_mode != :none
        println(io, "                device_path = ", repr(clause.device_path), ",")
        println(io, "                device_path_mode = ", repr(clause.device_path_mode), ",")
    end
    if !isempty(clause.parameter_conditions)
        println(io, "                parameter_conditions = [")
        for condition in clause.parameter_conditions
            println(io, "                    ", repr(condition.first), " => ", repr(condition.second), ",")
        end
        println(io, "                ],")
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
    println(io, "        ]),")
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
    println(io)
    println(io, "# Read this file: ", normpath(joinpath(@__DIR__, "..", "docs", "figure_scripts.md")))
    println(io)
    println(io, "root_path = ", repr(root_abs))
    println(io, "project = ", project_binding)
    println(io, "source_files = [")
    for source_file in normalized_source_files
        println(io, "    ", repr(source_file), ",")
    end
    println(io, "]")
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
