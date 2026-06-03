"""
FigureScripts.jl - Reusable helpers for publication script generation
"""

"""
    MeasurementFilterClause(; source_file=nothing, measurement_kind=nothing,
                             device_path=String[],
                             parameter_conditions=Pair{Symbol,Any}[])

One exact filter clause over logical measurements. All populated conditions must match.
"""
struct MeasurementFilterClause
    source_file::Union{Nothing,String}
    measurement_kind::Union{Nothing,Symbol}
    device_path::Vector{String}
    parameter_conditions::Vector{Pair{Symbol,Any}}
end

function MeasurementFilterClause(;
    source_file::Union{Nothing,AbstractString}=nothing,
    measurement_kind::Union{Nothing,Symbol}=nothing,
    device_path::AbstractVector{<:AbstractString}=String[],
    parameter_conditions::AbstractVector{<:Pair}=Pair{Symbol,Any}[],
)
    source = source_file === nothing ? nothing : String(source_file)
    path = String.(device_path)
    normalized_conditions = _normalize_parameter_conditions(parameter_conditions)
    isempty(path) || any(isempty, path) && error("device_path cannot contain empty segments")
    return MeasurementFilterClause(source, measurement_kind, path, normalized_conditions)
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
    return true
end

function _measurement_matches_filter(measurement::MeasurementInfo, parameters::Dict{Symbol,Any}, filter::MeasurementGroupFilter)
    for clause in filter.clauses
        _measurement_matches_clause(measurement, parameters, clause) && return true
    end
    return false
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
        _measurement_matches_filter(measurement, parameter_maps[index], group.filter) && push!(matched, measurement)
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
        indexed = index_source_file(source_file)
        items = interpret_measurements(project, indexed, meta; should_cancel=should_cancel)
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

function _build_figure_measurements(
    project::AbstractProject,
    measurements::Vector{MeasurementInfo};
    should_cancel::Union{Nothing,Function}=nothing,
)
    records = FigureMeasurement[]
    for measurement in measurements
        parameters = _measurement_parameters(measurement)
        df = only(data_of_measurements(project, [measurement]; should_cancel))
        loaded = _ruo2_plot_data(measurement, df)
        analyzed = _analyze_ruo2_file_plot(
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
        (isempty(clause.device_path) ? 0 : 1) +
        length(clause.parameter_conditions)
end

function _is_exact_source_selector_clause(clause::MeasurementFilterClause)
    clause.source_file === nothing && return 0
    clause.measurement_kind === nothing &&
        isempty(clause.device_path) &&
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

struct _FigureScriptFactIndex
    measurements::Vector{MeasurementInfo}
    parameter_maps::Vector{Dict{Symbol,Any}}
    facts_by_measurement::Vector{Vector{Any}}
    fact_matches::Dict{Any,BitVector}
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

function _solution_score(chosen::Vector{_CandidateClause})
    return (
        sum((candidate.condition_count for candidate in chosen); init=0),
        length(chosen),
        sum((candidate.exact_source_selector for candidate in chosen); init=0),
    )
end

function _score_lt(lhs::NTuple{3,Int}, rhs::NTuple{3,Int})
    lhs[1] != rhs[1] && return lhs[1] < rhs[1]
    lhs[2] != rhs[2] && return lhs[2] < rhs[2]
    return lhs[3] < rhs[3]
end

function _score_le(lhs::NTuple{3,Int}, rhs::NTuple{3,Int})
    lhs == rhs && return true
    return _score_lt(lhs, rhs)
end

function _filter_fact_list(measurement::MeasurementInfo, parameters::Dict{Symbol,Any})
    facts = Any[
        (:source_file, measurement.filepath),
        (:measurement_kind, measurement.measurement_kind),
    ]

    path = measurement.device_info.location
    for prefix_length in 1:max(length(path) - 1, 0)
        push!(facts, (:device_prefix, Tuple(path[1:prefix_length])))
    end
    !isempty(path) && push!(facts, (:device_exact, Tuple(path)))

    for condition in _measurement_parameter_conditions(parameters)
        push!(facts, (:parameter, condition.first, condition.second))
    end
    return facts
end

function _build_figure_script_fact_index(
    measurements::Vector{MeasurementInfo},
    profiler::Union{Nothing,_FigureScriptProfileAccumulator}=nothing,
)
    parameter_maps = _profile_section!(profiler, :collect_parameter_maps) do
        _measurement_parameter_maps(measurements)
    end
    facts_by_measurement = Vector{Vector{Any}}(undef, length(measurements))
    fact_matches = Dict{Any,BitVector}()
    _profile_section!(profiler, :build_fact_index) do
        for (index, measurement) in enumerate(measurements)
            facts = _filter_fact_list(measurement, parameter_maps[index])
            facts_by_measurement[index] = facts
            for fact in facts
                matches = get!(fact_matches, fact) do
                    falses(length(measurements))
                end
                matches[index] = true
            end
        end
    end
    return _FigureScriptFactIndex(measurements, parameter_maps, facts_by_measurement, fact_matches)
end

function _clause_key(clause::MeasurementFilterClause)
    return (
        clause.source_file,
        clause.measurement_kind,
        Tuple(clause.device_path),
        Tuple((condition.first, condition.second) for condition in clause.parameter_conditions),
    )
end

function _local_cover(mask::BitVector, selected_indices::Vector{Int})
    return BitVector(mask[index] for index in selected_indices)
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

function _clause_from_facts(facts::Vector{Any})
    source_file = nothing
    measurement_kind = nothing
    device_path = String[]
    parameter_conditions = Pair{Symbol,Any}[]

    for fact in facts
        tag = fact[1]
        if tag == :source_file
            source_file = String(fact[2])
        elseif tag == :measurement_kind
            measurement_kind = fact[2]
        elseif tag == :device_prefix
            device_path = collect(String, fact[2])
        elseif tag == :device_exact
            device_path = collect(String, fact[2])
        elseif tag == :parameter
            push!(parameter_conditions, Symbol(fact[2]) => fact[3])
        else
            error("Unknown figure-script fact '$tag'")
        end
    end

    return MeasurementFilterClause(
        source_file=source_file,
        measurement_kind=measurement_kind,
        device_path=device_path,
        parameter_conditions=parameter_conditions,
    )
end

function _exact_selector_clause(measurement::MeasurementInfo, parameters::Dict{Symbol,Any})
    return MeasurementFilterClause(
        source_file=measurement.filepath,
        measurement_kind=measurement.measurement_kind,
        device_path=measurement.device_info.location,
        parameter_conditions=_measurement_parameter_conditions(parameters),
    )
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

function _candidate_from_mask(
    selected_indices::Vector{Int},
    clause::MeasurementFilterClause,
    mask::BitVector,
)
    cover = _local_cover(mask, selected_indices)
    any(cover) || return nothing
    return _CandidateClause(
        clause,
        cover,
        _clause_condition_count(clause),
        _is_exact_source_selector_clause(clause),
    )
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
    return _candidate_from_mask(selected_indices, clause, mask)
end

function _collect_exact_candidates_by_fact_search(
    index::_FigureScriptFactIndex,
    selected_indices::Vector{Int},
    selected_mask::BitVector,
    profiler::Union{Nothing,_FigureScriptProfileAccumulator}=nothing;
    max_depth::Union{Nothing,Int}=nothing,
    seed_exact_selectors::Bool=false,
)
    candidates = Dict{Any,_CandidateClause}()

    if seed_exact_selectors
        for selected_index in selected_indices
            candidate = _candidate_from_clause(
                index.measurements,
                index.parameter_maps,
                selected_mask,
                selected_indices,
                _exact_selector_clause(index.measurements[selected_index], index.parameter_maps[selected_index]),
                profiler,
            )
            candidate === nothing && throw(FigureScriptResolutionError(
                "Could not build an exact selector for the current figure-script selection",
            ))
            candidates[_clause_key(candidate.clause)] = candidate
        end
    end

    for anchor_index in selected_indices
        _profile_counter!(profiler, :anchor_count)
        anchor_facts = index.facts_by_measurement[anchor_index]

        function search(start_index::Int, depth::Int, current_facts::Vector{Any}, current_mask::Union{Nothing,BitVector})
            for fact_index in start_index:length(anchor_facts)
                _profile_counter!(profiler, :candidate_evaluation_count)
                fact = anchor_facts[fact_index]
                fact_mask = index.fact_matches[fact]
                next_mask = current_mask === nothing ? fact_mask : (current_mask .& fact_mask)
                any(next_mask) || continue
                _count_common_true(next_mask, selected_mask) == 0 && continue

                push!(current_facts, fact)
                if !any(next_mask .& .!selected_mask)
                    clause = _clause_from_facts(current_facts)
                    key = _clause_key(clause)
                    haskey(candidates, key) || (candidates[key] = _candidate_from_mask(selected_indices, clause, next_mask))
                elseif max_depth === nothing || depth < max_depth
                    search(fact_index + 1, depth + 1, current_facts, next_mask)
                end
                pop!(current_facts)
            end
        end

        _profile_section!(profiler, :candidate_search) do
            search(1, 1, Any[], nothing)
        end
    end

    isempty(candidates) && throw(FigureScriptResolutionError(
        "Could not infer any exact figure-script filter clauses from the current selection",
    ))
    _profile_counter!(profiler, :candidate_count, length(candidates))
    return collect(values(candidates))
end

function _cover_superset(lhs::BitVector, rhs::BitVector)
    length(lhs) == length(rhs) || error("Coverage vectors must have equal length")
    for index in eachindex(lhs)
        rhs[index] && !lhs[index] && return false
    end
    return true
end

function _remove_dominated_candidates(candidates::Vector{_CandidateClause})
    keep = trues(length(candidates))
    for i in eachindex(candidates)
        keep[i] || continue
        for j in eachindex(candidates)
            i == j && continue
            keep[j] || continue
            _cover_superset(candidates[j].cover, candidates[i].cover) || continue
            score_i = (
                candidates[i].condition_count,
                candidates[i].exact_source_selector,
            )
            score_j = (
                candidates[j].condition_count,
                candidates[j].exact_source_selector,
            )
            if score_j < score_i || (score_j == score_i && candidates[j].cover != candidates[i].cover)
                keep[i] = false
                break
            end
        end
    end
    return [candidates[index] for index in eachindex(candidates) if keep[index]]
end

function _choose_group_candidates_exact_cover(candidates::Vector{_CandidateClause}, selected_count::Int)
    coverers = [Int[] for _ in 1:selected_count]
    for (candidate_index, candidate) in enumerate(candidates)
        for measurement_index in eachindex(candidate.cover)
            candidate.cover[measurement_index] && push!(coverers[measurement_index], candidate_index)
        end
    end

    any(isempty, coverers) && throw(FigureScriptResolutionError(
        "Could not infer an exact figure-script filter for the current selection",
    ))

    order = sortperm(candidates; by=candidate -> (
        -count(candidate.cover),
        candidate.condition_count,
        candidate.exact_source_selector,
    ))
    best_score = nothing
    best_choice = _CandidateClause[]
    uncovered = trues(selected_count)

    function search(chosen::Vector{_CandidateClause}, uncovered_mask::BitVector)
        if !any(uncovered_mask)
            score = _solution_score(chosen)
            if best_score === nothing || _score_lt(score, best_score)
                best_score = score
                best_choice = copy(chosen)
            end
            return
        end

        partial_score = _solution_score(chosen)
        best_score !== nothing && _score_le(best_score, partial_score) && return

        uncovered_indices = findall(uncovered_mask)
        anchor = uncovered_indices[argmin(length(coverers[index]) for index in uncovered_indices)]
        candidates_for_anchor = sort!(
            copy(coverers[anchor]);
            by=candidate_index -> findfirst(==(candidate_index), order),
        )

        for candidate_index in candidates_for_anchor
            candidate = candidates[candidate_index]
            next_uncovered = copy(uncovered_mask)
            next_uncovered .&= .!candidate.cover
            push!(chosen, candidate)
            search(chosen, next_uncovered)
            pop!(chosen)
        end
    end

    search(_CandidateClause[], uncovered)
    isempty(best_choice) && throw(FigureScriptResolutionError(
        "Could not choose an exact figure-script filter for the current selection",
    ))
    return MeasurementGroupFilter([candidate.clause for candidate in best_choice])
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

    index = _build_figure_script_fact_index(all_measurements, profiler)
    selected_mask = _selected_measurement_mask(selected_indices, length(all_measurements))
    return selected_ids, selected_indices, selected_mask, index
end

function _finalize_group_inference(
    name::AbstractString,
    selected_ids::Set{String},
    measurements::Vector{MeasurementInfo},
    parameter_maps::Vector{Dict{Symbol,Any}},
    filter::MeasurementGroupFilter,
    profiler::Union{Nothing,_FigureScriptProfileAccumulator},
)
    _profile_counter!(profiler, :final_clause_count, length(filter.clauses))
    group = NamedMeasurementGroup(name, filter)
    matched_ids = _profile_section!(profiler, :verify) do
        Set(measurement.unique_id for measurement in _matching_measurements(measurements, parameter_maps, group))
    end
    matched_ids == selected_ids || throw(FigureScriptResolutionError(
        "Inferred figure-script group '$(group.name)' does not reproduce the selected logical measurements exactly",
    ))
    return group
end

function _infer_measurement_group_fact_cover(
    name::AbstractString,
    selected_measurements::Vector{MeasurementInfo},
    all_measurements::Vector{MeasurementInfo},
    profiler::Union{Nothing,_FigureScriptProfileAccumulator},
)
    selected_ids, selected_indices, selected_mask, index = _prepare_group_inference(
        selected_measurements,
        all_measurements,
        profiler,
    )
    candidates = _profile_section!(profiler, :collect_candidates) do
        _collect_exact_candidates_by_fact_search(
            index,
            selected_indices,
            selected_mask,
            profiler;
            max_depth=3,
            seed_exact_selectors=true,
        )
    end
    candidates = _profile_section!(profiler, :prune_candidates) do
        _remove_dominated_candidates(candidates)
    end
    filter = _profile_section!(profiler, :choose_filter) do
        _choose_group_candidates_exact_cover(candidates, length(selected_indices))
    end
    return _finalize_group_inference(name, selected_ids, all_measurements, index.parameter_maps, filter, profiler)
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

function _infer_measurement_group_profiled(
    name::AbstractString,
    selected_measurements::Vector{MeasurementInfo},
    all_measurements::Vector{MeasurementInfo},
)
    profiler = _new_figure_script_profile_accumulator(name, length(selected_measurements), length(all_measurements))
    try
        group = _infer_measurement_group_fact_cover(
            name,
            selected_measurements,
            all_measurements,
            profiler,
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
