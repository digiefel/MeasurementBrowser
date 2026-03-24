"""
FigureScripts.jl - Reusable helpers for publication script generation
"""

"""
    NamedMeasurementGroup(name, measurement_ids)

Plain data container describing one named ordered dataset in a generated figure script.

- `name::String` becomes the lookup key in `data[name]`
- `measurement_ids::Vector{String}` stores exact `MeasurementInfo.id` values in script order
"""
struct NamedMeasurementGroup
    name::String
    measurement_ids::Vector{String}
end

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

function _validate_named_measurement_groups(groups::Vector{NamedMeasurementGroup})
    isempty(groups) && throw(FigureScriptValidationError("At least one measurement group is required"))

    seen_names = Set{String}()
    for group in groups
        name = strip(group.name)
        isempty(name) && throw(FigureScriptValidationError("Measurement group names cannot be empty"))
        name in seen_names && throw(FigureScriptValidationError("Duplicate measurement group name '$name'"))
        push!(seen_names, name)

        isempty(group.measurement_ids) && throw(FigureScriptValidationError("Measurement group '$name' is empty"))
        seen_ids = Set{String}()
        for measurement_id in group.measurement_ids
            id = strip(measurement_id)
            isempty(id) && throw(FigureScriptValidationError("Measurement group '$name' contains an empty measurement id"))
            id in seen_ids && throw(FigureScriptValidationError("Measurement group '$name' contains duplicate measurement id '$id'"))
            push!(seen_ids, id)
        end
    end
    return nothing
end

function _measurement_source_file_id(measurement_id::AbstractString)
    id = strip(String(measurement_id))
    isempty(id) && throw(FigureScriptResolutionError("Measurement id cannot be empty"))
    return first(split(id, '#'; limit=2))
end

function _path_within_root(root_path::AbstractString, path::AbstractString)
    root_abs = normpath(abspath(String(root_path)))
    path_abs = normpath(abspath(String(path)))
    root_parts = splitpath(root_abs)
    path_parts = splitpath(path_abs)
    length(path_parts) >= length(root_parts) || return false
    return path_parts[1:length(root_parts)] == root_parts
end

function _validate_source_file_id(root_path::AbstractString, measurement_id::AbstractString)
    source_file_id = _measurement_source_file_id(measurement_id)
    isabspath(source_file_id) || throw(FigureScriptResolutionError(
        "Measurement id '$measurement_id' must use an absolute source file path",
    ))
    _path_within_root(root_path, source_file_id) || throw(FigureScriptResolutionError(
        "Measurement id '$measurement_id' is outside root path '$root_path'",
    ))
    isfile(source_file_id) || throw(FigureScriptResolutionError(
        "Measurement source file does not exist: $source_file_id",
    ))
    return source_file_id
end

function _resolve_measurement_lookup(
    root_path::AbstractString,
    project::AbstractProject,
    measurement_ids::Vector{String};
    should_cancel::Union{Nothing,Function}=nothing,
)
    requested_ids = unique(copy(measurement_ids))
    meta = _load_scan_metadata(String(root_path))
    resolved = Dict{String,MeasurementInfo}()

    for source_file_id in unique(_validate_source_file_id(root_path, id) for id in requested_ids)
        indexed = index_csv_file(source_file_id)
        items = interpret_measurements(project, indexed, meta; should_cancel=should_cancel)
        isempty(items) && throw(FigureScriptResolutionError(
            "Measurement source file '$source_file_id' did not produce any measurements for project $(project_name(project))",
        ))
        for item in items
            measurement = _measurement_info_from_item(item)
            resolved[measurement.id] = measurement
        end
    end

    missing_ids = String[]
    for measurement_id in requested_ids
        haskey(resolved, measurement_id) || push!(missing_ids, measurement_id)
    end
    isempty(missing_ids) || throw(FigureScriptResolutionError(
        "Could not resolve measurement id(s): $(join(sort!(missing_ids), ", "))",
    ))
    return resolved
end

"""
    prepare_measurement_groups(root_path, project, groups; should_cancel=nothing)

Resolve, load, and analyze exact measurement ids for publication scripting.
Returns a dictionary keyed by group name with group order preserved.
"""
function prepare_measurement_groups(
    root_path::AbstractString,
    project::AbstractProject,
    groups::Vector{NamedMeasurementGroup};
    should_cancel::Union{Nothing,Function}=nothing,
)
    _validate_named_measurement_groups(groups)
    all_ids = String[]
    for group in groups
        append!(all_ids, group.measurement_ids)
    end

    measurement_lookup = _resolve_measurement_lookup(root_path, project, all_ids; should_cancel=should_cancel)
    data = Dict{String,Vector{NamedTuple}}()
    for group in groups
        data[group.name] = map(group.measurement_ids) do measurement_id
                measurement = measurement_lookup[measurement_id]
                plot_params = merge(
                    deepcopy(measurement.device_info.parameters),
                    deepcopy(measurement.parameters),
                )
                loaded = load_plot_for_file(
                    project,
                    measurement.filepath,
                    measurement.measurement_kind;
                    device_params=plot_params,
                    should_cancel=should_cancel,
                )
                analyzed = analyze_plot_for_file(
                    project,
                    measurement.measurement_kind,
                    loaded;
                    device_params=plot_params,
                    should_cancel=should_cancel,
                )
                (
                    measurement=measurement,
                    plot_params=plot_params,
                    loaded=loaded,
                    analyzed=analyzed,
                )
        end
    end
    return data
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

function _render_figure_script(
    root_path::AbstractString,
    project::AbstractProject,
    groups::Vector{NamedMeasurementGroup},
    measurement_lookup::AbstractDict{String,<:MeasurementInfo},
)
    _validate_named_measurement_groups(groups)
    root_abs = normpath(abspath(String(root_path)))
    project_binding = _script_project_binding(project)
    for group in groups
        for measurement_id in group.measurement_ids
            haskey(measurement_lookup, measurement_id) || throw(FigureScriptValidationError(
                "Measurement id '$measurement_id' is not available in the current scan",
            ))
        end
    end

    io = IOBuffer()
    println(io, "using MeasurementBrowser")
    println(io)
    println(io, "# Read this file: ", normpath(joinpath(@__DIR__, "..", "docs", "figure_scripts.md")))
    println(io)
    println(io, "root_path = ", repr(root_abs))
    println(io, "project = ", project_binding)
    println(io)
    println(io, "groups = [")
    for group in groups
        println(io, "    MeasurementBrowser.NamedMeasurementGroup(")
        println(io, "        ", repr(group.name), ",")
        println(io, "        [")
        for measurement_id in group.measurement_ids
            println(io, "            ", repr(measurement_id), ",")
        end
        println(io, "        ],")
        println(io, "    ),")
    end
    println(io, "]")
    println(io)
    println(io, "data = MeasurementBrowser.prepare_measurement_groups(root_path, project, groups)")
    return String(take!(io))
end

"""
    write_figure_script(output_directory, root_path, project, script_name, groups, measurement_lookup; overwrite=false)

Render and write a figure-preparation script.
"""
function write_figure_script(
    output_directory::AbstractString,
    root_path::AbstractString,
    project::AbstractProject,
    script_name::AbstractString,
    groups::Vector{NamedMeasurementGroup},
    measurement_lookup::AbstractDict{String,<:MeasurementInfo};
    overwrite::Bool=false,
)
    path = figure_script_path(output_directory, script_name)
    if isfile(path) && !overwrite
        throw(FigureScriptExistsError(path))
    end

    script_contents = _render_figure_script(root_path, project, groups, measurement_lookup)
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

function _append_group_measurements(group::NamedMeasurementGroup, measurement_ids::Vector{String})
    merged_ids = copy(group.measurement_ids)
    for measurement_id in measurement_ids
        measurement_id in merged_ids || push!(merged_ids, measurement_id)
    end
    return NamedMeasurementGroup(group.name, merged_ids)
end

function _remove_group_measurements(group::NamedMeasurementGroup, measurement_ids::Vector{String})
    to_remove = Set(measurement_ids)
    remaining_ids = [measurement_id for measurement_id in group.measurement_ids if measurement_id ∉ to_remove]
    return NamedMeasurementGroup(group.name, remaining_ids)
end
