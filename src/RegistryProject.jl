"""
RegistryProject.jl - the registration-based project API.

A project is built by mutation with `register_measurement!` / `register_device_stat!` /
`register_plot!`, each pointing at plain callbacks. `RegistryProject` then satisfies the existing
`AbstractProject` contract by dispatching to those callbacks, so all package scanning, caching,
workspace, and GUI machinery is reused unchanged.

Pipeline per measurement kind (fixed order, each fed the previous output):
    detect(file)             -> Bool
    read(file)               -> DataFrame            # whole file, parsed once
    measurements(file, data) -> Vector{MeasurementInfo}   # one entry per measurement (a single one is a vector of one)
    process(mi, data)        -> DataFrame            # optional; default passthrough
    stats(mi, processed)     -> Dict{Symbol,Any}     # optional
    label(mi)                -> String               # optional

`file` is a `SourceFile` (has `.filename`, `.filepath`, `.timestamp`).
"""

using DataFrames: DataFrame

"""One registered measurement recipe."""
mutable struct MeasurementRecipe
    kind::Symbol
    detect::Function
    read::Function
    measurements::Function
    process::Union{Nothing,Function}
    stats::Union{Nothing,Function}
    label::Union{Nothing,Function}
end

"""One registered cross-measurement (device-scoped) stat."""
struct DeviceStatRecipe
    measurement_kinds::Vector{Symbol}
    group_by::Function
    compute_stats::Function
end

"""One registered plot recipe for a measurement kind."""
struct PlotRecipe
    kind::Symbol
    label::String
    setup::Function
    draw::Function
end

"""A project assembled from registered recipes."""
mutable struct RegistryProject <: AbstractProject
    name::String
    description::String
    recipes::Vector{MeasurementRecipe}
    device_stats::Dict{Tuple{Vararg{Symbol}},DeviceStatRecipe}
    plots::Dict{Symbol,PlotRecipe}
end

"""Create an empty project. Register measurements/stats/plots into it afterwards."""
function define_project(name::AbstractString; description::AbstractString="")::RegistryProject
    return RegistryProject(
        String(name),
        String(description),
        MeasurementRecipe[],
        Dict{Tuple{Vararg{Symbol}},DeviceStatRecipe}(),
        Dict{Symbol,PlotRecipe}(),
    )
end

"""
Register (or replace) the recipe for one measurement kind.

`detect`, `read`, and `measurements` are required. `measurements` builds the `MeasurementInfo`s
(including their device identity); a single-measurement file just returns a vector of one. Re-calling
with the same `kind` replaces the recipe in place, so editing and re-running the line updates a live
project.
"""
function register_measurement!(
    project::RegistryProject,
    kind::Symbol;
    detect::Function,
    read::Function,
    measurements::Function,
    process::Union{Nothing,Function}=nothing,
    stats::Union{Nothing,Function}=nothing,
    label::Union{Nothing,Function}=nothing,
)::RegistryProject
    recipe = MeasurementRecipe(kind, detect, read, measurements, process, stats, label)
    index = findfirst(r -> r.kind === kind, project.recipes)
    index === nothing ? push!(project.recipes, recipe) : (project.recipes[index] = recipe)
    return project
end

"""
Register (or replace) a cross-measurement stat folded over each device's measurements.

`compute_stats` receives the group `Vector{MeasurementInfo}` and fills each `mi.stats` in place.
The engine groups by `group_by` (default: device path) and runs this after the full scan. Re-calling
with the same `measurement_kinds` replaces it.
"""
function register_device_stat!(
    project::RegistryProject;
    measurement_kinds::Vector{Symbol},
    compute_stats::Function,
    group_by::Function=_default_device_group,
)::RegistryProject
    key = Tuple(sort(measurement_kinds))
    project.device_stats[key] = DeviceStatRecipe(measurement_kinds, group_by, compute_stats)
    return project
end

"""Register (or replace) the plot recipe for one measurement kind."""
function register_plot!(
    project::RegistryProject,
    kind::Symbol;
    label::AbstractString,
    setup::Function,
    draw::Function,
)::RegistryProject
    project.plots[kind] = PlotRecipe(kind, String(label), setup, draw)
    return project
end

_default_device_group(mi::MeasurementInfo)::String = device_path_key(mi.device_info)

# ---------------------------------------------------------------------------
# Recipe lookup helpers
# ---------------------------------------------------------------------------

_recipe(project::RegistryProject, kind::Symbol)::Union{Nothing,MeasurementRecipe} =
    (i = findfirst(r -> r.kind === kind, project.recipes); i === nothing ? nothing : project.recipes[i])

function _detect_recipe(project::RegistryProject, file::SourceFile)::Union{Nothing,MeasurementRecipe}
    for recipe in project.recipes
        recipe.detect(file)::Bool && return recipe
    end
    return nothing
end

"""Engine-generated stable identity: file + kind + the params that distinguish siblings."""
function _mint_id(filepath::AbstractString, kind::Symbol, params::Dict{Symbol,Any})::String
    isempty(params) && return "$(filepath)#kind=$(kind)"
    ordered = sort!(collect(keys(params)))
    return "$(filepath)#kind=$(kind)," * join(("$(k)=$(params[k])" for k in ordered), ",")
end

# ---------------------------------------------------------------------------
# AbstractProject contract
# ---------------------------------------------------------------------------

project_name(project::RegistryProject)::String = project.name
project_description(project::RegistryProject)::String = project.description

function detect_kind(project::RegistryProject, filename::String)::Symbol
    recipe = _detect_recipe(project, index_source_file(filename))
    return recipe === nothing ? :unknown : recipe.kind
end

function parse_device_info(::RegistryProject, file::SourceFile)::DeviceInfo
    error("RegistryProject sets device identity inside `measurements`, not via parse_device_info ($(file.filepath))")
end

function kind_label(project::RegistryProject, kind::Symbol)::String
    plot = get(project.plots, kind, nothing)
    return plot === nothing ? string(kind) : plot.label
end

function display_label(project::RegistryProject, measurement::MeasurementInfo)::String
    recipe = _recipe(project, measurement.measurement_kind)
    (recipe === nothing || recipe.label === nothing) && return _default_label(measurement)
    return recipe.label(measurement)::String
end

function _default_label(measurement::MeasurementInfo)::String
    parts = Any[measurement.timestamp, string(measurement.measurement_kind)]
    return join(filter(!isnothing, parts), " ")
end

"""Interpret one file: detect its recipe, then enumerate its logical measurements."""
function interpret_file(project::RegistryProject, file::SourceFile)::Vector{MeasurementInfo}
    recipe = _detect_recipe(project, file)
    recipe === nothing && return MeasurementInfo[]
    data = recipe.read(file)::DataFrame
    enumerated = recipe.measurements(file, data)::Vector{MeasurementInfo}
    return [_stamp(recipe, mi) for mi in enumerated]
end

"""Re-stamp a user-built measurement with the recipe kind and an engine-minted id."""
_stamp(recipe::MeasurementRecipe, mi::MeasurementInfo)::MeasurementInfo = MeasurementInfo(
    mi;
    measurement_kind=recipe.kind,
    unique_id=_mint_id(mi.filepath, recipe.kind, mi.parameters),
)

"""Direct data for a measurement is the whole-file read; `process` slices it later."""
function load_source_data(
    project::RegistryProject,
    source_file::SourceFile;
    measurement::Union{Nothing,MeasurementInfo}=nothing,
)::DataFrame
    kind = measurement === nothing ?
        detect_kind(project, source_file.filename) : measurement.measurement_kind
    recipe = _recipe(project, kind)
    recipe === nothing && error("No registered measurement for kind $kind")
    return recipe.read(source_file)::DataFrame
end

"""Processed data: run the recipe's `process` over the cached whole-file read."""
function process_measurement_data(
    workspace::Workspace.Workspace{RegistryProject},
    measurement::MeasurementInfo,
)::DataFrame
    recipe = _recipe(workspace.project, measurement.measurement_kind)
    recipe === nothing && error("No registered measurement for kind $(measurement.measurement_kind)")
    raw = only(read_measurement_data(workspace, [measurement]))
    recipe.process === nothing && return raw
    return recipe.process(measurement, raw)::DataFrame
end

"""Per-measurement stats, then cross-measurement device stats. Failures carry file + step."""
function compute_and_add_measurement_stats!(
    project::RegistryProject,
    measurements::Vector{MeasurementInfo},
    files::Vector{SourceFile};
    on_progress::Union{Nothing,Function}=nothing,
)::Vector{MeasurementAnalysisFailure}
    failures = MeasurementAnalysisFailure[]
    file_by_path = Dict(file.filepath => file for file in files)
    read_cache = Dict{String,DataFrame}()
    total = count(mi -> _recipe(project, mi.measurement_kind) !== nothing, measurements)
    done = 0

    for measurement in measurements
        check_cancel()
        recipe = _recipe(project, measurement.measurement_kind)
        recipe === nothing && continue
        if recipe.stats !== nothing
            file = get(file_by_path, measurement.filepath) do
                index_source_file(measurement.filepath)
            end
            try
                data = get!(read_cache, measurement.filepath) do
                    recipe.read(file)::DataFrame
                end
                processed = recipe.process === nothing ? data : recipe.process(measurement, data)
                merge!(measurement.stats, recipe.stats(measurement, processed))
            catch err
                is_job_cancelled(err) && rethrow()
                push!(failures, _step_failure(measurement.filepath, measurement.unique_id, :stats, err))
            end
        end
        done += 1
        on_progress === nothing || on_progress((total=total, processed=done, current_path=measurement.filepath))
    end

    for recipe in values(project.device_stats)
        check_cancel()
        for group in values(_group_for_device_stat(measurements, recipe))
            isempty(group) && continue
            try
                recipe.compute_stats(group)
            catch err
                is_job_cancelled(err) && rethrow()
                push!(failures, _step_failure(group[1].filepath, group[1].unique_id, :compute_stats, err))
            end
        end
    end

    return failures
end

function _group_for_device_stat(
    measurements::Vector{MeasurementInfo},
    recipe::DeviceStatRecipe,
)::Dict{Any,Vector{MeasurementInfo}}
    kinds = Set(recipe.measurement_kinds)
    groups = Dict{Any,Vector{MeasurementInfo}}()
    for measurement in measurements
        measurement.measurement_kind in kinds || continue
        push!(get!(groups, recipe.group_by(measurement), MeasurementInfo[]), measurement)
    end
    return groups
end

function _step_failure(filepath::String, measurement_id::String, step::Symbol, err)::MeasurementAnalysisFailure
    @error("Measurement analysis failed", file=filepath, measurement=measurement_id, step, exception=(err, catch_backtrace()))
    return MeasurementAnalysisFailure(filepath, measurement_id, "step=$step: " * sprint(showerror, err))
end
