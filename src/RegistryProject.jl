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
    # Transient scan state: each source file is parsed once during interpretation and the result is
    # reused by the stats pass, so a file is never read twice in one scan. Cleared when stats finish.
    read_cache::Dict{String,DataFrame}
    read_lock::ReentrantLock
end

# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------

"""Name of a callback, or "λ" for an anonymous function."""
_callback_name(f)::String = (n = string(nameof(f)); startswith(n, "#") ? "λ" : n)

"""Count with a pluralized noun, e.g. `1 plot`, `3 plots`."""
_plural(n::Integer, word::AbstractString)::String = "$n $word" * (n == 1 ? "" : "s")

Base.show(io::IO, project::RegistryProject) =
    print(io, "RegistryProject(\"$(project.name)\", ", _plural(length(project.recipes), "measurement"), ")")

function Base.show(io::IO, ::MIME"text/plain", project::RegistryProject)
    print(io, "RegistryProject \"$(project.name)\"")
    isempty(project.description) || print(io, " · ", project.description)

    print(io, "\n  ", _plural(length(project.recipes), "measurement"))
    isempty(project.recipes) || print(io, " (detection order):")
    pad = isempty(project.recipes) ? 0 : maximum(length(string(r.kind)) for r in project.recipes)
    for recipe in project.recipes
        optional = [name for (name, f) in
            (("process", recipe.process), ("stats", recipe.stats), ("label", recipe.label))
            if f !== nothing]
        steps = isempty(optional) ? "" : " · " * join(optional, " · ")
        print(io, "\n    ", rstrip(rpad(string(recipe.kind), pad) * steps))
    end

    print(io, "\n  ", _plural(length(project.device_stats), "device stat"))
    isempty(project.device_stats) || print(io, ":")
    for recipe in values(project.device_stats)
        print(io, "\n    (", join(recipe.measurement_kinds, ", "), ") → ", _callback_name(recipe.compute_stats))
    end

    print(io, "\n  ", _plural(length(project.plots), "plot"))
    isempty(project.plots) || print(io, ":")
    for recipe in values(project.plots)
        print(io, "\n    ", recipe.kind, " → \"", recipe.label, "\"")
    end
end

"""Create an empty project. Register measurements/stats/plots into it afterwards."""
function define_project(name::AbstractString; description::AbstractString="")::RegistryProject
    return RegistryProject(
        String(name),
        String(description),
        MeasurementRecipe[],
        Dict{Tuple{Vararg{Symbol}},DeviceStatRecipe}(),
        Dict{Symbol,PlotRecipe}(),
        Dict{String,DataFrame}(),
        ReentrantLock(),
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
    # Keep the parse so the stats pass reuses it instead of re-reading (source files are slow to open).
    lock(project.read_lock) do
        project.read_cache[file.filepath] = data
    end
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
    # process/stats run on worker threads, so callbacks must be pure (no shared mutable state). Each
    # measurement only mutates its own `stats`, the file read-cache is filled once then read-only, and
    # failures are collected under a lock. The cancel callback lives in task-local storage, which spawned
    # tasks do not inherit, so it is captured here and re-established inside each task via `with_cancel`.
    failures = MeasurementAnalysisFailure[]
    failures_lock = ReentrantLock()
    add_failure! = (filepath, id, step, err) -> lock(failures_lock) do
        push!(failures, _step_failure(filepath, id, step, err))
    end
    file_by_path = Dict(file.filepath => file for file in files)
    cancel_requested = get(task_local_storage(), CANCEL_CALLBACK_KEY, nothing)

    todo = Tuple{MeasurementInfo,MeasurementRecipe}[]
    for measurement in measurements
        recipe = _recipe(project, measurement.measurement_kind)
        recipe === nothing && continue
        recipe.stats === nothing && continue
        push!(todo, (measurement, recipe))
    end

    # Reuse the parse from interpretation; only files missing from the cache (e.g. a single-file
    # refresh) are read here, once each, in parallel.
    missing_paths = unique(m.filepath for (m, _) in todo if !haskey(project.read_cache, m.filepath))
    @sync for path in missing_paths
        Base.Threads.@spawn with_cancel(cancel_requested) do
            check_cancel()
            recipe = _recipe(project, detect_kind(project, basename(path)))
            recipe === nothing && return
            file = get(file_by_path, path) do
                index_source_file(path)
            end
            try
                data = recipe.read(file)::DataFrame
                lock(project.read_lock) do; project.read_cache[path] = data end
            catch
                # leave it absent; the per-measurement loop records a read failure for it
            end
        end
    end

    # Run process + stats per measurement, in parallel over the now read-only file cache.
    total = length(todo)
    done = Base.Threads.Atomic{Int}(0)
    Base.Threads.@threads for index in eachindex(todo)
        measurement, recipe = todo[index]
        with_cancel(cancel_requested) do
            check_cancel()
            data = get(project.read_cache, measurement.filepath, nothing)
            if data === nothing
                add_failure!(measurement.filepath, measurement.unique_id, :read,
                    ErrorException("could not read source file"))
            else
                try
                    processed = recipe.process === nothing ? data : recipe.process(measurement, data)
                    merge!(measurement.stats, recipe.stats(measurement, processed))
                catch err
                    is_job_cancelled(err) && rethrow()
                    add_failure!(measurement.filepath, measurement.unique_id, :stats, err)
                end
            end
        end
        if on_progress !== nothing
            processed = Base.Threads.atomic_add!(done, 1) + 1
            (processed % 128 == 0 || processed == total) &&
                on_progress((total=total, processed=processed, current_path=measurement.filepath))
        end
    end

    # 3. Cross-measurement device-stat folds run after per-measurement stats are in place. Groups are
    #    disjoint, so they fold in parallel; separate stat recipes stay sequential.
    for recipe in values(project.device_stats)
        check_cancel()
        @sync for group in values(_group_for_device_stat(measurements, recipe))
            isempty(group) && continue
            Base.Threads.@spawn with_cancel(cancel_requested) do
                check_cancel()
                try
                    recipe.compute_stats(group)
                catch err
                    is_job_cancelled(err) && rethrow()
                    add_failure!(group[1].filepath, group[1].unique_id, :compute_stats, err)
                end
            end
        end
    end

    # Release the per-scan file parses now that every stat is computed.
    lock(project.read_lock) do
        empty!(project.read_cache)
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
