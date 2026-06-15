"""
Project.jl - methods of the registration-based project API.

The `Project` type and its recipe types are defined in the `Projects` module (early, so other
submodules can name the concrete type). This file defines the methods: `define_project`, the
`register_*` functions, and the project's implementation of the engine interface (interpretation,
data access, plotting, serialization).

A project is built by mutation with `register_measurement!` / `register_device_stat!` /
`register_plot!`, each pointing at plain callbacks, which the engine then drives:

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
import Serialization

# The project is reachable from a cached SourceScan (twice: source.project and hierarchy.project), so
# it is serialized into the HDF5 cache. Only the registered recipes are persisted; transient scan
# state (read cache, profiling, locks) is rebuilt empty on load. This keeps the cache format stable
# when transient fields change, so adding scan instrumentation never silently invalidates a cache.
# serialize_cycle/deserialize_cycle preserve object identity, so the two references still resolve to
# one project after load.
function Serialization.serialize(s::Serialization.AbstractSerializer, project::Project)
    Serialization.serialize_cycle(s, project) && return nothing
    Serialization.serialize_type(s, Project, true)
    Serialization.serialize(s, project.name)
    Serialization.serialize(s, project.description)
    Serialization.serialize(s, project.recipes)
    Serialization.serialize(s, project.device_stats)
    Serialization.serialize(s, project.plots)
    return nothing
end

function Serialization.deserialize(s::Serialization.AbstractSerializer, ::Type{Project})
    project = Project(
        "", "", MeasurementRecipe[],
        Dict{Tuple{Vararg{Symbol}},DeviceStatRecipe}(), Dict{Symbol,PlotRecipe}(),
        Tuple{String,String,String}[], ReentrantLock(),
        Dict{Symbol,KindProfile}(), ReentrantLock(),
    )
    Serialization.deserialize_cycle(s, project)
    project.name = Serialization.deserialize(s)
    project.description = Serialization.deserialize(s)
    project.recipes = Serialization.deserialize(s)
    project.device_stats = Serialization.deserialize(s)
    project.plots = Serialization.deserialize(s)
    return project
end

# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------

"""Name of a callback, or "λ" for an anonymous function."""
_callback_name(f)::String = (n = string(nameof(f)); startswith(n, "#") ? "λ" : n)

"""Count with a pluralized noun, e.g. `1 plot`, `3 plots`."""
_plural(n::Integer, word::AbstractString)::String = "$n $word" * (n == 1 ? "" : "s")

Base.show(io::IO, project::Project) =
    print(io, "Project(\"$(project.name)\", ", _plural(length(project.recipes), "measurement"), ")")

function Base.show(io::IO, ::MIME"text/plain", project::Project)
    print(io, "Project \"$(project.name)\"")
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
function define_project(name::AbstractString; description::AbstractString="")::Project
    return Project(
        String(name),
        String(description),
        MeasurementRecipe[],
        Dict{Tuple{Vararg{Symbol}},DeviceStatRecipe}(),
        Dict{Symbol,PlotRecipe}(),
        Tuple{String,String,String}[],
        ReentrantLock(),
        Dict{Symbol,KindProfile}(),
        ReentrantLock(),
    )
end

"""Accumulate one timed file read (and the measurements it expanded to) into the scan profile."""
function _record_read!(
    project::Project,
    kind::Symbol,
    seconds::Float64,
    measurements::Int,
)::Nothing
    lock(project.profile_lock) do
        entry = get!(KindProfile, project.scan_profile, kind)
        entry.files += 1
        entry.read_seconds += seconds
        entry.measurements += measurements
    end
    return nothing
end

"""Accumulate one timed process+stats computation into the per-kind scan profile."""
function _record_stats!(project::Project, kind::Symbol, seconds::Float64)::Nothing
    lock(project.profile_lock) do
        entry = get!(KindProfile, project.scan_profile, kind)
        entry.stats_seconds += seconds
    end
    return nothing
end

function reset_scan_profile!(project::Project)::Nothing
    lock(project.profile_lock) do
        empty!(project.scan_profile)
    end
    lock(project.stat_failures_lock) do
        empty!(project.stat_failures)
    end
    return nothing
end

"""Log and record one per-measurement analysis failure raised while interpreting a file."""
function _record_stat_failure!(
    project::Project,
    filepath::String,
    measurement_id::String,
    step::Symbol,
    err,
)::Nothing
    @error("Measurement analysis failed", file=filepath, measurement=measurement_id, step,
        exception=(err, catch_backtrace()))
    lock(project.stat_failures_lock) do
        push!(project.stat_failures, (filepath, measurement_id, "step=$step: " * sprint(showerror, err)))
    end
    return nothing
end

function scan_profile_summary(project::Project)::Vector{NamedTuple}
    rows = lock(project.profile_lock) do
        [(; kind, files=p.files, measurements=p.measurements,
            read_seconds=p.read_seconds, stats_seconds=p.stats_seconds)
         for (kind, p) in project.scan_profile]
    end
    sort!(rows; by=row -> row.read_seconds + row.stats_seconds, rev=true)
    return rows
end

"""
Register (or replace) the recipe for one measurement kind.

`detect`, `read`, and `measurements` are required. `measurements` builds the `MeasurementInfo`s
(including their device identity); a single-measurement file just returns a vector of one. Re-calling
with the same `kind` replaces the recipe in place, so editing and re-running the line updates a live
project.
"""
function register_measurement!(
    project::Project,
    kind::Symbol;
    detect::Function,
    read::Function,
    measurements::Function,
    process::Union{Nothing,Function}=nothing,
    stats::Union{Nothing,Function}=nothing,
    label::Union{Nothing,Function}=nothing,
)::Project
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
    project::Project;
    measurement_kinds::Vector{Symbol},
    compute_stats::Function,
    group_by::Function=_default_device_group,
)::Project
    key = Tuple(sort(measurement_kinds))
    project.device_stats[key] = DeviceStatRecipe(measurement_kinds, group_by, compute_stats)
    return project
end

"""
Register (or replace) the plot recipe for one measurement kind.

`setup(workspace, measurements, processed_data)` builds and returns the `Figure`;
`draw(workspace, measurements, processed_data, figure)` fills it in. The package resolves
`processed_data` (a `Vector{DataFrame}`, one per measurement and in the same order, through its
processed-data cache) before calling either callback, so neither calls `process_measurement_data`
itself.
"""
function register_plot!(
    project::Project,
    kind::Symbol;
    label::AbstractString,
    setup::Function,
    draw::Function,
)::Project
    project.plots[kind] = PlotRecipe(kind, String(label), setup, draw)
    return project
end

_default_device_group(mi::MeasurementInfo)::String = device_path_key(mi.device_info)

# ---------------------------------------------------------------------------
# Recipe lookup helpers
# ---------------------------------------------------------------------------

_recipe(project::Project, kind::Symbol)::Union{Nothing,MeasurementRecipe} =
    (i = findfirst(r -> r.kind === kind, project.recipes); i === nothing ? nothing : project.recipes[i])

function _detect_recipe(project::Project, file::SourceFile)::Union{Nothing,MeasurementRecipe}
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
# Engine interface implementation
# ---------------------------------------------------------------------------

project_name(project::Project)::String = project.name
project_description(project::Project)::String = project.description

function detect_kind(project::Project, filename::String)::Symbol
    recipe = _detect_recipe(project, index_source_file(filename))
    return recipe === nothing ? :unknown : recipe.kind
end

function parse_device_info(::Project, file::SourceFile)::DeviceInfo
    error("Project sets device identity inside `measurements`, not via parse_device_info ($(file.filepath))")
end

function kind_label(project::Project, kind::Symbol)::String
    plot = get(project.plots, kind, nothing)
    return plot === nothing ? string(kind) : plot.label
end

function display_label(project::Project, measurement::MeasurementInfo)::String
    recipe = _recipe(project, measurement.measurement_kind)
    (recipe === nothing || recipe.label === nothing) && return _default_label(measurement)
    return recipe.label(measurement)::String
end

function _default_label(measurement::MeasurementInfo)::String
    parts = Any[measurement.timestamp, string(measurement.measurement_kind)]
    return join(filter(!isnothing, parts), " ")
end

"""
Interpret one file: read it, enumerate its measurements, and compute their per-measurement stats
while the parse is still in hand.

Folding stats into this single pass frees the (often large) parsed table as soon as the file's
measurements are analyzed, instead of holding every file's parse until a separate stats pass — so
scan memory stays bounded to the files in flight rather than the whole tree. Per-measurement failures
are recorded on the project and drained by `compute_and_add_measurement_stats!`.
"""
function interpret_file(project::Project, file::SourceFile)::Vector{MeasurementInfo}
    recipe = _detect_recipe(project, file)
    recipe === nothing && return MeasurementInfo[]
    read_started = time_ns()
    data = recipe.read(file)::DataFrame
    read_seconds = (time_ns() - read_started) / 1e9
    enumerated = recipe.measurements(file, data)::Vector{MeasurementInfo}
    measurements = [_stamp(recipe, mi) for mi in enumerated]
    _record_read!(project, recipe.kind, read_seconds, length(measurements))
    if recipe.stats !== nothing
        for measurement in measurements
            try
                stats_started = time_ns()
                processed = recipe.process === nothing ? data : recipe.process(measurement, data)
                merge!(measurement.stats, recipe.stats(measurement, processed))
                _record_stats!(project, measurement.measurement_kind,
                    (time_ns() - stats_started) / 1e9)
            catch err
                is_job_cancelled(err) && rethrow()
                _record_stat_failure!(project, measurement.filepath, measurement.unique_id, :stats, err)
            end
        end
    end
    # `data` falls out of scope here, so the parse becomes collectable once this file is done.
    return measurements
end

"""Re-stamp a user-built measurement with the recipe kind and an engine-minted id."""
_stamp(recipe::MeasurementRecipe, mi::MeasurementInfo)::MeasurementInfo = MeasurementInfo(
    mi;
    measurement_kind=recipe.kind,
    unique_id=_mint_id(mi.filepath, recipe.kind, mi.parameters),
)

"""Direct data for a measurement is the whole-file read; `process` slices it later."""
function load_source_data(
    project::Project,
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
    workspace::Workspace.Workspace{Project},
    measurement::MeasurementInfo,
)::DataFrame
    recipe = _recipe(workspace.project, measurement.measurement_kind)
    recipe === nothing && error("No registered measurement for kind $(measurement.measurement_kind)")
    raw = only(read_measurement_data(workspace, [measurement]))
    recipe.process === nothing && return raw
    return recipe.process(measurement, raw)::DataFrame
end

"""
Run the cross-measurement device stats, then return every analysis failure for the scan.

Per-measurement process/stats already ran inside `interpret_file` (one pass over each file's parse),
so this only does the disjoint device-scoped folds, which read the measurements' computed `stats`
rather than any raw data. Per-measurement failures collected during interpretation are drained and
returned alongside any device-fold failures.
"""
function compute_and_add_measurement_stats!(
    project::Project,
    measurements::Vector{MeasurementInfo},
    files::Vector{SourceFile};
    on_progress::Union{Nothing,Function}=nothing,
)::Vector{MeasurementAnalysisFailure}
    # The cancel callback lives in task-local storage, which spawned tasks do not inherit, so capture
    # it here and re-establish it inside each task via `with_cancel`.
    failures = MeasurementAnalysisFailure[]
    failures_lock = ReentrantLock()
    add_failure! = (filepath, id, step, err) -> lock(failures_lock) do
        push!(failures, _step_failure(filepath, id, step, err))
    end
    cancel_requested = get(task_local_storage(), CANCEL_CALLBACK_KEY, nothing)

    # Cross-measurement device-stat folds. Groups within one recipe are disjoint, so they fold in
    # parallel; separate stat recipes stay sequential.
    recipes = collect(values(project.device_stats))
    total = length(recipes)
    done = Base.Threads.Atomic{Int}(0)
    for recipe in recipes
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
        if on_progress !== nothing
            processed = Base.Threads.atomic_add!(done, 1) + 1
            on_progress((total=total, processed=processed, current_path=""))
        end
    end

    # Drain the per-measurement failures gathered while interpreting files.
    lock(project.stat_failures_lock) do
        for (filepath, measurement_id, message) in project.stat_failures
            push!(failures, MeasurementAnalysisFailure(filepath, measurement_id, message))
        end
        empty!(project.stat_failures)
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
