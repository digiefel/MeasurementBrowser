module Projects

using DataFrames: DataFrame

# ---------------------------------------------------------------------------
# Project type
#
# A project is a value built by registration (define_project + register_*). It is the single project
# type in the package. The struct and its recipe types live here, before MeasurementIndex/Workspace,
# so those modules can name the concrete type directly. The methods that operate on it (registration,
# interpretation, plotting, serialization) are defined later, once MeasurementIndex types exist.
# ---------------------------------------------------------------------------

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

"""Accumulated read/stats timing for one measurement kind in the current scan."""
mutable struct KindProfile
    files::Int
    measurements::Int
    read_seconds::Float64
    stats_seconds::Float64
end
KindProfile() = KindProfile(0, 0, 0.0, 0.0)

"""
A measurement project assembled from registered recipes.

Source interpretation, data processing, and presentation are defined by the registered callbacks.
Package-owned cache, job, and browser state does not belong here.
"""
mutable struct Project
    name::String
    description::String
    recipes::Vector{MeasurementRecipe}
    device_stats::Dict{Tuple{Vararg{Symbol}},DeviceStatRecipe}
    plots::Dict{Symbol,PlotRecipe}
    # Transient scan state: each source file is parsed once during interpretation and the result is
    # reused by the stats pass, so a file is never read twice in one scan. Cleared when stats finish.
    read_cache::Dict{String,DataFrame}
    read_lock::ReentrantLock
    # Transient per-kind timing for the latest scan, surfaced in the performance window. Reset at the
    # start of each scan and replaced wholesale, so it stays bounded to one row per measurement kind.
    scan_profile::Dict{Symbol,KindProfile}
    profile_lock::ReentrantLock
end

const PROJECTS = Project[]
const DEFAULT_PROJECT = Ref{Union{Project,Nothing}}(nothing)

# ---------------------------------------------------------------------------
# Interface functions
#
# Declared here so the package's submodules can call them; the methods are defined later, after the
# types they depend on (MeasurementInfo, SourceFile, Workspace) are available.
# ---------------------------------------------------------------------------

"""Return the stable name used to identify a project."""
function project_name end

"""Return a short human-readable description of a project."""
function project_description end

"""Parse the device represented by a project source file."""
function parse_device_info end

"""Classify the measurement represented by a project source filename."""
function detect_kind end

"""Return the human-readable label for a project measurement kind."""
function kind_label end

"""Return the human-readable label for one logical measurement."""
function display_label end

"""Interpret one indexed source file into its logical measurements."""
function interpret_file end

"""Load the direct table belonging to one source file or logical measurement."""
function load_source_data end

"""Convert direct measurement data into the reusable processed table defined by the project."""
function process_measurement_data end

"""Compute project-specific measurement statistics after the complete scan is known."""
function compute_and_add_measurement_stats! end

"""Return the project-specific display label for one device path."""
function device_path_label end

"""Clear any per-scan timing a project accumulates. Called once at the start of every scan."""
function reset_scan_profile! end

"""
Per-measurement-kind timing for the most recent scan, newest scan replacing the last.

Returns one `(; kind, files, measurements, read_seconds, stats_seconds)` row per kind. Surfaced in
the performance window.
"""
function scan_profile_summary end

end
