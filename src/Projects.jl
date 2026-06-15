module Projects

using DataFrames: DataFrame

"""
Base type implemented by measurement projects.

Project values define source interpretation, data processing, and presentation. Package-owned
cache, job, and browser state does not belong in project implementations. A project method may
receive a workspace when it needs package-managed measurement data.
"""
abstract type AbstractProject end

const PROJECTS = AbstractProject[]
const DEFAULT_PROJECT = Ref{Union{AbstractProject,Nothing}}(nothing)

"""Return the stable name used to identify a project implementation."""
function project_name end

"""Return a short human-readable description of a project implementation."""
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
reset_scan_profile!(::AbstractProject)::Nothing = nothing

"""
Per-measurement-kind timing for the most recent scan, newest scan replacing the last.

Returns one `(; kind, files, measurements, read_seconds, stats_seconds)` row per kind, or an
empty vector for projects that do not profile. Surfaced in the performance window.
"""
scan_profile_summary(::AbstractProject)::Vector{NamedTuple} = NamedTuple[]

end
