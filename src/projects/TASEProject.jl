"""
TASEProject.jl - Project dispatch methods for TASE four-terminal IV measurements

Filename pattern:
  {Chip}_{Facet}_{DeviceType}_{DeviceID}_{YYYYMMDD}_{HHMMSS}_{T}K_FourTerminalIV.csv

Example:
  TASESNS1c1f_A_2TSNJunction_11_20260224_111623_298K_FourTerminalIV.csv
"""

using GLMakie
using DataPlotter: load_tase_four_terminal_iv

const REGEX_TASE = r"^([^_]+)_([^_]+)_([^_]+)_(\d+)_\d{8}_\d{6}_\d+K_FourTerminalIV\.csv"i

# ---------------------------------------------------------------------------
# Interface implementations
# ---------------------------------------------------------------------------

project_name(::TASEProject) = "TASE"
project_description(::TASEProject) = "GaN TASE four-terminal IV"

accepts_file(::TASEProject, filename::String) = match(REGEX_TASE, filename) !== nothing

function parse_device_info(::TASEProject, filename::String)
    m = match(REGEX_TASE, filename)
    m === nothing && error("Unrecognized TASE filename format: $filename")
    chip, facet, device_type, device_id = String.(m.captures)
    return DeviceInfo([chip, facet, device_type, device_id])
end

detect_kind(::TASEProject, filename::String)::Symbol =
    match(REGEX_TASE, filename) !== nothing ? :four_terminal_iv : :unknown

function kind_label(::TASEProject, kind::Symbol)::String
    kind === :four_terminal_iv && return "Four-Terminal I-V"
    return "Unknown"
end

function display_label(proj::TASEProject, meas::MeasurementInfo)
    label = kind_label(proj, meas.measurement_kind)
    temp = get(meas.parameters, :temperature_K, nothing)
    parts = Any[meas.timestamp, label]
    temp !== nothing && push!(parts, "$(temp)K")
    return join(parts, " ")
end

expand_measurement(::TASEProject, meas::MeasurementInfo) = [meas]

# ---------------------------------------------------------------------------
# Plot dispatch
# ---------------------------------------------------------------------------

function load_plot_input_for_file(::TASEProject, path::AbstractString, kind::Union{Symbol,Nothing}; kwargs...)
    kind === :four_terminal_iv || return nothing
    return load_tase_four_terminal_iv(path)
end

function draw_plot_from_input(::TASEProject, kind::Union{Symbol,Nothing}, loaded; kwargs...)
    kind === :four_terminal_iv || return nothing
    loaded === nothing && return nothing

    try
        return _tase_plot_four_terminal_iv(loaded.df, loaded.title)
    catch err
        @warn "draw_plot_from_input (TASE) failed" error = err
        return nothing
    end
end

function figure_for_file(proj::TASEProject, path::AbstractString, kind::Union{Symbol,Nothing}; kwargs...)
    try
        loaded = load_plot_input_for_file(proj, path, kind; kwargs...)
        return draw_plot_from_input(proj, kind, loaded; kwargs...)
    catch err
        @warn "figure_for_file (TASE) failed" path kind error = err
        return nothing
    end
end

figure_for_files(::TASEProject, paths, combined_kind; kwargs...) = nothing

combined_plot_types(::TASEProject) = [(nothing, "None", "No combined plot")]

compatible_kinds(::TASEProject, ::Symbol) = Symbol[]

# ---------------------------------------------------------------------------
# Plot implementation
# ---------------------------------------------------------------------------

function _tase_plot_four_terminal_iv(df, title_str="Four-Terminal I-V")
    nrow(df) == 0 && return nothing

    current_uA = df.current_source .* 1e6   # A → µA
    voltage_mV = df.voltage_drop .* 1e3      # V → mV

    fig = Figure(size=(700, 500))
    ax = Axis(fig[1, 1],
        xlabel="Current (µA)",
        ylabel="Voltage Drop (mV)",
        title=title_str)
    lines!(ax, current_uA, voltage_mV, linewidth=2)
    scatter!(ax, current_uA, voltage_mV, markersize=4)
    return fig
end

# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------

const TASE_PROJECT = TASEProject()
push!(KNOWN_PROJECTS, TASE_PROJECT)
