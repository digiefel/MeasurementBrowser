module DataPlotter

using GLMakie
using DataFrames
using Statistics

using DataAnalysis: analyze_breakdown, analyze_pund, extract_tlm_geometry_from_params, analyze_tlm_combined, calculate_sheet_resistance, analyze_pund_fatigue_combined
using DataLoader: read_iv_sweep, read_fe_pund, read_tlm_4p, read_wakeup, read_pund_fatigue

include("PUND.jl")
include("TLM.jl")

export plot_fe_pund, plot_wakeup, plot_pund_fatigue
export plot_tlm_4p, plot_tlm_combined, plot_tlm_temperature
export plot_iv_sweep_single
export load_fe_pund_single, load_iv_sweep_single, load_tlm_4p_single, load_wakeup_single
export load_tase_four_terminal_iv
export load_tlm_analysis_combined, draw_tlm_analysis_combined

# Extract temperature in Kelvin from device params or filename
function _extract_temperature_K(params::Dict{Symbol,Any}, filepath::String)
    if haskey(params, :temperature_K)
        return try
            Float64(params[:temperature_K])
        catch
            NaN
        end
    end

    base = basename(filepath)
    if (m = match(r"(\d+(?:\.\d+)?)K", base)) !== nothing
        return try
            parse(Float64, m.captures[1])
        catch
            NaN
        end
    end
    return NaN
end

# Extract a site/chip label from params or filename
function _extract_site_label(params::Dict{Symbol,Any}, filepath::String)
    for key in (:site, :chip, :device_id, :die, :location)
        if haskey(params, key)
            return string(params[key])
        end
    end
    # Fallback: basename up to first date/time or temp token
    name_part = replace(basename(filepath), r"\.(csv|txt)$" => "")
    tokens = split(name_part, '_')
    device_tokens = String[]
    for t in tokens
        # stop before date/time, temperature, or geometry token
        if occursin(r"^\d{6,}$", t) || occursin(r"^\d{4}-\d{2}-\d{2}$", t) || occursin(r"^\d+K$", t) || occursin(r"^L\d+W\d+$", t)
            break
        end
        push!(device_tokens, t)
    end
    return isempty(device_tokens) ? name_part : join(device_tokens, "_")
end

# Extract oxygen percentage from params or filename (optional)
function _extract_oxygen_percent(params::Dict{Symbol,Any}, filepath::String)
    for key in (:oxygen_percent, :o2_percent)
        if haskey(params, key)
            return try
                Float64(params[key])
            catch
                NaN
            end
        end
    end

    base = basename(filepath)
    if (m = match(r"(\d+(?:\.\d+)?)%O2", base)) !== nothing
        return try
            parse(Float64, m.captures[1])
        catch
            NaN
        end
    end
    return NaN
end

# Extract oxygen flow in sccm from params or filename (optional)
function _extract_oxygen_flow(params::Dict{Symbol,Any}, filepath::String)
    for key in (:oxygen_flow_sccm, :o2_flow_sccm, :oxygen_flow, :o2_flow)
        if haskey(params, key)
            return try
                Float64(params[key])
            catch
                NaN
            end
        end
    end

    base = basename(filepath)
    if (m = match(r"(\d+(?:\.\d+)?)\s*sccm", base)) !== nothing
        return try
            parse(Float64, m.captures[1])
        catch
            NaN
        end
    end
    return NaN
end

"""
Simple ordinary least squares to fit y = a + b*x

Returns (a, b, success)
"""
ols_ab(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}) = begin
    n = length(x)
    sx = sum(x)
    sy = sum(y)
    sxx = sum(x .^ 2)
    sxy = sum(x .* y)
    denom = n * sxx - sx^2
    if !isfinite(denom) || abs(denom) < 1e-12
        (NaN, NaN, false)
    else
        b = (n * sxy - sx * sy) / denom
        a = (sy - b * sx) / n
        (a, b, true)
    end
end

function _plot_title(path::AbstractString)
    return strip(replace(basename(String(path)), r"\.csv$" => ""))
end

function load_iv_sweep_single(path::AbstractString)
    isfile(path) || return nothing
    df = read_iv_sweep(basename(path), dirname(path))
    return (df=df, title=_plot_title(path))
end

function load_tlm_4p_single(path::AbstractString)
    isfile(path) || return nothing
    df = read_tlm_4p(basename(path), dirname(path))
    return (df=df, title=_plot_title(path))
end

function load_wakeup_single(path::AbstractString)
    isfile(path) || return nothing
    df = read_wakeup(basename(path), dirname(path))
    return (df=df, title=_plot_title(path))
end

function load_fe_pund_single(path::AbstractString; device_params=nothing)
    isfile(path) || return nothing
    params = device_params isa Dict ? device_params : Dict{Symbol,Any}()
    area_um2 = get(params, :area_um2, nothing)
    title = _plot_title(path)
    fatigue_cycle = get(params, :fatigue_cycle, nothing)

    if fatigue_cycle !== nothing
        all_cycles = read_pund_fatigue(basename(path), dirname(path))
        df = get(all_cycles, Int(fatigue_cycle), DataFrame())
        return (df=df, title=title * " cycle $fatigue_cycle (fatigue)", area_um2=area_um2)
    end

    df = read_fe_pund(basename(path), dirname(path))
    return (df=df, title=title, area_um2=area_um2)
end

function load_tase_four_terminal_iv(path::AbstractString)
    # TASE uses the same raw reader as TLM 4-point files.
    return load_tlm_4p_single(path)
end


"""
Plot I-V sweep data for a single DataFrame
"""
function plot_iv_sweep_single(df, title_str="I-V Sweep"; kwargs...)
    if nrow(df) == 0
        return nothing
    end

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], xlabel="Voltage (V)", ylabel="Current (A)", title=title_str)

    lines!(ax, df.v, abs.(df.i), color=df.i .> 0, colormap=:RdBu_3, linewidth=2)

    ax.yscale = log10

    return fig
end

# """
# Plot I-V sweep data for a vector of files, grouped by device.
# """
# function plot_iv_sweep(files::Vector{Tuple{String,DataFrame}})
#     device_data = Dict{String,Vector{Tuple{String,DataFrame}}}()
#
#     # Group files by device
#     for (file, df) in files
#         m = match(r"RuO2test_A2_([^_\[]+)_([^_\[\(]+)", file)
#         device = m !== nothing ? "$(m.captures[1])_$(m.captures[2])" : "Unknown Device"
#         if !haskey(device_data, device)
#             device_data[device] = Vector{Tuple{String,DataFrame}}()
#         end
#         push!(device_data[device], (file, df))
#     end
#
#     plot_iv_sweep_by_device(device_data)
# end
#
# """
# Plot I-V sweep data grouped by device in separate windows
# """
# function plot_iv_sweep_by_device(device_data::Dict{String,Vector{Tuple{String,DataFrame}}})
#     for (device, files) in sort(collect(device_data))
#         fig = Figure(size=(1000, 600))
#
#         ax = Axis(fig[1, 1], xlabel="Voltage (V)", ylabel="Current (A)",
#             title="I-V Sweep for $(device)")
#
#         # Create gradient of colors using the viridis colormap
#         nfiles = length(files)
#         color_gradient = cgrad(:matter, nfiles, categorical=true)
#
#         # Observable toggle for absolute/raw current
#         abs_mode = Observable(true)
#
#         # Plot each file's data
#         for (j, (file_label, df)) in enumerate(files)
#             label = clean_title(file_label)
#             lines!(ax, df.voltage, @lift(if $abs_mode
#                     abs.(df.current1)
#                 else
#                     df.current1
#                 end),
#                 color=color_gradient[j], label=label)
#         end
#
#         axislegend(ax, position=:rt)
#
#         # Add toggle button for absolute/raw current
#         gl = GridLayout(fig[2, 1], tellwidth=false)
#         Label(gl[1, 1], "Absolute Value")
#         toggle = Toggle(gl[1, 2], active=abs_mode[])
#         on(toggle.active) do val
#             abs_mode[] = val
#         end
#
#         # add toggle button for log scale
#         Label(gl[2, 1], "Log Scale")
#         toggle_log = Toggle(gl[2, 2], active=false)
#         on(toggle_log.active) do val
#             abs_mode[] = val ? true : abs_mode[]
#             ax.yscale = val ? log10 : identity
#         end
#
#         display_new_window(fig)
#     end
# end


end # module DataPlotter
