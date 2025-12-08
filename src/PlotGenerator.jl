module PlotGenerator

using GLMakie
using DataFrames
using Statistics
using PrecompileTools: @setup_workload, @compile_workload

include("DataLoader.jl")
using .DataLoader: read_iv_sweep, read_fe_pund, read_tlm_4p, read_wakeup
include("DataAnalysis.jl")
using .Analysis: analyze_breakdown, analyze_pund, extract_tlm_geometry_from_params, analyze_tlm_combined, calculate_sheet_resistance, analyze_pund_fatigue_combined

export figure_for_file, figure_for_files, plot_tlm_combined, plot_pund_fatigue, get_combined_plot_types

"""
    figure_for_file(path::AbstractString; kind::Union{Symbol,Nothing}=nothing) -> Union{Figure,Nothing}

Given a filepath to a measurement CSV, detect its measurement type
from the filename alone (no MeasurementInfo dependency), load the data
with the appropriate reader, and return a Makie Figure. Returns `nothing`
if unsupported or loading/plotting fails.
"""
function figure_for_file(path::AbstractString, kind::Union{Symbol,Nothing}; kwargs...)
    isfile(path) || return nothing
    fname = basename(path)
    dir = dirname(path)

    # Helper to derive a title (strip .csv)
    title = strip(replace(fname, r"\.csv$" => ""))

    df = nothing
    fig = nothing
    try
        if kind === :pund
            df = read_fe_pund(fname, dir)
            fig = plot_fe_pund(df, title; kwargs...)
        elseif kind === :iv
            df = read_iv_sweep(fname, dir)
            fig = plot_iv_sweep_single(df, title; kwargs...)
        elseif kind === :tlm4p
            df = read_tlm_4p(fname, dir)
            fig = plot_tlm_4p(df, title; kwargs...)
        elseif kind === :breakdown
            # Treat as breakdown I-V for now
            df = read_iv_sweep(fname, dir)
            fig = plot_iv_sweep_single(df, title * " (Breakdown)"; kwargs...)
        elseif kind === :wakeup
            df = read_wakeup(fname, dir)
            fig = plot_wakeup(df, title; kwargs...)
        else
            # Fallback attempt: try I-V sweep reader
            try
                @warn "figure_for_file: Invalid measurement type. Attempting I-V sweep reader."
                df = read_iv_sweep(fname, dir)
                fig = plot_iv_sweep_single(df, title; kwargs...)
            catch
                @warn "figure_for_file: Failed to plot."
                return nothing
            end
        end
    catch err
        @warn "figure_for_file failed" path error = err
        return nothing
    end

    return fig
end

"""
Available combined plot types - easily extensible
"""
function get_combined_plot_types()
    return [
        (nothing, "None", "No combined plot selected"),
        (:tlm_analysis, "TLM Analysis", "Width-normalized resistance vs length from multiple TLM 4-point measurements"),
        (:pund_fatigue, "PUND Fatigue", "P-E curve evolution or remnant polarization vs fatigue cycles from PUND measurements"),
    ]
end

"""
    figure_for_files(paths::Vector{String}, combined_kind::Symbol; device_params_list=[], kwargs...) -> Union{Figure,Nothing}

Given a vector of file paths, a combined plot type, and device parameters,
load the data and return a combined Makie Figure. Returns `nothing`
if unsupported or loading/plotting fails.
"""
function figure_for_files(paths::Vector{String}, combined_kind::Symbol; device_params_list::Vector{Dict{Symbol,Any}}=Dict{Symbol,Any}[], kwargs...)
    isempty(paths) && return nothing

    try
        if combined_kind === :tlm_analysis
            return plot_tlm_combined(paths; device_params_list=device_params_list, kwargs...)
        elseif combined_kind === :pund_fatigue
            return plot_pund_fatigue(paths; device_params_list=device_params_list, kwargs...)
        else
            @warn "figure_for_files: Unknown combined plot kind: $combined_kind"
            return nothing
        end
    catch err
        @warn "figure_for_files failed" paths combined_kind error = err
        return nothing
    end
end


# small OLS helper: y = a + b x
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

"""
Plot combined TLM analysis with a toggle:
- Global fit (one OLS over all widths) OR Per-width fits (one OLS per width).
Bottom axis always shows the mean of width-normalized resistance vs length.

Requires: calculate_sheet_resistance(::DataFrame) -> (R_sheet, R_cprime, rho_c, R2)
"""
function plot_tlm_combined(paths::Vector{String}; device_params_list::Vector{Dict{Symbol,Any}}=Dict{Symbol,Any}[], kwargs...)
    @info "plot_tlm_combined called with $(length(paths)) files"

    if isempty(paths)
        @warn "No files provided for TLM combined analysis"
        return nothing
    end

    # Load all TLM files
    files_data_params = Tuple{String,DataFrame,Dict{Symbol,Any}}[]
    for (i, path) in enumerate(paths)
        try
            filename = basename(path)
            dirname_path = dirname(path)
            df = read_tlm_4p(filename, dirname_path)
            device_params = (i <= length(device_params_list)) ? device_params_list[i] : Dict{Symbol,Any}()
            push!(files_data_params, (path, df, device_params))
        catch err
            @warn "Failed to load TLM file: $path" error = err
        end
    end

    if isempty(files_data_params)
        @warn "No valid TLM files could be loaded"
        return nothing
    end

    # Perform combined analysis
    analysis_df = analyze_tlm_combined(files_data_params)

    if nrow(analysis_df) == 0
        @warn "TLM combined analysis produced no data"
        return nothing
    end

    # Width-invariant fit → (R_sheet, R_c', rho_c, R2)
    R_sheet, R_cprime, rho_c, r_squared = calculate_sheet_resistance(analysis_df)

    # Create the plot with three rows: main, residuals (shorter), and secondary+info
    fig = Figure(size=(900, 600))

    # Main plot: Resistance vs Length/Width
    ax1 = Axis(fig[1, 1:2],
        xlabel="Length/Width Ratio (L/W)",
        ylabel="Resistance (kΩ)",
        title="TLM Analysis - Model from (R□, R_c')")

    # Residuals plot (underneath main)
    ax_res = Axis(fig[2, 1:2],
        xlabel="Length/Width Ratio (L/W)",
        ylabel="Residuals (%)",
        title="Residuals vs linear model")
    linkxaxes!(ax1, ax_res)

    # Secondary plot: Width-normalized resistance vs length
    ax2 = Axis(fig[3, 1],
        xlabel="Length (μm)",
        ylabel="Width-Normalized Resistance (Ω·μm)",
        title="Width-Normalized Resistance vs Length")

    info_ax = Axis(fig[3, 2], xlabel="", ylabel="", title="Analysis Results")
    hidedecorations!(info_ax)
    hidespines!(info_ax)

    # Make the residuals row shorter than the main plot
    rowsize!(fig.layout, 1, Relative(0.55))
    rowsize!(fig.layout, 2, Relative(0.20))
    rowsize!(fig.layout, 3, Relative(0.25))

    # Geometry-averaged points for R vs L/W
    geometry_groups = combine(groupby(analysis_df, [:length_um, :width_um]),
        :resistance_ohm => (x -> mean(filter(isfinite, x))) => :avg_resistance_ohm)

    if nrow(geometry_groups) == 0
        @warn "No geometry-averaged data available"
        return fig
    end

    widths = unique(geometry_groups.width_um)
    colors = cgrad(:tab10, length(widths), categorical=true)

    # Scatter points by width (main axis)
    for (i, w_um) in enumerate(widths)
        width_geom = filter(row -> row.width_um == w_um, geometry_groups)
        if nrow(width_geom) > 0
            lw_vals = width_geom.length_um ./ width_geom.width_um
            scatter!(ax1, lw_vals, width_geom.avg_resistance_ohm ./ 1e3,
                color=colors[i], label="$(w_um) μm", markersize=10)
        end
    end

    # Predicted lines from inferred (R_sheet, R_c') for each width
    if isfinite(R_sheet) && isfinite(R_cprime)
        for (i, w_um) in enumerate(widths)
            mask = geometry_groups.width_um .== w_um
            if any(mask)
                w_cm = w_um * 1e-4
                a_w = 2 * (R_cprime / w_cm)                 # Ω
                xspan = geometry_groups.length_um[mask] ./ geometry_groups.width_um[mask]
                x_max = maximum(xspan)
                lw_fit = range(0, max(x_max, 1.0), length=200)
                r_fit = (a_w .+ R_sheet .* lw_fit) ./ 1e3   # kΩ
                lines!(ax1, lw_fit, r_fit, color=colors[i], linewidth=2, linestyle=:dash)
            end
        end
    end

    # Residuals vs the same inferred model
    if isfinite(R_sheet) && isfinite(R_cprime)
        xvals = geometry_groups.length_um ./ geometry_groups.width_um
        W_cm = geometry_groups.width_um .* 1e-4
        pred_ohm = (2 .* (R_cprime ./ W_cm)) .+ R_sheet .* xvals
        resid_rel = 1e2 .* (geometry_groups.avg_resistance_ohm .- pred_ohm) ./ geometry_groups.avg_resistance_ohm

        xext = extrema(xvals)
        lines!(ax_res, [xext[1], xext[2]], [0.0, 0.0], color=:black, linewidth=1, linestyle=:dot)
        for (x, r) in zip(xvals, resid_rel)
            lines!(ax_res, [x, x], [0.0, r], color=:blue, linewidth=1)
        end
        scatter!(ax_res, xvals, resid_rel, color=:red, markersize=6)
    end

    # Bottom axis: keep original behavior (per-width points)
    for (i, w_um) in enumerate(widths)
        width_data = filter(row -> row.width_um == w_um, analysis_df)
        scatter!(ax2, width_data.length_um, width_data.resistance_normalized,
            color=colors[i], label="$(w_um) μm", markersize=6)
    end

    # Info panel (values inferred from the fit function)
    results_text = ""
    if isfinite(R_sheet)
        results_text *= "Sheet Resistance (R□):  $(round(R_sheet, digits=2)) Ω/□\n\n"
        if isfinite(R_cprime)
            results_text *= "Contact per width (R_c'):  $(round(R_cprime, sigdigits=4)) Ω·cm\n\n"
        end
        if isfinite(rho_c)
            results_text *= "Specific Contact (ρ_c):  $(round(rho_c, sigdigits=4)) Ω·cm²\n\n"
        end
        results_text *= "R² (width-invariant fit):  $(round(r_squared, digits=4))\n\n"
    else
        results_text *= "Width-invariant fit unavailable (need ≥2 geometries)\n\n"
    end
    results_text *= "Files analyzed: $(length(files_data_params))\n"
    results_text *= "Widths: $(length(widths))\n"
    results_text *= "Data points: $(nrow(analysis_df))"

    text!(info_ax, 0.05, 0.95, text=results_text, align=(:left, :top), fontsize=15, space=:relative)
    xlims!(info_ax, 0, 1)
    ylims!(info_ax, 0, 1)

    if length(widths) > 1
        axislegend(ax1, position=:lt)
    end

    return fig
end

"""
Plot PUND fatigue analysis with:
- Top: overlapped P–E (or Q–V) curves color-coded by cumulative fatigue cycles (log-scaled colorbar)
- Bottom: remnant polarization (Pr) vs fatigue cycles
Includes a "Last only" toggle button (inside the figure) to show only the last repetition per FE_PUND file.
"""
function plot_pund_fatigue(paths::Vector{String}; device_params_list::Vector{Dict{Symbol,Any}}=Dict{Symbol,Any}[], mode=:overlapped, kwargs...)
    isempty(paths) && return nothing

    # Prepare entries for analysis (chronological is enforced in analysis)
    n = length(paths)
    device_params_list = length(device_params_list) == n ? device_params_list : [Dict{Symbol,Any}() for _ in 1:n]
    entries = NamedTuple[]
    for i in 1:n
        path = paths[i]
        params = device_params_list[i]
        fname = lowercase(basename(path))
        dirpath = dirname(path)
        ts = stat(path).mtime
        if occursin("wakeup", fname)
            df_w = read_wakeup(basename(path), dirpath)
            push!(entries, (kind=:wakeup, df=df_w, params=params, timestamp=ts))
        else
            df_p = read_fe_pund(basename(path), dirpath)
            push!(entries, (kind=:pund, df=df_p, params=params, timestamp=ts))
        end
    end

    # Delegate analysis (per-repetition cycles, chronological accumulation, Pr extraction)
    res = analyze_pund_fatigue_combined(entries)
    traces = get(res, :traces, NamedTuple[])
    x_label = get(res, :x_label, "Voltage (V)")
    y_label = get(res, :y_label, "Polarization (μC/cm²)")
    pr_points = get(res, :pr_points, NamedTuple[])

    # Nothing to show
    isempty(traces) && isempty(pr_points) && return nothing

    # Figure layout: top plot + control area on the right, bottom plot spans width
    fig = Figure(size=(1200, 800))
    ax_top = Axis(fig[1, 1], xlabel=x_label, ylabel=y_label,
        title="Overlapped P–E curves (color-coded by fatigue cycles)")
    # Controls + colorbar column
    gl_ctrl = GridLayout(fig[1, 2])
    ax_bottom = Axis(fig[2, 1:2], xlabel="Fatigue cycles", ylabel="Remnant polarization (μC/cm²)",
        title="Remnant polarization vs fatigue cycles")

    # Toggle for "Last only"
    last_only = Observable(false)
    btn = Button(gl_ctrl[1, 1], label="Last only: OFF")
    on(btn.clicks) do _
        last_only[] = !last_only[]
        btn.label[] = last_only[] ? "Last only: ON" : "Last only: OFF"
    end

    # Map cycles -> color using log scale; replace legend with a colorbar
    cyc_vals = [t.cycles for t in traces if hasproperty(t, :cycles)]
    if isempty(cyc_vals)
        # Fallback text if there are no cycles
        text!(ax_top, 0.5, 0.5, text="No valid P–E traces", align=(:center, :center), color=:gray)
    else
        # Log10 scaling; ensure limits are valid
        logvals = log10.(Float64.(cyc_vals))
        vmin = minimum(logvals)
        vmax = maximum(logvals)
        # Nice tick labels for colorbar (1, 10, 100, ...)
        min_pow = floor(Int, vmin)
        max_pow = ceil(Int, vmax)
        tick_positions = collect(min_pow:max_pow)
        tick_labels = string.(Int.(10 .^ tick_positions))

        # Colorbar (log ticks)
        Colorbar(gl_ctrl[2, 1], colormap=:viridis, limits=(vmin, vmax), label="Fatigue cycles (log10)",
            ticks=(tick_positions, tick_labels))

        # Plot all traces once; color by log(cycles) and toggle visibility
        plotted_traces = NamedTuple{(:plt, :is_last)}[]
        for t in traces
            isempty(t.x) && continue
            isempty(t.y) && continue
            t_log = log10(float(t.cycles))
            plt = lines!(ax_top, t.x, t.y;
                color=t_log, colormap=:viridis, colorrange=(vmin, vmax), linewidth=2)
            is_last = hasproperty(t, :rep_index) && hasproperty(t, :rep_count) && (t.rep_index == t.rep_count)
            push!(plotted_traces, (plt=plt, is_last=is_last))
        end
        # Initialize visibility (show all)
        for p in plotted_traces
            p.plt.visible[] = true
        end
        on(last_only) do v
            for p in plotted_traces
                p.plt.visible[] = (!v) || p.is_last
            end
        end
    end

    # Bottom panel: Pr vs cycles (all vs last-only)
    if !isempty(pr_points)
        all_x = Float64[p.cycles for p in pr_points]
        all_y = Float64[p.Pr for p in pr_points]
        ord_all = sortperm(all_x)
        lines_all = lines!(ax_bottom, all_x[ord_all], all_y[ord_all], color=:black, linewidth=2, alpha=0.7)
        scatter_all = scatter!(ax_bottom, all_x, all_y, color=:red, markersize=8)

        last_x = Float64[]
        last_y = Float64[]
        for p in pr_points
            if hasproperty(p, :rep_index) && hasproperty(p, :rep_count) && (p.rep_index == p.rep_count)
                push!(last_x, p.cycles)
                push!(last_y, p.Pr)
            end
        end
        ord_last = sortperm(last_x)
        lines_last = isempty(last_x) ? nothing : lines!(ax_bottom, last_x[ord_last], last_y[ord_last], color=:black, linewidth=2, alpha=0.7)
        scatter_last = isempty(last_x) ? nothing : scatter!(ax_bottom, last_x, last_y, color=:red, markersize=8)

        # Start with "all" visible
        if lines_last !== nothing
            lines_last.visible[] = false
        end
        if scatter_last !== nothing
            scatter_last.visible[] = false
        end
        on(last_only) do v
            lines_all.visible[] = !v
            scatter_all.visible[] = !v
            if lines_last !== nothing
                lines_last.visible[] = v
            end
            if scatter_last !== nothing
                scatter_last.visible[] = v
            end
        end
    else
        text!(ax_bottom, 0.5, 0.5, text="No remnant polarization points", align=(:center, :center), color=:gray)
    end

    return fig
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

"""
Plot FE PUND data with comprehensive visualization
"""
function plot_fe_pund(df, title_str="FE PUND"; area_um2=nothing, DEBUG::Bool=false, kwargs...)
    if DEBUG
        return debug_fe_pund(df, title_str; area_um2, DEBUG=DEBUG, kwargs...)
    end
    if nrow(df) == 0
        return nothing
    end

    df = analyze_pund(df)

    time_us = df.time * 1e6
    fig = Figure(size=(1200, 800))

    # create axes
    ax1 = Axis(fig[1, 1:2], xlabel="Time (μs)", ylabel="Current (μA)", yticklabelcolor=:blue, title="$title_str - Area = $(isnothing(area_um2) ? "?" : round(area_um2, digits=2)) um²")
    ax1twin = Axis(fig[1, 1:2], yaxisposition=:right, ylabel="Voltage (V)", yticklabelcolor=:red)
    ax2 = Axis(fig[2, 1], xlabel="Voltage (V)", ylabel="Current (μA)", title="$title_str - I-V Characteristic")
    ax3 = Axis(fig[2, 2], xlabel="Voltage (V)", ylabel="Switching Charge (pC)", title="$title_str - Ferroelectric Switching Charge")
    linkxaxes!(ax1, ax1twin)

    # Combined I, V plot
    l1 = lines!(ax1, time_us, df.current * 1e6, color=:blue, linewidth=2)
    l2 = lines!(ax1twin, time_us, df.voltage, color=:red, linewidth=2, linestyle=:dash)
    l3 = lines!(ax1, time_us, df.I_FE * 1e6, color=:purple, linewidth=2)

    # Current vs Voltage (hysteresis loop)
    lines!(ax2, df.voltage, df.current * 1e6, color=:green, linewidth=2)
    lines!(ax2, df.voltage, df.I_FE * 1e6, color=:purple, linewidth=2)

    # Remant charge - align P and N pulses
    Q_FE = df.Q_FE .- mean(filter(!isnan, df.Q_FE))
    # Remnant polarization - divide by area (um² -> cm²)
    if !isnothing(area_um2)
        area_cm2 = area_um2 / 1e8  # convert from um² to cm²
        label = "Remnant Polarization (μC/cm²)"
        ax3.ylabel = label
        P_FE = Q_FE / area_cm2
    end

    # Plot each PUND repetition separately for legend
    for rep in 1:maximum(df.pulse_idx)÷5
        pulse_range = (rep-1)*5+1:rep*5
        mask = [p in pulse_range for p in df.pulse_idx]
        if any(mask)
            if !isnothing(area_um2)
                lines!(ax3, df.voltage[mask], P_FE[mask] * 1e6, linewidth=2, color=:purple, label="$rep")
            else
                lines!(ax3, df.voltage[mask], Q_FE[mask] * 1e12, linewidth=2, color=:purple, label="$rep")
            end
        end
    end

    # legends
    Legend(fig[1, 1], [l1, l2, l3], ["Current", "Voltage", "FE Current"], tellwidth=false, tellheight=false, halign=:left, valign=:top)
    axislegend(ax3)

    return fig
end

function debug_fe_pund(df, title_str="FE PUND"; area_um2=nothing, DEBUG::Bool=false, kwargs...)
    if nrow(df) == 0
        return nothing
    end

    # Run analysis with optional debug logging inside analyze_pund
    df = analyze_pund(df; DEBUG=DEBUG)

    time_us = df.time * 1e6

    # Debug view: emphasize time-domain only with very clear pulse boundaries
    fig = Figure(size=(1200, 800))

    # Time-domain I and V on shared x, dual y-axes
    ax1 = Axis(fig[1, 1:2], xlabel="Time (μs)", ylabel="Current (μA)",
        yticklabelcolor=:blue, title="$title_str - DEBUG (time domain)")
    ax1twin = Axis(fig[1, 1:2], yaxisposition=:right, ylabel="Voltage (V)", yticklabelcolor=:red)
    linkxaxes!(ax1, ax1twin)

    lI = lines!(ax1, time_us, df.current * 1e6, color=:blue, linewidth=2)
    lIFE = lines!(ax1, time_us, df.I_FE * 1e6, color=:purple, linewidth=2)
    lV = lines!(ax1twin, time_us, df.voltage, color=:red, linewidth=1, linestyle=:dash)

    # Q_FE vs time (only non-NaN points present for P and N)
    ax2 = Axis(fig[2, 1:2], xlabel="Time (μs)", ylabel="Switching Charge (pC)",
        title="$title_str - Q_FE(t)")
    lines!(ax2, time_us, df.Q_FE * 1e12, color=:orange, linewidth=2)

    # Helper to draw vertical boundaries where pulse index changes
    pid = df.pulse_idx
    if length(pid) > 1
        # Boundaries at any change in pulse_idx (including in/out of zero)
        boundaries = [i for i in 2:length(pid) if pid[i] != pid[i-1]]
        for b in boundaries
            t = time_us[b]
            if isfinite(t)
                vlines!(ax1, t, color=:black, linestyle=:dot, linewidth=1, alpha=0.5)
                vlines!(ax2, t, color=:black, linestyle=:dot, linewidth=1, alpha=0.5)
            end
        end

        # Label segments with P/U/N/D/Poling labels at segment centers
        function pulse_label(p::Int)
            r = p % 5
            r == 1 && return "Poling"
            r == 2 && return "P"
            r == 3 && return "U"
            r == 4 && return "N"
            r == 0 && return "D"
            return ""
        end
        # Build contiguous segments of constant pulse_idx
        segs = UnitRange{Int}[]
        s = 1
        for i in eachindex(pid)[2:end]
            if pid[i] != pid[i-1]
                push!(segs, s:(i-1))
                s = i
            end
        end
        push!(segs, s:length(pid))
        segs = [r for r in segs if pid[first(r)] > 0]  # only label actual pulses

        # Place labels near the top of current axis and charge axis (robust to NaNs)
        yI = df.current * 1e6
        finite_yI = filter(isfinite, yI)
        yImaxabs = isempty(finite_yI) ? 1.0 : maximum(abs.(finite_yI))
        yImaxabs = isfinite(yImaxabs) ? yImaxabs : 1.0
        yImin = isempty(finite_yI) ? -yImaxabs : minimum(finite_yI)
        yImaxv = isempty(finite_yI) ? yImaxabs : maximum(finite_yI)
        if !(yImin < yImaxv)
            yImin, yImaxv = -yImaxabs, yImaxabs
        end

        yQ = df.Q_FE * 1e12
        finite_yQ = filter(isfinite, yQ)
        yQmaxabs = isempty(finite_yQ) ? 1.0 : maximum(abs.(finite_yQ))
        yQmaxabs = isfinite(yQmaxabs) ? yQmaxabs : 1.0
        yQmin = isempty(finite_yQ) ? -yQmaxabs : minimum(finite_yQ)
        yQmaxv = isempty(finite_yQ) ? yQmaxabs : maximum(finite_yQ)
        if !(yQmin < yQmaxv)
            yQmin, yQmaxv = -yQmaxabs, yQmaxabs
        end

        ylims!(ax1, yImin - 0.1 * yImaxabs, yImaxv + 0.3 * yImaxabs)
        ylims!(ax2, yQmin - 0.1 * yQmaxabs, yQmaxv + 0.3 * yQmaxabs)

        for r in segs
            tmid = (time_us[first(r)] + time_us[last(r)]) / 2
            lab = pulse_label(pid[first(r)])
            if !isempty(lab) && isfinite(tmid)
                text!(ax1, tmid, yImaxv + 0.25 * yImaxabs; text=lab, align=(:center, :baseline), color=:black)
                text!(ax2, tmid, yQmaxv + 0.25 * yQmaxabs; text=lab, align=(:center, :baseline), color=:black)
            end
        end
    end

    Legend(fig[1, 1], [lI, lV, lIFE], ["Current", "Voltage", "FE Current"],
        tellwidth=false, tellheight=false, halign=:left, valign=:top)

    return fig
end



"""
Plot wakeup data showing pulse count and amplitude as text
"""
function plot_wakeup(df, title_str="Wakeup"; kwargs...)
    if nrow(df) == 0
        return nothing
    end

    pulse_count = df.pulse_count[1]
    amplitude = df.amplitude[1]

    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1], title=title_str)

    # Hide axis elements since we only want to show text
    hidedecorations!(ax)
    hidespines!(ax)

    # Display the pulse count and amplitude as text
    text_content = "$(pulse_count)× wakeup pulses\namplitude = $(amplitude) V"
    text!(ax, 0.5, 0.5, text=text_content, align=(:center, :center),
        fontsize=24, color=:black)

    # Set axis limits to center the text
    xlims!(ax, 0, 1)
    ylims!(ax, 0, 1)

    return fig
end

"""
Plot TLM 4-point data with detailed analysis
"""
function plot_tlm_4p(df, title_str="TLM 4-Point"; device_params=Dict{Symbol,Any}(), kwargs...)
    if nrow(df) == 0
        return nothing
    end

    # Extract data
    I = df.current_source
    V = df.v_gnd
    
    # Filter out NaNs and infinite values
    mask = isfinite.(I) .& isfinite.(V)
    I = I[mask]
    V = V[mask]
    
    if isempty(I)
        return nothing
    end

    # Linear Fit V = R*I + offset
    # Simple linear regression
    n = length(I)
    sx = sum(I)
    sy = sum(V)
    sxx = sum(I .^ 2)
    sxy = sum(I .* V)
    
    denom = n * sxx - sx^2
    
    R_fit = 0.0
    offset = 0.0
    
    if abs(denom) > 1e-20
        R_fit = (n * sxy - sx * sy) / denom
        offset = (sy - R_fit * sx) / n
    end
    
    # Calculate Resistance for plotting
    # Avoid division by zero
    R_meas = V ./ I
    
    # Extract Geometry
    L_um, W_um = extract_tlm_geometry_from_params(device_params, title_str)
    
    # Calculate Sheet Resistivity
    rho_sheet = NaN
    if !isnan(L_um) && !isnan(W_um) && L_um > 0
        rho_sheet = R_fit * W_um / L_um
    end
    
    # Parse Title Info
    # Expected format: Chip_Type_Subtype_Geometry_Date_Time_Temp_...
    # e.g. A11_VI_TLM_L100W2_20251205_151649_298K_FourTerminalIV
    parts = split(title_str, ['_', ' '])
    
    chip = length(parts) >= 1 ? parts[1] : "?"
    geometry_str = "?"
    temp = "?"
    
    # Try to find geometry in parts
    for p in parts
        if occursin(r"L\d+W\d+", p)
            geometry_str = p
        end
    end
    
    # Use temperature from device_params if available
    if haskey(device_params, :temperature_K)
        temp = "$(device_params[:temperature_K])K"
    else
        # Fallback to parsing title
        for p in parts
            if occursin(r"\d+K", p)
                temp = p
            end
        end
    end
    
    new_title = "$chip $geometry_str $temp"
    
    # Plotting
    fig = Figure(size=(1000, 500))
    
    # Scale units
    I_uA = I .* 1e6
    V_mV = V .* 1e3
    R_kOhm = R_meas ./ 1e3
    R_fit_kOhm = R_fit / 1e3
    
    # Panel 1: I-V
    ax1 = Axis(fig[1, 1], xlabel="Current (μA)", ylabel="Voltage (mV)", title="I-V Curve")
    scatter!(ax1, I_uA, V_mV, color=:blue, markersize=8, label="Data")
    
    # Plot fit line
    I_min, I_max = minimum(I_uA), maximum(I_uA)
    # Add some padding
    I_span = I_max - I_min
    if I_span == 0
        I_span = 1.0
    end
    I_range = [I_min - 0.1*I_span, I_max + 0.1*I_span]
    # V = R*I + offset => V_mV = (R_fit * I_uA*1e-6 + offset) * 1e3
    # V_mV = R_fit * I_uA * 1e-3 + offset * 1e3
    # V_mV = (R_fit/1e3) * I_uA + offset*1e3
    V_fit_mV = (R_fit / 1e3) .* I_range .+ (offset * 1e3)
    
    lines!(ax1, I_range, V_fit_mV, color=:red, linewidth=2, label="Fit: R=$(round(R_fit, digits=2)) Ω")
    
    # Add Rho to legend if available
    if !isnan(rho_sheet)
        # Add an invisible line for the legend entry
        lines!(ax1, [NaN], [NaN], color=:transparent, label="ρ_sq = $(round(rho_sheet, digits=2)) Ω/sq")
    end
    
    axislegend(ax1, position=:lt)
    
    # Panel 2: R vs I
    ax2 = Axis(fig[1, 2], xlabel="Current (μA)", ylabel="Resistance (kΩ)", title="Resistance vs Current")
    
    # Filter R_meas for plotting (remove outliers/infinities)
    valid_R = isfinite.(R_kOhm) .& (abs.(R_kOhm) .< 1e6) # Arbitrary large cutoff
    
    if any(valid_R)
        scatter!(ax2, I_uA[valid_R], R_kOhm[valid_R], color=:green, markersize=8, label="R = V/I")
    end
    hlines!(ax2, [R_fit_kOhm], color=:red, linestyle=:dash, linewidth=2, label="Fitted R")
    
    axislegend(ax2)
    
    Label(fig[0, :], new_title, fontsize=20, font=:bold)
    
    return fig
end

@setup_workload begin
    # In-memory synthetic data only (avoid any filesystem interaction)
    iv_df = DataFrame(
        voltage=[0.0, 0.5, 1.0],
        current1=[0.0, 5e-7, 1e-6],
        current2=[0.0, 5.5e-7, 1.1e-6],
    )
    tlm_df = DataFrame(
        current_source=[0.0, 1e-6, 2e-6],
        i1=[0.0, 9e-7, 1.8e-6],
        i2=[0.0, 0.0, 0.0],
        is=[0.0, 0.0, 0.0],
        v_gnd=[0.0, 1e-3, 2e-3],
    )
    fe_df = DataFrame(
        time=[0.0, 1e-6, 2e-6, 3e-6],
        current=[1e-6, 1.1e-6, 1.05e-6, 1.2e-6],
        voltage=[0.0, 0.1, 0.2, 0.0],
        current_time=[0.0, 0.0, 0.0, 0.0],
        voltage_time=[0.0, 0.0, 0.0, 0.0],
    )
    @compile_workload begin
        try
            plot_iv_sweep_single(iv_df, "PC IV")
        catch
        end
        try
            plot_tlm_4p(tlm_df, "PC TLM")
        catch
        end
        try
            plot_fe_pund(fe_df, "PC PUND")
        catch
        end
        # Dispatcher requires a filesystem path; skip to avoid creating files.
    end
end

end # module
