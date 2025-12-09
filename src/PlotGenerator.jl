module PlotGenerator

using GLMakie
using DataFrames
using Statistics
using PrecompileTools: @setup_workload, @compile_workload

include("DataLoader.jl")
using .DataLoader: read_iv_sweep, read_fe_pund, read_tlm_4p, read_wakeup
include("DataAnalysis.jl")
using .Analysis: analyze_breakdown, analyze_pund, extract_tlm_geometry_from_params, analyze_tlm_combined, calculate_sheet_resistance, analyze_pund_fatigue_combined

export figure_for_file, figure_for_files, plot_tlm_combined, plot_tlm_temperature, plot_pund_fatigue, get_combined_plot_types

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
        (:tlm_temperature, "TLM vs Temperature", "Sheet resistance vs temperature (groups by site/chip)"),
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
        elseif combined_kind === :tlm_temperature
            return plot_tlm_temperature(paths; device_params_list=device_params_list, kwargs...)
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

"""
Plot sheet resistance vs temperature for TLM measurements.

Groups files by site/chip and temperature, computes sheet/contact parameters
via `calculate_sheet_resistance(analyze_tlm_combined(...))`, then plots
R□ (or thickness-normalized resistivity when thickness is present for all
files at that site/temperature) vs temperature with one series per site/chip.
Adds a linear fit per site to report TCR.
"""
function plot_tlm_temperature(paths::Vector{String}; device_params_list::Vector{Dict{Symbol,Any}}=Dict{Symbol,Any}[], kwargs...)
    isempty(paths) && return nothing

    entries = NamedTuple{(:path, :df, :params, :tempK, :site, :oxygen_percent)}[]

    for (i, path) in enumerate(paths)
        try
            df = read_tlm_4p(basename(path), dirname(path))
            params = i <= length(device_params_list) ? device_params_list[i] : Dict{Symbol,Any}()
            tempK = _extract_temperature_K(params, path)
            if !isfinite(tempK)
                @warn "Skipping TLM file without temperature" path
                continue
            end
            site = _extract_site_label(params, path)
            o2 = _extract_oxygen_percent(params, path)
            push!(entries, (path=path, df=df, params=params, tempK=tempK, site=site, oxygen_percent=o2))
        catch err
            @warn "Failed to load TLM file for temperature plot" path error = err
        end
    end

    isempty(entries) && return nothing

    results = DataFrame(site=String[], temperature_K=Float64[], R_val=Float64[], R_cprime=Float64[], rho_c=Float64[], r_squared=Float64[], n_files=Int[], thickness_cm=Float64[], oxygen_percent=Float64[])

    by_site = Dict{String, Vector{NamedTuple}}()
    for e in entries
        vec = get!(by_site, e.site, Vector{NamedTuple}())
        push!(vec, e)
    end

    for (site, vec) in by_site
        temps = unique([e.tempK for e in vec])
        for temp in temps
            subset = filter(e -> e.tempK == temp, vec)
            files_data_params = [(e.path, e.df, e.params) for e in subset]
            combined_df = analyze_tlm_combined(files_data_params)
            if nrow(combined_df) == 0
                @warn "No valid TLM data for site/temp" site tempK = temp
                continue
            end
            R_sheet, R_cprime, rho_c, r2 = calculate_sheet_resistance(combined_df)
            if !isfinite(R_sheet)
                @warn "Insufficient geometries for sheet resistance" site tempK = temp
                continue
            end

            thickness_values = Float64[]
            for s in subset
                if haskey(s.params, :t_RuO2_nm)
                    t_nm = try
                        Float64(s.params[:t_RuO2_nm])
                    catch
                        NaN
                    end
                    if isfinite(t_nm) && t_nm > 0
                        push!(thickness_values, t_nm)
                    end
                end
            end

            thickness_cm = NaN
            if length(thickness_values) == length(subset) && !isempty(thickness_values)
                thickness_cm = mean(thickness_values) * 1e-7
                if maximum(thickness_values) - minimum(thickness_values) > 1e-3 * maximum(thickness_values)
                    @warn "Thickness varies across files; using mean for normalization" site tempK = temp
                end
            elseif !isempty(thickness_values)
                @warn "Thickness missing for some files; skipping normalization" site tempK = temp
            end

            R_val = isfinite(thickness_cm) ? R_sheet * thickness_cm : R_sheet
            o2 = begin
                vals = filter(isfinite, [s.oxygen_percent for s in subset])
                isempty(vals) ? NaN : mean(vals)
            end
            push!(results, (site, temp, R_val, R_cprime, rho_c, r2, length(subset), thickness_cm, o2))
        end
    end

    nrow(results) == 0 && return nothing

    fig = Figure(size=(1100, 700))
    any_thickness = any(isfinite, results.thickness_cm)
    y_label = any_thickness ? "Resistivity (mΩ·cm)" : "Sheet Resistance (Ω/□)"
    title_str = any_thickness ? "TLM Resistivity vs Temperature" : "TLM Sheet Resistance vs Temperature"
    ax = Axis(fig[1, 1:2], xlabel="Temperature (K)", ylabel=y_label, title=title_str)

    sites = unique(results.site)
    site_o2 = Dict{String,Float64}()
    for s in sites
        vals = filter(isfinite, results.oxygen_percent[results.site .== s])
        site_o2[s] = isempty(vals) ? NaN : mean(vals)
    end
    sites = sort(sites; by = s -> begin
        o2 = site_o2[s]
        isfinite(o2) ? o2 : Inf
    end)

    o2_vals = filter(isfinite, collect(values(site_o2)))
    use_o2_colors = !isempty(o2_vals)
    cmap = to_colormap(Reverse(:seaborn_flare_gradient))
    o2_min, o2_max = use_o2_colors ? (minimum(o2_vals), maximum(o2_vals)) : (0.0, 1.0)
    base_colors = to_colormap(:tab10)

    rt_points = NamedTuple{(:site, :o2, :resistivity)}[]
    tcr_points = NamedTuple{(:site, :o2, :tcr_ppm)}[]

    for (i, site) in enumerate(sites)
        sub = sort(filter(row -> row.site == site, results), :temperature_K)
        label = site
        if any(isfinite, sub.oxygen_percent)
            o2 = mean(filter(isfinite, sub.oxygen_percent))
            label = string(site, " (", round(o2, digits=2), "% O2)")
        end
        color = use_o2_colors && isfinite(site_o2[site]) && o2_max > o2_min ?
            cmap[clamp(Int(round(1 + (length(cmap)-1) * (site_o2[site]-o2_min)/(o2_max-o2_min))), 1, length(cmap))] :
            base_colors[mod1(i, length(base_colors))]

        plot_vals = [isfinite(sub.thickness_cm[j]) ? sub.R_val[j] * 1e3 : sub.R_val[j] for j in eachindex(sub.R_val)]

        scatterlines!(ax, sub.temperature_K, plot_vals;
            color=color, marker=:circle, markersize=10, linewidth=2, label=label)

        # R(T) = R0*(1 + alpha*T) => alpha = R0*alpha / R0 = b/a
        # units of R don't matter since we divide by R0
        a, b, ok = ols_ab(sub.temperature_K, plot_vals)
        if ok && isfinite(a) && a != 0
            alpha = b / a
            tmin, tmax = extrema(sub.temperature_K)
            tspan = range(tmin, tmax; length=50)
            fit_vals = a .+ b .* tspan
            lines!(ax, tspan, fit_vals; color=color, linestyle=:dash, linewidth=1.5)
            t_annot = maximum(sub.temperature_K)
            R_annot = a + b * t_annot
            text!(ax, t_annot, R_annot; text="TCR=$(round(alpha*1e3, digits=2)) ppm/K", align=(:left, :center), color=color, offset=(8, 0))

            if isfinite(site_o2[site])
                push!(tcr_points, (site=site, o2=site_o2[site], tcr_ppm=alpha*1e3))
            end
        end

        if isfinite(site_o2[site]) && !isempty(sub.temperature_K)
            idx = argmin(abs.(sub.temperature_K .- 298.0))
            rt_val = plot_vals[idx]
            push!(rt_points, (site=site, o2=site_o2[site], resistivity=rt_val))
        end
    end

    axislegend(ax, position=:rt)

    ax_rt = Axis(fig[2, 1], xlabel="O2 (%)", ylabel=any_thickness ? "Resistivity(~298K) (mΩ·cm)" : "Sheet R(~298K) (Ω/□)", title="Resistivity vs O2%")
    ax_tcr = Axis(fig[2, 2], xlabel="O2 (%)", ylabel="TCR (ppm/K)", title="TCR vs O2")

    site_color = Dict{String, Any}()
    for (i, s) in enumerate(sites)
        color = use_o2_colors && isfinite(site_o2[s]) && o2_max > o2_min ?
            cmap[clamp(Int(round(1 + (length(cmap)-1) * (site_o2[s]-o2_min)/(o2_max-o2_min))), 1, length(cmap))] :
            base_colors[mod1(i, length(base_colors))]
        site_color[s] = color
    end

    for p in rt_points
        scatter!(ax_rt, [p.o2], [p.resistivity]; color=site_color[p.site], marker=:circle, markersize=9)
    end
    for p in tcr_points
        scatter!(ax_tcr, [p.o2], [p.tcr_ppm]; color=site_color[p.site], marker=:diamond, markersize=9)
    end

    rowsize!(fig.layout, 1, Relative(0.6))
    rowsize!(fig.layout, 2, Relative(0.4))

    return fig
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

    # Create the plot with two rows: main and secondary+info
    fig = Figure(size=(900, 600))

    # Main plot: Resistance vs Length/Width
    ax1 = Axis(fig[1, 1:2],
        xlabel="Length/Width Ratio (L/W)",
        ylabel="Resistance (kΩ)",
        title="TLM Analysis - Model from (R□, R_c')")

    # Secondary plot: Width-normalized resistance vs length
    ax2 = Axis(fig[2, 1],
        xlabel="Length (μm)",
        ylabel="Width-Normalized Resistance (Ω·μm)",
        title="Width-Normalized Resistance vs Length")

    info_ax = Axis(fig[2, 2], xlabel="", ylabel="", title="Analysis Results")
    hidedecorations!(info_ax)
    hidespines!(info_ax)

    # Balance the two rows (slightly larger main plot)
    rowsize!(fig.layout, 1, Relative(0.65))
    rowsize!(fig.layout, 2, Relative(0.35))

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
    V = df.voltage_drop
    
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

end # module
