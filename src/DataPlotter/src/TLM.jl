export plot_tlm_4p, plot_tlm_combined, plot_tlm_temperature

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
    L_um, W_um = extract_tlm_geometry_from_params(device_params)
    
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