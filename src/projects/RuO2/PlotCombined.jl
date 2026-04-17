function load_plot_for_files(::RuO2Project, paths::Vector{String}, combined_kind::Symbol;
                             device_params_list::Vector{Dict{Symbol,Any}}=Dict{Symbol,Any}[],
                             should_cancel::Union{Nothing,Function}=nothing, kwargs...)
    isempty(paths) && return nothing

    if combined_kind === :tlm_analysis
        files_data_params = Tuple{String,DataFrame,Dict{Symbol,Any}}[]
        for (i, path) in enumerate(paths)
            _check_plot_cancel(should_cancel)
            filename = basename(path)
            dirname_path = dirname(path)
            df = read_tlm_4p(filename, dirname_path)
            device_params = i <= length(device_params_list) ? device_params_list[i] : Dict{Symbol,Any}()
            push!(files_data_params, (path, df, device_params))
        end
        return (files_data_params=files_data_params,)
    elseif combined_kind === :tlm_temperature
        entries = NamedTuple{(:path, :df, :params, :tempK, :site, :oxygen_percent, :oxygen_flow_sccm)}[]
        for (i, path) in enumerate(paths)
            _check_plot_cancel(should_cancel)
            df = read_tlm_4p(basename(path), dirname(path))
            params = i <= length(device_params_list) ? device_params_list[i] : Dict{Symbol,Any}()
            tempK = _extract_temperature_K(params, path)
            isfinite(tempK) || continue
            site = _extract_site_label(params, path)
            o2 = _extract_oxygen_percent(params, path)
            o2_flow = _extract_oxygen_flow(params, path)
            push!(entries, (path=path, df=df, params=params, tempK=tempK, site=site,
                            oxygen_percent=o2, oxygen_flow_sccm=o2_flow))
        end
        return (entries=entries, use_o2_flow=get(kwargs, :use_o2_flow, false))
    elseif combined_kind === :pund_fatigue
        entries = NamedTuple[]
        n = length(paths)
        params_list = length(device_params_list) == n ? device_params_list : [Dict{Symbol,Any}() for _ in 1:n]
        for i in 1:n
            _check_plot_cancel(should_cancel)
            path = paths[i]
            params = params_list[i]
            fname = lowercase(basename(path))
            dirpath = dirname(path)
            ts = stat(path).mtime
            if occursin("wakeup", fname)
                df_w = read_wakeup(basename(path), dirpath)
                push!(entries, (kind=:wakeup, df=df_w, params=params, timestamp=ts))
            elseif haskey(params, :fatigue_cycle)
                df_p = read_pund_fatigue_cycle(basename(path), dirpath, Int(params[:fatigue_cycle]))
                push!(entries, (kind=:pund, df=df_p, params=params, timestamp=ts))
            else
                df_p = read_fe_pund(basename(path), dirpath)
                push!(entries, (kind=:pund, df=df_p, params=params, timestamp=ts))
            end
        end
        return (entries=entries,)
    end

    return nothing
end

function analyze_plot_for_files(::RuO2Project, combined_kind::Symbol, loaded; kwargs...)
    loaded === nothing && return nothing
    should_cancel = get(kwargs, :should_cancel, nothing)
    _check_plot_cancel(should_cancel)

    if combined_kind === :tlm_analysis
        analysis_df = analyze_tlm_combined(loaded.files_data_params)
        nrow(analysis_df) == 0 && return nothing
        R_sheet, R_cprime, rho_c, r_squared = calculate_sheet_resistance(analysis_df)
        geometry_groups = combine(groupby(analysis_df, [:length_um, :width_um]),
            :resistance_ohm => (x -> mean(filter(isfinite, x))) => :avg_resistance_ohm)
        widths = unique(geometry_groups.width_um)
        info_text = ""
        if isfinite(R_sheet)
            info_text *= "Sheet Resistance (R□):  $(round(R_sheet, digits=2)) Ω/□\n\n"
            isfinite(R_cprime) && (info_text *= "Contact per width (R_c'):  $(round(R_cprime, sigdigits=4)) Ω·cm\n\n")
            isfinite(rho_c) && (info_text *= "Specific Contact (ρ_c):  $(round(rho_c, sigdigits=4)) Ω·cm²\n\n")
        else
            info_text *= "Width-invariant fit unavailable (need ≥2 geometries)\n\n"
        end
        info_text *= "Files analyzed: $(length(loaded.files_data_params))\nWidths: $(length(widths))\nData points: $(nrow(analysis_df))"
        return (
            analysis_df=analysis_df,
            geometry_groups=geometry_groups,
            widths=widths,
            fit=(R_sheet=R_sheet, R_cprime=R_cprime, rho_c=rho_c, r_squared=r_squared),
            info_text=info_text,
        )
    elseif combined_kind === :tlm_temperature
        entries = loaded.entries
        isempty(entries) && return nothing
        results = DataFrame(site=String[], temperature_K=Float64[], R_val=Float64[], R_cprime=Float64[],
                            rho_c=Float64[], r_squared=Float64[], n_files=Int[], thickness_cm=Float64[],
                            oxygen_percent=Float64[], oxygen_flow_sccm=Float64[])
        by_site = Dict{String,Vector{NamedTuple}}()
        for entry in entries
            _check_plot_cancel(should_cancel)
            push!(get!(by_site, entry.site, NamedTuple[]), entry)
        end
        for (site, vec) in by_site
            _check_plot_cancel(should_cancel)
            temps = unique([entry.tempK for entry in vec])
            for temp in temps
                _check_plot_cancel(should_cancel)
                subset = filter(entry -> entry.tempK == temp, vec)
                files_data_params = [(entry.path, entry.df, entry.params) for entry in subset]
                combined_df = analyze_tlm_combined(files_data_params)
                nrow(combined_df) == 0 && continue
                R_sheet, R_cprime, rho_c, r2 = calculate_sheet_resistance(combined_df)
                isfinite(R_sheet) || continue
                thickness_values = Float64[]
                for sample in subset
                    if haskey(sample.params, :t_RuO2_nm)
                        t_nm = try
                            Float64(sample.params[:t_RuO2_nm])
                        catch
                            NaN
                        end
                        isfinite(t_nm) && t_nm > 0 && push!(thickness_values, t_nm)
                    end
                end
                thickness_cm = NaN
                if length(thickness_values) == length(subset) && !isempty(thickness_values)
                    thickness_cm = mean(thickness_values) * 1e-7
                end
                R_val = isfinite(thickness_cm) ? R_sheet * thickness_cm : R_sheet
                o2 = begin
                    vals = filter(isfinite, [sample.oxygen_percent for sample in subset])
                    isempty(vals) ? NaN : mean(vals)
                end
                o2_flow = begin
                    vals = filter(isfinite, [sample.oxygen_flow_sccm for sample in subset])
                    isempty(vals) ? NaN : mean(vals)
                end
                push!(results, (site, temp, R_val, R_cprime, rho_c, r2, length(subset), thickness_cm, o2, o2_flow))
            end
        end
        nrow(results) == 0 && return nothing

        sites = sort(unique(results.site))
        site_summary = DataFrame(
            site=String[],
            color_value=Float64[],
            rt_val=Float64[],
            tcr_ppm=Float64[],
        )
        for site in sites
            sub = sort(filter(row -> row.site == site, results), :temperature_K)
            plot_vals = [isfinite(sub.thickness_cm[j]) ? sub.R_val[j] * 1e3 : sub.R_val[j] for j in eachindex(sub.R_val)]
            a, b, ok = _ols_ab(sub.temperature_K, plot_vals)
            alpha_ppm = (ok && isfinite(a) && a != 0) ? (b / a) * 1e6 : NaN
            idx = isempty(plot_vals) ? 1 : argmin(abs.(sub.temperature_K .- 298.0))
            rt_val = isempty(plot_vals) ? NaN : plot_vals[idx]
            o2_vals = filter(isfinite, results.oxygen_percent[results.site .== site])
            color_value = isempty(o2_vals) ? NaN : mean(o2_vals)
            push!(site_summary, (site, color_value, rt_val, alpha_ppm))
        end
        return (results=results, site_summary=site_summary, use_o2_flow=loaded.use_o2_flow)
    elseif combined_kind === :pund_fatigue
        result = analyze_pund_fatigue_combined(loaded.entries)
        traces = get(result, :traces, NamedTuple[])
        pr_points = get(result, :pr_points, NamedTuple[])
        traces_df = DataFrame(
            cycles=Float64[],
            x=Vector{Float64}[],
            y=Vector{Float64}[],
            rep_index=Int[],
            rep_count=Int[],
        )
        for trace in traces
            _check_plot_cancel(should_cancel)
            push!(traces_df, (
                float(trace.cycles),
                collect(trace.x),
                collect(trace.y),
                Int(trace.rep_index),
                Int(trace.rep_count),
            ))
        end
        pr_df = DataFrame(
            cycles=Float64[],
            pr=Float64[],
            rep_index=Int[],
            rep_count=Int[],
        )
        for point in pr_points
            _check_plot_cancel(should_cancel)
            push!(pr_df, (
                float(point.cycles),
                float(point.Pr),
                Int(point.rep_index),
                Int(point.rep_count),
            ))
        end
        return (
            traces_df=traces_df,
            pr_df=pr_df,
            x_label=get(result, :x_label, "Voltage (V)"),
            y_label=get(result, :y_label, "Polarization (μC/cm²)"),
        )
    end

    return loaded
end

function draw_plot_for_files(::RuO2Project, combined_kind::Symbol, analyzed; kwargs...)
    analyzed === nothing && return nothing

    if combined_kind === :tlm_analysis
        analysis_df = analyzed.analysis_df
        geometry_groups = analyzed.geometry_groups
        fit = analyzed.fit
        fig = Figure(size=(900, 600))
        ax1 = Axis(fig[1, 1:2], xlabel="Length/Width Ratio (L/W)", ylabel="Resistance (kΩ)",
            title="TLM Analysis - Model from (R□, R_c')")
        ax2 = Axis(fig[2, 1], xlabel="Length (μm)", ylabel="Width-Normalized Resistance (Ω·μm)",
            title="Width-Normalized Resistance vs Length")
        info_ax = Axis(fig[2, 2], xlabel="", ylabel="", title="Analysis Results")
        hidedecorations!(info_ax)
        hidespines!(info_ax)
        rowsize!(fig.layout, 1, Relative(0.65))
        rowsize!(fig.layout, 2, Relative(0.35))
        widths = analyzed.widths
        colors = cgrad(:tab10, length(widths), categorical=true)
        for (i, width_um) in enumerate(widths)
            width_geom = filter(row -> row.width_um == width_um, geometry_groups)
            if nrow(width_geom) > 0
                scatter!(ax1, width_geom.length_um ./ width_geom.width_um, width_geom.avg_resistance_ohm ./ 1e3,
                    color=colors[i], label="$(width_um) μm", markersize=10)
            end
        end
        if isfinite(fit.R_sheet) && isfinite(fit.R_cprime)
            for (i, width_um) in enumerate(widths)
                mask = geometry_groups.width_um .== width_um
                if any(mask)
                    w_cm = width_um * 1e-4
                    a_w = 2 * (fit.R_cprime / w_cm)
                    xspan = geometry_groups.length_um[mask] ./ geometry_groups.width_um[mask]
                    lw_fit = range(0, max(maximum(xspan), 1.0), length=200)
                    r_fit = (a_w .+ fit.R_sheet .* lw_fit) ./ 1e3
                    lines!(ax1, lw_fit, r_fit, color=colors[i], linewidth=2, linestyle=:dash)
                end
            end
        end
        for (i, width_um) in enumerate(widths)
            width_data = filter(row -> row.width_um == width_um, analysis_df)
            scatter!(ax2, width_data.length_um, width_data.resistance_normalized,
                color=colors[i], label="$(width_um) μm", markersize=6)
        end
        text!(info_ax, 0.05, 0.95, text=analyzed.info_text, align=(:left, :top), fontsize=15, space=:relative)
        xlims!(info_ax, 0, 1)
        ylims!(info_ax, 0, 1)
        length(widths) > 1 && axislegend(ax1, position=:lt)
        return fig
    elseif combined_kind === :tlm_temperature
        results = analyzed.results
        site_summary = analyzed.site_summary
        fig = Figure(size=(1100, 720))
        any_thickness = any(isfinite, results.thickness_cm)
        y_label = any_thickness ? "Resistivity (mΩ·cm)" : "Sheet Resistance (Ω/□)"
        title_str = any_thickness ? "TLM Resistivity vs Temperature" : "TLM Sheet Resistance vs Temperature"
        controls = GridLayout(fig[0, 1:2], tellwidth=false)
        Label(controls[1, 1], "O2 metric")
        toggle_flow = Toggle(controls[1, 2], active=analyzed.use_o2_flow, width=80)
        Label(controls[1, 3], "Use O2 flow (sccm)")
        ax = Axis(fig[1, 1:2], xlabel="Temperature (K)", ylabel=y_label, title=title_str)
        ax_rt = Axis(fig[2, 1], xlabel="O2 (%)",
            ylabel=any_thickness ? "Resistivity(~298K) (mΩ·cm)" : "Sheet R(~298K) (Ω/□)",
            title="Resistivity vs O2%")
        ax_tcr = Axis(fig[2, 2], xlabel="O2 (%)", ylabel="TCR (ppm/K)", title="TCR vs O2")
        sites = sort(unique(results.site))
        cmap = to_colormap(Reverse(:seaborn_flare_gradient))
        base_colors = to_colormap(:tab10)
        finite_color_vals = filter(isfinite, site_summary.color_value)
        color_min, color_max = isempty(finite_color_vals) ? (0.0, 1.0) : (minimum(finite_color_vals), maximum(finite_color_vals))
        colors = [
            (isfinite(site_summary.color_value[i]) && color_max > color_min) ?
                cmap[clamp(Int(round(1 + (length(cmap) - 1) * (site_summary.color_value[i] - color_min) / (color_max - color_min))), 1, length(cmap))] :
                base_colors[mod1(i, length(base_colors))]
            for i in 1:nrow(site_summary)
        ]
        rt_x = [Observable([NaN]) for _ in sites]
        tcr_x = [Observable([NaN]) for _ in sites]
        data_lines = Any[]
        for (i, site) in enumerate(sites)
            sub = sort(filter(row -> row.site == site, results), :temperature_K)
            plot_vals = [isfinite(sub.thickness_cm[j]) ? sub.R_val[j] * 1e3 : sub.R_val[j] for j in eachindex(sub.R_val)]
            push!(data_lines, scatterlines!(ax, sub.temperature_K, plot_vals; color=colors[i], marker=:circle, markersize=10, linewidth=2, label=site))
            summary_row = filter(row -> row.site == site, site_summary)
            nrow(summary_row) == 1 || continue
            rt_val = summary_row.rt_val[1]
            tcr_ppm = summary_row.tcr_ppm[1]
            scatter!(ax_rt, rt_x[i], [rt_val]; color=colors[i], marker=:circle, markersize=9)
            scatter!(ax_tcr, tcr_x[i], [tcr_ppm]; color=colors[i], marker=:diamond, markersize=9)
            if isfinite(tcr_ppm)
                a, b, ok = _ols_ab(sub.temperature_K, plot_vals)
                if ok
                    tspan = range(minimum(sub.temperature_K), maximum(sub.temperature_K); length=50)
                    fit_vals = a .+ b .* tspan
                    lines!(ax, tspan, fit_vals; color=colors[i], linestyle=:dash, linewidth=1.5)
                    text!(ax, maximum(sub.temperature_K), a + b * maximum(sub.temperature_K);
                        text="TCR=$(round(tcr_ppm)) ppm/K", align=(:left, :center), color=colors[i], offset=(8, 0))
                end
            end
        end
        axislegend(ax, position=:rt)
        function update_o2_metric!()
            use_flow = toggle_flow.active[]
            field = use_flow ? :oxygen_flow_sccm : :oxygen_percent
            suffix = use_flow ? " sccm O2" : "% O2"
            ax_rt.xlabel[] = use_flow ? "O2 Flow (sccm)" : "O2 (%)"
            ax_tcr.xlabel[] = ax_rt.xlabel[]
            ax_rt.title[] = use_flow ? "Resistivity vs O2 flow" : "Resistivity vs O2%"
            ax_tcr.title[] = use_flow ? "TCR vs O2 flow" : "TCR vs O2"
            for (i, site) in enumerate(sites)
                vals = filter(isfinite, results[results.site .== site, field])
                metric_val = isempty(vals) ? NaN : mean(vals)
                data_lines[i].label[] = isfinite(metric_val) ? string(site, " (", round(metric_val, digits=2), suffix, ")") : site
                rt_x[i][] = [metric_val]
                tcr_x[i][] = [metric_val]
            end
        end
        update_o2_metric!()
        on(toggle_flow.active) do _
            update_o2_metric!()
            reset_limits!(ax_rt)
            reset_limits!(ax_tcr)
        end
        rowsize!(fig.layout, 1, Relative(0.6))
        rowsize!(fig.layout, 2, Relative(0.4))
        return fig
    elseif combined_kind === :pund_fatigue
        traces_df = analyzed.traces_df
        pr_df = analyzed.pr_df
        (nrow(traces_df) == 0 && nrow(pr_df) == 0) && return nothing
        fig = Figure(size=(1200, 800))
        ax_top = Axis(fig[1, 1], xlabel=analyzed.x_label, ylabel=analyzed.y_label,
            title="Overlapped P-E curves (color-coded by fatigue cycles)")
        gl_ctrl = GridLayout(fig[1, 2])
        ax_bottom = Axis(fig[2, 1:2], xlabel="Fatigue cycles", ylabel="Remnant polarization (μC/cm²)",
            title="Remnant polarization vs fatigue cycles")
        last_only = Observable(false)
        btn = Button(gl_ctrl[1, 1], label="Last only: OFF")
        on(btn.clicks) do _
            last_only[] = !last_only[]
            btn.label[] = last_only[] ? "Last only: ON" : "Last only: OFF"
        end
        cyc_vals = traces_df.cycles
        if isempty(cyc_vals)
            text!(ax_top, 0.5, 0.5, text="No valid P-E traces", align=(:center, :center), color=:gray)
        else
            logvals = log10.(Float64.(cyc_vals))
            vmin = minimum(logvals)
            vmax = maximum(logvals)
            tick_positions = collect(floor(Int, vmin):ceil(Int, vmax))
            tick_labels = string.(Int.(10 .^ tick_positions))
            Colorbar(gl_ctrl[2, 1], colormap=:viridis, limits=(vmin, vmax), label="Fatigue cycles (log10)",
                ticks=(tick_positions, tick_labels))
            plotted_traces = NamedTuple{(:plt, :is_last)}[]
            for row in eachrow(traces_df)
                isempty(row.x) && continue
                isempty(row.y) && continue
                plt = lines!(ax_top, row.x, row.y; color=log10(float(row.cycles)), colormap=:viridis,
                    colorrange=(vmin, vmax), linewidth=2)
                push!(plotted_traces, (plt=plt, is_last=row.rep_index == row.rep_count))
            end
            on(last_only) do value
                for entry in plotted_traces
                    entry.plt.visible[] = (!value) || entry.is_last
                end
            end
        end
        if nrow(pr_df) > 0
            ord_all = sortperm(pr_df.cycles)
            all_x = pr_df.cycles[ord_all]
            all_y = pr_df.pr[ord_all]
            lines_all = lines!(ax_bottom, all_x, all_y, color=:black, linewidth=2, alpha=0.7)
            scatter_all = scatter!(ax_bottom, pr_df.cycles, pr_df.pr, color=:red, markersize=8)
            last_mask = pr_df.rep_index .== pr_df.rep_count
            last_lines = nothing
            last_scatter = nothing
            if any(last_mask)
                ord_last = sortperm(pr_df.cycles[last_mask])
                last_x = pr_df.cycles[last_mask][ord_last]
                last_y = pr_df.pr[last_mask][ord_last]
                last_lines = lines!(ax_bottom, last_x, last_y, color=:black, linewidth=2, alpha=0.7)
                last_scatter = scatter!(ax_bottom, pr_df.cycles[last_mask], pr_df.pr[last_mask], color=:red, markersize=8)
                last_lines.visible[] = false
                last_scatter.visible[] = false
            end
            on(last_only) do value
                lines_all.visible[] = !value
                scatter_all.visible[] = !value
                if last_lines !== nothing
                    last_lines.visible[] = value
                    last_scatter.visible[] = value
                end
            end
        end
        return fig
    end

    return nothing
end
