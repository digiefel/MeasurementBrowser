function load_plot_for_file(::RuO2Project, path::AbstractString, kind::Union{Symbol,Nothing};
                            device_params=Dict{Symbol,Any}(), should_cancel::Union{Nothing,Function}=nothing, kwargs...)
    _check_plot_cancel(should_cancel)
    if kind === :pund
        params = device_params isa Dict ? device_params : Dict{Symbol,Any}()
        area_um2 = get(params, :area_um2, nothing)
        title = _plot_title(path)
        fatigue_cycle = get(params, :fatigue_cycle, nothing)
        if fatigue_cycle !== nothing
            df = read_pund_fatigue_cycle(basename(path), dirname(path), Int(fatigue_cycle))
            return (df=df, title=title * " cycle $fatigue_cycle (fatigue)", area_um2=area_um2, debug=get(kwargs, :DEBUG, false))
        end
        df = read_fe_pund(basename(path), dirname(path))
        return (df=df, title=title, area_um2=area_um2, debug=get(kwargs, :DEBUG, false))
    elseif kind === :iv
        df = read_iv_sweep(basename(path), dirname(path))
        return (df=df, title=_plot_title(path))
    elseif kind === :breakdown
        df = read_iv_sweep(basename(path), dirname(path))
        return (df=df, title=_plot_title(path) * " (Breakdown)")
    elseif kind === :tlm4p
        df = read_tlm_4p(basename(path), dirname(path))
        return (df=df, title=_plot_title(path), device_params=device_params isa Dict ? device_params : Dict{Symbol,Any}())
    elseif kind === :wakeup
        df = read_wakeup(basename(path), dirname(path))
        return (df=df, title=_plot_title(path))
    elseif kind === :cvsweep
        loaded = read_cv_sweep(basename(path), dirname(path))
        return (df=loaded.df, title=_plot_title(path), secondary_kind=loaded.secondary_kind)
    end
    @warn "Unsupported RuO2 single-file plot kind" kind path
    return nothing
end

function analyze_plot_for_file(::RuO2Project, kind::Union{Symbol,Nothing}, loaded; kwargs...)
    loaded === nothing && return nothing
    should_cancel = get(kwargs, :should_cancel, nothing)
    _check_plot_cancel(should_cancel)

    if kind === :pund
        debug_mode = get(kwargs, :DEBUG, getproperty(loaded, :debug))
        df = debug_mode ? analyze_pund(loaded.df; DEBUG=true) : analyze_pund(loaded.df)
        analyzed_df = DataFrame(df)
        analyzed_df.time_us = analyzed_df.time .* 1e6
        analyzed_df.current_uA = analyzed_df.current .* 1e6
        analyzed_df.i_fe_uA = analyzed_df.I_FE .* 1e6
        analyzed_df.q_fe_pC = analyzed_df.Q_FE .* 1e12

        finite_q = filter(isfinite, analyzed_df.Q_FE)
        q_mean = isempty(finite_q) ? 0.0 : mean(finite_q)
        q_centered = analyzed_df.Q_FE .- q_mean
        analyzed_df.q_centered_pC = q_centered .* 1e12
        y_label = "Switching Charge (pC)"
        if loaded.area_um2 !== nothing
            area_cm2 = loaded.area_um2 / 1e8
            analyzed_df.p_fe_uC_cm2 = (q_centered ./ area_cm2) .* 1e6
            y_label = "Remnant Polarization (μC/cm²)"
        end

        pulse_groups = Tuple{Int,BitVector}[]
        max_pulse = maximum(analyzed_df.pulse_idx)
        for rep in 1:(max_pulse ÷ 5)
            pulse_range = (rep - 1) * 5 + 1:rep * 5
            mask = BitVector([pulse in pulse_range for pulse in analyzed_df.pulse_idx])
            any(mask) && push!(pulse_groups, (rep, mask))
        end

        debug_labels = NamedTuple[]
        if debug_mode
            pid = analyzed_df.pulse_idx
            time_us = analyzed_df.time_us
            boundaries = Float64[]
            for i in 2:length(pid)
                pid[i] != pid[i - 1] && push!(boundaries, time_us[i])
            end

            function pulse_label(p::Int)
                r = p % 5
                r == 1 && return "Poling"
                r == 2 && return "P"
                r == 3 && return "U"
                r == 4 && return "N"
                r == 0 && return "D"
                return ""
            end

            seg_start = 1
            for i in 2:length(pid)
                if pid[i] != pid[i - 1]
                    if pid[seg_start] > 0
                        label = pulse_label(pid[seg_start])
                        if !isempty(label)
                            tmid = (time_us[seg_start] + time_us[i - 1]) / 2
                            push!(debug_labels, (time_us=tmid, label=label))
                        end
                    end
                    seg_start = i
                end
            end
            if !isempty(pid) && pid[seg_start] > 0
                label = pulse_label(pid[seg_start])
                if !isempty(label)
                    tmid = (time_us[seg_start] + time_us[end]) / 2
                    push!(debug_labels, (time_us=tmid, label=label))
                end
            end
            return (df=analyzed_df, title=loaded.title, area_um2=loaded.area_um2, pulse_groups=pulse_groups,
                    remnant_y_label=y_label, debug=true, debug_boundaries=boundaries, debug_labels=debug_labels)
        end

        return (df=analyzed_df, title=loaded.title, area_um2=loaded.area_um2, pulse_groups=pulse_groups,
                remnant_y_label=y_label, debug=false)
    elseif kind === :iv || kind === :breakdown || kind === :unknown || kind === nothing
        df = DataFrame(v=loaded.df.v, i_abs=abs.(loaded.df.i), i_positive=loaded.df.i .> 0)
        return (df=df, title=loaded.title)
    elseif kind === :tlm4p
        df = loaded.df
        mask = isfinite.(df.current_source) .& isfinite.(df.voltage_drop)
        I = df.current_source[mask]
        V = df.voltage_drop[mask]
        isempty(I) && return nothing

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

        L_um, W_um = extract_tlm_geometry_from_params(loaded.device_params)
        rho_sheet = (!isnan(L_um) && !isnan(W_um) && L_um > 0) ? R_fit * W_um / L_um : NaN

        parts = split(loaded.title, ['_', ' '])
        chip = !isempty(parts) ? parts[1] : "?"
        geometry_str = "?"
        temp = "?"
        for part in parts
            occursin(r"L\d+W\d+", part) && (geometry_str = part)
            occursin(r"\d+K", part) && (temp = part)
        end
        if haskey(loaded.device_params, :temperature_K)
            temp = "$(loaded.device_params[:temperature_K])K"
        end

        i_uA = I .* 1e6
        v_mV = V .* 1e3
        r_kohm = (V ./ I) ./ 1e3
        i_min, i_max = minimum(i_uA), maximum(i_uA)
        i_span = i_max - i_min
        i_span == 0 && (i_span = 1.0)
        i_fit = [i_min - 0.1 * i_span, i_max + 0.1 * i_span]
        v_fit_mV = (R_fit / 1e3) .* i_fit .+ (offset * 1e3)

        analyzed_df = DataFrame(
            current_uA=i_uA,
            voltage_mV=v_mV,
            resistance_kohm=r_kohm,
            valid_resistance=isfinite.(r_kohm) .& (abs.(r_kohm) .< 1e6),
        )
        return (
            df=analyzed_df,
            title="$chip $geometry_str $temp",
            fit_resistance_ohm=R_fit,
            fit_resistance_kohm=R_fit / 1e3,
            fit_current_uA=i_fit,
            fit_voltage_mV=v_fit_mV,
            rho_sheet=rho_sheet,
        )
    elseif kind === :wakeup
        pulse_count = loaded.df.pulse_count[1]
        amplitude = loaded.df.amplitude[1]
        return (df=loaded.df, title=loaded.title, pulse_count=pulse_count, amplitude=amplitude,
                text_content="$(pulse_count)× wakeup pulses\namplitude = $(amplitude) V")
    elseif kind === :cvsweep
        df = loaded.df
        isempty(df) && return nothing
        secondary_kind = loaded.secondary_kind
        secondary_label = secondary_kind === :conductance ? "G (nS)" : "Rp (MΩ)"
        return (
            df=df,
            title=loaded.title,
            frequencies_Hz=sort(unique(df.frequency_Hz)),
            secondary_kind=secondary_kind,
            secondary_label=secondary_label,
        )
    end

    return loaded
end

function draw_plot_for_file(::RuO2Project, kind::Union{Symbol,Nothing}, analyzed; kwargs...)
    analyzed === nothing && return nothing

    if kind === :pund
        df = analyzed.df
        nrow(df) == 0 && return nothing
        if analyzed.debug
            fig = Figure(size=(1200, 800))
            ax1 = Axis(fig[1, 1:2], xlabel="Time (μs)", ylabel="Current (μA)",
                yticklabelcolor=:blue, title="$(analyzed.title) - DEBUG (time domain)")
            ax1twin = Axis(fig[1, 1:2], yaxisposition=:right, ylabel="Voltage (V)", yticklabelcolor=:red)
            linkxaxes!(ax1, ax1twin)
            lI = lines!(ax1, df.time_us, df.current_uA, color=:blue, linewidth=2)
            lIFE = lines!(ax1, df.time_us, df.i_fe_uA, color=:purple, linewidth=2)
            lV = lines!(ax1twin, df.time_us, df.voltage, color=:red, linewidth=1, linestyle=:dash)
            ax2 = Axis(fig[2, 1:2], xlabel="Time (μs)", ylabel="Switching Charge (pC)", title="$(analyzed.title) - Q_FE(t)")
            lines!(ax2, df.time_us, df.q_fe_pC, color=:orange, linewidth=2)
            for t in analyzed.debug_boundaries
                vlines!(ax1, t, color=:black, linestyle=:dot, linewidth=1, alpha=0.5)
                vlines!(ax2, t, color=:black, linestyle=:dot, linewidth=1, alpha=0.5)
            end
            yI = df.current_uA
            yImaxabs = isempty(yI) ? 1.0 : maximum(abs.(filter(isfinite, yI)))
            yImaxabs = isfinite(yImaxabs) && yImaxabs > 0 ? yImaxabs : 1.0
            yImaxv = isempty(yI) ? yImaxabs : maximum(filter(isfinite, yI))
            yQ = df.q_fe_pC
            yQmaxabs = isempty(yQ) ? 1.0 : maximum(abs.(filter(isfinite, yQ)))
            yQmaxabs = isfinite(yQmaxabs) && yQmaxabs > 0 ? yQmaxabs : 1.0
            yQmaxv = isempty(yQ) ? yQmaxabs : maximum(filter(isfinite, yQ))
            for entry in analyzed.debug_labels
                text!(ax1, entry.time_us, yImaxv + 0.25 * yImaxabs; text=entry.label, align=(:center, :baseline), color=:black)
                text!(ax2, entry.time_us, yQmaxv + 0.25 * yQmaxabs; text=entry.label, align=(:center, :baseline), color=:black)
            end
            Legend(fig[1, 1], [lI, lV, lIFE], ["Current", "Voltage", "FE Current"],
                tellwidth=false, tellheight=false, halign=:left, valign=:top)
            return fig
        end

        fig = Figure(size=(1200, 800))
        area_str = isnothing(analyzed.area_um2) ? "?" : round(analyzed.area_um2, digits=2)
        ax1 = Axis(fig[1, 1:2], xlabel="Time (μs)", ylabel="Current (μA)", yticklabelcolor=:blue,
            title="$(analyzed.title) - Area = $area_str um²")
        ax1twin = Axis(fig[1, 1:2], yaxisposition=:right, ylabel="Voltage (V)", yticklabelcolor=:red)
        ax2 = Axis(fig[2, 1], xlabel="Voltage (V)", ylabel="Current (μA)", title="$(analyzed.title) - I-V Characteristic")
        ax3 = Axis(fig[2, 2], xlabel="Voltage (V)", ylabel=analyzed.remnant_y_label, title="$(analyzed.title) - Ferroelectric Switching Charge")
        linkxaxes!(ax1, ax1twin)
        l1 = lines!(ax1, df.time_us, df.current_uA, color=:blue, linewidth=2)
        l2 = lines!(ax1twin, df.time_us, df.voltage, color=:red, linewidth=2, linestyle=:dash)
        l3 = lines!(ax1, df.time_us, df.i_fe_uA, color=:purple, linewidth=2)
        lines!(ax2, df.voltage, df.current_uA, color=:green, linewidth=2)
        lines!(ax2, df.voltage, df.i_fe_uA, color=:purple, linewidth=2)
        for (rep, mask) in analyzed.pulse_groups
            if hasproperty(df, :p_fe_uC_cm2)
                lines!(ax3, df.voltage[mask], df.p_fe_uC_cm2[mask], linewidth=2, color=:purple, label="$rep")
            else
                lines!(ax3, df.voltage[mask], df.q_centered_pC[mask], linewidth=2, color=:purple, label="$rep")
            end
        end
        Legend(fig[1, 1], [l1, l2, l3], ["Current", "Voltage", "FE Current"],
            tellwidth=false, tellheight=false, halign=:left, valign=:top)
        axislegend(ax3)
        return fig
    elseif kind === :iv || kind === :breakdown || kind === :unknown || kind === nothing
        nrow(analyzed.df) == 0 && return nothing
        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1], xlabel="Voltage (V)", ylabel="Current (A)", title=analyzed.title)
        lines!(ax, analyzed.df.v, analyzed.df.i_abs, color=analyzed.df.i_positive, colormap=:RdBu_3, linewidth=2)
        ax.yscale = log10
        return fig
    elseif kind === :tlm4p
        df = analyzed.df
        nrow(df) == 0 && return nothing
        fig = Figure(size=(1000, 500))
        ax1 = Axis(fig[1, 1], xlabel="Current (μA)", ylabel="Voltage (mV)", title="I-V Curve")
        scatter!(ax1, df.current_uA, df.voltage_mV, color=:blue, markersize=8, label="Data")
        lines!(ax1, analyzed.fit_current_uA, analyzed.fit_voltage_mV, color=:red, linewidth=2,
            label="Fit: R=$(round(analyzed.fit_resistance_ohm, digits=2)) Ω")
        if isfinite(analyzed.rho_sheet)
            lines!(ax1, [NaN], [NaN], color=:transparent, label="ρ_sq = $(round(analyzed.rho_sheet, digits=2)) Ω/sq")
        end
        axislegend(ax1, position=:lt)
        ax2 = Axis(fig[1, 2], xlabel="Current (μA)", ylabel="Resistance (kΩ)", title="Resistance vs Current")
        valid = df.valid_resistance
        any(valid) && scatter!(ax2, df.current_uA[valid], df.resistance_kohm[valid], color=:green, markersize=8, label="R = V/I")
        hlines!(ax2, [analyzed.fit_resistance_kohm], color=:red, linestyle=:dash, linewidth=2, label="Fitted R")
        axislegend(ax2)
        Label(fig[0, :], analyzed.title, fontsize=20, font=:bold)
        return fig
    elseif kind === :wakeup
        nrow(analyzed.df) == 0 && return nothing
        fig = Figure(size=(600, 400))
        ax = Axis(fig[1, 1], title=analyzed.title)
        hidedecorations!(ax)
        hidespines!(ax)
        text!(ax, 0.5, 0.5, text=analyzed.text_content, align=(:center, :center), fontsize=24, color=:black)
        xlims!(ax, 0, 1)
        ylims!(ax, 0, 1)
        return fig
    elseif kind === :cvsweep
        df = analyzed.df
        nrow(df) == 0 && return nothing

        fig = Figure(size=(1100, 500))
        ax_cp = Axis(fig[1, 1], xlabel="Bias (V)", ylabel="Cp (pF)", title="Capacitance")
        ax_secondary = Axis(fig[1, 2], xlabel="Bias (V)", ylabel=analyzed.secondary_label, title="Secondary Channel")

        colors = cgrad(:viridis, length(analyzed.frequencies_Hz); categorical=true)
        for (idx, freq_Hz) in enumerate(analyzed.frequencies_Hz)
            mask = df.frequency_Hz .== freq_Hz
            color = colors[idx]
            label = _format_frequency_label(freq_Hz)
            lines!(ax_cp, df.bias_V[mask], df.cp_F[mask] .* 1e12, color=color, linewidth=2, label=label)

            if analyzed.secondary_kind === :conductance
                lines!(ax_secondary, df.bias_V[mask], df.secondary_value[mask] .* 1e9, color=color, linewidth=2)
            else
                lines!(ax_secondary, df.bias_V[mask], df.secondary_value[mask] ./ 1e6, color=color, linewidth=2)
            end
        end

        Label(fig[0, :], analyzed.title, fontsize=20, font=:bold)
        axislegend(ax_cp, position=:lt)
        return fig
    end

    return nothing
end
