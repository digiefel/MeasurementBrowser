function _ruo2_plot_data(
    measurement::MeasurementInfo,
    df::DataFrame;
    debug::Bool=false,
)
    params = merge(measurement.device_info.parameters, measurement.parameters)
    kind = measurement.measurement_kind
    title = _plot_title(measurement.filepath)
    if kind === :pund || kind === :pn
        segment = kind === :pn ? :pn : nothing
        if is_pund_fatigue_file(measurement.filepath)
            fatigue_count = Int(measurement.parameters[:fatigue_idx])
            title *= " cycle $fatigue_count (fatigue)"
        elseif segment === :pn
            title *= " PN"
        end
        return (df=df, title=title, area_um2=get(params, :area_um2, nothing), debug=debug, segment=segment)
    elseif kind === :wakeup_pn || kind === :wakeup_pund
        segment = kind === :wakeup_pn ? :pn : :pund
        return (df=df, title=title, area_um2=get(params, :area_um2, nothing), debug=debug, segment=segment)
    elseif kind === :iv
        return (df=df, title=title)
    elseif kind === :breakdown
        return (df=df, title=title * " (Breakdown)")
    elseif kind === :tlm4p
        return (df=df, title=title, device_params=params)
    end
    error("Unsupported RuO2 single-measurement plot kind '$kind'")
end

function load_source_data(
    ::RuO2Project,
    source_file::SourceFile;
    measurement::Union{Nothing,MeasurementInfo}=nothing,
    should_cancel::Union{Nothing,Function}=nothing,
)::DataFrame
    _check_plot_cancel(should_cancel)
    kind = measurement === nothing ? detect_kind(RUO2_PROJECT, source_file.filename) : measurement.measurement_kind
    if kind === :pund || kind === :pn
        if is_pund_fatigue_file(source_file.filepath)
            measurement === nothing && error("PUND fatigue data requires a logical measurement")
            fatigue_count = Int(measurement.parameters[:fatigue_idx])
            fatigue_df = read_pund_fatigue_file(source_file.filepath; should_cancel)
            return _select_pund_fatigue_cycle(fatigue_df, fatigue_count)
        end
        return read_pund_file(source_file.filename, dirname(source_file.filepath))
    elseif kind === :wakeup_pn || kind === :wakeup_pund
        measurement === nothing && error("PUND wakeup data requires a logical measurement")
        wakeup_V = Float64(measurement.parameters[:wakeup_V])
        return _select_pund_wakeup_readout(source_file.filepath, wakeup_V)
    elseif kind === :iv || kind === :breakdown
        return read_iv_sweep(source_file.filename, dirname(source_file.filepath))
    elseif kind === :tlm4p
        return read_tlm_4p(source_file.filename, dirname(source_file.filepath))
    elseif kind === :cvsweep
        return read_cv_sweep(source_file.filename, dirname(source_file.filepath))
    end
    error("RuO2 source data API is not implemented for $kind")
end

# this should be moved to PUNDAnalysis or similar once PlotKind is in
function debug_plot(::RuO2Project, measurements::Vector{MeasurementInfo}, loaded; kwargs...)
    measurement = only(measurements)
    measurement.measurement_kind in (:pund, :wakeup_pund) ||
        error("RuO2 debug plots are currently only implemented for PUND measurements")

    df = loaded.df
    time_us = df.time .* 1e6
    current_uA = df.current .* 1e6
    voltage = df.voltage

    fig = Figure(size=(1250, 900))
    controls = GridLayout(fig[0, 1:2], tellwidth=false)
    ax_wave = Axis(fig[1, 1:2], xlabel="Time (μs)", ylabel="Current (μA)",
        title="$(loaded.title) - PUND detector debug")
    ax_voltage = Axis(fig[1, 1:2], yaxisposition=:right, ylabel="Voltage (V)")
    ax_dv = Axis(fig[2, 1:2], xlabel="Time (μs)", ylabel="dV", title="Smoothed derivative and pulse boundaries")
    ax_iv = Axis(fig[3, 1], xlabel="Voltage (V)", ylabel="Current (μA)", title="I-V")
    ax_qv = Axis(fig[3, 2], xlabel="Voltage (V)", ylabel="Q_FE (pC)", title="Q-V")
    linkxaxes!(ax_wave, ax_voltage)

    lines!(ax_wave, time_us, current_uA; color=:blue, linewidth=1)
    lines!(ax_voltage, time_us, voltage; color=:red, linewidth=1)

    smooth = Slider(controls[1, 1], range=3:2:51, startvalue=9)
    dV_fraction = Slider(controls[1, 3], range=0.05:0.05:0.95, startvalue=0.5)
    min_len = Slider(controls[2, 1], range=5:5:100, startvalue=20)
    min_V = Slider(controls[2, 3], range=0.05:0.05:0.8, startvalue=0.25)
    Label(controls[1, 2], lift(v -> "smooth window: $(Int(v))", smooth.value))
    Label(controls[1, 4], lift(v -> "dV threshold: $(round(v, digits=2))", dV_fraction.value))
    Label(controls[2, 2], lift(v -> "min pulse len: $(Int(v))", min_len.value))
    Label(controls[2, 4], lift(v -> "min V amp: $(round(v, digits=2))", min_V.value))

    debug_state = lift(smooth.value, dV_fraction.value, min_len.value, min_V.value) do sw, dvf, ml, mv
        detector_kwargs = (
            smooth_window=Int(sw),
            min_pulse_len=Int(ml),
            dV_threshold_fraction=Float64(dvf),
            min_V_amp_fraction=Float64(mv),
        )
        result = analyze_pund_and_pn(DataFrame(df); DEBUG=true, detector_kwargs...)

        detector = result.pund_df === nothing ? nothing :
            detect_pund_pulses(result.pund_df.time, result.pund_df.voltage, result.pund_df.current; detector_kwargs...)
        return (result=result, detector=detector)
    end

    dV_points = lift(debug_state) do state
        state.detector === nothing && return Point2f[]
        return Point2f.(state.result.pund_df.time .* 1e6, state.detector.dV)
    end
    dV_pos_points = lift(debug_state) do state
        state.detector === nothing && return Point2f[]
        t = state.result.pund_df.time
        isempty(t) && return Point2f[]
        return [
            Point2f(first(t) * 1e6, state.detector.dV_threshold),
            Point2f(last(t) * 1e6, state.detector.dV_threshold),
        ]
    end
    dV_neg_points = lift(debug_state) do state
        state.detector === nothing && return Point2f[]
        t = state.result.pund_df.time
        isempty(t) && return Point2f[]
        return [
            Point2f(first(t) * 1e6, -state.detector.dV_threshold),
            Point2f(last(t) * 1e6, -state.detector.dV_threshold),
        ]
    end
    boundaries = lift(debug_state) do state
        state.detector === nothing && return Float64[]
        t = state.result.pund_df.time .* 1e6
        return sort!([t[i] for pulse in state.detector.pulses for i in (first(pulse), last(pulse))])
    end
    iv_points = lift(debug_state) do state
        analyzed = state.result.pund
        analyzed === nothing && return Point2f.(voltage, current_uA)
        return Point2f.(analyzed.voltage, analyzed.current .* 1e6)
    end
    qv_points = lift(debug_state) do state
        analyzed = state.result.pund
        analyzed === nothing && return Point2f[]
        return Point2f.(analyzed.voltage, analyzed.Q_FE .* 1e12)
    end
    status = lift(debug_state) do state
        state.detector === nothing && return "No PUND segment detected"
        n_pulses = length(state.detector.pulses)
        group_size = state.result.pund_group_size
        threshold = round(state.detector.dV_threshold, sigdigits=4)
        return "$n_pulses pulses, group size $group_size, threshold $threshold"
    end

    lines!(ax_dv, dV_points; color=:black, linewidth=1)
    lines!(ax_dv, dV_pos_points; color=:orange, linestyle=:dash)
    lines!(ax_dv, dV_neg_points; color=:orange, linestyle=:dash)
    vlines!(ax_dv, boundaries; color=:black, linestyle=:dot, linewidth=1)
    lines!(ax_iv, iv_points; color=:green, linewidth=1)
    lines!(ax_qv, qv_points; color=:purple, linewidth=1)
    Label(controls[3, 1:4], status)

    return fig
end

function _ruo2_pulse_groups(df::DataFrame, group_size::Int)
    pulse_groups = Tuple{Int,BitVector}[]
    hasproperty(df, :pulse_idx) || return pulse_groups
    isempty(df.pulse_idx) && return pulse_groups
    max_pulse = maximum(df.pulse_idx)
    for rep in 1:(max_pulse ÷ group_size)
        pulse_range = (rep - 1) * group_size + 1:rep * group_size
        mask = BitVector([pulse in pulse_range for pulse in df.pulse_idx])
        any(mask) && push!(pulse_groups, (rep, mask))
    end
    return pulse_groups
end

function _ruo2_pulse_label(pulse_idx::Int, group_size::Int)
    r = pulse_idx % group_size
    if group_size == 1
        return "PN"
    elseif group_size == 2
        r == 1 && return "P"
        r == 0 && return "N"
    elseif group_size == 4
        r == 1 && return "P"
        r == 2 && return "U"
        r == 3 && return "N"
        r == 0 && return "D"
    else
        r == 1 && return "Poling"
        r == 2 && return "P"
        r == 3 && return "U"
        r == 4 && return "N"
        r == 0 && return "D"
    end
    return ""
end

function _analyze_ruo2_file_plot(::RuO2Project, kind::Union{Symbol,Nothing}, loaded; kwargs...)
    loaded === nothing && return nothing
    should_cancel = get(kwargs, :should_cancel, nothing)
    _check_plot_cancel(should_cancel)

    if kind === :pund || kind === :pn || kind === :wakeup_pund || kind === :wakeup_pn
        debug_mode = get(kwargs, :DEBUG, getproperty(loaded, :debug))
        segment = hasproperty(loaded, :segment) ? loaded.segment : (kind === :pn ? :pn : nothing)
        group_size = 5
        if segment === :pn
            result = analyze_pund_and_pn(loaded.df; DEBUG=debug_mode)
            df = result.pn
            group_size = result.pn_group_size
        elseif segment === :pund
            result = analyze_pund_and_pn(loaded.df; DEBUG=debug_mode)
            df = result.pund
            group_size = result.pund_group_size
        else
            df = debug_mode ? analyze_pund(loaded.df; DEBUG=true) : analyze_pund(loaded.df)
        end
        df === nothing && return nothing
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

        pulse_groups = _ruo2_pulse_groups(analyzed_df, group_size)

        debug_labels = NamedTuple[]
        if debug_mode
            pid = analyzed_df.pulse_idx
            time_us = analyzed_df.time_us
            boundaries = Float64[]
            for i in 2:length(pid)
                pid[i] != pid[i - 1] && push!(boundaries, time_us[i])
            end

            seg_start = 1
            for i in 2:length(pid)
                if pid[i] != pid[i - 1]
                    if pid[seg_start] > 0
                        label = _ruo2_pulse_label(pid[seg_start], group_size)
                        if !isempty(label)
                            tmid = (time_us[seg_start] + time_us[i - 1]) / 2
                            push!(debug_labels, (time_us=tmid, label=label))
                        end
                    end
                    seg_start = i
                end
            end
            if !isempty(pid) && pid[seg_start] > 0
                label = _ruo2_pulse_label(pid[seg_start], group_size)
                if !isempty(label)
                    tmid = (time_us[seg_start] + time_us[end]) / 2
                    push!(debug_labels, (time_us=tmid, label=label))
                end
            end
            return (df=analyzed_df, title=loaded.title, area_um2=loaded.area_um2, pulse_groups=pulse_groups,
                    pulse_group_size=group_size, remnant_y_label=y_label, debug=true,
                    debug_boundaries=boundaries, debug_labels=debug_labels)
        end

        return (df=analyzed_df, title=loaded.title, area_um2=loaded.area_um2, pulse_groups=pulse_groups,
                pulse_group_size=group_size, remnant_y_label=y_label, debug=false)
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
    end

    return loaded
end

function setup_plot(
    ::RuO2Project,
    plot_kind::Symbol,
    measurements::Vector{MeasurementInfo},
)::Figure
    plot_kind === :cvsweep || error("Unsupported RuO2 plot kind '$plot_kind'")
    isempty(measurements) && error("RuO2 CVSweep plot requires at least one measurement")
    title = length(measurements) == 1 ? only(measurements).clean_title : "RuO2 CVSweep overlay"
    fig = Figure(size=(1100, 500))
    Axis(fig[1, 1], xlabel="Bias (V)", ylabel="C (pF)", title="Capacitance")
    Axis(fig[1, 2], xlabel="Bias (V)", ylabel="|Z| (MΩ)", title="Impedance")
    Label(fig[0, :], title, fontsize=20, font=:bold)
    return fig
end

function plot_data!(
    project::RuO2Project,
    plot_kind::Symbol,
    measurements::Vector{MeasurementInfo},
    figure::Figure,
)::Nothing
    plot_kind === :cvsweep || error("Unsupported RuO2 plot kind '$plot_kind'")
    axes_c = contents(figure[1, 1])
    axes_z = contents(figure[1, 2])
    isempty(axes_c) && error("RuO2 CVSweep plot figure has no capacitance axis")
    isempty(axes_z) && error("RuO2 CVSweep plot figure has no impedance axis")
    ax_c = only(axes_c)
    ax_z = only(axes_z)

    data = data_of_measurements(project, measurements)
    for (measurement, df) in zip(measurements, data)
        nrow(df) == 0 && continue
        frequencies_Hz = sort(unique(df.frequency_Hz))
        colors = cgrad(:viridis, length(frequencies_Hz); categorical=true)
        for (idx, freq_Hz) in enumerate(frequencies_Hz)
            mask = df.frequency_Hz .== freq_Hz
            color = colors[idx]
            label = length(measurements) == 1 ?
                _format_frequency_label(freq_Hz) :
                "$(measurement.clean_title) $(_format_frequency_label(freq_Hz))"
            lines!(ax_c, df.bias_V[mask], df.Cp_F[mask] .* 1e12; color, linewidth=2, label)
            lines!(ax_z, df.bias_V[mask], df.Z_Ohm[mask] ./ 1e6; color, linewidth=2)
        end
    end
    axislegend(ax_c, position=:lt)
    return nothing
end

function _draw_ruo2_file_plot(::RuO2Project, kind::Union{Symbol,Nothing}, analyzed; kwargs...)
    analyzed === nothing && return nothing

    if kind === :pund || kind === :pn || kind === :wakeup_pund || kind === :wakeup_pn
        df = analyzed.df
        nrow(df) == 0 && return nothing
        if hasproperty(analyzed, :debug) && analyzed.debug
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
            yI_finite = filter(isfinite, df.current_uA)
            yImaxabs = isempty(yI_finite) ? 1.0 : maximum(abs.(yI_finite))
            yImaxabs = isfinite(yImaxabs) && yImaxabs > 0 ? yImaxabs : 1.0
            yImaxv = isempty(yI_finite) ? yImaxabs : maximum(yI_finite)
            yQ_finite = filter(isfinite, df.q_fe_pC)
            yQmaxabs = isempty(yQ_finite) ? 1.0 : maximum(abs.(yQ_finite))
            yQmaxabs = isfinite(yQmaxabs) && yQmaxabs > 0 ? yQmaxabs : 1.0
            yQmaxv = isempty(yQ_finite) ? yQmaxabs : maximum(yQ_finite)
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
    end

    return nothing
end
