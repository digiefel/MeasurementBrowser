using DataFrames: DataFrame, nrow
using GLMakie

function setup_pund_plot(_workspace, measurements)::Figure
    isempty(measurements) && error("PUND plot requires at least one measurement")
    title = length(measurements) == 1 ? item_label(only(measurements)) : "PUND overlay"
    figure = Figure(size=(1100, 720))
    Axis(figure[1, 1:2], xlabel="Time (μs)", ylabel="Current (μA)", title=title)
    Axis(figure[1, 1:2], yaxisposition=:right, ylabel="Voltage (V)")
    Axis(figure[2, 1], xlabel="Voltage (V)", ylabel="Current (μA)", title="I-V")
    Axis(figure[2, 2], xlabel="Voltage (V)", ylabel="Switching charge (pC)", title="Switching")
    return figure
end

function draw_pund_plot(_workspace, measurements, figure::Figure)::Nothing
    processed = item_data.(measurements)
    time_axes = contents(figure[1, 1:2])
    length(time_axes) < 2 && error("PUND plot figure has no voltage axis")
    ax_time, ax_voltage = time_axes[1], time_axes[2]
    ax_iv = only(contents(figure[2, 1]))
    ax_q = only(contents(figure[2, 2]))
    linkxaxes!(ax_time, ax_voltage)

    colors = cgrad(:tab10, max(length(measurements), 1); categorical=true)
    single = length(measurements) == 1
    for (index, (measurement, raw)) in enumerate(zip(measurements, processed))
        df = with_plot_units(raw, measurement)
        nrow(df) == 0 && continue

        color = colors[index]
        label = single ? nothing : item_label(measurement)
        current_color = single ? :blue : color
        voltage_color = single ? :red : color
        fe_color = single ? :purple : color

        current_line = lines!(ax_time, df.time_us, df.current_uA; color=current_color, linewidth=2, label)
        voltage_line = lines!(ax_voltage, df.time_us, df.voltage; color=voltage_color, linewidth=2, linestyle=:dash)
        fe_line = lines!(ax_time, df.time_us, df.i_fe_uA; color=fe_color, linewidth=2, linestyle=:dot)
        lines!(ax_iv, df.voltage, df.current_uA; color=single ? :green : color, linewidth=2)
        lines!(ax_iv, df.voltage, df.i_fe_uA; color=fe_color, linewidth=2, linestyle=:dot)

        y_column = hasproperty(df, :p_fe_uC_cm2) ? :p_fe_uC_cm2 : :q_centered_pC
        ax_q.ylabel = y_column === :p_fe_uC_cm2 ? "Polarization (μC/cm²)" : "Switching charge (pC)"
        for (rep, mask) in pulse_groups(df)
            lines!(ax_q, df.voltage[mask], df[mask, y_column]; color, linewidth=2, label=single ? string(rep) : nothing)
        end

        single && Legend(
            figure[1, 1],
            [current_line, voltage_line, fe_line],
            ["Current", "Voltage", "FE current"];
            tellwidth=false,
            tellheight=false,
            halign=:left,
            valign=:top,
        )
    end
    single ? axislegend(ax_q) : axislegend(ax_time)
    return nothing
end

function setup_iv_plot(_workspace, measurements)::Figure
    isempty(measurements) && error("IV plot requires at least one measurement")
    title = length(measurements) == 1 ? item_label(only(measurements)) : "IV overlay"
    figure = Figure(size=(820, 580))
    Axis(figure[1, 1], xlabel="Voltage (V)", ylabel="|Current| (A)", title=title, yscale=log10)
    return figure
end

function draw_iv_plot(_workspace, measurements, figure::Figure)::Nothing
    processed = item_data.(measurements)
    axis = only(contents(figure[1, 1]))
    for (measurement, df) in zip(measurements, processed)
        nrow(df) == 0 && continue
        label = length(measurements) == 1 ? nothing : item_label(measurement)
        lines!(axis, df.v, abs.(df.i); color=df.i .> 0, colormap=:RdBu_3, linewidth=2, label)
    end
    length(measurements) > 1 && axislegend(axis)
    return nothing
end

function fit_line(df::DataFrame)::Tuple{Vector{Float64},Vector{Float64}}
    i_min, i_max = extrema(df.i)
    span = max(i_max - i_min, 1e-6)
    fit_i = [i_min - 0.1 * span, i_max + 0.1 * span]
    return fit_i, first(df.fit_resistance_ohm) .* fit_i .+ first(df.fit_offset_v)
end

function setup_tlm_plot(_workspace, measurements)::Figure
    isempty(measurements) && error("TLM plot requires at least one measurement")
    title = length(measurements) == 1 ? item_label(only(measurements)) : "TLM 4-point overlay"
    figure = Figure(size=(1000, 560))
    Label(figure[0, :], title, fontsize=18)
    Axis(figure[1, 1], xlabel="Current (μA)", ylabel="Voltage (mV)", title="Linear fit")
    Axis(figure[1, 2], xlabel="Current (μA)", ylabel="Resistance (kΩ)", title="Point resistance")
    return figure
end

function draw_tlm_plot(_workspace, measurements, figure::Figure)::Nothing
    processed = item_data.(measurements)
    ax_iv = only(contents(figure[1, 1]))
    ax_r = only(contents(figure[1, 2]))
    data = processed
    colors = cgrad(:tab10, max(length(measurements), 1); categorical=true)

    for (index, (measurement, df)) in enumerate(zip(measurements, data))
        nrow(df) == 0 && continue
        color = colors[index]
        label = length(measurements) == 1 ? "Data" : item_label(measurement)
        scatter!(ax_iv, df.i .* 1e6, df.v .* 1e3; color, markersize=8, label)

        fit_i, fit_v = fit_line(df)
        fit_label = length(measurements) == 1 ? "Fit: R=$(round(first(df.fit_resistance_ohm); digits=2)) Ω" : nothing
        lines!(ax_iv, fit_i .* 1e6, fit_v .* 1e3; color=length(measurements) == 1 ? :red : color, linewidth=2, linestyle=:dash, label=fit_label)

        valid = df.valid_resistance
        any(valid) && scatter!(
            ax_r,
            df.i[valid] .* 1e6,
            df.resistance_ohm[valid] ./ 1e3;
            color,
            markersize=8,
            label=length(measurements) == 1 ? "V/I" : item_label(measurement),
        )
        hlines!(ax_r, [first(df.fit_resistance_ohm) / 1e3]; color=length(measurements) == 1 ? :red : color, linestyle=:dash)
    end

    axislegend(ax_iv, position=:lt)
    axislegend(ax_r)
    draw_tlm_summary!(figure, measurements, data)
    return nothing
end

function setup_cv_plot(_workspace, measurements)::Figure
    isempty(measurements) && error("C-V plot requires at least one measurement")
    title = length(measurements) == 1 ? item_label(only(measurements)) : "C-V overlay"
    figure = Figure(size=(1050, 500))
    Label(figure[0, :], title, fontsize=18)
    Axis(figure[1, 1], xlabel="Bias (V)", ylabel="C (pF)", title="Capacitance")
    Axis(figure[1, 2], xlabel="Bias (V)", ylabel="|Z| (MΩ)", title="Impedance")
    return figure
end

function draw_cv_plot(_workspace, measurements, figure::Figure)::Nothing
    processed = item_data.(measurements)
    ax_c = only(contents(figure[1, 1]))
    ax_z = only(contents(figure[1, 2]))

    for (measurement, df) in zip(measurements, processed)
        nrow(df) == 0 && continue
        frequencies = sort(unique(df.frequency_Hz))
        colors = cgrad(:viridis, max(length(frequencies), 1); categorical=true)
        for (index, freq_Hz) in enumerate(frequencies)
            mask = df.frequency_Hz .== freq_Hz
            label = length(measurements) == 1 ?
                engineering_label(freq_Hz, "Hz") :
                "$(item_label(measurement)) $(engineering_label(freq_Hz, "Hz"))"
            lines!(ax_c, df.bias_V[mask], df.Cp_F[mask] .* 1e12; color=colors[index], linewidth=2, label)
            lines!(ax_z, df.bias_V[mask], df.Z_Ohm[mask] ./ 1e6; color=colors[index], linewidth=2)
        end
    end

    axislegend(ax_c, position=:lt)
    return nothing
end

function setup_wgfmu_sampling_plot(_workspace, measurements)::Figure
    isempty(measurements) && error("WGFMU sampling plot requires at least one measurement")
    title = length(measurements) == 1 ? item_label(only(measurements)) : "WGFMU Sampling overlay"
    figure = Figure(size=(900, 520))
    Axis(figure[1, 1], xlabel="Time (s)", ylabel="Current (A)", title=title)
    return figure
end

function draw_wgfmu_sampling_plot(_workspace, measurements, figure::Figure)::Nothing
    axis = only(contents(figure[1, 1]))
    for measurement in measurements
        df = item_data(measurement)
        nrow(df) == 0 && continue
        label = length(measurements) == 1 ? nothing : item_label(measurement)
        lines!(axis, df.Time_s, df.Current_1_A; linewidth=1.5, label)
        lines!(axis, df.Time_s, df.Current_2_A; linewidth=1.5, linestyle=:dash)
    end
    length(measurements) > 1 && axislegend(axis)
    return nothing
end
