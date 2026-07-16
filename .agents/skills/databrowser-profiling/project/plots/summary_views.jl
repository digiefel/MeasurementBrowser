using DataFrames: DataFrame, nrow
using DataFrames: combine, groupby
using Dates: DateTime
using GLMakie
using Statistics: mean

function setup_tlm_analysis_plot(_workspace, measurements)::Figure
    isempty(measurements) && error("TLM analysis plot requires at least one measurement")
    figure = Figure(size=(900, 620))
    Axis(
        figure[1, 1:2],
        xlabel="Length/width ratio (L/W)",
        ylabel="Resistance (kΩ)",
        title="TLM model",
    )
    Axis(
        figure[2, 1],
        xlabel="Length (μm)",
        ylabel="Width-normalized resistance (Ω·μm)",
        title="Width-normalized resistance",
    )
    info_axis = Axis(figure[2, 2], xlabel="", ylabel="", title="Fit")
    hidedecorations!(info_axis)
    hidespines!(info_axis)
    rowsize!(figure.layout, 1, Relative(0.62))
    rowsize!(figure.layout, 2, Relative(0.38))
    return figure
end

function draw_tlm_analysis_plot(_workspace, measurements, figure::Figure)::Nothing
    processed = item_data.(measurements)
    analysis = tlm_combined_table(measurements, processed)
    nrow(analysis) == 0 && return nothing

    ax_model = only(contents(figure[1, 1:2]))
    ax_norm = only(contents(figure[2, 1]))
    ax_info = only(contents(figure[2, 2]))
    fit = fit_tlm_sheet_resistance(analysis)
    geometry = combine(
        groupby(analysis, [:length_um, :width_um]),
        :resistance_ohm => (values -> mean(filter(isfinite, values))) => :R_ohm,
    )
    widths = sort(unique(geometry.width_um))
    colors = cgrad(:tab10, max(length(widths), 1); categorical=true)

    for (index, width_um) in enumerate(widths)
        rows = geometry[geometry.width_um .== width_um, :]
        nrow(rows) == 0 && continue
        scatter!(
            ax_model,
            rows.length_um ./ rows.width_um,
            rows.R_ohm ./ 1e3;
            color=colors[index],
            markersize=10,
            label="$(width_um) μm",
        )

        if isfinite(fit.R_sheet) && isfinite(fit.R_cprime)
            width_cm = width_um * 1e-4
            contact_term = 2 * fit.R_cprime / width_cm
            x_max = max(maximum(rows.length_um ./ rows.width_um), 1.0)
            x_fit = range(0, x_max; length=200)
            lines!(
                ax_model,
                x_fit,
                (contact_term .+ fit.R_sheet .* x_fit) ./ 1e3;
                color=colors[index],
                linewidth=2,
                linestyle=:dash,
            )
        end

        raw_rows = analysis[analysis.width_um .== width_um, :]
        scatter!(
            ax_norm,
            raw_rows.length_um,
            raw_rows.resistance_normalized;
            color=colors[index],
            markersize=6,
            label="$(width_um) μm",
        )
    end

    text!(
        ax_info,
        0.04,
        0.94;
        text=tlm_fit_text(fit, nrow(analysis), length(widths)),
        align=(:left, :top),
        fontsize=14,
        space=:relative,
    )
    xlims!(ax_info, 0, 1)
    ylims!(ax_info, 0, 1)
    length(widths) > 1 && axislegend(ax_model, position=:lt)
    return nothing
end

function tlm_fit_text(fit::NamedTuple, rows::Integer, width_count::Integer)::String
    lines = String[]
    if isfinite(fit.R_sheet)
        push!(lines, "R□: $(round(fit.R_sheet; digits=2)) Ω/□")
        isfinite(fit.R_cprime) && push!(lines, "R_c': $(round(fit.R_cprime; sigdigits=4)) Ω·cm")
        isfinite(fit.rho_c) && push!(lines, "ρ_c: $(round(fit.rho_c; sigdigits=4)) Ω·cm²")
        isfinite(fit.r_squared) && push!(lines, "R²: $(round(fit.r_squared; digits=3))")
    else
        push!(lines, "Fit unavailable")
        push!(lines, "Need at least two valid geometries")
    end
    append!(lines, ["", "Widths: $width_count", "Data points: $rows"])
    return join(lines, "\n")
end

function draw_tlm_summary!(figure::Figure, measurements, data)::Nothing
    table = tlm_combined_table(measurements, data)
    nrow(table) == 0 && return nothing

    fit = fit_tlm_sheet_resistance(table)
    isfinite(fit.R_sheet) || return nothing

    info_axis = Axis(figure[2, 1:2], xlabel="", ylabel="", title="Combined TLM fit")
    hidedecorations!(info_axis)
    hidespines!(info_axis)
    text!(
        info_axis,
        0.02,
        0.9;
        text="R□: $(round(fit.R_sheet; digits=2)) Ω/□\nR_c': $(round(fit.R_cprime; sigdigits=4)) Ω·cm\nρ_c: $(round(fit.rho_c; sigdigits=4)) Ω·cm²\nR²: $(round(fit.r_squared; digits=3))",
        align=(:left, :top),
        fontsize=14,
        space=:relative,
    )
    xlims!(info_axis, 0, 1)
    ylims!(info_axis, 0, 1)
    rowsize!(figure.layout, 2, Fixed(135))
    return nothing
end

function fatigue_entries(measurements, processed)::Vector{NamedTuple}
    entries = NamedTuple[]
    for (measurement, df) in zip(measurements, processed)
        params = metadata(measurement)
        cycle = get(params, :cycle, get(params, :fatigue_count, 0))
        params = merge(params, Dict(:fatigue_count => cycle))
        timestamp = something(get(params, :timestamp, nothing), DateTime(1))
        push!(entries, (kind=:pund, df=df, params=params, timestamp=timestamp))
    end
    return entries
end

function setup_fatigue_plot(_workspace, measurements)::Figure
    isempty(measurements) && error("PUND fatigue plot requires at least one measurement")
    figure = Figure(size=(1120, 760))
    Axis(figure[1, 1], title="P-E curves by fatigue cycle")
    Axis(figure[2, 1], xlabel="Fatigue cycles", ylabel="Remnant polarization (μC/cm²)", title="Pr vs fatigue")
    return figure
end

function draw_fatigue_plot(_workspace, measurements, figure::Figure)::Nothing
    processed = item_data.(measurements)
    result = analyze_pund_fatigue_combined(fatigue_entries(measurements, processed))
    traces = get(result, :traces, NamedTuple[])
    points = get(result, :pr_points, NamedTuple[])
    isempty(traces) && isempty(points) && return nothing

    ax_traces = only(contents(figure[1, 1]))
    ax_pr = only(contents(figure[2, 1]))
    ax_traces.xlabel = get(result, :x_label, "Voltage (V)")
    ax_traces.ylabel = get(result, :y_label, "Polarization (μC/cm²)")

    if !isempty(traces)
        log_cycles = log10.(max.(Float64.(getproperty.(traces, :cycles)), 1.0))
        limits = extrema(log_cycles)
        for trace in traces
            color_value = log10(max(Float64(trace.cycles), 1.0))
            lines!(ax_traces, trace.x, trace.y; color=color_value, colormap=:viridis, colorrange=limits, linewidth=2)
        end
        Colorbar(figure[1, 2], colormap=:viridis, limits=limits, label="log10(cycles)")
    end

    if !isempty(points)
        ordered = sort(points; by=point -> point.cycles)
        cycles = getproperty.(ordered, :cycles)
        pr = getproperty.(ordered, :Pr)
        lines!(ax_pr, cycles, pr; color=:black, linewidth=2)
        scatter!(ax_pr, cycles, pr; color=:red, markersize=8)
    end

    return nothing
end
