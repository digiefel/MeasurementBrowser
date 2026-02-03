export plot_fe_pund, plot_pund_fatigue, plot_wakeup

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
        elseif haskey(params, :fatigue_cycle)
            df_p = read_pund_fatigue_cycle(basename(path), dirpath, Int(params[:fatigue_cycle]))
            push!(entries, (kind=:pund, df=df_p, params=params, timestamp=ts))
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
