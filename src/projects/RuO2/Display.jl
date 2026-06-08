function kind_label(::RuO2Project, kind::Symbol)::String
    kind === :pund && return "FE PUND"
    kind === :pn && return "PN"
    kind === :pund_wakeup && return "PUND Wakeup"
    kind === :pund_fatigue && return "PUND Fatigue"
    kind === :cvsweep && return "C-V Sweep"
    kind === :iv && return "I-V Sweep"
    kind === :tlm4p && return "TLM 4-Point"
    kind === :breakdown && return "Breakdown"
    kind === :wakeup_pn   && return "Wakeup PN"
    kind === :wakeup_pund && return "Wakeup PUND"
    return "Unknown"
end

function display_label(proj::RuO2Project, meas::MeasurementInfo)
    label = kind_label(proj, meas.measurement_kind)
    temp = get(meas.parameters, :temperature_K, nothing)
    parts = Any[meas.timestamp, label]
    voltage_label = _ruo2_voltage_label(meas.stats)

    if meas.measurement_kind === :wakeup_pn || meas.measurement_kind === :wakeup_pund
        voltage_label !== nothing && push!(parts, voltage_label)
        temp !== nothing && push!(parts, "$(temp)K")
        return join(parts, " ")
    end

    if meas.measurement_kind === :pund || meas.measurement_kind === :pn
        voltage_label !== nothing && push!(parts, voltage_label)
        freq = get(meas.stats, :fatigue_f, NaN)
        if isfinite(freq)
            push!(parts, freq >= 1000 ? "$(round(freq/1000, digits=1)) kHz" : "$(round(freq, digits=1)) Hz")
        end
        fatigue_count = get(meas.stats, :fatigue_count, 0)
        fatigue_count > 0 && push!(parts, "cycle $fatigue_count (fatigue)")
    end

    temp !== nothing && push!(parts, "$(temp)K")
    return join(parts, " ")
end

function _ruo2_voltage_label(stats::Dict{Symbol,Any})::Union{Nothing,String}
    haskey(stats, :V_base) && haskey(stats, :V_amp) || return nothing
    return "$(stats[:V_base]) ± $(stats[:V_amp]) V"
end
