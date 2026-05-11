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

    if meas.measurement_kind === :wakeup_pn || meas.measurement_kind === :wakeup_pund
        amp = get(meas.parameters, :amplitude_V, nothing)
        amp !== nothing && push!(parts, "$(amp)V")
        temp !== nothing && push!(parts, "$(temp)K")
        return join(parts, " ")
    end

    if meas.measurement_kind === :pund || meas.measurement_kind === :pn
        m = match(r"(\d+(?:\.\d+)?)V", meas.filename)
        voltage = m !== nothing ? tryparse(Float64, m.captures[1]) : get(meas.parameters, :voltage_V, nothing)
        voltage !== nothing && push!(parts, "$(voltage)V")
        freq = get(meas.parameters, :frequency_Hz, nothing)
        if freq !== nothing
            push!(parts, freq >= 1000 ? "$(round(freq/1000, digits=1)) kHz" : "$(round(freq, digits=1)) Hz")
        end
        haskey(meas.parameters, :fatigue_cycle) && push!(parts, "cycle $(meas.parameters[:fatigue_cycle]) (fatigue)")
    end

    temp !== nothing && push!(parts, "$(temp)K")
    return join(parts, " ")
end
