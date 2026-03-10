function kind_label(::RuO2Project, kind::Symbol)::String
    kind === :pund && return "FE PUND"
    kind === :pund_fatigue && return "PUND Fatigue"
    kind === :iv && return "I-V Sweep"
    kind === :tlm4p && return "TLM 4-Point"
    kind === :breakdown && return "Breakdown"
    kind === :wakeup && return "Wakeup"
    return "Unknown"
end

function display_label(proj::RuO2Project, meas::MeasurementInfo)
    label = kind_label(proj, meas.measurement_kind)
    temp = get(meas.parameters, :temperature_K, nothing)
    parts = Any[meas.timestamp, label]

    if meas.measurement_kind == :wakeup
        meas.wakeup_pulse_count !== nothing && push!(parts, "$(meas.wakeup_pulse_count)×")
    elseif meas.measurement_kind == :pund
        m = match(r"(\d+(?:\.\d+)?)V", meas.filename)
        voltage = m !== nothing ? tryparse(Float64, m.captures[1]) : get(meas.parameters, :voltage_V, nothing)
        voltage !== nothing && push!(parts, "$(voltage)V")
        haskey(meas.parameters, :fatigue_cycle) && push!(parts, "cycle $(meas.parameters[:fatigue_cycle]) (fatigue)")
    end

    temp !== nothing && push!(parts, "$(temp)K")
    return join(parts, " ")
end
