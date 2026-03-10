function kind_label(::TASEProject, kind::Symbol)::String
    kind === :four_terminal_iv && return "Four-Terminal I-V"
    return "Unknown"
end

function display_label(proj::TASEProject, meas::MeasurementInfo)
    label = kind_label(proj, meas.measurement_kind)
    temp = get(meas.parameters, :temperature_K, nothing)
    parts = Any[meas.timestamp, label]
    temp !== nothing && push!(parts, "$(temp)K")
    return join(parts, " ")
end
