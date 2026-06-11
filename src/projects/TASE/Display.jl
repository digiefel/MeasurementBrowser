"""Return the browser label for one TASE measurement kind."""
function kind_label(::TASEProject, kind::Symbol)::String
    kind === :four_terminal_iv && return "Four-Terminal I-V"
    return "Unknown"
end

"""Build the measurement label shown by the browser."""
function display_label(project::TASEProject, measurement::MeasurementInfo)::String
    label = kind_label(project, measurement.measurement_kind)
    temp = get(measurement.parameters, :temperature_K, nothing)
    parts = Any[measurement.timestamp, label]
    temp !== nothing && push!(parts, "$(temp)K")
    return join(parts, " ")
end
