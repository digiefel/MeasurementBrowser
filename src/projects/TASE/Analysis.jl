using DataFrames
using DataLoader: read_tlm_4p

const TASE_CONDUCTANCE_QUANTUM_S = 77.480_917_30e-6
const TASE_WIRE_COUNT_MATRIX = [
    1 1 1 1;
    1 1 1 1;
    4 4 4 4;
    4 4 4 4;
    16 16 16 16;
    16 16 0 0;
]
const TASE_WIRE_DIAMETERS_NM = [32.0, 36.0, 40.0, 48.0, 0.0]

function available_analyses(::TASEProject, measurements)
    tase_measurements = _tase_four_terminal_measurements(measurements)
    isempty(tase_measurements) && return NamedTuple[]
    return [
        (
            key=:iv_trace_table,
            label="IV Trace Table",
            description="One row per IV datapoint with TASE device metadata.",
            row_kind=:iv_point,
            parameters=NamedTuple[],
        ),
        (
            key=:iv_device_summary,
            label="IV Device Summary",
            description="One row per IV measurement with fitted resistance and conductance metrics.",
            row_kind=:device_measurement,
            parameters=NamedTuple[],
        ),
    ]
end

function run_analysis(::TASEProject, key::Symbol, measurements; should_cancel::Union{Nothing,Function}=nothing, kwargs...)
    tase_measurements = _tase_four_terminal_measurements(measurements)
    isempty(tase_measurements) && return nothing
    if key === :iv_trace_table
        return _tase_iv_trace_table_result(tase_measurements; should_cancel)
    elseif key === :iv_device_summary
        return _tase_iv_device_summary_result(tase_measurements; should_cancel)
    end
    return nothing
end

function _tase_four_terminal_measurements(measurements)
    return [
        meas for meas in measurements
        if meas isa MeasurementInfo && meas.measurement_kind === :four_terminal_iv
    ]
end

function _tase_iv_trace_table_result(measurements::Vector{MeasurementInfo}; should_cancel::Union{Nothing,Function}=nothing)
    rows = NamedTuple[]
    for meas in measurements
        _check_plot_cancel(should_cancel)
        chip, facet, device_type, device_id = _tase_measurement_parts(meas)
        num_wires, wire_diameter_nm = _tase_wire_layout(device_id)
        df = read_tlm_4p(basename(meas.filepath), dirname(meas.filepath))
        for row in eachrow(df)
            push!(rows, (
                measurement_id=meas.filepath,
                filepath=meas.filepath,
                chip=chip,
                facet=facet,
                device_type=device_type,
                device_id=device_id,
                num_wires=num_wires,
                wire_diameter_nm=wire_diameter_nm,
                current_A=row.current_source,
                current_uA=row.current_source * 1e6,
                voltage_V=row.voltage_drop,
                voltage_mV=row.voltage_drop * 1e3,
                timestamp=meas.timestamp,
            ))
        end
    end
    table = DataFrame(rows)
    return AnalysisResult(
        :iv_trace_table,
        "TASE IV Trace Table",
        :iv_point,
        table,
        [
            (kind=:line, x=:voltage_mV, y=:current_uA, group=:measurement_id, color=:wire_diameter_nm, title="IV traces"),
            (kind=:facet_line, x=:voltage_mV, y=:current_uA, group=:measurement_id, facet=:num_wires, color=:wire_diameter_nm, title="IV traces by wire count"),
        ],
        Dict{Symbol,Any}(:project => :TASE, :measurement_count => length(measurements)),
    )
end

function _tase_iv_device_summary_result(measurements::Vector{MeasurementInfo}; should_cancel::Union{Nothing,Function}=nothing)
    rows = NamedTuple[]
    for meas in measurements
        _check_plot_cancel(should_cancel)
        chip, facet, device_type, device_id = _tase_measurement_parts(meas)
        num_wires, wire_diameter_nm = _tase_wire_layout(device_id)
        df = read_tlm_4p(basename(meas.filepath), dirname(meas.filepath))
        resistance_ohm, ci_low_ohm, ci_high_ohm = _fit_resistance_through_origin(df)
        conductance_s = isfinite(resistance_ohm) && resistance_ohm != 0 ? 1 / resistance_ohm : NaN
        g_over_g0 = conductance_s / TASE_CONDUCTANCE_QUANTUM_S
        g_per_wire_g0 = num_wires > 0 ? g_over_g0 / num_wires : NaN
        push!(rows, (
            measurement_id=meas.filepath,
            filepath=meas.filepath,
            chip=chip,
            facet=facet,
            device_type=device_type,
            device_id=device_id,
            num_wires=num_wires,
            wire_diameter_nm=wire_diameter_nm,
            R_ohm=resistance_ohm,
            R_ci_low_ohm=ci_low_ohm,
            R_ci_high_ohm=ci_high_ohm,
            G_S=conductance_s,
            G_over_G0=g_over_g0,
            G_per_wire_G0=g_per_wire_g0,
            timestamp=meas.timestamp,
        ))
    end
    table = DataFrame(rows)
    return AnalysisResult(
        :iv_device_summary,
        "TASE IV Device Summary",
        :device_measurement,
        table,
        [
            (kind=:scatter, x=:num_wires, y=:G_over_G0, color=:wire_diameter_nm, title="G/G0 vs number of wires"),
            (kind=:scatter, x=:num_wires, y=:G_per_wire_G0, color=:wire_diameter_nm, title="G per wire vs number of wires"),
            (kind=:scatter, x=:wire_diameter_nm, y=:G_per_wire_G0, color=:num_wires, title="G per wire vs diameter"),
        ],
        Dict{Symbol,Any}(:project => :TASE, :measurement_count => length(measurements)),
    )
end

function _tase_measurement_parts(meas::MeasurementInfo)
    parts = meas.device_info.location
    length(parts) == 4 || error("Expected 4 TASE device path parts, got $(length(parts)) for $(meas.filepath)")
    return parts[1], parts[2], parts[3], parts[4]
end

function _tase_wire_layout(device_id::AbstractString)
    length(device_id) >= 2 || error("Invalid TASE device id: $device_id")
    row = parse(Int, device_id[1:1])
    column = parse(Int, device_id[2:2])
    num_wires = TASE_WIRE_COUNT_MATRIX[row, column]
    diameter_nm = num_wires == 0 ? 0.0 : TASE_WIRE_DIAMETERS_NM[column]
    return num_wires, diameter_nm
end

function _fit_resistance_through_origin(df::DataFrame)
    mask = isfinite.(df.current_source) .& isfinite.(df.voltage_drop)
    currents = collect(df.current_source[mask])
    voltages = collect(df.voltage_drop[mask])
    isempty(currents) && return (NaN, NaN, NaN)
    denom = sum(abs2, currents)
    denom == 0 && return (NaN, NaN, NaN)
    resistance_ohm = sum(currents .* voltages) / denom
    length(currents) == 1 && return (resistance_ohm, resistance_ohm, resistance_ohm)

    residuals = voltages .- resistance_ohm .* currents
    sigma2 = sum(abs2, residuals) / max(length(currents) - 1, 1)
    stderr = sqrt(sigma2 / denom)
    delta = 1.96 * stderr
    return (resistance_ohm, resistance_ohm - delta, resistance_ohm + delta)
end
