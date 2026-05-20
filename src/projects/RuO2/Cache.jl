using DataFrames

function project_cache_semantic_fields(::RuO2Project)
    return Dict{Symbol,Vector{Symbol}}(
        :measurement => [:measurement_id, :measurement_kind, :timestamp, :source_file],
        :device => [:device_key, :device_path, :area_um2, :x, :y],
        :signal => [:time, :voltage, :current, :cycle, :frequency],
        :selection => [:wakeup_V, :fatigue_idx],
        :summary => [:wakeup_count, :wakeup_f, :wakeup_V,
                     :fatigue_count, :fatigue_f, :fatigue_V,
                     :V_base, :V_min, :V_max, :V_amp],
    )
end

function project_cache_write_file_payload!(
    file_group,
    ::RuO2Project,
    file::SourceFile,
    measurements::Vector{MeasurementInfo};
    should_cancel::Union{Nothing,Function}=nothing,
)
    if !isempty(measurements) && is_pund_fatigue_file(file.filepath)
        _ruo2_write_cached_pund_fatigue_file!(file_group, file, measurements; should_cancel)
        return nothing
    end
    invoke(
        project_cache_write_file_payload!,
        Tuple{Any,AbstractProject,SourceFile,Vector{MeasurementInfo}},
        file_group,
        RUO2_PROJECT,
        file,
        measurements;
        should_cancel,
    )
    return nothing
end

function project_cache_write_measurement_payload!(
    measurement_group,
    ::RuO2Project,
    measurement::MeasurementInfo;
    should_cancel::Union{Nothing,Function}=nothing,
)
    params = _measurement_parameters(measurement)
    loaded = load_plot_for_file(
        RUO2_PROJECT,
        measurement.filepath,
        measurement.measurement_kind;
        device_params=params,
        should_cancel,
    )
    analyzed = analyze_plot_for_file(
        RUO2_PROJECT,
        measurement.measurement_kind,
        loaded;
        device_params=params,
        should_cancel,
    )
    analyzed === nothing && error("RuO2 cache transform produced no plot payload for $(measurement.unique_id)")
    _ruo2_write_cached_analyzed_plot!(measurement_group, measurement.measurement_kind, analyzed)
    return nothing
end

function project_cache_read_plot_payload(
    ::RuO2Project,
    measurement::MeasurementInfo,
    file_group,
    measurement_group,
)
    if haskey(file_group, "status") && read(file_group["status"]) == "error"
        msg = haskey(file_group, "error_message") ? read(file_group["error_message"]) : "unknown cache error"
        error("Cached file failed during cache build: $msg")
    end
    haskey(measurement_group, "plot") ||
        throw(ProjectCacheInvalidError("", "cached measurement '$(measurement.unique_id)' is missing plot payload"))
    plot_group = measurement_group["plot"]
    source = read(plot_group["source"])
    if source == "analyzed"
        return _ruo2_read_cached_analyzed_plot(plot_group, measurement.measurement_kind, "")
    elseif source == "file_pund_fatigue"
        return _ruo2_read_cached_pund_fatigue_cycle(file_group, measurement)
    end
    error("Unsupported RuO2 cached plot source '$source' for $(measurement.unique_id)")
end

function _ruo2_write_cached_analyzed_plot!(measurement_group, kind::Symbol, analyzed)
    plot_group = _replace_group(measurement_group, "plot")
    _write_dataset!(plot_group, "source", "analyzed"; compress=false)
    _write_dataset!(plot_group, "measurement_kind", String(kind); compress=false)
    hasproperty(analyzed, :df) || error("RuO2 analyzed payload for '$kind' is missing df")
    _write_dataframe!(plot_group, "df", analyzed.df)

    scalars = Dict{Symbol,Any}()
    arrays = Dict{Symbol,Any}()
    for name in propertynames(analyzed)
        name === :df && continue
        value = getproperty(analyzed, name)
        if name === :pulse_groups || name === :debug_boundaries || name === :debug_labels
            continue
        elseif _cache_scalar_supported(value)
            scalars[name] = value
        elseif value isa AbstractVector && _cache_vector_supported(value)
            arrays[name] = collect(value)
        else
            error("Unsupported RuO2 cached plot field '$name' of type $(typeof(value))")
        end
    end
    _write_parameters!(plot_group, "scalars", scalars)
    arrays_group = _replace_group(plot_group, "arrays")
    _write_symbol_vector!(arrays_group, "names", collect(keys(arrays)))
    for (name, values) in arrays
        _write_dataset!(arrays_group, String(name), values)
    end
    return nothing
end

function _cache_scalar_supported(value)
    return value === nothing ||
        value isa Bool ||
        value isa Integer ||
        value isa AbstractFloat ||
        value isa AbstractString ||
        value isa Symbol ||
        value isa Date ||
        value isa DateTime
end

function _cache_vector_supported(value::AbstractVector)
    isempty(value) && return true
    T = eltype(value)
    return T <: Bool || T <: Integer || T <: AbstractFloat || T <: AbstractString
end

function _ruo2_read_cached_analyzed_plot(plot_group, kind::Symbol, cache_path::AbstractString)
    df = _read_dataframe(plot_group, "df", cache_path)
    scalars = _read_parameters(plot_group, "scalars", cache_path)
    arrays = Dict{Symbol,Any}()
    if haskey(plot_group, "arrays")
        arrays_group = plot_group["arrays"]
        for name in _read_symbol_vector(arrays_group, "names", cache_path)
            arrays[name] = read(arrays_group[String(name)])
        end
    end

    title = String(get(scalars, :title, ""))
    if kind === :pund || kind === :pn || kind === :wakeup_pn || kind === :wakeup_pund
        default_group_size = (kind === :pn || kind === :wakeup_pn) ? 1 : 5
        group_size = Int(get(scalars, :pulse_group_size, default_group_size))
        return (
            df=df,
            title=title,
            area_um2=get(scalars, :area_um2, nothing),
            pulse_groups=_ruo2_cached_pulse_groups(df, group_size),
            pulse_group_size=group_size,
            remnant_y_label=String(get(scalars, :remnant_y_label, "Switching Charge (pC)")),
            debug=Bool(get(scalars, :debug, false)),
        )
    elseif kind === :iv || kind === :breakdown || kind === :unknown || kind === nothing
        return (df=df, title=title)
    elseif kind === :tlm4p
        return (
            df=df,
            title=title,
            fit_resistance_ohm=Float64(get(scalars, :fit_resistance_ohm, NaN)),
            fit_resistance_kohm=Float64(get(scalars, :fit_resistance_kohm, NaN)),
            fit_current_uA=get(arrays, :fit_current_uA, Float64[]),
            fit_voltage_mV=get(arrays, :fit_voltage_mV, Float64[]),
            rho_sheet=Float64(get(scalars, :rho_sheet, NaN)),
        )
    elseif kind === :cvsweep
        return (
            df=df,
            title=title,
            frequencies_Hz=get(arrays, :frequencies_Hz, Float64[]),
        )
    end
    return (df=df, title=title)
end

function _ruo2_cached_pulse_groups(df::DataFrame, group_size::Int)
    hasproperty(df, :pulse_idx) || return Tuple{Int,BitVector}[]
    isempty(df.pulse_idx) && return Tuple{Int,BitVector}[]
    max_pulse = maximum(df.pulse_idx)
    pulse_groups = Tuple{Int,BitVector}[]
    for rep in 1:(max_pulse ÷ group_size)
        pulse_range = (rep - 1) * group_size + 1:rep * group_size
        mask = BitVector([pulse in pulse_range for pulse in df.pulse_idx])
        any(mask) && push!(pulse_groups, (rep, mask))
    end
    return pulse_groups
end

function _ruo2_write_cached_pund_fatigue_file!(
    file_group,
    file::SourceFile,
    measurements::Vector{MeasurementInfo};
    should_cancel::Union{Nothing,Function}=nothing,
)
    df = _load_ruo2_pund_fatigue_file(file.filepath; should_cancel=should_cancel)
    signals_group = _ensure_group(file_group, "signals")
    _write_dataframe!(signals_group, "pund_fatigue", df)

    measurements_group = file_group["measurements"]
    for measurement in measurements
        measurement_group = measurements_group[_measurement_group_key(measurement.unique_id)]
        plot_group = _replace_group(measurement_group, "plot")
        _write_dataset!(plot_group, "source", "file_pund_fatigue"; compress=false)
        _write_dataset!(plot_group, "measurement_kind", String(:pund); compress=false)
    end
    return nothing
end

function _ruo2_read_cached_pund_fatigue_cycle(file_group, measurement::MeasurementInfo)
    haskey(file_group, "signals") && haskey(file_group["signals"], "pund_fatigue") ||
        throw(ProjectCacheInvalidError("", "cached fatigue file is missing pund_fatigue signals"))
    full_df = _read_dataframe(file_group["signals"], "pund_fatigue", "")
    cycle = Int(measurement.parameters[:fatigue_idx])
    df = _select_pund_fatigue_cycle(full_df, cycle)
    params = _measurement_parameters(measurement)
    loaded = (
        df=df,
        title=measurement.clean_title,
        area_um2=get(params, :area_um2, nothing),
        debug=false,
    )
    return analyze_plot_for_file(RUO2_PROJECT, :pund, loaded; device_params=params)
end
