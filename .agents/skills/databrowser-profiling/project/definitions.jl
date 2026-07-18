#  RuO2 measurement project.
#
#  Launch:
#      julia> using Revise
#      julia> includet("browser.jl")          # scans and opens the browser window
#  Edit any analysis function and Revise reloads it live; re-run a register_item! line to
#  change a recipe.

# using Revise  # not wanted in profiling runs
using DataBrowser
import DataBrowser: collection, id, item_data, item_label, kind, metadata
using CSV
using DataFrames
using Statistics: mean

include("data/pund.jl")       # PUND readers and table normalization
include("data/cv.jl")         # C-V reader
include("data/iv.jl")         # IV/TLM/breakdown reader

contains_ci(s::AbstractString, sub::AbstractString)::Bool = endswith(s, ".csv") && occursin(lowercase(sub), lowercase(s))

const RUO2_KIND_LABELS = Dict(
    :pund => "FE PUND",
    :wakeup => "PUND Wakeup",
    :pund_fatigue => "PUND Fatigue",
    :cvsweep => "C-V Sweep",
    :iv => "I-V Sweep",
    :tlm4p => "TLM 4-Point",
    :breakdown => "Breakdown",
    :wgfmu_sampling => "WGFMU Sampling",
)

is_iv_file(file)::Bool =
    (contains_ci(file.filename, "IVSweep") || contains_ci(file.filename, "FourTerminal")) &&
    !contains_ci(file.filename, "_TLM_") && !contains_ci(file.filename, "Breakdown") &&
    !contains_ci(file.filename, "FeCapBD") && !contains_ci(file.filename, "FeCap_BD")

function ruo2_location(filename::AbstractString)::Vector{String}
    m = match(r"^(RuO2.+?)_\d{8}_\d{6}_", String(filename))
    m === nothing && error("Unrecognized RuO2 filename: $filename")
    identifier = m.captures[1]
    structured = match(r"^(RuO2test_?[A-Z0-9]+)_([XVI]+)_(.+)_([A-Z][0-9]+(?:[A-Z][0-9]+)*)$", identifier)
    structured !== nothing && return [
        String(structured[1]),
        String(structured[2]),
        replace(String(structured[3]), "_" => ""),
        String(structured[4]),
    ]
    parts = split(identifier, '_')
    length(parts) >= 4 || error("Invalid RuO2 identifier: $identifier")
    return [String(join(parts[1:end-3], "_")), String(parts[end-2]), String(parts[end-1]), String(parts[end])]
end

const RUO2_LOCATION_KEYS = (:chip, :site, :subsite, :device)

function ruo2_file_metadata(file)::Dict{Symbol,Any}
    location = ruo2_location(file.filename)
    md = Dict{Symbol,Any}(
        :filepath => file.filepath,
        :filename => file.filename,
    )
    for (key, value) in zip(RUO2_LOCATION_KEYS, location)
        md[key] = value
    end
    temp = match(r"_(\d+(?:\.\d+)?)K(?:_|\.|$)", file.filename)
    temp === nothing || (md[:temperature_K] = parse(Float64, temp.captures[1]))
    file.timestamp === nothing || (md[:timestamp] = file.timestamp)
    return md
end

function ruo2_entry(data, metadata; params...)
    additions = Dict{Symbol,Any}(params)

    # Labels are built during indexing, before DataBrowser runs `process` and `analyze`.
    # Voltage is measurement identity here, so keep the raw envelope on the item metadata.
    if !haskey(metadata, :voltage) && !haskey(additions, :voltage) &&
            !haskey(metadata, :V_base) && !haskey(additions, :V_base) &&
            hasproperty(data, :voltage) && !isempty(data.voltage)
        V = data.voltage
        n_pre = min(9, length(V))
        round_one(x) = (y = round(x; digits=1); iszero(y) ? 0.0 : y)
        V_base = round_one(mean(V[1:n_pre]))
        additions[:V_base] = V_base
        additions[:V_amp] = round_one(maximum(abs.(V .- V_base)))
    end
    return (data=data, metadata=additions)
end

function ruo2_read(data, metadata)
    entry = ruo2_entry(data, metadata)
    return (data=data, metadata=merge(Dict{Symbol,Any}(metadata), entry.metadata))
end

ruo2_collection(_data, metadata)::Vector{String} =
    String[String(metadata[key]) for key in RUO2_LOCATION_KEYS]

function ruo2_id(_data, metadata)::String
    id_parts = String[]
    for key in (:cycle, :device, :segment, :voltage)
        haskey(metadata, key) && push!(id_parts, "$key=$(metadata[key])")
    end
    return isempty(id_parts) ? "item" : join(id_parts, ",")
end

strcat(parts...) = join(filter(!isempty, string.(parts)), " ")
ruo2_timestamp_label(metadata) = (t = get(metadata, :timestamp, nothing); t === nothing ? "" : string(t))
ruo2_voltage_label(metadata) =
    haskey(metadata, :voltage) ? "$(metadata[:voltage]) V" :
    haskey(metadata, :V_base) && haskey(metadata, :V_amp) ?
        "$(metadata[:V_base]) ± $(metadata[:V_amp]) V" : ""
ruo2_frequency_label(freq_Hz) = isfinite(freq_Hz) ? (freq_Hz >= 1e6 ? "$(round(freq_Hz / 1e6; digits=3)) MHz" : freq_Hz >= 1e3 ? "$(round(freq_Hz / 1e3; digits=3)) kHz" : "$(round(freq_Hz; digits=3)) Hz") : ""
ruo2_temperature_label(metadata) = haskey(metadata, :temperature_K) ? "$(metadata[:temperature_K])K" : ""

include("analysis/pund.jl")   # PUND, PN, and fatigue analysis
include("analysis/stats.jl")  # stats callbacks and device history
include("analysis/tlm.jl")    # combined TLM fitting helpers
include("plots/common.jl")    # shared plot formatting and unit conversion
include("plots/summary_views.jl")
include("plots/measurement_views.jl")

const PROJECT = define_project(
    "DataBrowserProfilingRuO2";
    description="Ferroelectric RuO2test profiling workload",
)

# PUND fatigue — one file expands to one PUND per fatigue cycle.
register_item!(PROJECT, :pund_fatigue;
    detect = file -> contains_ci(file.filename, "PUND_Fatigue"),
    read   = file -> (
        data=read_pund_fatigue_file(file.filepath),
        metadata=ruo2_file_metadata(file),
    ),
    entries = (data, metadata) -> [
        let cycle_data = select_pund_fatigue_cycle(data, Int(c))
            ruo2_entry(
                cycle_data,
                metadata;
                cycle=Int(c),
                frequency_kHz=pund_frequency_kHz(cycle_data),
            )
        end
        for c in sort(unique(data.cycle))
    ],
    process = (data, _metadata) -> analyze_pund(data),
    analyze = pund_stats,
    label   = (_data, metadata) -> strcat(ruo2_timestamp_label(metadata), "PUND", ruo2_voltage_label(metadata), ruo2_frequency_label(get(metadata, :frequency_kHz, NaN) * 1e3), "cycle $(get(metadata, :cycle, "?"))"),
    collection = ruo2_collection,
    id = ruo2_id,
)

# PUND wakeup — one file expands to one readout per (voltage, PN/PUND segment).
register_item!(PROJECT, :wakeup;
    detect = file -> contains_ci(file.filename, "PUND_Wakeup"),
    read   = read_wakeup,
    entries = wakeup_measurements,
    process = process_wakeup,
    analyze = pund_stats,
    label   = (_data, metadata) -> strcat(ruo2_timestamp_label(metadata), get(metadata, :segment, :pund) === :pn ? "Wakeup PN" : "Wakeup PUND", ruo2_voltage_label(metadata), ruo2_temperature_label(metadata)),
    collection = ruo2_collection,
    id = ruo2_id,
)

# Single PUND.
register_item!(PROJECT, :pund;
    detect = file -> contains_ci(file.filename, "pund"),
    read   = file -> ruo2_read(
        read_pund_file(file.filename, dirname(file.filepath)),
        ruo2_file_metadata(file),
    ),
    process = (data, _metadata) -> analyze_pund(data),
    analyze = pund_stats,
    label   = (_data, metadata) -> strcat(ruo2_timestamp_label(metadata), "PUND", ruo2_voltage_label(metadata), ruo2_temperature_label(metadata)),
    collection = ruo2_collection,
    id = ruo2_id,
)

# Breakdown — paired-device files split into one measurement per device.
register_item!(PROJECT, :breakdown;
    detect = file -> contains_ci(file.filename, "Breakdown") || contains_ci(file.filename, "break") ||
        contains_ci(file.filename, "FeCapBD") || contains_ci(file.filename, "FeCap_BD"),
    read   = file -> (
        data=read_iv_sweep(file.filepath),
        metadata=ruo2_file_metadata(file),
    ),
    entries = breakdown_measurements,
    process = (data, _metadata) -> iv_resistance(data),
    analyze = iv_stats,
    label   = (_data, metadata) -> strcat(ruo2_timestamp_label(metadata), "Breakdown", ruo2_temperature_label(metadata)),
    collection = ruo2_collection,
    id = ruo2_id,
)

# C-V sweep.
register_item!(PROJECT, :cvsweep;
    detect = file -> contains_ci(file.filename, "CVSweep"),
    read   = file -> (
        data=read_cv_sweep(file.filename, dirname(file.filepath)),
        metadata=ruo2_file_metadata(file),
    ),
    label  = (_data, metadata) -> strcat(ruo2_timestamp_label(metadata), "C-V Sweep", ruo2_temperature_label(metadata)),
    collection = ruo2_collection,
    id = ruo2_id,
)

# WGFMU sampling — raw two-channel current trace.
register_item!(PROJECT, :wgfmu_sampling;
    detect = file -> contains_ci(file.filename, "WGFMU_Sampling"),
    read   = file -> (
        data=CSV.read(file.filepath, DataFrame; comment="#", ntasks=1),
        metadata=ruo2_file_metadata(file),
    ),
    label   = (_data, metadata) -> strcat(ruo2_timestamp_label(metadata), "WGFMU Sampling", ruo2_temperature_label(metadata)),
    collection = ruo2_collection,
    id = ruo2_id,
)

# TLM 4-point — IVSweep-format files in a TLM subsite.
register_item!(PROJECT, :tlm4p;
    detect = file -> contains_ci(file.filename, "_TLM_") && !contains_ci(file.filename, "WGFMU_Sampling"),
    read   = file -> (
        data=read_iv_sweep(file.filepath),
        metadata=ruo2_file_metadata(file),
    ),
    process = (data, _metadata) -> iv_resistance(data),
    analyze = iv_stats,
    label   = (_data, metadata) -> strcat(ruo2_timestamp_label(metadata), "TLM 4-Point", ruo2_temperature_label(metadata)),
    collection = ruo2_collection,
    id = ruo2_id,
)

# Current-voltage sweep (FeCap / four-terminal).
register_item!(PROJECT, :iv;
    detect = is_iv_file,
    read   = file -> (
        data=read_iv_sweep(file.filepath),
        metadata=ruo2_file_metadata(file),
    ),
    process = (data, _metadata) -> iv_resistance(data),
    analyze = iv_stats,
    label   = (_data, metadata) -> strcat(ruo2_timestamp_label(metadata), "I-V Sweep", ruo2_temperature_label(metadata)),
    collection = ruo2_collection,
    id = ruo2_id,
)

register_plot!(PROJECT, :pund;
    label = "PUND",
    setup = setup_pund_plot,
    draw = draw_pund_plot,
)

register_plot!(PROJECT, :wakeup;
    label = "Wakeup PUND/PN",
    setup = setup_pund_plot,
    draw = draw_pund_plot,
)

register_plot!(PROJECT, :pund_fatigue;
    label = "PUND",
    setup = setup_pund_plot,
    draw = draw_pund_plot,
)

register_plot!(PROJECT, :pund_fatigue;
    label = "PUND Fatigue",
    setup = setup_fatigue_plot,
    draw = draw_fatigue_plot,
)

register_plot!(PROJECT, :iv;
    label = "I-V",
    setup = setup_iv_plot,
    draw = draw_iv_plot,
)

register_plot!(PROJECT, :breakdown;
    label = "Breakdown I-V",
    setup = setup_iv_plot,
    draw = draw_iv_plot,
)

register_plot!(PROJECT, :tlm4p;
    label = "TLM 4-point",
    setup = setup_tlm_plot,
    draw = draw_tlm_plot,
)

register_plot!(PROJECT, :tlm4p;
    label = "TLM analysis",
    setup = setup_tlm_analysis_plot,
    draw = draw_tlm_analysis_plot,
)

register_plot!(PROJECT, :cvsweep;
    label = "C-V Sweep",
    setup = setup_cv_plot,
    draw = draw_cv_plot,
)

register_plot!(PROJECT, :wgfmu_sampling;
    label = "WGFMU Sampling",
    setup = setup_wgfmu_sampling_plot,
    draw = draw_wgfmu_sampling_plot,
)
