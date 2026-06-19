Moving the code to the RuO2 project highlighted that the reworked API isn't particularly good. A good summary of its problems is the following snippet:

```
"""MeasurementBrowser project for the RuO2 experiment."""
struct RuO2Project <: AbstractProject end

include("ruo2/measurement_labels.jl")
include("ruo2/measurements.jl")
include("ruo2/measurement_analysis.jl")

"""Read the direct table represented by one RuO2 measurement."""
function load_source_data(
    ::RuO2Project,
    source_file::SourceFile;
    measurement::Union{Nothing,MeasurementInfo}=nothing,
)::DataFrame
    kind = measurement === nothing ?
        detect_kind(RUO2_PROJECT, source_file.filename) :
        measurement.measurement_kind
    if kind === :pund || kind === :pn
        if is_pund_fatigue_file(source_file.filepath)
            measurement === nothing && error("PUND fatigue data requires a logical measurement")
            fatigue_count = Int(measurement.parameters[:fatigue_idx])
            fatigue_df = read_pund_fatigue_file(source_file.filepath)
            return select_pund_fatigue_cycle(fatigue_df, fatigue_count)
        end
        return read_pund_file(source_file.filename, dirname(source_file.filepath))
    elseif kind === :wakeup_pn || kind === :wakeup_pund
        measurement === nothing && error("PUND wakeup data requires a logical measurement")
        wakeup_V = Float64(measurement.parameters[:wakeup_V])
        return _select_pund_wakeup_readout(source_file.filepath, wakeup_V)
    elseif kind === :iv || kind === :breakdown
        return read_ruo2_iv_sweep(source_file.filename, dirname(source_file.filepath))
    elseif kind === :tlm4p
        return read_ruo2_tlm_4p(source_file.filename, dirname(source_file.filepath))
    elseif kind === :cvsweep
        return read_cv_sweep(source_file.filename, dirname(source_file.filepath))
    end
    error("RuO2 source data API is not implemented for $kind")
end
```

The random struct definition with no fields (AbstractProject), and the inclusion of many files which all do stuff like the following:
```
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
```
These constant massive kind switches are the product of spaghetti code accumulated through bad AI slop accumulating for a while, but are also endemic to the architecture requiring many functions that are all called regardless of the type of measurement, and relying on symbols and hacks like filenames.

A new API, initially proposed by Fable, could look like the following:
```
using MeasurementBrowser
using Revise # should be leveraged/supported as much as possible

include("analysis_code.jl") # if not defined in the file itself, of course

project = define_project("RuO2", {other kwargs as needed})

workspace = open_workspace(project, rootpath)
browser = open_browser(workspace; wait=false) # non-blocking REPL

register_measurement!(project, :pund;
    detect = file -> occursin("PUND", file.name)::Bool,
    read = file -> read_pund(file)::DataFrame, # usually just file -> CSV.read(file, DataFrame)
    measurements = (file, data) -> pund_measurements(file, data)::Vector{MeasurementInfo},
    label = mi -> "FE PUND $(mi.params[:device]) $(mi.params[:temperature])",
    process = mi -> analyze_pund(mi)::DataFrame,                 # reads mi.data
    stats = mi -> pund_stats(mi)::Dict{Symbol,Any},              # reads processed mi.data
)
# re-calling it updates everything on the fly whether the browser is open or not,
# but using Revise to re-define the individual functions is also enough to update information on its own.
# the package must be smart enough to detect when that has happened, and whether a new rescan/cache rebuild is needed.

register_measurement!(project, :iv;
...) # repeat as many times as measurement types

register_device_stat!(project;
    measurement_kinds = [:pund, :wakeup_pn, :wakeup_pund],
    group_by = mi -> device_key(mi),   # OPTIONAL, defaults to the device leaf
    compute_stats = pund_history!,     # (ms::Vector{MeasurementInfo}) -> ms; fills each mi.stats in place
)

# OPTIONAL: global stats, e.g. how many devices do we have? Or anything like that. Might be useful?
register_global_stat!(project; ...)

register_plot!(project, :pund; 
    label = "PUND plot"
    setup = setup_pund_plot()::Figure, # or something else? GridLayout, or completely left to the user even?
    draw = draw_pund_plot(figure, measurements; ...) # supports multiple selections, and so on
    ...
)
```

`register_measurement!` defines a recipe on how to process a certain file, with a pipeline that tries to be logical, intuitive, and cache-friendly. The package then executes that recipe as cleverly as possible optimizing for responsivity and memory usage through caching and lazy on-demand off-thread work.

This new API should be superior in various ways. The metrics we care about are:
b) how much boilerplate is there to write? How much scaffolding is needed vs actual project scripting code that a user wants to write?
p) how clear is the purpose of that scaffolding? How much does it read like something that is intuitive to someone for whose scripting is just a means to an ends, who has no CS experience?
e) how extensible is it? Not in terms of how generic the verbs sound, but how easy does this public API translate to extensibility and generality?
r) how well does this scale with on-the-fly Revise.jl recompilation without restarting the browser, with a dynamic set of measurements and analysis that can get recomputed on the fly?

All of these should be fully satisfied to call the proposal complete.

Settled design:

Pipeline (per measurement kind). The callbacks run in a fixed order, each fed the previous output:
- per file, once (cached by file fingerprint):
    detect(file) -> Bool
    read(file) -> data                                   # the whole file, parsed once
    measurements(file, data) -> Vector{MeasurementInfo}  # one by default; many for expansion (fatigue cycles,
                                                          # wakeup voltages, paired breakdown devices). Omit when
                                                          # one file is one measurement.
- per measurement (cached by id):
    process(mi) -> processed                             # reads mi.data
    stats(mi) -> Dict                                    # reads processed mi.data
    label(mi) -> String
Because `read` runs once per file and its output is threaded into `measurements` and every `process` call, a
file is parsed once and sliced many times. That is the old multiple-reads-during-expansion problem solved by
construction.

Identity. `id` is engine-generated, not user-facing: derived from (filepath, kind, params).
`measurements` returns the existing `MeasurementInfo` (it already carries `parameters` and `stats`); the project
sets kind + params, the engine sets the id. Sibling measurements from one file must differ in params (they must,
or they are the same measurement).

No split of register_measurement!. The package can see when a function changes (and when it doesn't), so it
doesn't need the split to decide what to re-run.

register_device_stat! — cross-measurement stats that one measurement cannot compute alone (e.g. a cumulative
cycle count across a device's measurements in time order). Replaces the old six-dict
`compute_and_add_measurement_stats!`. The engine groups (by device, or `group_by`), hands the whole group to
`compute_stats`, and caches the result keyed by the group's file fingerprints (recompute a device only when one
of its files changes; the fold reads cached scalars, no disk). All accumulator state is plain local variables in
`compute_stats` — no init/step split. Worked example:
```
function pund_history!(ms::Vector{MeasurementInfo})
    sort!(ms; by = mi -> something(mi.timestamp, typemax(DateTime)))
    wakeup_counts = Dict{Tuple{Float64,Float64},Float64}()
    seen = Set{Tuple{String,Float64,Float64}}()
    wakeup_count = 0.0; fatigue_count = 0
    for mi in ms
        if mi.measurement_kind in (:wakeup_pn, :wakeup_pund)
            wV, wf = mi.parameters[:wakeup_V], mi.parameters[:wakeup_f]
            ek = (mi.filepath, wV, wf)
            if !(ek in seen)
                wakeup_counts[(wV,wf)] = get(wakeup_counts,(wV,wf),0.0) + mi.parameters[:wakeup_count]
                push!(seen, ek)
            end
            wakeup_count = wakeup_counts[(wV,wf)]
        elseif mi.measurement_kind === :pund && haskey(mi.parameters, :fatigue_idx)
            fatigue_count = mi.parameters[:fatigue_idx]
        end
        mi.stats[:wakeup_count] = wakeup_count
        mi.stats[:fatigue_count] = fatigue_count
    end
    return ms
end
```

Failures. Errors are handled and propagated, never swallowed. Every callback invocation is wrapped so a failure
records: the source file that produced the data, the measurement (when applicable), and the exact step that
failed (detect / read / measurements / process / stats / compute_stats / label / plot), plus the original error
and backtrace. A failure in one file or device group must not stop the others; failures are collected and
surfaced so any error is traceable to its file and its step.

Still open / deferred:
- register_global_stat! (project-wide aggregates, e.g. device counts) — shape deferred.
- group ordering for register_device_stat! — could take a `sort` argument like DataFrames `groupby`; irrelevant
  for now, `compute_stats` sorts internally.
- whether MeasurementInfo gains per-kind subtypes — pinned; pure author ergonomics, does not block the engine.
