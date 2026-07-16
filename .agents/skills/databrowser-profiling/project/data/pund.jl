using CSV
using DataFrames

const _PUND_COLUMNS_V1 = ("Time", "MeasResult1_value", "MeasResult2_value")
const _PUND_COLUMNS_V2 = ("Time_s", "Current_A", "Voltage_V")
const _PUND_WAKEUP_COLUMNS = ("Vmax_V", "Time_s", "Voltage_V", "Current_A")
const _PUND_FATIGUE_COLUMNS = ("Cycle", "Time_s", "Voltage_V", "Current_A")

# we need this right now because the old pund files have extra non-comment
# and non-table-header lines at the top that we want to skip
# if we completely drop those, we can just use CSV's header argument directly
# and avoid this custom logic
function _find_table_header_row(filepath::AbstractString, schemas)
    row = open(filepath, "r") do io
        for (row, line) in enumerate(eachline(io))
            columns = Set(String.(strip.(split(line, ','))))
            for schema in schemas
                issubset(schema, columns) && return row
            end
        end
        return nothing
    end
    row !== nothing && return row
    error("Could not find a supported PUND data header in $filepath")
end

"""
Old and new PUND exports use different column names, but downstream analysis expects the
same four time/current/voltage columns.
"""
function read_pund_file(filename, workdir=".")
    filepath = joinpath(workdir, filename)
    header_row = _find_table_header_row(filepath, (_PUND_COLUMNS_V1, _PUND_COLUMNS_V2))
    source = CSV.read(filepath, DataFrame; header=header_row, ntasks=1)

    if all(column -> column in names(source), _PUND_COLUMNS_V1)
        # current_time = "MeasResult1_time" in names(source) ? source.MeasResult1_time : source.Time
        # voltage_time = "MeasResult2_time" in names(source) ? source.MeasResult2_time : source.Time
        rename!(source, 
            [
                :Time => :time,
                :MeasResult1_value => :current,
                :MeasResult2_value => :voltage,
                # :MeasResult1_time => :current_time,
                # :MeasResult2_time => :voltage_time
            ])
    elseif all(column -> column in names(source), _PUND_COLUMNS_V2)
        rename!(source,
            [
                :Time_s => :time,
                :Current_A => :current,
                :Voltage_V => :voltage,
                # :Time_s => :current_time,
                # :Time_s => :voltage_time
            ])
    end
    return source[!, [:time, :current, :voltage]]
end

"""
Wakeup files contain one readout table with one block per configured wakeup amplitude.
"""
function read_pund_wakeup_file(filepath::AbstractString)
    # header_row = _find_table_header_row(filepath, (_PUND_WAKEUP_COLUMNS,))
    source = CSV.read(
        filepath,
        DataFrame;
        comment="#",
        # header=header_row,
        select=collect(Symbol.(_PUND_WAKEUP_COLUMNS)),
        types=Dict(
            :Vmax_V => Float64,
            :Time_s => Float64,
            :Voltage_V => Float64,
            :Current_A => Float64,
        ),
        ntasks=1,
    )
    rename!(source, [:Time_s => :time, :Voltage_V => :voltage, :Current_A => :current, :Vmax_V => :wakeup_V])
    return source
end

function read_pund_fatigue_file(filepath::AbstractString)
    # header_row = _find_table_header_row(filepath, (_PUND_FATIGUE_COLUMNS,))
    source = CSV.read(
        filepath,
        DataFrame;
        comment="#",
        # header=header_row,
        select=collect(Symbol.(_PUND_FATIGUE_COLUMNS)),
        types=Dict(
            :Cycle => Int,
            :Time_s => Float64,
            :Voltage_V => Float64,
            :Current_A => Float64,
        ),
        # The scan reads files in parallel, so parse single-threaded to avoid CPU oversubscription
        # (CSV multithreads each file otherwise; with 705 large fatigue files that just thrashes).
        ntasks=1,
    )
    rename!(source, [:Time_s => :time, :Voltage_V => :voltage, :Current_A => :current, :Cycle => :cycle])
    return source
end

"""Return one fatigue readout as a standard PUND waveform."""
function select_pund_fatigue_cycle(fatigue_df::DataFrame, cycle::Integer)
    mask = fatigue_df.cycle .== Int(cycle)
    any(mask) || error("PUND fatigue data has no rows for cycle $cycle")
    return @view fatigue_df[mask, :]
end
