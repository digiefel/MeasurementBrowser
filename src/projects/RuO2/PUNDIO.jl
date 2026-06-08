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
    source = CSV.read(filepath, DataFrame; header=header_row)
    if all(column -> column in names(source), _PUND_COLUMNS_V1)
        current_time = "MeasResult1_time" in names(source) ? source.MeasResult1_time : source.Time
        voltage_time = "MeasResult2_time" in names(source) ? source.MeasResult2_time : source.Time
        return DataFrame(
            time=source.Time,
            current=source.MeasResult1_value,
            voltage=source.MeasResult2_value,
            current_time=current_time,
            voltage_time=voltage_time,
        )
    end

    return DataFrame(
        time=source.Time_s,
        current=source.Current_A,
        voltage=source.Voltage_V,
        current_time=source.Time_s,
        voltage_time=source.Time_s,
    )
end

"""
Wakeup files contain one readout table with one block per configured wakeup amplitude.
"""
function read_pund_wakeup_file(filepath::AbstractString)
    header_row = _find_table_header_row(filepath, (_PUND_WAKEUP_COLUMNS,))
    source = CSV.read(
        filepath,
        DataFrame;
        header=header_row,
        select=collect(Symbol.(_PUND_WAKEUP_COLUMNS)),
        types=Dict(
            :Vmax_V => Float64,
            :Time_s => Float64,
            :Voltage_V => Float64,
            :Current_A => Float64,
        ),
    )
    return DataFrame(
        wakeup_V=source.Vmax_V,
        time=source.Time_s,
        voltage=source.Voltage_V,
        current=source.Current_A,
    )
end

function read_pund_fatigue_file(filepath::AbstractString)
    header_row = _find_table_header_row(filepath, (_PUND_FATIGUE_COLUMNS,))
    _check_cancel()
    source = CSV.read(
        filepath,
        DataFrame;
        header=header_row,
        select=[:Cycle, :Time_s, :Voltage_V, :Current_A],
        types=Dict(
            :Cycle => Int,
            :Time_s => Float64,
            :Voltage_V => Float64,
            :Current_A => Float64,
        ),
    )
    _check_cancel()
    return DataFrame(
        cycle=source.Cycle,
        time=source.Time_s,
        voltage=source.Voltage_V,
        current=source.Current_A,
    )
end

function _select_pund_fatigue_cycle(fatigue_df::DataFrame, cycle::Integer)
    mask = fatigue_df.cycle .== Int(cycle)
    any(mask) || error("PUND fatigue data has no rows for cycle $cycle")
    return DataFrame(
        time=fatigue_df.time[mask],
        current=fatigue_df.current[mask],
        voltage=fatigue_df.voltage[mask],
    )
end
