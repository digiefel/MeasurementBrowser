using CSV
using DataFrames

function _load_ruo2_pund_fatigue_file(
    filepath::AbstractString;
    should_cancel::Union{Nothing,Function}=nothing,
)
    header_row = nothing
    header = nothing
    open(filepath, "r") do io
        for (row, raw_line) in enumerate(eachline(io))
            _check_cancel(should_cancel)
            line = strip(raw_line)
            isempty(line) && continue
            if startswith(line, "Cycle,")
                header_row = row
                header = line
                break
            end
        end
    end
    header_row !== nothing || error(
        "Fatigue file '$filepath' is missing a data header with columns Cycle, Time_s, Voltage_V, Current_A",
    )
    _ruo2_fatigue_columns(filepath, header)

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
    _check_cancel(should_cancel)
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
