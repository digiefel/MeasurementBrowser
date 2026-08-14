module FerroelectricProject

using CSV
using DataBrowser
using DataFrames: DataFrame

export define_ferroelectric_project

function filename_parts(filename::AbstractString)
    stem = splitext(basename(filename))[1]
    parts = split(stem, '_')
    length(parts) >= 2 || error("Expected CHIP_DEVICE_... filename, got $filename")
    return (chip=parts[1], device=parts[2])
end

function clean_iv(table, metadata::Dict)::DataFrame
    table = copy(table)
    sort!(table, :voltage_v)
    return table
end

function analyze_iv(table::DataFrame, metadata::Dict)::Dict{Symbol,Any}
    return Dict{Symbol,Any}(
        :points => size(table, 1),
        :maximum_current_a => maximum(abs, table.current_a),
    )
end

function define_ferroelectric_project()::Project
    project = define_project("Semiconductor and ferroelectric characterization")

    register_item!(project, :iv;
        # `detect` and `read` receive a source item, whatever kind of source produced it. Reach it
        # through the contract — `label`, `source_item_path` — rather than any one source's fields.
        detect = source_item -> endswith(label(source_item), "_iv.csv"),
        read = source_item -> begin
            parts = filename_parts(label(source_item))
            (
                data=CSV.read(source_item_path(source_item), DataFrame),
                metadata=Dict(:chip => parts.chip, :device => parts.device),
            )
        end,
        label = (table, metadata::Dict) -> "$(metadata[:device]) I-V",
        collection = (table, metadata::Dict) -> [metadata[:chip], metadata[:device]],
        process = clean_iv,
        analyze = analyze_iv,
    )

    register_item!(project, :pund;
        detect = source_item -> endswith(label(source_item), "_pund.csv"),
        read = source_item -> begin
            table = CSV.read(source_item_path(source_item), DataFrame)
            parts = filename_parts(label(source_item))
            (
                data=table,
                metadata=Dict(:chip => parts.chip, :device => parts.device),
            )
        end,
        entries = (table::DataFrame, metadata::Dict) -> [
            (
                data=view(table, findall(==(pulse), table.pulse), :),
                metadata=Dict(:pulse => String(pulse)),
            )
            for pulse in unique(table.pulse)
        ],
        id = (table, metadata::Dict) -> metadata[:pulse],
        label = (table, metadata::Dict) ->
            "$(metadata[:device]) $(metadata[:pulse])",
        collection = (table, metadata::Dict) -> [metadata[:chip], metadata[:device]],
        process = (table, metadata) -> DataFrame(table),
        analyze = (table, metadata) -> Dict{Symbol,Any}(:points => size(table, 1)),
    )

    register_collection_analysis!(project, :pund;
        analyze = (data, metadata) -> Dict{Symbol,Any}(:pulses => length(data)),
    )

    return project
end

end
