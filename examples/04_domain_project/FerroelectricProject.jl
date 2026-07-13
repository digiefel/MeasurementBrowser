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

function clean_iv(table, metadata::AbstractDict)::DataFrame
    table = copy(table)
    sort!(table, :voltage_v)
    return table
end

function analyze_iv(table::DataFrame, metadata::AbstractDict)::Dict{Symbol,Any}
    return Dict{Symbol,Any}(
        :points => size(table, 1),
        :maximum_current_a => maximum(abs, table.current_a),
    )
end

function define_ferroelectric_project()::Project
    project = define_project("Semiconductor and ferroelectric characterization")

    register_item!(project, :iv;
        detect = (file::SourceFile) -> endswith(file.filename, "_iv.csv"),
        read = (file::SourceFile) -> begin
            parts = filename_parts(file.filename)
            DataItem(CSV.read(file.filepath, DataFrame);
                metadata=(chip=parts.chip, device=parts.device),
                label="$(parts.device) I-V",
                collection=[parts.chip, parts.device],
            )
        end,
        process = clean_iv,
        analyze = analyze_iv,
    )

    register_item!(project, :pund;
        detect = (file::SourceFile) -> endswith(file.filename, "_pund.csv"),
        read = (file::SourceFile) -> begin
            table = CSV.read(file.filepath, DataFrame)
            parts = filename_parts(file.filename)
            [
                DataItem(view(table, findall(==(pulse), table.pulse), :);
                    metadata=(chip=parts.chip, device=parts.device, pulse=String(pulse)),
                    id=String(pulse),
                    label="$(parts.device) $(String(pulse))",
                    collection=[parts.chip, parts.device],
                )
                for pulse in unique(table.pulse)
            ]
        end,
        process = (table, metadata) -> DataFrame(table),
        analyze = (table, metadata) -> Dict{Symbol,Any}(:points => size(table, 1)),
    )

    register_collection_analysis!(project, :pund;
        analyze = (data, metadata) -> Dict{Symbol,Any}(:pulses => length(data)),
    )

    return project
end

end
