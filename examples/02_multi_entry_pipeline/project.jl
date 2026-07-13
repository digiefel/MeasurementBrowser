using CSV
using DataBrowser
using DataFrames: DataFrame
using GLMakie: Axis, Figure, contents, lines!

length(ARGS) == 1 || error("Usage: julia --project project.jl DATA_DIRECTORY")

project = define_project("Fatigue cycles")

function clean_cycle(table, metadata::Dict)::DataFrame
    table = DataFrame(table)
    sort!(table, :time_s)
    return table
end

function analyze_cycle(table::DataFrame, metadata::Dict)::Dict{Symbol,Any}
    return Dict{Symbol,Any}(
        :points => size(table, 1),
        :maximum_voltage_v => maximum(abs, table.voltage_v),
        :maximum_current_a => maximum(abs, table.current_a),
    )
end

register_item!(project, :cycles;
    detect = (file::SourceFile) -> endswith(file.filename, "_fatigue.csv"),
    read = (file::SourceFile) -> CSV.read(file.filepath, DataFrame),
    entries = (table::DataFrame, metadata::Dict) -> [
        (
            data=view(table, findall(==(cycle), table.cycle), :),
            metadata=Dict(:cycle => Int(cycle)),
        )
        for cycle in unique(table.cycle)
    ],
    id = (table, metadata::Dict) -> metadata[:cycle],
    label = (table, metadata::Dict) -> "Cycle $(metadata[:cycle])",
    collection = (table, metadata::Dict) -> ["Fatigue"],
    process = clean_cycle,
    analyze = analyze_cycle,
)

register_collection_analysis!(project, :cycles;
    analyze = (data, metadata) -> Dict{Symbol,Any}(
        :cycles => length(data),
        :last_cycle => maximum(row[:cycle] for row in metadata),
    ),
)

register_plot!(project, :cycles;
    label = "Current",
    setup = (data, metadata) -> begin
        figure = Figure()
        Axis(figure[1, 1]; xlabel="Time (s)", ylabel="Current (A)")
        figure
    end,
    draw = (figure, data, metadata) -> begin
        axis = only(contents(figure[1, 1]))
        for (table, values) in zip(data, metadata)
            lines!(axis, table.time_s, table.current_a; label="Cycle $(values[:cycle])")
        end
        nothing
    end,
)

workspace = open_workspace(project, only(ARGS))
open_browser(workspace)
