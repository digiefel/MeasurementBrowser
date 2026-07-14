using DataBrowser
using DataBrowser: items_for_file, setup_plot, plot_data!
using Test
using DataFrames: DataFrame, nrow
using GLMakie: Figure, Axis, lines!, contents

const MB = DataBrowser

# A realistic registry project (TASE four-terminal IV) exercised against the bundled fixtures. TASE
# itself now lives as an external registry project; this keeps coverage of a real read + plot path.
const REGEX_TASE =
    r"^([^_]+)_([^_]+)_([^_]+)_(\d+)_\d{8}_\d{6}_\d+K_FourTerminalIV\.csv"i

function _read_tase(filepath::AbstractString)::DataFrame
    lines = readlines(filepath)
    data_start = 1
    has_header = false
    for (i, line) in enumerate(lines)
        if occursin("Current_A,VoltageHigh_V,VoltageLow_V", line)
            data_start = i + 1
            has_header = true
            break
        end
        if occursin(r"^-?\d+\.?\d*[eE]?-?\d*,(-?\d+\.?\d*[eE]?-?\d*,){2,}", line)
            data_start = i
            break
        end
    end
    current = Float64[]
    voltage = Float64[]
    for line in lines[data_start:end]
        isempty(strip(line)) && continue
        parts = split(line, ',')
        try
            push!(current, parse(Float64, parts[1]))
            push!(voltage, has_header ?
                parse(Float64, parts[2]) - parse(Float64, parts[3]) :
                parse(Float64, filter(p -> !isempty(strip(p)), parts)[2]))
        catch
            continue
        end
    end
    return DataFrame(i=current, v=voltage)
end

function _tase_project()
    project = MB.define_project("TASE")
    MB.register_item!(
        project,
        :four_terminal_iv;
        detect=file -> match(REGEX_TASE, file.filename) !== nothing,
        read=file -> _read_tase(file.filepath),
        collection=function (_data, metadata)
            m = match(REGEX_TASE, metadata[:filename])
            location = String.(collect(m.captures))
            return location
        end,
        label=(_data, metadata) -> join(
            String.(collect(match(REGEX_TASE, metadata[:filename]).captures)), "_"),
    )
    MB.register_plot!(
        project,
        :four_terminal_iv;
        label="Four-Terminal I-V",
        setup=(_ws, _items) -> (figure = Figure(); Axis(figure[1, 1]); figure),
        draw=function (_workspace, items, figure)
            axis = only(contents(figure[1, 1]))
            for item in items
                df = item_data(item)
                nrow(df) == 0 && continue
                lines!(axis, df.i, df.v)
            end
            nothing
        end,
    )
    return project
end

@testset "TASE registry project" begin
    fixture1 = joinpath(@__DIR__, "fixtures", "TASE", "TASESNS1c1f_A_2TSNJunction_11_20260224_111623_298K_FourTerminalIV.csv")
    fixture2 = joinpath(@__DIR__, "fixtures", "TASE", "TASESNS1c1f_A_2TSNJunction_31_20260224_111700_298K_FourTerminalIV.csv")
    project = _tase_project()
    items = [
        only(items_for_file(project, fixture1)),
        only(items_for_file(project, fixture2)),
    ]

    interpretation = DataBrowserCore.interpret_source_item(
        project,
        DirectorySource(dirname(fixture1)),
        DataBrowserSources.index_source_file(fixture1),
    )
    @test DataBrowser.collection(only(interpretation.interpreted_items)) ==
        ["TASESNS1c1f", "A", "2TSNJunction", "11"]

    @testset "plot data api" begin
        workspace = DataBrowserCore.Workspace.Workspace(project, DirectorySource(dirname(fixture1)))
        try
            data = item_data.(items)
            @test length(data) == 2
            @test all(nrow(df) == 3 for df in data)
            @test all(names(df) == ["i", "v"] for df in data)

            plot_kind = RegisteredPlot{:four_terminal_iv,Symbol("Four-Terminal I-V")}
            @test plot_kind === RegisteredPlot{:four_terminal_iv,Symbol("Four-Terminal I-V")}
            figure = setup_plot(workspace, plot_kind, items)
            @test plot_data!(workspace, plot_kind, items, figure) === nothing
        finally
            DataBrowserCore.Workspace.close_workspace!(workspace)
        end
    end
end
