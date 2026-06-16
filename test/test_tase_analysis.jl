using MeasurementBrowser
using Test
using DataFrames: DataFrame, nrow
using GLMakie: Figure, Axis, lines!, contents

const MB = MeasurementBrowser

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
    MB.register_measurement!(
        project,
        :four_terminal_iv;
        detect=file -> match(REGEX_TASE, file.filename) !== nothing,
        read=file -> _read_tase(file.filepath),
        measurements=function (file, _data)
            m = match(REGEX_TASE, file.filename)
            location = String.(collect(m.captures))
            return [MB.ItemRecord(
                filepath=file.filepath,
                kind=:four_terminal_iv,
                collection=location,
                timestamp=file.timestamp,
                clean_title=join(location, "_"),
            )]
        end,
    )
    MB.register_plot!(
        project,
        :four_terminal_iv;
        label="Four-Terminal I-V",
        setup=(_ws, _items) -> (figure = Figure(); Axis(figure[1, 1]); figure),
        draw=function (_workspace, items, figure)
            axis = only(contents(figure[1, 1]))
            for item in items
                df = item.data
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
    measurements = [
        only(items_for_file(project, fixture1)),
        only(items_for_file(project, fixture2)),
    ]

    @test measurements[1].collection == ["TASESNS1c1f", "A", "2TSNJunction", "11"]

    @testset "plot data api" begin
        workspace = MB.Workspace.Workspace(project, dirname(fixture1))
        data = read_item_data(workspace, measurements)
        @test length(data) == 2
        @test all(nrow(df) == 3 for df in data)
        @test all(names(df) == ["i", "v"] for df in data)

        plot_kind = MB.RegistryPlot{:four_terminal_iv,Symbol("Four-Terminal I-V")}
        @test plot_kind === MB.RegistryPlot{:four_terminal_iv,Symbol("Four-Terminal I-V")}
        figure = setup_plot(workspace, plot_kind, measurements)
        @test plot_data!(workspace, plot_kind, measurements, figure) === nothing
    end
end
