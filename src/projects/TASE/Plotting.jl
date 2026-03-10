using GLMakie
using DataFrames

using DataLoader: read_tlm_4p

function load_plot_for_file(::TASEProject, path::AbstractString, kind::Union{Symbol,Nothing};
                            should_cancel::Union{Nothing,Function}=nothing, kwargs...)
    _check_plot_cancel(should_cancel)
    kind === :four_terminal_iv || return nothing
    df = read_tlm_4p(basename(path), dirname(path))
    title = strip(replace(basename(String(path)), r"\.csv$" => ""))
    return (df=df, title=title)
end

function analyze_plot_for_file(::TASEProject, kind::Union{Symbol,Nothing}, loaded; kwargs...)
    kind === :four_terminal_iv || return nothing
    loaded === nothing && return nothing
    df = loaded.df
    return (
        df=df,
        title=loaded.title,
        current_uA=df.current_source .* 1e6,
        voltage_mV=df.voltage_drop .* 1e3,
    )
end

function draw_plot_for_file(::TASEProject, kind::Union{Symbol,Nothing}, analyzed; kwargs...)
    kind === :four_terminal_iv || return nothing
    analyzed === nothing && return nothing
    nrow(analyzed.df) == 0 && return nothing
    fig = Figure(size=(700, 500))
    ax = Axis(fig[1, 1], xlabel="Current (µA)", ylabel="Voltage Drop (mV)", title=analyzed.title)
    lines!(ax, analyzed.current_uA, analyzed.voltage_mV, linewidth=2)
    scatter!(ax, analyzed.current_uA, analyzed.voltage_mV, markersize=4)
    return fig
end
