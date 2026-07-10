import DataBrowserGUI
const Browser = DataBrowserGUI.Browser

using DataBrowserAPI:
    PlotKind,
    plot_kind_from_name,
    plot_kind_name

mutable struct PlotsExtension <: Browser.GuiExtension
    plots::PlotState
    profile_next_plot::Bool
    table_plot::TablePlotState
end

PlotsExtension() = PlotsExtension(PlotState(), false, TablePlotState())

function _persisted_plot_view_from_dict(data::AbstractDict)::PersistedPlotView
    items = get(data, "items", Any[])
    return PersistedPlotView(
        id=String(get(data, "id", "")),
        title=String(get(data, "title", "")),
        plot_kind=String(get(data, "plot_kind", "")),
        live=Bool(get(data, "live", false)),
        items=items isa AbstractVector ? String.(items) : String[],
    )
end

function _plot_view_from_persisted(plot_view::PersistedPlotView)::PlotViewState
    return PlotViewState(
        id=plot_view.id,
        title=isempty(plot_view.title) ? "Plot" : plot_view.title,
        live=plot_view.live,
        item_ids=copy(plot_view.items),
        plot_kind=isempty(plot_view.plot_kind) ?
            nothing :
            plot_kind_from_name(plot_view.plot_kind),
    )
end

function Browser.init!(::PlotsExtension, ::Browser.BrowserState)
    return nothing
end

function Browser.menu!(ext::PlotsExtension, state::Browser.BrowserState)
    render_table_plot_menu!(state, ext.table_plot)
    return nothing
end

function Browser.draw!(ext::PlotsExtension, state::Browser.BrowserState)
    render_plot_window(state)
    render_additional_plot_windows(state)
    render_table_plot_window!(state, ext.table_plot)
    return nothing
end

function Browser.reset!(ext::PlotsExtension, ::Browser.BrowserState)
    ext.plots = PlotState()
    ext.profile_next_plot = false
    ext.table_plot = TablePlotState()
    return nothing
end

function Browser.shutdown!(ext::PlotsExtension, ::Browser.BrowserState)
    clear_plot_views!(ext.plots)
    ext.table_plot.figure = nothing
    return nothing
end

function Browser.is_ready(ext::PlotsExtension, state::Browser.BrowserState)
    ensure_plot_runtime_warmed!(state)
    return ext.plots.runtime_warmed
end

function Browser.save_view(ext::PlotsExtension, ::Browser.BrowserState)
    plots = ext.plots
    return Dict{String,Any}(
        "plot_kinds" => Dict(
            String(kind) => plot_kind_name(plot_kind)
            for (kind, plot_kind) in plots.kind_by_item
        ),
        "main_plot" => Dict{String,Any}(
            "id" => plots.main.id,
            "title" => plots.main.title,
            "plot_kind" => plots.main.plot_kind === nothing ?
                "" : plot_kind_name(plots.main.plot_kind),
            "live" => plots.main.live,
            "items" => copy(plots.main.item_ids),
        ),
        "plot_windows" => [
            Dict{String,Any}(
                "id" => view.id,
                "title" => view.title,
                "plot_kind" => view.plot_kind === nothing ? "" : plot_kind_name(view.plot_kind),
                "live" => view.live,
                "items" => copy(view.item_ids),
            )
            for view in plots.windows
        ],
    )
end

function Browser.load_view!(ext::PlotsExtension, ::Browser.BrowserState, view::Dict{String,Any})
    plots = ext.plots
    empty!(plots.kind_by_item)
    plot_kinds = get(view, "plot_kinds", Dict{String,Any}())
    if plot_kinds isa AbstractDict
        for (kind, plot_kind_name_value) in plot_kinds
            isempty(plot_kind_name_value) && continue
            plot_kind = plot_kind_from_name(String(plot_kind_name_value))
            plot_kind === nothing && continue
            plots.kind_by_item[Symbol(String(kind))] = plot_kind
        end
    end

    main_data = get(view, "main_plot", Dict{String,Any}())
    plots.main = _plot_view_from_persisted(
        main_data isa AbstractDict ?
            _persisted_plot_view_from_dict(main_data) :
            PersistedPlotView(id="main", title="Plot Area", live=true),
    )

    window_data = get(view, "plot_windows", Any[])
    plots.windows = PlotViewState[]
    if window_data isa AbstractVector
        for entry in window_data
            entry isa AbstractDict || continue
            push!(plots.windows, _plot_view_from_persisted(_persisted_plot_view_from_dict(entry)))
        end
    end

    counters = [
        parse(Int, only(match_result.captures))
        for plot_view in plots.windows
        for match_result in (match(r"^plot_(\d+)$", plot_view.id),)
        if match_result !== nothing
    ]
    plots.next_window_id = isempty(counters) ? 0 : maximum(counters)
    return nothing
end

function Browser.open_detached_plot!(
    ext::PlotsExtension,
    ::Browser.BrowserState;
    item_ids::Vector{String},
    title::String,
    kind::Union{Nothing,Symbol},
)::Bool
    plots = ext.plots
    plots.next_window_id += 1
    push!(
        plots.windows,
        PlotViewState(
            id="plot_$(plots.next_window_id)",
            title=title,
            live=false,
            item_ids=copy(item_ids),
            plot_kind=kind === nothing ? nothing : get(plots.kind_by_item, kind, nothing),
        ),
    )
    return true
end

function plots_extension(state::Browser.BrowserState)::PlotsExtension
    for ext in state.extensions
        ext isa PlotsExtension && return ext
    end
    error("PlotsExtension is not loaded")
end

function request_plot_profile!(state::Browser.BrowserState)::Nothing
    ext = plots_extension(state)
    ext.profile_next_plot = true
    for view in (ext.plots.main, ext.plots.windows...)
        view.last_key = nothing
    end
    return nothing
end

export PlotsExtension, plots_extension, request_plot_profile!
