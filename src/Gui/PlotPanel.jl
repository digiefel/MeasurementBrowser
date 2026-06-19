using GLMakie: Axis, Figure, lines!
import GLMakie.Makie as Makie
import CImGui as ig
import CImGui.CSyntax: @c
using NativeFileDialog: save_file

using ..Projects:
    source_label
using ..ItemIndex: ItemRecord
import ..Workspace
using ..Visualization:
    PlotKind,
    debug_plot,
    plot_data!,
    plot_kind_label,
    plot_kind_name,
    registered_plot_kinds,
    setup_plot
using .MakieImguiIntegration: MakieFigure

const PLOT_HELP_TEXT = "Live follows the browser selection.\nDetach opens an independent plot window.\nExport saves the current figure.\nScroll zooms, right-drag pans, Ctrl-click resets limits."

"""Discard rendered figures so every open plot is redrawn."""
function clear_plot_views!(plots::PlotState)::Nothing
    for view in (plots.main, plots.windows...)
        view.figure = nothing
        view.last_key = nothing
        view.error = ""
        view.export_error = ""
    end
    return nothing
end

"""
Draw one plot and keep a short UI error when drawing fails.

The console receives the complete exception together with the source, plot type, and source items.
"""
function draw_plot_view!(
    state::BrowserState,
    view::PlotViewState,
    plot_key::Tuple,
    records::Vector{ItemRecord},
    plot_kind::Type{<:PlotKind},
)::Nothing
    workspace = state.workspace::Workspace.Workspace
    source = workspace.source
    plots = state.plots
    started_ns = time_ns()
    draw_alloc = 0
    try
        figure = nothing
        draw_alloc = @allocated begin
            figure = if plots.debug
                debug_plot(workspace, Workspace.materialize_items(workspace, records); plot_kind)
            else
                result = setup_plot(workspace, plot_kind, records)
                plot_data!(workspace, plot_kind, records, result)
                result
            end
        end
        figure === nothing && error("Plot renderer returned no figure.")
        _append_perf_sample!(
            state,
            :plot_draw,
            (time_ns() - started_ns) / 1e6,
            draw_alloc,
        )
        view.figure = figure
        view.last_key = plot_key
        view.debug = plots.debug
        view.error = ""
    catch err
        bt = catch_backtrace()
        _append_perf_sample!(
            state,
            :plot_draw,
            (time_ns() - started_ns) / 1e6,
            draw_alloc,
        )
        view.figure = nothing
        view.last_key = plot_key
        view.debug = plots.debug
        summary = first(split(sprint(showerror, err), '\n'; limit=2))
        view.error = "Plot failed: $summary. See the console for full details."
        item_context = join(
            [
                "$(record.item_label) ($(record.kind))\n" *
                "  $(record.source_item_id)"
                for record in records
            ],
            "\n",
        )
        @error(
            "Plot drawing failed\n" *
            "Source: $(source_label(source))\n" *
            "Plot: $(plot_kind_name(plot_kind))\n" *
            "Items: $(length(records))\n$item_context",
            exception=(err, bt),
        )
    end
    return nothing
end

"""Render the controls shared by the main and detached plot windows."""
function render_plot_toolbar!(
    state::BrowserState,
    view::PlotViewState,
    records::Vector{ItemRecord},
    selected_records::Vector{ItemRecord},
    available::Vector{Type{<:PlotKind}},
)::Nothing
    plots = state.plots
    project = (state.workspace::Workspace.Workspace).project
    current = view.plot_kind
    if ig.BeginCombo(
        "##plot_kind_$(view.id)",
        current === nothing ? "Choose plot kind" : plot_kind_label(project, current),
    )
        for candidate in available
            if ig.Selectable(plot_kind_label(project, candidate), current === candidate)
                view.plot_kind = candidate
                view.last_key = nothing
                view.error = ""
                for record in records
                    plots.kind_by_item[record.kind] = candidate
                end
            end
        end
        ig.EndCombo()
    end

    live = view.live
    ig.SameLine()
    if @c ig.Checkbox("Live##live_$(view.id)", &live)
        view.live = live
        view.last_key = nothing
        if !live
            view.item_ids = [record.id for record in selected_records]
        end
    end

    if view === plots.main
        can_detach = current !== nothing && !isempty(records)
        ig.SameLine()
        !can_detach && ig.BeginDisabled()
        if ig.Button("Detach") && can_detach
            plots.next_window_id += 1
            push!(
                plots.windows,
                PlotViewState(
                    id="plot_$(plots.next_window_id)",
                    title=length(records) == 1 ?
                          only(records).item_label :
                          "$(length(records)) items",
                    live=false,
                    item_ids=[record.id for record in records],
                    plot_kind=current,
                ),
            )
        end
        !can_detach && ig.EndDisabled()
    end

    can_export = view.figure !== nothing
    ig.SameLine()
    !can_export && ig.BeginDisabled()
    if ig.Button("Export##export_$(view.id)") && can_export
        name = current === nothing ? "plot" : plot_kind_name(current)
        default_name = length(records) == 1 ?
            "$(only(records).id)-$name.png" :
            "$(length(records))-items-$name.png"
        path = save_file(default_name; filterlist="png,jpg,jpeg,svg,pdf")
        if !isempty(path)
            isempty(splitext(path)[2]) && (path *= ".png")
            try
                Makie.save(path, view.figure)
                view.export_error = ""
            catch err
                bt = catch_backtrace()
                summary = first(split(sprint(showerror, err), '\n'; limit=2))
                view.export_error =
                    "Export failed: $summary. See the console for full details."
                @error("Plot export failed\nOutput: $path", exception=(err, bt))
            end
        end
    end
    !can_export && ig.EndDisabled()

    ig.SameLine()
    if ig.Button("?##plot_help_$(view.id)")
        ig.OpenPopup("plot_help_popup_$(view.id)")
    end
    if ig.BeginItemTooltip()
        ig.TextUnformatted("Help")
        ig.EndTooltip()
    end
    if ig.BeginPopup("plot_help_popup_$(view.id)")
        ig.PushTextWrapPos(ig.GetFontSize() * 38.0)
        ig.TextUnformatted(PLOT_HELP_TEXT)
        ig.PopTextWrapPos()
        ig.EndPopup()
    end
    isempty(view.export_error) || ig.TextWrapped(view.export_error)
    return nothing
end

"""Render a figure or the reason no figure is available."""
function render_plot_body!(
    state::BrowserState,
    view::PlotViewState,
    status::Symbol,
)::Nothing
    plots = state.plots
    if plots.debug
        ig.TextColored((0.2, 0.8, 0.2, 1.0), "Debug Plot Mode")
        ig.SameLine()
        _helpmarker("Debug mode bypasses cache and reads source files directly.")
    end

    if view.figure !== nothing
        makie_id = view === plots.main ?
            "item_plot" :
            "item_plot_$(replace(view.id, '/' => '_'))"
        _time!(state, :makie_fig) do
            MakieFigure(
                makie_id,
                view.figure;
                auto_resize_x=true,
                auto_resize_y=true,
            )
        end
    elseif status == :needs_kind
        ig.TextDisabled("No plot kind selected")
    elseif status == :empty
        ig.TextDisabled("No item selected")
    else
        ig.TextColored((1.0, 0.4, 0.4, 1.0), "Plot generation failed")
        isempty(view.error) || ig.TextWrapped(view.error)
    end
    return nothing
end

"""
Update and render one plot window.

Live views use the current browser selection. Static views resolve their saved item ids against the
current workspace index.
"""
function render_plot_view!(
    state::BrowserState,
    view::PlotViewState,
    selected_records::Vector{ItemRecord},
)::Nothing
    plots = state.plots

    if view.live
        view.item_ids = [record.id for record in selected_records]
    end
    records = view.live ?
        selected_records :
        _items_for_ids(state, view.item_ids)

    workspace = state.workspace::Workspace.Workspace
    source = workspace.source
    project = workspace.project
    kind = if isempty(records)
        nothing
    else
        kind = first(records).kind
        all(record -> record.kind == kind, records) ? kind : nothing
    end
    available = kind === nothing ?
        Type{<:PlotKind}[] :
        registered_plot_kinds(project, kind)
    if view === plots.main && view.live
        if view.kind !== kind
            view.kind = kind
            if kind === nothing
                view.plot_kind = nothing
            else
                remembered = get(plots.kind_by_item, kind, nothing)
                view.plot_kind =
                    remembered !== nothing && remembered in available ? remembered :
                    isempty(available) ? nothing : first(available)
            end
            view.last_key = nothing
        end
    end

    if view.plot_kind === nothing && !isempty(available)
        view.plot_kind = first(available)
        view.last_key = nothing
    elseif view.plot_kind !== nothing && !(view.plot_kind in available)
        view.plot_kind = nothing
        view.figure = nothing
        view.last_key = nothing
    end

    status = isempty(records) ?
        :empty :
        view.plot_kind === nothing ? :needs_kind : :ready
    if status == :ready
        plot_key = (
            source_label(workspace.source),
            view.id,
            plot_kind_name(view.plot_kind),
            [record.id for record in records],
            plots.debug,
        )
        view.last_key == plot_key ||
            draw_plot_view!(
                state,
                view,
                plot_key,
                records,
                view.plot_kind,
            )
    elseif view.live
        view.figure = nothing
        view.last_key = nothing
    end

    render_plot_toolbar!(state, view, records, selected_records, available)
    ig.Separator()
    render_plot_body!(state, view, status)
    return nothing
end

"""Warm GLMakie once in a hidden window before the first visible plot."""
function ensure_plot_runtime_warmed!(state::BrowserState)::Nothing
    plots = state.plots
    plots.runtime_warmed && return nothing
    flags = ig.ImGuiWindowFlags_NoDecoration | ig.ImGuiWindowFlags_NoInputs |
            ig.ImGuiWindowFlags_NoBackground | ig.ImGuiWindowFlags_NoSavedSettings |
            ig.ImGuiWindowFlags_NoNav | ig.ImGuiWindowFlags_NoFocusOnAppearing
    ig.SetNextWindowPos((-10_000.0, -10_000.0), ig.ImGuiCond_Always)
    ig.SetNextWindowSize((8.0, 8.0), ig.ImGuiCond_Always)
    if ig.Begin("###plot_runtime_warmup", C_NULL, flags)
        if plots.warmup_figure === nothing
            figure = Figure(size=(64, 48))
            axis = Axis(figure[1, 1], xlabel="x", ylabel="y", title="warm")
            lines!(axis, [0.0, 1.0], [0.0, 1.0], color=:blue)
            plots.warmup_figure = figure
        end
        MakieFigure(
            "_plot_runtime_warmup",
            plots.warmup_figure;
            auto_resize_x=false,
            auto_resize_y=false,
        )
        plots.runtime_warmed = true
    end
    ig.End()
    return nothing
end

"""Render the main plot window."""
function render_plot_window(state::BrowserState)::Nothing
    _, selected, _ = _project_visible_selection(state)
    main = state.plots.main
    if ig.Begin(main.title)
        render_plot_view!(state, main, selected)
    end
    ig.End()
    return nothing
end

"""Render every detached plot window and remove windows closed by the user."""
function render_additional_plot_windows(state::BrowserState)::Nothing
    plots = state.plots
    _, selected, _ = _project_visible_selection(state)
    kept = PlotViewState[]

    for view in plots.windows
        open = Ref(true)
        window_id = replace(view.id, '/' => '_')
        if ig.Begin("Plot: $(view.title)###plot_window_$window_id", open)
            view.debug == plots.debug || (view.last_key = nothing)
            render_plot_view!(state, view, selected)
        end
        ig.End()
        open[] && push!(kept, view)
    end
    plots.windows = kept
    return nothing
end
