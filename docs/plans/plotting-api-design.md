# Plotting API Design

## Purpose

Project code describes how one plot is laid out and how processed items are drawn. The package owns
selection, materialization, plot windows, composition, and export; the GUI shell owns browser
lifecycle and hosts the registered plot windows.

```julia
register_plot!(project, :iv;
    label="I-V",
    setup=(workspace, items) -> Figure(),
    draw=(workspace, items, figure) -> draw_iv!(figure, items),
)
```

`setup` returns the Makie `Figure`. `draw` fills that figure and returns `nothing`. Both callbacks
receive processed `AbstractDataItem`s whose data is available through `item_data(item)` or
`item.data` for `DataItem`.

## Ownership

| Package owns | Plot callback owns |
|---|---|
| resolving stable selection ids | figure layout required by that plot |
| loading or producing processed items | axes, legends, labels, and series |
| materializing the selection once | interpreting the selected item's domain data |
| main and detached plot windows | drawing into the supplied figure |
| Live selection behavior | project-specific visual meaning |
| export and error presentation | clear callback failures |

The browser never passes internal `ItemRecord`s to project callbacks. One render materializes its
selection once and passes the same item objects to `setup` and `draw`. This prevents duplicate DuckDB
reads or processing and gives both callbacks one consistent snapshot.

## Dispatch

Plots are registered by `(kind, label)`. The browser restores the most recently selected label for an
item kind and resolves it to the corresponding `RegisteredPlot` type. `Visualization.setup_plot` and
`Visualization.plot_data!` are the internal bridge; project code uses only `register_plot!`.

The main plot follows browser selection while Live is enabled. Detached plots start with Live off and
keep their stable item ids until the user enables it. Both use the same rendering path.

## Data Rules

- Plot callbacks always receive processed data.
- Cache policy changes how the data is obtained, not the plotting contract.
- Callbacks treat item data as read-only.
- A plot failure is attached to that plot state and does not discard valid processed data.
- Returning no figure from `setup` is an error.

There is no debug-plot mode. A future diagnostic tool should be designed after the processing and plot
lifecycle are stable; it should not add a second source-loading or processing path.

## Composition

Today one registered plot produces one figure. Future composition should combine figures or plot
descriptions in `DataBrowserPlots` without making project callbacks manage browser windows, selection,
or workspace state. The current `setup`/`draw` contract remains the small project-facing unit.
