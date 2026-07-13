# Figure Composition and Workflows

The public plot callback contracts are defined in [../api.md](../api.md). A registration plot
receives processed data and effective metadata; a type-based plot receives concrete
`AbstractDataItem` values. This plan owns composition above those units.

## Ownership

DataBrowser owns selection, materialization, plot windows, live updates, composition, persistence,
and export. A plot callback owns the axes, legends, labels, and series required to present its data.

One render materializes its selection once and passes one consistent snapshot to setup and drawing.
Plot callbacks treat their inputs as read-only. A failure belongs to that plot state and leaves
processed data available to other views.

## Composition

A registered or type-based plot produces one figure. Figure composition combines those figures or
their plot descriptions without moving browser selection, workspace state, or source access into
project callbacks.

```julia
figure = Figure(workspace)
plot!(figure, items, X(:voltage_v), Y(:current_a); visualizer=LinePlot)
fit = linear_fit!(figure, items, X(:voltage_v), Y(:current_a))
annotate!(figure, Arrow(fit; label="linear region"))
```

The package owns:

- live selections and matching rules;
- editable visualizer state;
- multi-panel layout and linked axes;
- figure annotations;
- data snapshots for self-contained figures;
- export and workflow persistence.

## Workflows

A workflow records actions in public package terms: select data, create a visualizer, map columns or
dimensions, change style, fit, annotate, and export. Replaying it restores every visualizer window
and its independent state. Concrete selections and live matching rules use the same action model.
