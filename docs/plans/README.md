# Design Plans

These docs hold focused architecture and remaining work. Each document states whether it describes
current or intended behavior. Use [../ARCHITECTURE.md](../ARCHITECTURE.md) as the general entry point.

## Reading Order

1. [workspace-vision.md](workspace-vision.md) — product and architecture north star. Start here for
   scriptable API, GUI/API parity, visualizers, workflows, figure annotations, self-contained figures, and
   non-blocking REPL use.
2. [source-cache-processing.md](source-cache-processing.md) — current source access, interpreted and
   processed caches, processing queue, interactive priority, invalidation, and task ownership.
3. [data-model-generalization.md](data-model-generalization.md) — generalizing the core item/collection
   types toward DataBrowser: two extension APIs, `Any` payloads, the cache trait, and the `collect`
   model (the tree as a derived view via `collection`).
4. [plotting-api-design.md](plotting-api-design.md) — project-facing plot API (`register_plot!`,
   `setup`/`draw`) and GUI-owned composition.
5. [refactor-gap-report.md](refactor-gap-report.md) — current branch gap report and immediate
   cleanup order. This is the most time-sensitive plan doc.
6. [spatial-browser.md](spatial-browser.md) — spatial browser and annotations.
7. [measurement-parameters-and-stats.md](measurement-parameters-and-stats.md) — metadata bucket
   contract for parsed parameters and computed stats.
8. [streaming-status-contract.md](streaming-status-contract.md) — one streaming status contract
   (`WorkspaceStatus` + event stream) between the cache/scan engine and the UI: live feedback, the
   single status/progress model, isolated live rescans, and removal of HDF5-era status shapes.

## Plan Boundaries

The plans are deliberately split by ownership boundary:

| Doc | Owns | Does not own |
|---|---|---|
| `workspace-vision.md` | The end-state user model and shared vocabulary. | Specific signatures, migration steps, or current-state claims. |
| `source-cache-processing.md` | Source access, interpreted and processed cache levels, processing priority, invalidation, and safe task ownership. | Cache table layout, plotting signatures, or widget presentation. |
| `data-model-generalization.md` | The core item/collection type system, two extension APIs, payload/cache trait, and the `collect` model (`collection`, canonical vs. view). | Visualizers, figures, workflows, annotation storage. |
| `plotting-api-design.md` | The project plot API (`register_plot!`, `setup`/`draw`) and GUI-owned composition. | The item/data model, generic built-in visualizers, figure/workflow persistence. |
| `refactor-gap-report.md` | Near-term cleanup list for the branch at the time it was written. | Long-term product vision. |
| `spatial-browser.md` | Spatial navigation plus device and measurement annotations. | Figure annotations and workflow persistence. |
| `measurement-parameters-and-stats.md` | Public meaning of metadata fields. | UI layout, cache schema, plotting signatures. |
| `streaming-status-contract.md` | The status/progress/error boundary between engine and watchers (`WorkspaceStatus` + event stream), live-feedback rules, and removal of legacy status shapes. | The metadata field meanings, the item/collection type system, plotting signatures. |

## Shared Direction

All plans should preserve the same core model:

```text
project code
  -> identifies measurements
  -> loads or computes measurement data
  -> optionally defines domain-specific visualizers

package code
  -> owns opened-project state, cache, background work, selection, workflow history,
     generic visualizers, plot windows, annotations, figure annotations, and persistence
```

When two docs overlap, prefer the one that owns the boundary in the table above. Update the other
doc with a short cross-reference instead of duplicating the design.

Older project/cache/plotting rewrite sketches were removed because they duplicated the workspace
vision and the gap report. Use git history if a deleted sketch is needed for context.
