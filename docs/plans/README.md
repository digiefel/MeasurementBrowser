# Design Plans

Focused design docs for work not yet built. The product and architecture north star is
[../vision.md](../vision.md); the staged plan to get there is [roadmap.md](roadmap.md). Use
[../ARCHITECTURE.md](../ARCHITECTURE.md) for the current-state entry point.

## Reading Order

1. [../vision.md](../vision.md) — product and architecture north star: the DataBrowser family,
   workspace/project/package/command vocabulary, scriptable API, generic plotter, figures and
   workflows, caching, and the long-term arc.
2. [roadmap.md](roadmap.md) — the ordered plan from `MeasurementBrowser` to the DataBrowser family,
   starting with the package split.
3. [package-split.md](package-split.md) — the concrete execution of Stage 1: file-by-file mapping,
   extraction order, and the keep-it-green procedure for cutting the package family.
4. [data-model-generalization.md](data-model-generalization.md) — generalizing the core item/collection
   types toward DataBrowser: two extension APIs, `Any` payloads, the cache trait, and the `collect`
   model (the tree as a derived view via `collection`).
5. [project-persistence.md](project-persistence.md) — project files, immutable config,
   script-started projects, data-source config, and export boundaries.
6. [plotting-api-design.md](plotting-api-design.md) — project-facing plot API (`register_plot!`,
   `setup`/`draw`) and GUI-owned composition.
7. [spatial-browser.md](spatial-browser.md) — spatial browser, device/measurement annotations.

## Plan Boundaries

The plans are deliberately split by ownership boundary:

| Doc | Owns | Does not own |
|---|---|---|
| `../vision.md` | The product/architecture north star and shared vocabulary. | Migration steps, current-state claims, specific signatures. |
| `roadmap.md` | The ordered stages from today to the vision. | File-by-file split mechanics, type-system shape, plot signatures, cache schema, behavior changes. |
| `package-split.md` | The concrete Stage 1 execution: which code lands in which `DataBrowser*` package, extraction order, keep-it-green procedure. | Why/when (the roadmap owns that), the type-system generalization, runtime behavior changes. |
| `data-model-generalization.md` | The core item/collection type system, two extension APIs, payload/cache trait, and the `collect` model (`collection`, canonical vs. view). | Visualizers, figures, workflows, annotation storage. |
| `project-persistence.md` | Project files, immutable config, script-started projects, data-source config, and export boundaries. | Item interfaces, plot API, package split mechanics. |
| `plotting-api-design.md` | The project plot API (`register_plot!`, `setup`/`draw`) and GUI-owned composition. | The item/data model, generic built-in visualizers, figure/workflow persistence. |
| `spatial-browser.md` | Spatial navigation plus device and measurement annotations. | Figure annotations and workflow persistence. |

When two docs overlap, prefer the one that owns the boundary in the table above. Update the other
doc with a short cross-reference instead of duplicating the design.
