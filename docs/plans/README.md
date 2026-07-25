# Design Plans

Focused plans for product work beyond the public contracts in [../api.md](../api.md) and the current
architecture in [../ARCHITECTURE.md](../ARCHITECTURE.md).

## Reading order

1. [../vision.md](../vision.md) — product and architecture north star.
2. [roadmap.md](roadmap.md) — version targets, current progress, and longer-term goals.
3. [typed-pipeline.md](typed-pipeline.md) — type-first read, entries, item, and collection stages,
   with registration as an adapter.
4. [project-persistence.md](project-persistence.md) — project files, immutable configuration,
   source configuration, and export boundaries.
5. [diagnostics.md](diagnostics.md) — lightweight workspace health and explicit diagnostic runs.
6. [plotting-api-design.md](plotting-api-design.md) — figure composition and workflow persistence
   above the public plot callback contract.
7. [spatial-browser.md](spatial-browser.md) — spatial navigation, device layout, tags, and notes.

## Ownership

| Document | Owns |
|---|---|
| `../api.md` | Public registration, source, item, processing, and plot callback contracts. |
| `../data-model.md` | Public identity, metadata, collection, and pipeline semantics. |
| `../vision.md` | Product direction, shared vocabulary, and package family. |
| `roadmap.md` | Version focus, feature targets, and current progress. |
| `typed-pipeline.md` | Type-first pipeline stages and registration-adapter architecture. |
| `project-persistence.md` | Project files, immutable configuration, source configuration, and archives. |
| `diagnostics.md` | Performance-window scope and diagnostic capture. |
| `plotting-api-design.md` | Figure composition and workflow persistence. |
| `spatial-browser.md` | Spatial navigation and item annotations. |
