# Future Plans

These docs describe intended architecture and work order. They are not current-state reference
docs; use [../ARCHITECTURE.md](../ARCHITECTURE.md) when you need to know how the app works today.

## Reading Order

1. [workspace-vision.md](workspace-vision.md) — product and architecture north star. Start here for
   scriptable API, GUI/API parity, visualizers, workflows, figure annotations, self-contained figures, and
   non-blocking REPL use.
2. [refactor-gap-report.md](refactor-gap-report.md) — current branch gap report and immediate
   cleanup order. This is the most time-sensitive plan doc.
3. [spatial-browser.md](spatial-browser.md) — spatial browser and annotations.
4. [measurement-parameters-and-stats.md](measurement-parameters-and-stats.md) — metadata bucket
   contract for parsed parameters and computed stats.

## Plan Boundaries

The plans are deliberately split by ownership boundary:

| Doc | Owns | Does not own |
|---|---|---|
| `workspace-vision.md` | The end-state user model and shared vocabulary. | Specific signatures, migration steps, or current-state claims. |
| `refactor-gap-report.md` | Near-term cleanup list for the branch at the time it was written. | Long-term product vision. |
| `spatial-browser.md` | Spatial navigation plus device and measurement annotations. | Figure annotations and workflow persistence. |
| `measurement-parameters-and-stats.md` | Public meaning of metadata fields. | UI layout, cache schema, plotting signatures. |

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
