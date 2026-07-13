# Example Projects

These projects show the public API in increasing levels of specialization.

| Example | Style | Main idea |
|---|---|---|
| [`01_simple_csv`](01_simple_csv/project.jl) | Registration | One reader turns each CSV file into one browsable table. |
| [`02_multi_entry_pipeline`](02_multi_entry_pipeline/project.jl) | Registration | One file expands into several items with identity, metadata, processing, analysis, and a plot. |
| [`03_typed_micrographs`](03_typed_micrographs/project.jl) | Type API | A custom source and concrete `AbstractDataItem` use multiple dispatch. |
| [`04_domain_project`](04_domain_project/project.jl) | Packaged registration API | A domain module exposes its own project constructor over multiple registered formats. |

Every example takes a data directory as its first command-line argument:

```bash
julia --project examples/01_simple_csv/project.jl /path/to/data
```

The scripts contain the complete project definitions. Domain projects can move the same definitions
into a Julia package and expose a narrower constructor or launch function to their users.
