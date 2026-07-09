# Project Persistence

## Purpose

Make a project reopenable without turning it into a bundle format.

A project is one editable text file. That file records the durable configuration for one analysis:
which code defines the project, which data sources it opens, and where generated state such as cache
may live. The Julia script still defines the callbacks and plots. The project file is the stable
thing both the GUI and scripts can load and save.

## Core Model

There are only three moving pieces:

- **project file**: the persistent text file for this analysis;
- **project config**: the immutable value loaded from that file;
- **workspace**: one live opening of the project config plus the registered Julia callbacks.

Data-source configuration is core project state. It is not transient like the current selection.
Changing data-source config is a meaningful edit to the project file.

Selections, and live GUI focus are workspace/browser state. They should be restored as
convenience state, but they should not be part of the project config.
Open windows are a more fundamental part of the analysis (domain-specific analysis),
so they may be part of the config.

## Script Flow

Current use starts from a Julia script. The persistent version should keep that shape:

```julia
# non-definitive, pseudo-api inspiration
config = ProjectConfig(
    sources = [DirectorySource("data")],
)

project = define_project("ruo2.toml"; config)
register_item!(project, :iv; ...)
register_plot!(project, :iv; ...)

ws = open_workspace(project)
```

`define_project(project_file; config=nothing)` is the boundary:

- if the file does not exist and `config === nothing`, our reasonable default is written to file;
- if the file does not exist, `config` may be constructed and explicitly passed by the user;
- if the file exists and `config === nothing`, the config is loaded;
- if the file exists and `config` matches the file, loading succeeds;
- if the file exists, `config` is explicitly passed, and it differs, loading fails.

That makes rerunning a script safe. It reopens the same project file; it does not create a second
project by accident. If the script and file disagree, the error is explicit instead of silently
choosing one. The user may then decide to change project name or delete the old file.

The current `define_project(name)` plus `open_workspace(project, root)` shape remains the in-memory
development path until persistence exists. It maps naturally to the project-file version: the name
becomes a file path, and the `root` argument moves into the saved config.

## What The File Contains

The project file should stay small and hand-editable:

```toml
name = "RuO2"
entry = "browser.jl"

[[sources]]
kind = "directory"
path = "../data/electrical"
metadata_file = "device_info.txt"

[cache]
path = ".databrowser/cache.duckdb"
```

Exact keys are deferred. The important shape is:

- code entry: script or package entry needed to register callbacks;
- data-source configs: enough to construct each source;
- cache location: where generated state for this project is stored.

The file does not contain data, figures, cache payloads, or a copied Julia environment. Those belong
to export/archive workflows, not normal project persistence.

## GUI Edits

The GUI can edit data-source config because it is project config.

Adding folder B edits the project file from "A" to "A + B". Replacing A with B edits the project file
from "A" to "B". The GUI should save that change explicitly, just like editing a configuration file.

Changing selection, moving focus, or opening a transient plot does not make the project file dirty.
Those can be remembered separately as app/workspace convenience state if needed.

API functions may also edit project settings through a live workspace. That is fine: they are commands
against an opened project, not the primary way to define one. The project file and immutable loaded
config remain the definition boundary.

## Cache And Removed Sources

The first implementation should stay close to the current source/cache model.

If a source is removed from the project file, its items are removed from the workspace view. Cache
rows for that source can be deleted or left for later cleanup; either way, the user-facing behavior is
simple: removed sources are no longer browsable.

Keeping detached cached items visible after their source is removed would require a larger model:
items would need an explicit missing-source state, clicks would need clear failure behavior when
processed data is absent, and the UI would need tools to distinguish archived, missing, trashed, and
live items. That is not the default project-persistence design.

A trash/bin behavior can be built later as an item annotation/tag workflow. That is separate from
removing a data source from the project file.

## Export

Export is not normal open/save.

An exported project is a packaged Julia environment around a project file. Export options may choose
to copy data into the package and may choose to include cache, but those choices are source-dependent
and may be impossible for some sources. Export is therefore a later workflow, not part of the core
project-file format.

Normal sharing can be as simple as putting the project file, scripts, environment files, and data in a
repo or directory under the user's own conventions.

## Deferred

- exact project-file extension;
- exact TOML schema;
- how the GUI discovers/runs the `entry` script when opening a project file directly;
- whether cache is deleted immediately or left for cleanup when a source is removed;
- export/archive format and source-specific data-copy support.
