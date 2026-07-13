# Workspaces

A project explains data. A workspace applies that explanation to one source and owns the resulting
live state.

```julia
project = define_project("My project")
workspace = open_workspace(project, source)
open_browser(workspace)
```

For a directory-backed project, `source` can be a root path:

```julia
workspace = open_workspace(project, "/path/to/data")
```

## Project, source, and workspace

These values have separate responsibilities:

| Value | Responsibility |
|---|---|
| project | interpretation, processing, analysis, and description rules |
| source | discovery, change reporting, and access to source items |
| workspace | scanning, scheduling, persistent results, queries, selections, and browser state |

Project callbacks do not manage cache files, worker tasks, or GUI state. The same project can be
opened on different sources, and the same source can be reopened with its persistent workspace.

## Opening and building

When a workspace opens, DataBrowser discovers source items and compares their fingerprints with the
persistent index. Unchanged results are reused. New or changed sources enter the pipeline.

Items are published as soon as their interpretation and cheap description are available. Deferred
processing and analysis continue in the background. Selecting an unfinished item raises the
priority of the work needed to use it.

## Live changes

Directory sources watch for added, changed, and removed files. Custom sources can implement the same
update contract for databases, remote stores, or streams.

A changed source invalidates the work that depends on its fingerprint. Other items remain attached
to their cached results, selections, annotations, and views.

Project code can also be rerun during development. Re-registering the same name replaces that
registration in place, allowing a live project definition to evolve predictably with Revise.

## Closing

```julia
close_workspace!(workspace)
```

Closing releases source resources, watchers, background tasks, cache connections, and browser state
owned by the live workspace. Persistent results remain available for the next open.
