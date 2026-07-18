# Profile navigation

## Shared artifacts

Each CPU or wall-time capture writes three views of the same `Profile.retrieve()` data and line
dictionary:

- `*.jls`: portable raw sample/frame records retained because native JLPROF discards some
  system-specific frame-category information;
- `*.jlprof`: ProfileView-native saved profile;
- `*.pb.gz`: Google's pprof format and the agent's direct CLI input.

Do not collect a new profile merely to change a filter or visualization. `cpu_profile.txt` is not
emitted because a flat dump loses the useful call-path structure.

Sources:

- [ProfileView README](https://github.com/timholy/ProfileView.jl/blob/master/README.md)
- [FlameGraphs documentation](https://timholy.github.io/FlameGraphs.jl/stable/)
- [PProf.jl README](https://github.com/JuliaPerf/PProf.jl)
- [Google pprof documentation](https://github.com/google/pprof/blob/main/doc/README.md)
- [Julia profiling manual](https://docs.julialang.org/en/v1/manual/profile/)

## Agent CLI queries

PProf.jl bundles Google's pprof executable; it does not need to be installed on `PATH`. Resolve it
from the benchmark environment once, then query the exact `.pb.gz` printed by `profile_live!`:

```bash
PPROF="$(julia --project=bench -e 'using PProf; PProf.pprof_jll.pprof() do path; print(path); end')"
PROFILE="/absolute/path/to/cpu_profile.pb.gz"

"$PPROF" -top -nodecount=25 "$PROFILE"
```

The first table must always be that unfiltered raw top. Useful next queries are:

```bash
# Cumulative project-focused top: find the widest project-owned branch.
"$PPROF" -top -cum -nodecount=25 \
  -focus='DataBrowser|Workspace|ProjectCache|WorkGraph' "$PROFILE"

# Callers and callees around a named hotspot.
"$PPROF" -peek='interpret_source_item|read' -nodecount=25 "$PROFILE"

# Per-source-line samples for one function.
"$PPROF" -list='interpret_source_item' "$PROFILE"
```

In a top table, **flat** is samples attributed to that frame itself; **cum** is samples whose stack
passes through that frame, including descendants. Sort with `-cum` when choosing the widest branch.
Report the widest project frame with its cumulative count, then the widest shown child with the
child's own count. Do not attach the ancestor's count to an entire call chain.

Filtering semantics:

- `-focus=REGEX` keeps samples whose stack passes through a matching frame;
- `-hide=REGEX` removes matching frames from the displayed path while preserving their descendants
  and sample total;
- `-ignore=REGEX` removes paths/samples passing through a match and therefore changes the total;
- `-peek=REGEX` shows callers and callees around matching functions;
- `-list=REGEX` annotates matching source lines;
- `-nodecount=N` bounds the displayed table or graph.

Use `hide`, not `ignore`, to get render/runtime wrappers out of the way while retaining DataBrowser
work beneath them. For example:

```bash
"$PPROF" -top -cum -nodecount=25 \
  -hide='renderloop|#render#|try_yieldto|poptask|trypoptask|multiq_deletemin' \
  "$PROFILE"
```

Print every filter expression with its resulting table. Return to the unfiltered total before
accepting a conclusion.

## Human graphical views

The human—not the agent—uses ProfileView for Julia-aware flamegraph navigation, source opening,
clicked-frame `code_warntype`, and task/thread tabs:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.refresh_profileview!()'
```

The human can also start or refresh pprof's web UI from the latest `.pb.gz`:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.refresh_pprof!()'
```

Open `http://localhost:57599`. **Flame Graph** gives the broad call-path view; **Top** shows flat and
cumulative rankings; **Graph**, **Peek**, and **Source** narrow a named branch. The Config menu can
apply focus/hide filters. Refreshing after the next capture keeps the URL and switches it to the new
`.pb.gz`.
