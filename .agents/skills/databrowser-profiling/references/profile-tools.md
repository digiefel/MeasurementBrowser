# Profile navigation

Each profile produces three representations of the same samples:

- `*.jls`: the most complete saved Julia sample and frame data;
- `*.jlprof`: the file ProfileView opens;
- `*.pb.gz`: the file pprof reads.

Changing the pprof view or filter does not require another profile.

Julia's profiler checks the running program many times during the measurement and saves the chain
of functions running at each check. Use these samples to find where time went. Use
`debug_timings.txt` from `profile_scan!` for the durations of marked DataBrowser calls.

## Save inspectable pprof reports

PProf.jl includes the pprof executable. From the repository root, set the run directory printed by
`profile_live!` and save each query as a text file beside the other artifacts:

```bash
PPROF="$(julia --project=bench -e 'using PProf; PProf.pprof_jll.pprof() do path; print(path); end')"
RUN="/absolute/path/to/results/run"
PROFILE="$RUN/cpu_profile.pb.gz"

"$PPROF" -top -nodecount=25 "$PROFILE" > "$RUN/pprof-top.txt"
```

`pprof-top.txt` is the unfiltered reference. Julia profiles often sample threads while they are
waiting. Common entries include:

- `__psynch_cvwait`, `kevent`, and `mach_msg2_trap`: operating-system waits;
- `ijl_task_get_next` and `multiq_deletemin`: Julia waiting for runnable tasks;
- `jl_apply` and `ijl_apply_generic`: Julia runtime call machinery;
- `[unknown function]`: a frame pprof could not name.

These names are context, not automatically things to optimize. If waiting dominates a profile that
was meant to capture active scanning, the profile probably missed the scan work. If the app was
intentionally idle except for the GUI, a large waiting share is normal.

Save a second report that removes waiting samples and keeps stacks that pass through DataBrowser:

```bash
"$PPROF" -top -cum -nodecount=25 \
  -ignore='__psynch_cvwait|kevent|mach_msg2_trap|ijl_task_get_next|multiq_deletemin' \
  -focus='DataBrowser|Workspace|ProjectCache|WorkGraph' \
  "$PROFILE" > "$RUN/pprof-databrowser.txt"
```

This filtered file has a different total because `-ignore` removes waiting samples. Keep
`pprof-top.txt` beside it so the user can see what was removed.

In the table:

- **flat** is samples spent directly in a function;
- **cum** includes that function and everything it called.

Look for the first wide function whose name or source belongs to DataBrowser or the copied RuO2
project. Then ask pprof to show its callers and callees, and save the result:

```bash
"$PPROF" -peek='FUNCTION_NAME' -nodecount=25 "$PROFILE" > "$RUN/pprof-peek.txt"
"$PPROF" -list='FUNCTION_NAME' "$PROFILE" > "$RUN/pprof-source.txt"
```

`pprof-peek.txt` shows how execution reaches the function and where it goes next.
`pprof-source.txt` assigns samples to source lines when source information is available.

Give the user direct links to these text files. Explain the finding using the surrounding
DataBrowser operation—for example source discovery, measurement reading, item processing, cache
loading, tree construction, or GUI drawing. Do not present runtime frames such as `jl_apply` as the
application-level finding.

## pprof web view

Start or refresh the server from the Julia session:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.refresh_pprof!()'
```

Open `http://localhost:57599` and use:

- **Flame Graph** to see call paths and relative widths;
- **Top** to sort functions by direct or cumulative samples;
- **Refine → Focus** with `DataBrowser` to retain stacks containing DataBrowser code;
- **Refine → Ignore** with the wait-function expression above to remove idle samples;
- **Source** after narrowing to one recognizable function.

In a flame graph, start at a recognizable DataBrowser frame and follow boxes toward the work it
calls. A wide operating-system wait box means the sampled thread was asleep, not that the wait
function is consuming CPU.

## ProfileView

ProfileView reads the same capture but preserves Julia-specific navigation:

```bash
jmux --project /Users/davide/code/Julia/DataBrowser/bench \
  'DataBrowserInteractiveProfile.refresh_profileview!()'
```

The user can click frames, open their source, and inspect type inference. The agent should use the
saved pprof text reports and leave graphical inspection to the user.

Sources:

- [ProfileView README](https://github.com/timholy/ProfileView.jl/blob/master/README.md)
- [PProf.jl README](https://github.com/JuliaPerf/PProf.jl)
- [Google pprof documentation](https://github.com/google/pprof/blob/main/doc/README.md)
- [Julia profiling manual](https://docs.julialang.org/en/v1/manual/profile/)
