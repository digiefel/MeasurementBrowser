# Profiling

Profiling belongs to one open workspace. It separates cheap metrics, detailed traces, sampled CPU
profiles, and crash forensics instead of making one mechanism serve incompatible costs.

## Normal Diagnostics

Project callback timings, scan summaries, plot timings, failures, process memory, rendering state,
and aggregate cache-write metrics remain available normally. `BuildMetrics` stores only atomic totals
and call counts used by status, the Performance window, and benchmarks. These diagnostics do not
retain engine events or serialize files.

## Internal Capture

Open a workspace with `profile_internal=true` to reveal the Internal Performance tab and permit
buffered tracing. `profile_cpu=true` adds Julia CPU sampling to each manual capture.
`profile_output="trace.json"` exports Perfetto JSON when capture stops; later captures use numeric
suffixes. The matching environment flags are `MB_PROFILE_INTERNAL`, `MB_PROFILE_CPU`, and
`MB_PROFILE_OUTPUT`. Boolean environment values accept only `0` or `1`, and explicit keywords take
precedence.

A session moves through disabled, idle, preparing, recording, complete, and error states. Starting
from idle records subsequent work immediately. Starting during a build cancels scanning and analysis,
cancels queued processing while currently executing project callbacks finish, waits for complete
idleness, and then starts one clean profiled rebuild. Stopping ends capture without canceling current
application work.

The structured trace records completed spans for source discovery and callbacks, interpretation,
processing and statistics, queueing, cache reads and reconstruction, writer acquisition,
transactions, collection analysis, and plot materialization/setup/drawing. Fixed attributes may
include stage, kind, source/item identity, row/item/byte counts, batch size, cache size, and writer
wait/service time. Measurement values and metadata payload contents are never recorded.

Process-wide RSS, GC totals and pauses, queue progress, and writer counters are sampled every 100 ms
on counter tracks. They are not attributed to overlapping spans. A capture retains at most 500,000
events and reports any dropped events explicitly. When internal profiling is disabled or idle, span
attributes are not evaluated.

CPU sampling uses Julia's standard profiler only during a manual capture. Starting fails if another
Julia profiler is active. Stopping stores a reduced source-line hotspot table and clears the raw
sampling buffers.

## Crash Forensics

`crash_trace="trace.jsonl"` or `MB_CRASH_TRACE` synchronously writes and flushes span start/end records
from workspace construction until close. It is independent of buffered profiling, so a terminated
process leaves an unmatched start for the active operation. Failure to open the requested file aborts
workspace creation. Use it only for a reproducible crash because synchronous flushing is intentionally
expensive.

Start, stop, reset, and export functions remain module-qualified internals used by the GUI,
benchmarks, and tests. They are not part of the project-facing API.
