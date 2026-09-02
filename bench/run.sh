#!/bin/sh
# Peak RSS of the performance snapshot, in kilobytes. `ps -o rss=` is kB on macOS and Linux.
# After a successful run that wrote status.txt, appends peak_rss_kb and prints the file.
set -eu
BENCH=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
export MB_BENCH_SCALE="${MB_BENCH_SCALE:-0.1}" # default 0.1 if unspecified

julia --project="$BENCH" --threads=auto "$BENCH/run.jl" &
pid=$!
peak=0
while kill -0 "$pid" 2>/dev/null; do
    rss=$(ps -o rss= -p "$pid" | awk '{print $1; exit}')
    if [ -n "$rss" ] && [ "$rss" -gt "$peak" ]; then
        peak=$rss
    fi
    sleep 0.5
done
set +e
wait "$pid"
status=$?
set -e
printf 'peak rss: %s kB\n' "$peak"
if [ -s "$BENCH/status.txt" ]; then
    grep -v '^peak_rss_kb' "$BENCH/status.txt" > "$BENCH/status.txt.tmp"
    printf '%-40s  %-24s  # %s\n' peak_rss_kb "$peak" "High-water mark of the Julia process RSS during the whole run, in kilobytes, from ps -o rss= about twice per second (kB on macOS and Linux)." >> "$BENCH/status.txt.tmp"
    mv "$BENCH/status.txt.tmp" "$BENCH/status.txt"
    cat "$BENCH/status.txt"
fi
exit "$status"
