#!/bin/sh
# Peak RSS of a Julia bench, in kilobytes. `ps -o rss=` is kB on macOS and Linux.
#   bench/run.sh                         # combined harness (writes status.txt)
#   bench/run.sh run.jl
#   bench/run.sh realistic_browse.jl [scale]
#   bench/run.sh scaling.jl [n1,n2,...]
# After a successful run that wrote status.txt, appends peak_rss_kb.
set -eu
BENCH=$(CDPATH= cd -- "$(dirname "$0")" && pwd)

if [ $# -eq 0 ]; then
    script=run.jl
else
    script=$1
    shift
fi

julia --project="$BENCH" --threads=auto "$BENCH/$script" "$@" &
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
fi
exit "$status"
