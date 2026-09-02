#!/bin/sh
# Peak RSS of a Julia bench, in kilobytes. `ps -o rss=` is kB on macOS and Linux.
#   bench/run.sh realistic_browse.jl [scale]
#   bench/run.sh scaling.jl [n1,n2,...]
set -eu
BENCH=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
stamp=$(date +%Y%m%d-%H%M%S)
outdir=${MB_BENCH_OUTDIR:-"$BENCH/results/run-$stamp"}
mkdir -p "$outdir"
export MB_BENCH_OUTDIR="$outdir"

script=$1
shift
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
printf '%s\n' "$peak" > "$outdir/peak_rss_kb.txt"
printf 'peak rss: %s kB\n' "$peak"
exit "$status"
