#!/bin/sh
set -eu
BENCH=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
exec julia --project="$BENCH" --threads=auto "$BENCH/../scripts/check.jl" bench "$@"
