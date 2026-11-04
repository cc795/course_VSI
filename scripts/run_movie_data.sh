#!/usr/bin/env bash
# Reproduce missing early integer-orbit snapshots without changing the main run.
set -euo pipefail
project=$(cd "$(dirname "$0")/.." && pwd)
run=${1:-"$project/results/movie_early"}
if [[ -e "$run/vsi.hst" ]]; then
  echo 'Output directory already has data; use a new directory to avoid overwriting.' >&2
  exit 1
fi
mkdir -p "$run"
export OMP_NUM_THREADS=4
"$project/bin/athena" -i "$project/inputs/athinput.figure2" -d "$run" \
  time/tlim=69.11503837897544 output1/dt=6.283185307179586
