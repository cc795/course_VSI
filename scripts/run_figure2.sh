#!/usr/bin/env bash
# Foreground runner; can be launched with nohup for a persistent calculation.
set -euo pipefail
project=$(cd "$(dirname "$0")/.." && pwd)
run=${1:-"$project/results/figure2"}
mkdir -p "$run"
run=$(cd "$run" && pwd)
exec 9>"$run/.runner.lock"
flock -n 9 || { echo 'A runner already owns this output directory.' >&2; exit 1; }
export OMP_NUM_THREADS=4
export MPLCONFIGDIR=/tmp/vsi-mpl
if [[ -f "$run/vsi.final.rst" ]]; then
  restart="$run/vsi.final.rst"
else
  restart=$(find "$run" -maxdepth 1 -name 'vsi.[0-9][0-9][0-9][0-9][0-9].rst' -print | sort | tail -n1)
fi
if [[ -n "$restart" ]]; then
  output_anchor=$(python3 "$project/scripts/restart_next_output.py" "$restart")
  "$project/bin/athena" -r "$restart" -d "$run" time/nlim=-1 time/tlim=3141.592653589793 \
    output1/dt=6.283185307179586 output1/next_time="$output_anchor"
else
  "$project/bin/athena" -i "$project/inputs/athinput.figure2" -d "$run"
fi
python3 "$project/scripts/plot_figure2.py" "$run" --require-target \
  --output "$run/figure2.png" > "$run/diagnostics.json"
printf 'Completed target calculation and plot.\n' > "$run/COMPLETE"
