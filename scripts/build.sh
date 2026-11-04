#!/usr/bin/env bash
set -euo pipefail
project=$(cd "$(dirname "$0")/.." && pwd)
athena=${ATHENA_DIR:-/workspace/athena}
if [[ ! -d "$athena/.git" ]]; then
  git clone https://github.com/PrincetonUniversity/athena.git "$athena"
  git -C "$athena" checkout 053bc738d5c85dac717326b4338b34c7db4e715d
fi
if [[ $(git -C "$athena" rev-parse HEAD) != 053bc738d5c85dac717326b4338b34c7db4e715d ]]; then
  echo 'Use Athena++ commit 053bc738d5c85dac717326b4338b34c7db4e715d (separate checkout if needed).' >&2
  exit 1
fi
cp "$project/src/svanberg.cpp" "$athena/src/pgen/svanberg.cpp"
cd "$athena"
python3 configure.py --prob=svanberg --coord=spherical_polar --eos=adiabatic --flux=hllc -omp
make -j4
mkdir -p "$project/bin"
cp bin/athena "$project/bin/athena"
git rev-parse HEAD > "$project/bin/athena-commit.txt"
