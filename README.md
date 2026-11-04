# Svanberg, Cui & Latter (2022), Figure 2: Athena++ simulation

An Athena++ implementation of the axisymmetric, two-dimensional hydrodynamic disc described in Section 3 of *Wavelike nature of the vertical shear instability in global protoplanetary discs*, MNRAS 514, 4581–4587. DOI: [10.1093/mnras/stac1598](https://doi.org/10.1093/mnras/stac1598).

The target is the Figure 2 snapshot at **t = 500 P0**. The paper's later wave analysis extends to 3000 P0.

## Physical and numerical setup

- Spherical polar coordinates `(r, theta, phi)`, axisymmetry (`nx3=1`), no magnetic field, and no explicit viscosity.
- `GM = R0 = rho0 = 1`, `P0 = 2π`, `gamma = 7/5`, and `h = 0.05`.
- `cs² = 0.0025/R` and `rho_mid = R^-2`. The vertically stratified equilibrium density and rotation follow equations (24) and (25). Gravity uses Athena++'s built-in point-mass source term.
- `r ∈ [1,30]` and `theta ∈ [π/2−0.25, π/2+0.25]`, with a `704×448` grid and adjacent radial cell-width ratio `x1rat=1.005`.
- Ghost cells at all four active boundaries are reset to the equilibrium solution. No additional damping zones are applied.
- Each velocity component receives independent uniform perturbations in `[-0.01 cs, 0.01 cs]`, with seed `202211`. Coordinate-based hashing makes the perturbations independent of the block decomposition.
- After each complete timestep, internal energy is relaxed using the exact exponential update in equation (27), with `tau = 10^-10 P0`. Kinetic energy is preserved, and both conserved energy and primitive pressure are updated.
- The solver uses the adiabatic energy equation, VL2 integration, second-order PLM reconstruction, HLLC fluxes, `CFL=0.3`, and double precision. The supplied input uses four OpenMP threads.

The paper does not specify every numerical option or the random seed. These are documented implementation choices, and pixel-for-pixel agreement is not guaranteed.

The paper's resolution estimates are internally inconsistent: 448 cells across a `10H` angular domain imply approximately 44.8 cells/H, rather than the stated approximately 12 cells/H. This implementation follows the explicit grid dimensions and domain. A strictly logarithmic grid over `[1,30]` with 704 cells would have width ratio `30^(1/704) ≈ 1.004843`; here, the stated ratio `1.005` is used, with Athena++ preserving the endpoints and geometrically increasing cell widths.

## Repository contents

| File | Purpose |
| --- | --- |
| [`src/svanberg.cpp`](src/svanberg.cpp) | Equilibrium disc, velocity perturbations, boundary conditions, and thermal relaxation. |
| [`inputs/athinput.figure2`](inputs/athinput.figure2) | Physical parameters, mesh, integration settings, and output intervals. |
| [`scripts/build.sh`](scripts/build.sh) | Build the problem generator with the pinned Athena++ revision. |
| [`scripts/run_figure2.sh`](scripts/run_figure2.sh) | Run or resume the simulation and generate the target static plot. |
| [`scripts/plot_figure2.py`](scripts/plot_figure2.py) | Generate the paper-style two-panel plot. |
| [`scripts/make_movie.py`](scripts/make_movie.py) | Assemble consecutive snapshots into an MP4. |
| [`scripts/watch_movie.py`](scripts/watch_movie.py) | Monitor new snapshots and refresh the movie. |
| [`scripts/run_movie_data.sh`](scripts/run_movie_data.sh) | Generate supplemental early-time snapshots. |
| [`scripts/validate.py`](scripts/validate.py) | Check equilibrium initialization, cooling, and restart consistency. |
| [`scripts/plot_style.py`](scripts/plot_style.py) | Configure Times-style typography. |

## Run locally on macOS

Run the following commands in macOS Terminal, in order, using the same terminal session. Homebrew must be installed. These commands use Homebrew GCC for OpenMP support and launch Athena++ directly, avoiding the Linux-specific `flock` dependency in the shell runner.

### 1. Install dependencies

```bash
brew install gcc python ffmpeg
```

### 2. Download the repositories and prepare Python

```bash
mkdir -p ~/vsi
cd ~/vsi

git clone https://github.com/cc795/course_VSI.git
git clone https://github.com/PrincetonUniversity/athena.git

python3 -m venv vsi-venv
source vsi-venv/bin/activate
python -m pip install numpy matplotlib

export ATHENA_DIR="$HOME/vsi/athena"
export PROJECT_DIR="$HOME/vsi/course_VSI"
```

If either repository already exists, reuse its checkout instead of cloning it again.

### 3. Compile Athena++

```bash
cd "$ATHENA_DIR"

git checkout 053bc738d5c85dac717326b4338b34c7db4e715d
cp "$PROJECT_DIR/src/svanberg.cpp" src/pgen/

GXX="$(find "$(brew --prefix gcc)/bin" -name 'g++-[0-9]*' | head -n 1)"

python configure.py \
  --prob=svanberg \
  --coord=spherical_polar \
  --eos=adiabatic \
  --flux=hllc \
  --cxx=g++ \
  --ccmd="$GXX" \
  -omp

make -j4

mkdir -p "$PROJECT_DIR/bin"
cp bin/athena "$PROJECT_DIR/bin/"
```

### 4. Run a 200-step smoke test

```bash
cd "$PROJECT_DIR"
mkdir -p results/mac_test

OMP_NUM_THREADS=4 bin/athena \
  -i inputs/athinput.figure2 \
  -d results/mac_test \
  time/nlim=200
```

A successful run ends with `Terminating on cycle limit`.

### 5. Run the full simulation in the background

Start this command once, using an output directory that is not being used by another simulation:

```bash
mkdir -p "$PROJECT_DIR/results/figure2"

OMP_NUM_THREADS=4 nohup caffeinate -i \
  "$PROJECT_DIR/bin/athena" \
  -i "$PROJECT_DIR/inputs/athinput.figure2" \
  -d "$PROJECT_DIR/results/figure2" \
  > "$PROJECT_DIR/results/figure2/run.log" 2>&1 &
```

This runs to **500 P0** with four OpenMP threads. Full flow snapshots are saved every 1 P0. The `caffeinate -i` command prevents idle sleep while the simulation runs; keep a laptop connected to power with its lid open.

### 6. Monitor progress

```bash
tail -f ~/vsi/course_VSI/results/figure2/run.log
```

Press `Control+C` to stop viewing the log; the background simulation continues. Divide the logged `time` by `2π` to obtain the elapsed simulation time in P0.
