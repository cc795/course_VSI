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

## Build and run

Requirements: g++, make, Python 3, NumPy, Matplotlib, and a Times-family font. Movie encoding additionally requires FFmpeg.

The upstream Athena++ revision is pinned to `053bc738d5c85dac717326b4338b34c7db4e715d`. The build script copies the standalone problem generator into the Athena++ checkout without changing `disk.cpp` or the solver.

```bash
git clone https://github.com/cc795/course_VSI.git
cd course_VSI
ATHENA_DIR="$PWD/../athena" bash scripts/build.sh
bash scripts/run_figure2.sh
```

Set `ATHENA_DIR` to the desired Athena++ checkout location. The build script creates the checkout if it does not exist and verifies the pinned revision.

The runner prevents duplicate runs in the same output directory and supports resuming from a restart file. The input targets 500 P0, with full VTK snapshots every 1 P0, restart files every 10 P0, and history records every 0.1 P0. Reported simulation time is in code units; divide by `2π` to express it in P0.

## Plotting and movies

The static Figure 2 plot shows a signed logarithmic map of `vz/cs` and a meridional velocity-vector panel. The coordinate transformation uses spherical volume-averaged cell centers:

```text
vz = vr cos(theta) − vtheta sin(theta)
vR = vr sin(theta) + vtheta cos(theta)
```

The movie uses a single flow panel with a dark background and a cyan–true-red signed logarithmic colormap. Color limits and coordinate ranges remain fixed across frames. Titles display the nominal integer orbit; actual snapshot times are retained in the metadata. Regular snapshots are written on the first timestep that reaches or crosses each scheduled time.

Each snapshot is displayed for **0.25 seconds**. Encoding is H.264 at `1440×960` and 24 fps. Frames follow the contiguous sequence `0,1,2,…`, stopping at the first missing orbit. The movie uses simulation data without temporal interpolation. Rendered frames and clips are cached for incremental updates.

Titles, tick labels, and mathematical expressions use Times-style typography, preferring Times New Roman and falling back to the Times-compatible Nimbus Roman font.

- [Download the uploaded MP4](media/figure2_evolution.mp4)
- [View the style preview](media/figure2_preview.png)
- [Inspect frame times and data provenance](media/movie_metadata.json)

The frame metadata describes the time range included in the uploaded movie. Generated simulation data and build products are excluded from version control.

## Validation

The validation script checks the equilibrium density profile, initial perturbation amplitude, temperature relaxation, finite and positive fluid variables, and agreement between continuous and restarted integrations. VTK output is single precision, while the solver uses double precision; comparisons account for output rounding.
