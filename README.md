# Svanberg, Cui & Latter (2022)

## Figure 2: Athena++ simulation

An Athena++ implementation of the axisymmetric, two-dimensional hydrodynamic disc described in Section 3 of *Wavelike nature of the vertical shear instability in global protoplanetary discs*, MNRAS 514, 4581–4587. DOI: [10.1093/mnras/stac1598](https://doi.org/10.1093/mnras/stac1598).

## Athena++ code

Use this problem generator with the official [Athena++ code](https://github.com/PrincetonUniversity/athena). Copy [`src/svanberg.cpp`](src/svanberg.cpp) into Athena++'s `src/pgen/` directory, configure Athena++ with `--prob=svanberg --coord=spherical_polar --eos=adiabatic --flux=hllc -omp`, and compile it. Run the resulting executable with [`inputs/athinput.figure2`](inputs/athinput.figure2):

```bash
/path/to/athena/bin/athena -i /path/to/course_VSI/inputs/athinput.figure2
```

[Read the Section 3 implementation guide (PDF, Chinese)](readme.pdf).

## Numerical setup

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

## Movie

[Watch or download the movie](media/figure2_evolution.mp4).
