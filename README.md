# Svanberg, Cui & Latter (2022)

## Figure 2: Athena++ simulation

An Athena++ implementation of the axisymmetric, two-dimensional hydrodynamic disc described in Section 3 of *Wavelike nature of the vertical shear instability in global protoplanetary discs*, MNRAS 514, 4581–4587. DOI: [10.1093/mnras/stac1598](https://doi.org/10.1093/mnras/stac1598).

## Athena++ code

Use this problem generator with the official [Athena++ code](https://github.com/PrincetonUniversity/athena). Copy [`src/svanberg.cpp`](src/svanberg.cpp) into Athena++'s `src/pgen/` directory, configure Athena++ with `--prob=svanberg --coord=spherical_polar --eos=adiabatic --flux=hllc -omp`, and compile it. Run the resulting executable with [`inputs/athinput.figure2`](inputs/athinput.figure2):

```bash
/path/to/athena/bin/athena -i /path/to/course_VSI/inputs/athinput.figure2
```

[Read the Section 3 implementation guide (PDF, Chinese)](readme.pdf).

## Physical setup

The disc is vertically isothermal, with cylindrical radius $R=r\sin\theta$, height $z=r\cos\theta$, and spherical radius $r=\sqrt{R^2+z^2}$. The prescribed mid-plane density and local isothermal sound speed are

$$
\rho_{\mathrm{mid}}(R)=\rho_0\left(\frac{R}{R_0}\right)^{-q_D},
\qquad
c_s^2(R)=c_0^2\left(\frac{R}{R_0}\right)^{-q_T}.
$$

We use $q_D=2$, $q_T=1$, $h=H/R=0.05$, and $GM=R_0=\rho_0=1$, giving $c_0=0.05$ and $c_s^2(R)=0.0025/R$.

### Equilibrium density — equation (24)

Vertical hydrostatic balance in the stellar potential gives

$$
\rho(R,z)=\rho_0\left(\frac{R}{R_0}\right)^{-q_D}
\exp\left[\frac{GM}{c_s^2(R)}\left(\frac{1}{r}-\frac{1}{R}\right)\right].
\tag{24}
$$

This profile is implemented by `DenProfileCyl()` in [`src/svanberg.cpp`](src/svanberg.cpp#L206), using `PoverR()` for $c_s^2(R)$. The equilibrium pressure is $P=\rho c_s^2(R)$.

### Equilibrium rotation — equation (25)

The angular velocity includes radial pressure support and vertical shear:

$$
\Omega(R,z)=\Omega_K(R)
\left[-(q_D+q_T)h^2+(1-q_T)+q_T\frac{R}{r}\right]^{1/2},
\qquad \Omega_K(R)=\sqrt{\frac{GM}{R^3}}.
\tag{25}
$$

`VelProfileCyl()` in [`src/svanberg.cpp`](src/svanberg.cpp#L228) returns the azimuthal linear velocity $v_\phi=R\Omega$. For the adopted parameters, $v_\phi=\sqrt{GM/R}\sqrt{-3h^2+R/r}$. The unperturbed meridional velocities are $v_r=v_\theta=0$.

### Thermal relaxation — equation (27)

The gas temperature relaxes toward the prescribed local equilibrium temperature $T_0(R)$:

$$
\Delta T=\left[T_0(R)-T\right]
\left[1-\exp\left(-\frac{\Delta t}{\tau}\right)\right],
\qquad \tau=10^{-10}P_0.
\tag{27}
$$

`MeshBlock::UserWorkInLoop()` in [`src/svanberg.cpp`](src/svanberg.cpp#L515) implements this exponential update through the internal energy, whose target is $\rho c_s^2(R)/(\gamma-1)$ with $\gamma=7/5$. Density and momentum remain unchanged during the relaxation step, and the original kinetic energy is retained when reconstructing the total energy. The very short relaxation time maintains the locally isothermal temperature profile.

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
