# I1.2 — Light-off in a monolithic catalytic converter

An ignition front sweeping down a monolith channel, and the wrong-way
temperature rise that appears when the same converter is suddenly cooled.

- **Structures:** `S4` (transient 1-D convection–reaction), `S7` (two-phase
  gas/solid with interfacial transfer)
- **Reference:** Oh, S. H.; Cavendish, J. C. (1982), *Ind. Eng. Chem. Prod. Res.
  Dev.* 21(1) 29–37, doi:10.1021/i300005a006
- **Runtime:** ~2 min 15 s

## Status

`in-progress` — staged in `queue_cases/I1.2/page/`, not yet promoted to
`pages/`. See `queue_cases/I1.2.yaml`.

## Agreement

**All six numbers of Table III, with nothing fitted.** 0.63 % mean absolute
deviation in the peak wall temperature (990.4 / 986.1 / 966.9 K against
996 / 991 / 975), and at worst 1.25 cells of the paper's own 1/80 grid in its
location. Both trends reproduced: raising the solid conductivity lowers the peak
and moves it downstream. The transient rise is 236 °C against their 241, and the
peak is reached at 9.75 s against 10.1 s.

Free checks that cost nothing and test different things: Table II reproduces
from the cell density and the stated wall thickness alone to 0.09 %; the Lewis
number of H₂ comes out 3.89 against the ~4 the paper states, which tests the
diffusivities, the gas conductivity correlation and the assumed molecular weight
together. The species balance closes to 3.7e-8 of the species fed and the energy
balance to 1.3e-11 of the heat released, at every one of 3001 time levels.

Light-off (a prediction, not a comparison — see below): 47.1 s to 50 % CO
conversion, moving by 0.45 s (1.0 %) over a fourfold change in grid and a
fourfold change in time step.

**The Table III residual is a bias, not scatter, and it is not the grid.** All
three temperatures are low by roughly the same fraction; the value is converged
to 1.5 K between 80 and 320 cells while the gap to the paper is 5 K. Something
small and systematic separates this reimplementation from theirs, most plausibly
one of the four unprinted quantities below. The page does not claim to know
which.

## Provenance: tier 6, not experimental

The paper contains **no experimental data of any kind**. Its two experimental
references — Hegedus (1975) for the steady-state overshoot above adiabatic, and
Mondt (1981) for the absence of melting during a deceleration — are qualitative
appeals to other people's work with no numbers attached. Everything checked here
is the authors' own computed result, from a Galerkin discretisation on 81 grid
points integrated with GEAR. Do not describe this page as validated against
experiment.

The asymmetry between the two cases is the thing an integrator should know:
**case 2 (step decrease) has six published numbers; case 1 (light-off) has
none.** The light-off time on the page is a prediction, and the only claims made
about it are that it is converged and that the qualitative sequence the paper
describes in words comes out.

## Four things the paper needs and does not print

Each is taken from a source the paper itself names, and each is checked on the
page rather than assumed:

1. **`S`**, geometric area per unit volume — follows from the channel geometry,
   `S = 2*eps/R_h`. Checked by reproducing all six numbers of Table II from the
   cell density and the stated 0.0254 cm wall thickness.
2. **`D_i`** — the Slattery–Bird formula the paper cites (Bird et al. 1960,
   p. 505). Checked against the paper's own statement that H₂ has a Lewis number
   near 4.
3. **`Nu_inf`, `Sh_inf`** — the paper says only "constant wall heat flux, Shah
   and London (1971)", which for a square duct is three different numbers
   (3.608 / 3.091 / 2.976). Table III is used as the arbiter; the page prints
   what each alternative gives.
4. **`(-dH)_i`** — standard formation enthalpies, water as vapour.

## Two readings of the printed equations, and how they were settled

Equation 12 is printed as `d(C_ps*T_s)/dt`, not `C_ps*dT_s/dt`, and the
Nomenclature confirms it with
`Psi = (1-eps)*rho_s*[C_ps + T_s*dC_ps/dT_s]`. With the tabulated `C_ps` the two
readings differ by about a factor of two in effective thermal inertia at 300 K.
The page runs both and lets Table III choose — the time of the peak discriminates
much more sharply than the peak temperature.

`R_h` is defined in the paper as `2*(area/perimeter)`, i.e. **half** the usual
hydraulic diameter, which is why eqs 17–18 divide by `2*R_h`. Table II confirms
the reading.

## Four pymrm traps recorded here

1. **The Newton system must be scaled on both sides.** Twelve fields spanning
   750 K to 5e-5 mole fraction, residuals spanning eight orders of magnitude:
   without the constant diagonal row/column scaling in `_build_scaling`, the
   first sparse solve returns "matrix is exactly singular".
2. **The sign of the old-time term.** The accumulation contribution of the
   previous step enters the residual as `+ enthalpy(T_old)/dt` and is
   *subtracted*; getting that sign backwards produces a Newton iteration that
   converges happily onto the temperature lower bound at every step and never
   reports an error.
3. **A conservation check needs the outlet FACE value, not the last cell.** With
   a zero-gradient outlet the face value is extrapolated from the interior and
   differs from the last cell centre by 3 parts in 10⁴. Written on cell values
   the species balance closes to 1 %, which is small enough to look like solver
   noise and large enough to hide a real error; written on
   `compute_boundary_values(..., bc=self.bc_conv[1], bound_id=1)` it closes to
   1e-12.
4. **Nothing that touches a DataFrame belongs inside the residual.** Reading
   parameters with a pandas mask made the notebook take 33 minutes instead of
   2 minutes — twelve times the cost of the linear algebra it wrapped. Resolve
   the table to plain floats before the time loop.

## Rebuilding the page

```bash
python build_page.py                    # regenerate index.ipynb from the builder
python - <<'EOF'
import nbformat; from nbclient import NotebookClient
nb=nbformat.read("index.ipynb",as_version=4)
NotebookClient(nb,timeout=1800,kernel_name="python3",resources={"metadata":{"path":"."}}).execute()
nbformat.write(nb,"index.ipynb")
EOF
```

No PDF is needed to rebuild: all three datasets are committed CSVs.
