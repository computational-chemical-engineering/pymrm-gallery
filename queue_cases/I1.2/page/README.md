# I1.2 — Light-off in a monolithic catalytic converter

An ignition front sweeping down a monolith channel, and the wrong-way
temperature rise that appears when the same converter is suddenly cooled.

- **Structures:** `S4` (transient 1-D convection–reaction), `S7` (two-phase
  gas/solid with interfacial transfer)
- **Reference:** Oh, S. H.; Cavendish, J. C. (1982), *Ind. Eng. Chem. Prod. Res.
  Dev.* 21(1) 29–37, doi:10.1021/i300005a006
- **Runtime:** see `meta.yaml`

## Status

`in-progress` — staged in `queue_cases/I1.2/page/`, not yet promoted to
`pages/`. See `queue_cases/I1.2.yaml`.

## Agreement

PLACEHOLDER — filled after the run.

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

## Two pymrm traps recorded here

1. **The Newton system must be scaled on both sides.** Twelve fields spanning
   750 K to 5e-5 mole fraction, residuals spanning eight orders of magnitude:
   without the constant diagonal row/column scaling in `_build_scaling`, the
   first sparse solve returns "matrix is exactly singular".
2. **The sign of the old-time term.** The accumulation contribution of the
   previous step enters the residual as `+ enthalpy(T_old)/dt` and is
   *subtracted*; getting that sign backwards produces a Newton iteration that
   converges happily onto the temperature lower bound at every step and never
   reports an error.

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
