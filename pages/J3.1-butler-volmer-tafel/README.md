# J3.1 — Butler–Volmer kinetics and its linear and Tafel limits

The interfacial charge-transfer law in the form a reactor model uses it, and the
overpotential at which each of its two standard approximations goes wrong by a
stated amount — checked against the overpotentials the published `J3.4` cell
actually runs at.

- **Structures:** `S1` (time-marched IVP, for the cell used as an operating
  envelope), `S10` (porous electrode with electroneutrality and a distributed
  interfacial reaction)
- **Reference:** Doyle, Fuller & Newman (1993), *J. Electrochem. Soc.* **140**(6)
  1526–1533, doi:10.1149/1.2221597 — Eqs. 6, 8, 16, 17, 18, 28, 29, 30
- **Also:** Marquis, Sulzer, Timms, Please & Chapman (2019), *J. Electrochem.
  Soc.* **166**(15) A3693–A3706, doi:10.1149/2.0341915jes — the same kinetics in
  hyperbolic-sine form (Eqs. 1g–1h), and its printed closed-form inversion
  (Eq. 18)
- **Runtime:** ~75 s

## The source situation — read this first

The Butler–Volmer equation's historical sources — Tafel (1905), J. A. V. Butler
(1924, 1932), Erdey-Grúz & Volmer (1930) — **were not consulted**. They are not
on disk, are not open access, and **no DOI is asserted for them**; they are
recorded as conventional attribution only.

The source that *was* read is Doyle, Fuller & Newman (1993), which prints the law
twice: as the two-exponential kinetics of a flat electrode (Eq. 6 with Eq. 8) and
as a mass-action expression for an insertion reaction (Eq. 17 with Eq. 30). Every
equation used here was read off 600 dpi renders of journal pages 1527, 1528 and
1530, because the scan's text layer mangles them (Eq. 17 comes out as
`i = Fk2(emax - c)~~ ... ~R~01 - U'))`).

The page states this distinction in its Background section. Nothing on it is a
restatement of a source that was not read.

## Provenance: tier 6, not experimental

Nothing here is compared with a measurement. Doyle's Figure 2 contains no
experimental points — confirmed by the maintainer on review of the original — and
every check is an algebraic identity or a number the paper prints. Do **not**
describe this page as validated against experiment.

## Why this is not covered by `J3.4`

`J3.4` uses Eq. 17 inside a full P2D cell and never examines the law: it does not
derive the exchange current density (it calls Eq. 30 only from `kappa_from_nu`,
never in the kinetics), does not ask what happens when the transfer coefficients
fail to sum to one, and does not touch either limit — its `arcsinh` is the exact
inversion of the symmetric law, not an approximation. This page is about the
constitutive law, and the two pages check each other — the law's limits say what
the cell may safely assume, the cell says which part of the law's domain it lives
in.

Two qualifications, both stated on the page. `J3.5`'s `bv_invert` silently
*contains* the reduction derived here (at zero current its root is Eq. 16, and
its prefactor has Eq. 30's structure) without ever saying so or stating the
α_a + α_c = 1 condition it needs. And `J3.5`'s `C_r_doyle` is the same
interfacial/superficial ratio, printed there as ~1e−4 against 1.0e−3 here — the
same quantity on a different *i*₀ scale, differing by exactly √(u₀(1−u₀)) =
0.0995, which the page computes so the two entries cannot be read as
contradicting each other.

## Agreement

| what | result |
|---|---|
| Eq. 17 reduced with Eq. 16 vs the printed Eq. 30 | 4.6e−12 relative, 4000 states × 4 α pairs |
| Eq. 17's zero-current potential vs Eq. 16 | 6.6e−16 V (and fails by O(10 mV) when α_a+α_c ≠ 1) |
| *i* = 0 at zero surface overpotential | 7.8e−15 |
| anodic and cathodic partial currents balance at equilibrium | 1.5e−14 — **structural**, and reported as such |
| that common value = Eq. 30's *i*₀ | 5.3e−15 over four α pairs *including asymmetric ones* |
| *R*_ct = *RT*/(*nFi*₀), *n* = α_a+α_c | 1.8e−9 vs a central difference on Eq. 17 |
| arcsinh inversion fed forward through Eq. 17 | 9.7e−11 over five decades of current |
| closed-form Tafel threshold vs root-finding | 3.1e−15 |
| Doyle's ν²/δ vs Eq. 30 (κ cancels) | −15.7 %, unexplained — *and the same number as `J3.4`'s 19 % κ spread, not a second test* |
| standard form vs Eq. 17 as printed, in the cell | 2.2e−12 mV |
| Tafel-substituted first step vs prediction | 0.46 % |
| grid spread over an 8× refinement | 0.149 mV — **spatial only** |
| temporal error at the production schedule | **0.273 mV**, observed order 0.95 — 1.8× the grid spread |
| total salt conservation | 8.9e−16 — structural; a 0.1 % collector leak takes it to 2.2e−2 |

### The validity thresholds (α_a + α_c = 1)

| deviation | Tafel, 298 K | Tafel, 373 K | linear (α=0.5), 298 K | linear, 373 K |
|---|---|---|---|---|
| 10 % | 61.6 mV | 77.1 mV | 41.3 mV | 51.7 mV |
| 5 % | 78.2 mV | 97.9 mV | 28.7 mV | 35.9 mV |
| 1 % | 118.6 mV | 148.4 mV | 12.6 mV | 15.8 mV |

The Tafel column has a closed form, η = *RT* ln(1 + 1/ε)/((α_a+α_c)*F*), and
depends on the transfer coefficients **only through their sum**. The linear
column does depend on the split.

### The published cell, at *I* = 10 A/m²

The verdict depends on the tolerance, so the page states it, and both are printed.

| | *i*/*i*₀ | η_s | × the 5 % linear threshold | × the 1 % one |
|---|---|---|---|---|
| porous cathode | ~1e−3 | 0.032 mV mean, 0.74 mV peak¹ | 0.001 (mean) | 0.002 (mean) |
| lithium foil | 0.79 | 24.8 mV | 0.69 — inside | 1.57 — outside |
| lithium foil at *I* = 20 | 1.59 | 46.6 mV | 1.30 — outside | 2.95 — outside |

¹ The peak is a *first-step* value and is converged in neither the grid nor the
time step — 1.519 mV with both refined eightfold. The sustained value past the
initial transient is 0.146 mV. See "A fifth" below.

At *I* = 10 both electrodes are inside the **5 %** linear window; they are in
genuinely different regimes only at the 1 % tolerance, or at *I* = 20. What is
tolerance-independent is the **2.9 decades** between them.

Substituting the linear law in the cathode moves the cell potential by 1.5e−10
mV. Substituting Tafel moves it by 519 mV — and of the wrong sign, because Tafel
with *i* ≪ *i*₀ has no solution near η_s = 0.

### What the identity checks can and cannot catch

Measured by injecting five mis-transcriptions one at a time, not asserted:

- at Doyle's own α_a = α_c = 0.5, **every** exponent swap is a no-op — the whole
  discriminating power of the Eq. 17 vs Eq. 30 check comes from the sweep's
  hypothetical asymmetric rows, where each single-equation defect lifts the
  residual from 4.5e−12 to between 0.95 and 1.0;
- swapping the `(c_max−c)^α_c c^α_a` prefactor in **both** equations at once is
  undetectable even there, because it is the same code expression on both sides
  and cancels; only the `(c_T−c_s)^α_a c_s^α_c` half is actually tested;
- the partial-current balance never moves under any defect — it is an exact
  algebraic consequence of α_a + α_c = 1.

## Three defects the checks caught, and inspection did not

1. A bracketing root-find on the linear deviation returned a *later* crossing:
   the α_a = 0.3 threshold came out as 203 mV instead of 8.3 mV. The deviation is
   non-monotonic for asymmetric α and the threshold is the **first** crossing.
2. Comparing two kinetic laws under independently adaptive time stepping put
   0.44 mV of time-discretisation error into a quantity whose true value is
   1e−10 mV. The reference step sequence is now replayed exactly.
3. The kinetic switch switched *both* electrodes, so the "cathode" comparison was
   reading the anode's error — and produced almost the same 0.44 mV. Two
   different defects, one plausible-looking answer.

## A fourth, found by cross-page audit (2026-08-02)

The *I* = 20 A/m² envelope numbers were read at the **march endpoint** instead of
at the 1.7 V cutoff the loaded `J3.4` parameter file carries. `Cell.march`
inherits `v_stop = 1.55` V, 150 mV below Doyle's own cutoff, so its endpoint is a
stepper artefact. The page said "collapses at *u* = 0.40" with a cathode peak of
3.7 mV, "the largest kinetic overpotential anywhere". Read at the cutoff:

| | value |
|---|---|
| *u* at 1.9 V, *I* = 20 | 0.2638 — `J3.4`'s 0.264 |
| *u* at the 1.7 V cutoff | **0.3707** (Doyle states "about 30 %"; the CSV carries 0.30) |
| march endpoint, *V* = 1.473 V | 0.4033 — the old "0.40" |
| cathode peak inside the cutoff | **1.679 mV**, against 3.705 mV over the full run |
| same run, *I* = 10, at the cutoff | 0.8312 — `J3.4`'s 0.831, Doyle's 0.84 |

The model was right; only the readout was. The superlative was wrong
independently — the lithium foil at *I* = 20 runs at 46.6 mV — and has been
removed rather than restated. No agreement metric moved: the *I* = 5 and *I* = 10
readouts use a `u < 0.80` mask stricter than their cutoffs.

## A fifth, from the guard-structure sweep (2026-08-05)

`Cell` is `J3.4`'s, and `Cell.march` sets `dt` and `dt_max` from the current
alone — from `I` and the reference discharge time, and from **neither** `n_s` nor
`n_c`. So the grid ladder refined space at a *constant* temporal error, and the
three `grid_*` metrics this page publishes were being read as if they bounded the
whole discretisation error. `J3.5` had already fixed this in the same helper by
adding a `dt_scale` that multiplies `dt` **and** `dt_max` — the only knob that
refines a schedule growing ×1.3 per step — and `J3.4` took the fix on 2026-08-05.
It had never reached here.

Ported, with `dt_scale = 1` as the production schedule, so **all 47 previously
published metrics are bit-identical**; 22 were added. What the missing half turns
out to be worth:

| quantity | 8× grid | 8× dt | ratio |
|---|---|---|---|
| cell potential spread | 0.1488 mV | 0.2377 mV (0.273 mV extrapolated) | 1.6× (1.8×) |
| cell-averaged \|η_s\| spread | 0.403 % | 0.949 % | 2.4× |
| peak \|η_s\|, coarsest → finest | 1.90× | 1.34× | not converged in either |

Observed order 0.95 — first order, as implicit Euler gives. **On both quantities
that converge, the larger of the two numerical errors was the one nobody had
measured.** `J3.4`, the same cell and the same schedule, reports 0.273 mV against
a 0.138 mV grid spread; the agreement is not a coincidence, it is the shared
helper.

The peak is the exception, and refining time is what shows why: **it sits in the
first time step**, 0.39 s into a discharge of 1.7×10⁴ s. It is a boundary layer
in *t* as much as in *x*, and it grows under both knobs:

| peak \|η_s\| | value | margin to the 1 % threshold |
|---|---|---|
| production grid, production dt | 0.738 mV | 21.4× |
| 8× grid only | 1.008 mV | 15.7× |
| 8× dt only | 0.992 mV | 15.9× |
| **both together** | **1.519 mV** | **10.4×** |

The classification survives — an order of magnitude of margin — but 10.4× is the
honest number and it replaces the 15.7× taken from the grid ladder alone. The
*sustained* peak, past the initial transient, is 0.146 mV and moves only 10.6 %
over the whole refinement box; that is the converged number.

Three further things did not travel with the helper, all latent and all recorded
on the page: `march` here has no `record`/`t_stop` snapshot interface, which is
what `J3.4`'s two checks with real discriminating power are built on (this page
inherited the structural salt identity instead); the charge branch of the
stopping rule, `(I < 0 and v > v_top)` with `v_top = 3.2` V, is absent; and
`kappa`, `eps` and `delta_c` are fixed constants rather than constructor
arguments.

## Every metric has a row, and eleven are invisible to CI

Applying `AGENTS.md`'s break-table rule to the whole page: of **69** published
metrics, **63 carry a row that moves them** and **6 are labelled structural**,
each stating what it cannot detect. Rows were built where they were missing —
the α-offset closed form, the zero-current check, the *R*_ct identity, the
arcsinh inversion, the Tafel threshold's log form, the mirror identity, the
general *i*₀, the linear threshold's first-crossing rule, and salt conservation.

**Eleven of the 69 sit below `check_agreement.py`'s `ABS_FLOOR = 1e-12` and are
therefore not compared by CI at all** — unprotected, not proven. Every break row
added is far above the floor, so the evidence is CI-visible even where the
residual it guards is not. Both counts are published as metrics, so a later edit
that pushes a metric under the floor or relabels one structural fails CI.

## Data

No dataset of its own, and **no figure was digitised**. Two files are loaded
cross-page from `J3.4`: its reviewed parameter set and its stated-results table.
`J3.4` must stay published for this page to execute.
