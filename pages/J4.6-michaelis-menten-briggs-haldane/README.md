# J4.6 — Michaelis–Menten and Briggs–Haldane

**One rate law, two constants, and what data can decide between them.**

Section 1.5.1 of Froment, De Wilde & Bischoff derives the enzyme rate law twice.
Assume the complex decomposition is rate determining and the binding step
equilibrates, and you get eq. (1.5.1-8) with **K_M = k₋₁/k₁** — the line after it
names this "the Michaelis-Menten equation". Assume only a pseudo steady state on
the complex — "the second approach, formulated by Briggs and Haldane" — and you
get the *identical* expression with **K_M = (k₋₁ + k₂)/k₁**, eq. (1.5.1-14). The
two constants differ by the factor 1 + k₂/k₋₁, and nothing in the rate law bounds
it.

This page asks what an experiment can do about that.

## What it finds

**No steady-state rate measurement can separate them, and this is exact.** The
result is algebraic and is **proved that way**: substituting k₂ = r_m/C_E⁰,
k₋₁ = k₂/ρ, k₁ = (k₋₁+k₂)/K_M into eq. (1.5.1-13) returns eq. (1.5.1-15)
identically and ∂r/∂ρ ≡ 0, in sympy, in the notebook's symbolic cell under
*The published model*. The
numerical sweep is a check on the code, not the evidence: hold r_m and K_M fixed,
sweep ρ = k₂/k₋₁ over **eight decades**, evaluate the rate law from the
elementary constants over four decades of concentration, and the curves are
identical to double precision — maximum relative difference **0.0** on the nine
swept ρ, and at most **4.4×10⁻¹⁶** over 2000 random (ρ, C_E⁰) draws, where the
float round trip through the elementary constants is no longer exact. That exact
zero is a property of the values swept; the degeneracy is a property of the
algebra, and both are reported. Meanwhile K_M(BH)/K_M(MM) = 1 + ρ runs from
1.0001 to 10⁴. The map from four elementary constants to the two observable ones
has a two-dimensional fibre and ρ runs along it.

**A batch transient sees ρ, but almost entirely where nobody samples.** With the
observables fixed, C_A(t) spreads by **0.9773 ε** across those eight decades,
ε = C_E⁰/(c_A0 + K_M) — and the maximum sits at **t = 2.0×10⁻⁷ min, inside the
induction layer**, whose duration 1/[k₁(c_A0+K_M)] runs from 2.1×10⁻⁸ min at
ρ = 10⁻⁴ to 2.1×10⁻⁴ min at ρ = 10⁴. A uniform 181-point scan of the 18-minute
run reads that maximum **66.3 % low**, which is why it is root-found on a
logarithmic scan. Froment's one sentence — *"only transient experimentation
(stopped flow, or relaxation …) can help"* — is exactly right, and this is the
number behind it.

**Over the window a slow experiment can see, the requirement is brutal.** For
t ≥ 3 min (the source table's first sample after zero) only **0.242 ε** survives.
Root-found against that table's 5×10⁻⁴ cmol/L print resolution: ε = 2.07×10⁻³ —
and that holds every other constant fixed. Let r_m, K_M and C_E⁰ float and it
becomes **ε = 5.80×10⁻²**. At that loading the quasi-steady-state approximation,
the thing that makes the rate law Michaelis-Menten at all, is itself off by
**88.9 times** the resolution the discrimination is read at (that is against the
free substrate C_A; against total substrate C_A + C_A‑E it is 38.8 times — both
far above 1, which is the claim).

**How much of that is compensation.** Those two ε's are not the same
measurement: 2.07×10⁻³ is a **maximum over t ∈ [3,18] min** of the ρ-spread,
5.80×10⁻² is an **RMS over a 25-point grid starting at t = 0** of a *refitted*
misfit. Measured like for like — same statistic, same grid, nothing free to
compensate — the threshold is **ε = 3.18×10⁻³**. So of the raw factor **28.0**
between the two published numbers, **1.54 is the change of statistic and window**
and only **18.2 is the price of compensation**. 18.2 is the number this page
means by that price.

**On the actual data, ρ is free.** Rawlings & Ekerdt's Exercise 9.15 gives
r_m = 0.162583 cmol/(L·min) and K_M = 0.503055 cmol/L at 4.02×10⁻³ RMS, against
3.90×10⁻² for the best first-order null. Pin ρ anywhere in [10⁻⁴, 10⁴] and let
the rest re-optimise: the achievable sum of squares moves **8.58 %**,
monotonically, no interior optimum. Adding both elementary-step parameters gives
**F = 1.90 against F(2,3) = 9.55** at 95 %.

## What it does not claim

It does **not** validate "the Michaelis-Menten mechanism". What the seven points
support is that a two-parameter saturating law beats the best one-parameter
simple order by 9.71× in RMS — a claim about *saturation*, compatible with every
mechanism in the two-parameter family constructed here, which has both of the
book's derivations as end members.

## Two printed defects, reported and not repaired

Both are settled from the book's own equations, and both were read at digit scale
on 300 ppi crops.

- **p. 24**, verbatim: *"the rate levels off and becomes zero order with respect
  to the reactant, r_A = k₁C_E⁰."* [sic] The limit of the book's own
  eq. (1.5.1-8) is k₂C_E⁰; its eq. (1.5.1-15) calls k₂C_E⁰ "the maximum possible
  rate"; its Problem 2.4(b) on p. 144 asks the student to show that limit; and
  k₁C_E⁰ has units of s⁻¹, which a rate cannot have.
- **p. 25**, verbatim: *"the Lineweaver-Burke plot of 1/r versus C_A."* [sic]
  1/r is exactly affine in 1/C_A (9.67×10⁻¹⁶ off a straight line) and strictly
  convex in C_A, d²(1/r)/dC_A² = 2K_M/(r_m C_A³) > 0. The convexity settles the
  defect and is convention-free. So is the *size*, once the right line is named:
  the statistic is a **maximum**, so the line that minimises it is the Chebyshev
  (minimax) line, and **no straight line comes closer than 31.7 %** of the
  curve's own range over 0.05–1.0 cmol/L — closed form, with three-point
  equioscillation as the optimality certificate. The distance from the best
  **least-squares** line is larger and *is* a convention: **50.3 %** under the
  page's log-weighted convention, **66.9 %** under a uniform one, both in closed
  form with the extremum root-found. Corroborated independently by Levenspiel p. 615
  at 600 ppi native — where, incidentally, **no surname is printed at all**: the
  string "Burk" occurs nowhere in that book.

And one that is *not* a defect: the clean, born-digital text layer silently drops
the leading term of eq. (1.5.1-17) on p. 26, making a correct equation look wrong.
The render shows it is fine.

## Sources

**Read from** — Froment, G. F., De Wilde, J. & Bischoff, K. B., *Chemical Reactor
Analysis and Design*, 3rd edn, Wiley (2011), ISBN 978-0-470-56541-4, section
1.5.1, book pp. 23–26, eqs. (1.5.1-1) to (1.5.1-15). The filename says
"Froment_Bischoff"; this is the **third** edition and De Wilde is a full author.

**Origins, cited and NOT consulted** — Michaelis, L. & Menten, M. L., *Biochem.
Z.* **49**, 333–369 (1913); Briggs, G. E. & Haldane, J. B. S., *Biochem. J.*
**19**, 338–339 (1925). Neither is on disk; neither was read. Everything the page
says about them is bounded by the two clauses the book prints. **The book itself
prints no citation for either** — a full-text search of all 902 PDF pages for
`michaelis|menten|briggs|haldane` returns fourteen lines (five in §1.5.1–1.5.2,
three in Problem 2.4, three in later running text, three in the Subject Index),
the same fourteen in all three `pdftotext` modes, none a bibliographic entry; the
Chapter 1 reference list was read
in full at the alphabetical positions the entries would occupy, and Graef &
Andrews, Monod and Williams — the other named results of section 1.5 — *are*
there.

**Data** — Rawlings, J. B. & Ekerdt, J. G., *Chemical Reactor Analysis and Design
Fundamentals*, 2nd edn (2025 printing), Exercise 9.15: seven c_S-versus-time
pairs. **Provenance unstated.** The exercise says only *"The following
measurements of cS versus time were taken in your laboratory"* — no citation, no
enzyme, no loading, no temperature, no pH. The page therefore does not call it a
laboratory record. It is not an exactly-rounded model evaluation either (maximum
residual 13.87× the print resolution), and the fitted constants sit within 0.61 %
of the round values K = 2 L/cmol, K_M = 0.5 cmol/L. Nothing the page concludes
depends on it.

**Levenspiel**, *Chemical Reaction Engineering*, 3rd edn (1999), p. 615, is used
for exactly two printed things — the reciprocal plot's axes, and his claim that
the direct C_A-versus-τ fit "is direct, is less prone to fiddling, and is more
reliable" (measured here: the reciprocal construction run as instructed returns
K_M **46.7 % low** and r_m **34.1 % low**, of which +14.0 points is the
central-difference rate estimate alone). He is **not** a second source for the
derivation and nothing is adjudicated between the two books.

## Files

| file | what it is |
|---|---|
| `index.ipynb` | the page; nine sections, executed clean, 195 s |
| `build_page.py` | regenerates `index.ipynb` |
| `data/rawlings-ekerdt-ex9.15-batch.csv` | the seven printed points, with a sidecar that leads on the unstated provenance |
| `data/froment-1.5.1-printed.csv` | the first fifteen of the section's eighteen numbered equations, and nine prose claims, verbatim, defects flagged and not repaired |
| `agreement.json` | 53 metrics; 7 sit below CI's `ABS_FLOOR` and are named structural with above-floor companions |

## Validation summary

- ten symbolic identities close eqs. (1.5.1-2) → (1.5.1-15), both routes;
- two conservation balances in the pymrm mass-action marcher (1.2×10⁻¹⁷,
  1.7×10⁻¹⁵) on a model that integrates all four fields independently;
- backward Euler at observed order **0.9980**, Richardson vs Radau
  **4.19×10⁻⁷**, LSODA vs Radau **3.85×10⁻¹¹**;
- **two second, independent computations**: the fit re-done through the pymrm
  marcher instead of the Lambert-W closed form (agreeing to 5.9×10⁻⁸ / 5.6×10⁻⁹,
  sharing no line of code), and ρ* from an analytic sensitivity with its extremum
  root-found (2.1934×10⁻³ vs 2.1951×10⁻³, **0.080 %** apart);
- a defect-injection table of **52 rows: 51 defect injections, all 51 of which
  move their metric** by more than 1×10⁻⁶ relative, covering 42 of the 53
  metrics, **plus one robustness row** — written to show a metric does *not*
  move — which is labelled as such and is **not** counted as coverage anywhere.
  The coverage map prints how hard each metric's strongest row hits it; the
  other eleven metrics are named individually, and **both** generic
  not-coverage labels — the `UNCOVERED` fallback and the "robustness rows only"
  branch — are asserted never to be used;
- an assertion cell checking **every number quoted in prose, in `meta.yaml` and
  in this README** against the live computation (68 checks), including the four
  values that exist only as break-table entries, which are looked up from the
  table rather than retyped. The page fails to execute if any drifts.

## Reuse

Fit eq. (1.5.1-15) — or its integral, for a batch run — directly and report
(r_m, K_M). Do **not** report k₋₁/k₁ or (k₋₁+k₂)/k₁ as though your data
distinguished them. If you are linearising, the abscissa is 1/C_A; and consider
not linearising. If you are checking the quasi-steady state, the group is
C_E⁰/(c_A0 + K_M) — comparing C_E⁰ with the substrate alone mispredicts the error
by up to a factor **95** once K_M is comparable with c_A0. To put this rate law
inside transport, take `BatchReduced.source` to `J4.7` or `B1.1` unchanged.
