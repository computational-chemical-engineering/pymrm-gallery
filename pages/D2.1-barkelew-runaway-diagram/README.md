# D2.1 — Barkelew's runaway diagram

The 1959 chart that collapses the runaway boundary of a cooled tubular reactor
onto two dimensionless groups — rebuilt from the paper that reprints it, and
tested where a similarity reduction is supposed to fail.

- **Structures:** `S2` (plug flow with reaction)
- **Reference (the origin of the result, NOT consulted):** Barkelew (1959),
  Chem. Engng Prog. Symp. Ser. 55 37
- **Read from:** Van Welsenaere & Froment (1970), Chem Eng Sci 25(10) 1503–1516,
  doi:10.1016/0009-2509(70)85073-4 — their §5 is a comparison with Barkelew and
  reprints his diagram as their Fig. 9
- **Runtime:** ~26 s

## The reprint route, and what it does and does not licence

Barkelew's paper is pre-DOI, is not on disk and has no open-access route. It was
**not** consulted. Everything on this page comes from Van Welsenaere & Froment's
Section 5, which names his method, reproduces his diagram, says what his axes
mean, states the one numerical property of his critical locus, and applies its
own criterion in his notation. That satisfies the `AGENTS.md` test: the reprint
both *names* and *uses* the result.

What it does not supply is a definition. Their notation list gives "*N*
dimensionless group defined by Barkelew" and nothing more, so **the
dimensionless temperature is a reconstruction**, labelled as such on the page.
Most of it is forced by print, and one step of it is not:

- **Forced.** The *N/S* sentence fixes the cooling term as *Nτ* and the
  generation as *S y* e^τ, hence τ ∝ (T − T_w) with τ(T_w) = 0. The *S*
  sentence — the one that names the wall temperature, by feeding the reactants
  at it — plus their Eq. 32 gives S = a(ΔT)_ad/T_w². And their own printed
  τ_m = 1 for criterion 1 forces the exponent to be τ itself and not a multiple
  of it, because the maxima curve y ∝ τ e^(−λτ) peaks at 1/λ.
- **Not forced.** The temperature at which the Arrhenius exponent is
  linearised. a(T−T_w)/T_w² is the tangent at the wall, but the tangent at T_M
  and the chord from T_w to T_M satisfy every printed sentence equally. That is
  a Frank–Kamenetskii convention, not a constraint from the paper — so the page
  says so rather than calling the whole definition pinned.

The unforced step is nevertheless **testable**, and the page tests it twice: as
a deliberate defect in the break table (1.77 % → 7.22 %), and against the
blow-up positions of their Fig. 10 (see below).

## Agreement

**The tangency hot spot against their stated 1.275.** Computed: 1.2879, 1.2763,
1.2547, 1.2293 at S = 4, 8, 16, 32 — mean deviation **1.575 %** over the four
curves Barkelew drew, **3.512 %** over S = 3 to 200, where it falls to 1.1622.
Nothing was fitted and 1.275 enters no computation.

**The tangency is converged, not merely computed.** The envelope is located from
a difference in ln S over a step h, which carries an O(h²) truncation error; the
page Richardson-extrapolates the pair (0.01, 0.005) and a second pair
(0.005, 0.0025) confirms it to 3e-5 %. This matters for the conclusion, not just
the digits: the bias grows with S, in the *same* direction as the drift the page
reports, so a plausible fixed h = 0.04 would have read the fall from S = 4 to
S = 200 as −10.28 % instead of −9.76 %, inflating the page's own result by 5 %.

**Barkelew's criterion against Van Welsenaere & Froment's.** Their critical
inlet pressures agree to **1.774 % mean, 4.080 % worst** over T_w = 600–700 K,
which is what their *"agree extremely well"* is worth. The two cross at
**645.1 K** (printed on the page from a bracketed root-find): Barkelew is more
conservative below and less above. At 625 K the criterion-1 back-integration
reproduces their printed 0.01651 atm to **0.01 %**.

**And the disagreement factorises.** §5's comparison mixes two approximations.
Computing the intermediate case separates them exactly:

| T_w | rate law alone | criterion alone | total |
|---|---|---|---|
| 600 K | −5.06 % | +2.81 % | −2.39 % |
| 650 K | −6.17 % | +6.96 % | +0.36 % |
| 700 K | −7.69 % | +12.76 % | +4.08 % |

Opposite signs at all 21 wall temperatures, and the total smaller in magnitude
than both components at all of them. The two errors partly cancel, which is why
the methods agree better than either approximation does alone. Neither paper
says this.

**pymrm against an independent route.** The finite-volume reactor and an
adaptive phase-plane quadrature that never forms a grid — and that transcribes
the reduced model separately rather than calling into it — agree to **0.260 %**
at n_zeta = 1600, converging at order **0.99** (first-order upwind, as expected).

## Provenance: tier 6, not experimental

Barkelew's chart correlates his own numerical integrations; every Van Welsenaere
& Froment number is likewise their own computed value. Neither paper contains an
experiment. **No figure was digitised into the committed data** — the only
values in the CSVs from Figs. 9 and 10 are labels printed inside them.

One measurement was taken off a curve, and it is quarantined. The four positions
at which the Fig. 10 profiles leave the printed 800 K axis are the only printed
object that can discriminate the unforced step of the reconstruction. They were
measured by `queue_cases/D2.1/review/extract_figure10.py` (0.4082, 0.4725,
0.5448, 0.7555 m), the page computes 0.4090, 0.4613, 0.5398, 0.7511 m — 0.2 to
2.4 % — and every alternative linearisation fails qualitatively: the T_M tangent
never reaches 800 K in three of the four cases, the chord in two, and dropping
the square runs away in none of them. **This is awaiting the maintainer's visual
review** of the numbered overlay in `review/`; until then the page presents it
as an indication, the case file carries a non-blocking follow-up, and nothing
else on the page depends on it. The reconstruction stands on the printed
sentences either way.

## What the page adds

Barkelew's criterion is a tangency between a curve and the envelope of a family,
which in 1959 had to be drawn and judged by eye. Computed instead — and
extrapolated to zero step size — the tangency hot spot becomes a quantity with a
value at every S, and the collapse is seen to drift monotonically downward past
S ≈ 4, so 1.275 is a good summary of his four curves and a steadily worse one
outside them. The page also separates the two approximations §5's comparison
mixes, his tangency criterion and his exponential rate law, and finds their
errors partly cancel.

## What the checks cannot see — measured, not asserted

The notebook prints a deliberate-break table and, separately, a table of the
defects the checks are *blind* to. Four results matter to anyone reading a
headline number:

1. **The 1.275 check is blind to the dimensional groups.** Linearising the rate
   about T_M instead of T_w moves the critical-pressure comparison from 1.774 %
   to 7.221 %, and leaves the 1.275 check the *identical computation* — the
   reduced system contains no reference temperature.
2. **It is nearly blind to the reaction order too.** Second order in the
   reactant gives 1.339 %, better than the baseline. Over orders 0.5 to 2 the
   tangency hot spot stays within 1.211–1.315, so their 1.275 characterises the
   exponential temperature dependence and almost nothing else. The abscissa of
   the tangency does move — by a factor 2.24 at S = 4 — so the chart itself is
   order-specific even though that one number is not.
3. **The critical-pressure check does not move at all under a rescaling of B.**
   That is structural: p⁰ ∝ 1/B on both sides, so it says nothing whatever about
   (−ΔH), ρ_b or p_B⁰. It is also nearly blind to C — U × 1.2 moves it by 6 %,
   and the c_p × ρ_g trap below by 14 %, both inside anything a reader would
   call agreement. The check that catches that trap is the parameter reading,
   which goes from 0.090 % to 12.4 %.
4. **One break row is printed because it cannot fail.** Dropping the square in
   τ = a(T−T_w)/T_w² multiplies p⁰ by exactly 1/T_w and nothing else, so its
   99.845 % is arithmetic — the page reproduces the same number from the
   baseline table without re-running an ODE. It tests nothing about the
   reconstruction.

The check on the pymrm route, by contrast, has more power than a reader might
assume: the two routes do **not** share a `reaction` function, and injecting a
flipped cooling sign or a 10 % error in `BarkelewTube.reaction` alone sends it
from 0.260 % to `inf`.

## Two traps, both inherited from `D2.2`

1. **`c_p` is printed as kcal/m³·°C — it is already volumetric.** Multiplying by
   `rho_g` changes `C = 2U/(c_p R)` by 1.293 and moves the boundary. It is
   caught by the parameter reading (check 5, ln K −2.0568 → −1.7999 against the
   printed −2.055), not by the critical-pressure comparison, which barely
   notices.
2. **The Elsevier full-text API is useless for numbers here.** It discards the
   1970 mid-dot decimal separator: `1275` for 1.275, `001651` for 0.01651. It is
   excellent for prose — it is how §5 was located — and every number was read
   from a 600 dpi page render.

## Rebuilding the page

```bash
python build_page.py                    # regenerate index.ipynb from the builder
python ../../scripts/run_pages.py       # execute it
python ../../scripts/check_agreement.py # metric regression check
```

No PDF is needed to rebuild: both datasets are committed CSVs. The PDF *is*
needed to re-run `review/extract_figure10.py`, which is review material and not
part of the page build.

## See also

[`D2.2`](../D2.2-van-welsenaere-froment-runaway/) is the other half of the same
paper — the explicit closed-form criteria in dimensional variables. The two
pages are two answers to one question and disagree by a few percent in a
direction that depends on the wall temperature.
