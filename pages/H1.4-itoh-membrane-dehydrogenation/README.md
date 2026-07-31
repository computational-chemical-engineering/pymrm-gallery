# H1.4 — Itoh Pd membrane reactor: shifting a dehydrogenation equilibrium

Reproduces Itoh (1987), *A membrane reactor using palladium*, AIChE Journal
**33**(9) 1576–1578, [doi:10.1002/aic.690330921](https://doi.org/10.1002/aic.690330921).

Cyclohexane dehydrogenates over Pt/Al₂O₃ inside a 200 µm palladium tube;
hydrogen — and only hydrogen — leaves through the wall into an argon purge, so
a reaction whose equilibrium conversion is 18.7 % runs to a measured 99.7 %.
The page solves the paper's two-channel plug-flow model (Eqs. 1–8) as a pymrm
finite-volume/Newton system, lands on the measured conversion with nothing
fitted, and then measures how much that agreement can actually resolve — which
is the part worth reading.

## Validation route

Stated numerical results plus internal identities — no figure digitisation:

1. the stated measured conversion (0.997) at the stated flows — the one
   experimental datum, modelled as 0.9983 (+0.13 %). **Read the resolving power
   with it:** at this operating point the model sits exactly on its algebraic
   fast-permeation/equilibrium ceiling, so varying α_H over a factor 33 or
   `k·V_r` over a factor 100 does not move the number at all, and K_p may be
   divided by 3 or multiplied by 10 with the deviation staying inside 0.3 %. On
   the quantity that does vary, `1 − X`, the model is 45 % low. The comparison
   establishes that the reactor reached its co-current fast-permeation
   asymptote and that the model puts that asymptote in the right place; it is
   not evidence about the kinetics or the permeance;
2. the printed permeation constant α_H recomputed from Eq. 1's own inputs
   (0.03 %) — a genuine internal identity that propagates to nothing;
3. K_p(473 K) reconstructed from the paper's stated equilibrium conversion
   (18.7 %), cross-checked against an independent van 't Hoff estimate
   (ratio 1.58);
4. FV/Newton vs an independent stiff-IVP integration, first-order grid
   convergence (ratio 2.0), hydrogen conservation reported separately for the
   co-current case (exact by construction — the sums telescope) and the
   counter-current case (a genuine O(h) error, with its ladder), and a
   kinetics-free fast-permeation ceiling that bounds the whole operating map.

## Two misprints in the paper, and one printed value confirmed

The scan is a 300 dpi bilevel CCITT image, so a 600 dpi render is a 2×
upsample: every marginal glyph was read on the native bitmap.

- **Eq. 5** prints `v_H = u_H - 3(u_C0 - u_C)`. Substituting the paper's own
  hydrogen balance `u_H + v_H = 3(u_C0 - u_C)` into it gives `v_H = -v_H`, i.e.
  the printed equation forces `v_H ≡ 0` and contradicts Eq. 3's permeation term.
  Sign flip. No convention escape: Eq. 8, Figure 3 and Figure 1 all establish
  co-current flow, and the sentence below Eq. 6 makes `v_i` the separation-side
  flow rate, positive in `+L`.
- **K_P = 4.89e35 exp(+3190/T) Pa³** contradicts the paper's own 18.7 %
  equilibrium conversion (it implies ~100 %). A units re-reading leaves it 27
  orders out, and the exponent's sign is independently wrong for an endothermic
  reaction. Decision table in the notebook. `K_p(473) = 2.357e11 Pa³` is
  reconstructed from the stated equilibrium conversion — every input printed.
- **`v_A0 = 11.8e-5` mol/s is NOT a misprint.** The superscript reads `-5`
  cleanly on the native bitmap. The page confirms it rather than correcting it,
  and the confirmation is kinetics-free: under a `10⁻⁶` reading the
  thermodynamic ceiling — infinite kinetics *and* infinite permeance — is
  X = 0.716, so the measured 0.997 would be impossible. Figure 4's
  hand-lettered abscissa exponent *is* illegible on the scan; nothing here
  depends on reading it.

Native-resolution verification crops are in `../review/` (git-ignored: they are
page images).

## Files

- `build_page.py` — writes `index.ipynb` (run from this directory)
- `index.ipynb` — the executed page
- `data/itoh-1987-stated-values.csv` + `.meta.yaml` — the stated values
  (tier 2: numbers printed in the paper) with full provenance
- `agreement.json` — CI regression metrics, including the sensitivity spans, so
  a later reader cannot mistake +0.13 % for evidence about the kinetics

## Honesty notes

The single measured point sits exactly on the model's thermodynamic ceiling.
Check 1 of the notebook quantifies what that costs and the page's summary lines
are written to match it. Figure 4's experimental markers at other purge rates
were deliberately not digitised (ranked validation policy). The counter-current
comparison in "What pymrm adds" is a pure prediction, computed on a converged
grid (n_z = 3200, ladder shown to 6400) with its hydrogen closure reported.

A non-blocking follow-up is recorded in `queue_cases/H1.4.yaml`: a trial
digitisation of Figure 4's *calculated* curves, which would add a check that
does discriminate K_p, is staged unreviewed in `../review/` and is deliberately
not used anywhere on this page.
