# PDFs wanted

Papers needed to build planned gallery pages. Publishers block scripted
downloads (ScienceDirect returns 403 to non-browser requests regardless of
institutional IP), so these have to be fetched through a browser.

**Where to put them:** anywhere outside the repository — the session scratchpad
is fine. Never in the repo: `.gitignore` blocks `*.pdf` and
`scripts/check_metadata.py` errors on committed PDFs, deliberately.

**A better long-term fix.** Elsevier's TDM API (`dev.elsevier.com`) issues
institutional API keys for exactly this purpose — programmatic full-text access
aligned with the Art. 15n research exception, and not bot-blocked. Wiley and
Springer have equivalents. One key would unblock most of this list permanently.

Ordered by value. "Need" says what the page actually requires, which is often
less than the whole paper.

## Priority 1 — unblocks the next pages

| # | Paper | DOI / locator | Publisher | Need |
|---|---|---|---|---|
| 1 | **Xu & Froment (1989)**, *Methane steam reforming, methanation and water-gas shift: I. Intrinsic kinetics*, AIChE J 35(1) 88–96 | `10.1002/aic.690350109` | Wiley | Rate-parameter table (k, K_ads, activation energies) and the rate-vs-conversion data or parity plot. Kinetics papers usually tabulate — if so no digitising is needed at all. Page `C2.1`. |
| 2 | **Krishna & Ellenberger (1996)**, *Gas holdup in bubble column reactors operating in the churn-turbulent flow regime*, AIChE J 42(9) 2627–2634 | `10.1002/aic.690420923` | Wiley | Holdup-vs-superficial-velocity data for the small/large bubble split, plus column diameters. Page `F1.4`. |
| 3 | **Weisz & Hicks (1962)**, *The behaviour of porous catalyst particles in view of internal mass and heat diffusion effects*, Chem Eng Sci 17(4) 265–275 | `10.1016/0009-2509(62)85005-2` | Elsevier | Their computed η(φ) curves for γ = 20. Page `B1.1` is already published and validated against exact solutions; digitising these would add a third, independent comparison and let me check my β and γ conventions match theirs. **Nice-to-have, not blocking.** |

## Priority 2 — strong pages, ready to build once available

| # | Paper | DOI / locator | Publisher | Need |
|---|---|---|---|---|
| 4 | **Van Welsenaere & Froment (1970)**, *Parametric sensitivity and runaway in fixed bed catalytic reactors*, Chem Eng Sci 25(10) 1503–1516 | `10.1016/0009-2509(70)80070-1` | Elsevier | The two runaway criteria and the worked o-xylene case (kinetics, bed and coolant conditions). Page `D2.2`. |
| 5 | **Wakao & Funazkri (1978)**, *Effect of fluid dispersion coefficients on particle-to-fluid mass transfer coefficients in packed beds*, Chem Eng Sci 33(10) 1375–1384 | `10.1016/0009-2509(78)85120-3` | Elsevier | The Sh-vs-Re dataset they correlated (a compilation of many studies — likely a large scatter plot worth digitising). Page `A3.4`. |
| 6 | **Itoh (1987)**, *A membrane reactor using palladium*, AIChE J 33(9) 1576–1578 | `10.1002/aic.690330921` | Wiley | Conversion-vs-parameter data showing the equilibrium shift, plus membrane permeability and reactor dimensions. Page `H1.4`, paired with `H1.1`. |
| 7 | **Oh & Cavendish (1982)**, *Transients of monolithic catalytic converters*, I&EC Prod Res Dev 21(1) 29–37 | `10.1021/i300005a006` | ACS | Light-off temperature curves and the Voltz kinetics parameters. Page `I1.2`. |
| 8 | **Kunii & Levenspiel (1968)**, *Bubbling bed model*, I&EC Proc Des Dev 7(4) 481–492 | `10.1021/i260028a001` | ACS | The interchange-coefficient definitions and any conversion data. Page `E2.1`. Note the 1991 book *Fluidization Engineering* is the fuller source if it is easier to reach. |

## Priority 3 — useful later

| # | Paper | DOI / locator | Publisher | Need |
|---|---|---|---|---|
| 9 | **Graaf et al. (1988)**, *Kinetics of low-pressure methanol synthesis*, Chem Eng Sci 43(12) 3185–3195 | `10.1016/0009-2509(88)85127-3` | Elsevier | Rate parameters + data. Page `C2.4`, to be paired with `C2.5`. |
| 10 | **Vanden Bussche & Froment (1996)**, *A steady-state kinetic model for methanol synthesis…*, J Catal 161(1) 1–10 | `10.1006/jcat.1996.0156` | Elsevier | Rate parameters + parity data. Page `C2.5`; the direct comparison with `C2.4` is the point. |
| 11 | **van Deemter, Zuiderweg & Klinkenberg (1956)**, *Longitudinal diffusion and resistance to mass transfer…*, Chem Eng Sci 5(6) 271–289 | `10.1016/0009-2509(56)80003-1` | Elsevier | HETP-vs-velocity data. Page `J1.10`. |
| 12 | **Ergun (1952)**, *Fluid flow through packed columns*, Chem Eng Prog 48 89–94 | not indexed by DOI | — | The original Δp dataset. **Likely hard to obtain** — *Chemical Engineering Progress* from 1952 is poorly digitised. If it proves awkward, skip it: modern open Δp datasets exist and page `A1.1` can be built against one of those instead. |

## Already obtained

- Duncan & Toor (1962), AIChE J 8(1) 38–41 — used for page `A4.9`. Figure 2
  digitised, all parameters verified.
- **All three priority-1 papers**, delivered 2026-07-26. Inspected but not yet
  used — see [`pdf-findings.md`](pdf-findings.md) for what is in each, the OCR
  obstacle on the Xu & Froment scan, and the recommended order of work. Items 1,
  2 and 3 in the priority-1 table above can be considered supplied.

## What I do with each

For every paper I check the parameter table first and the data second — wrong
parameters make the simulation wrong, whereas missing data only costs the
overlay. If measurements are tabulated I transcribe them directly (provenance
tier 2, no digitising error). If they are only plotted, I digitise
programmatically and validate the extraction by closure, monotonicity and any
available conservation law, exactly as recorded in
`pages/A4.9-duncan-toor/data/duncan-toor-1962-run1.meta.yaml`.
