# Design: published numbers that hold on any machine

Status: REVISED 2026-10-29 after review (`think`, verdict "adopt with changes"; the
changes are folded in below). Part of step 0 of [`refresh-plan.md`](refresh-plan.md).

## 1. The problem, measured

Every page was run on one machine (Threadripper 2950X, no AVX-512) under three
OpenBLAS kernels (`OPENBLAS_CORETYPE` default Haswell, Sandybridge, Nehalem), which
reproduces the cross-runner failures CI has shown since at least 2026-09-14.

**Pages.** Default kernel: 82 of 82 pass. Sandybridge: 6 fail (`D1.1 F3.2 J4.1 J4.6
J4.8 J5.3`). Nehalem: 5 fail (the same without `J4.8`). Causes: a prose check that
asserts a round-off-level number (`J4.6`, `J5.3`, `F3.2`); the set of metrics decided
by comparing a value with `ABS_FLOOR` (`J4.1`, `J5.3`); a division by a quantity that
is exactly zero on one kernel (`J4.8`); a bisection that tracks where Newton converges
and which branch it reaches, instead of the physical threshold (`D1.1`).

**Metrics.** 3037 metrics in 82 `agreement.json` files. Relative spread across the
three kernels:

| spread | metrics |
|---|---|
| any | 713 |
| > 1e-9 | 216 |
| > 1e-6 | 128 (on 27 pages) |
| > 5e-4 | 96 |
| > 5 % (fails CI today) | 64 |

The size of a value does not predict its spread. Of the 286 metrics below `ABS_FLOOR =
1e-12`, 33 move by more than 1e-6; of the 2403 above 1e-4, 36 do: optimiser and
root-finder results at their tolerance (`F3.2`, about 1e-5), convergence orders
(`A4.3` 1e-4, `B1.7` 2.3 %), a nonuniformity ratio (`A4.9` 4.4 %), and quantities that
change by factors (`B2.3` optimiser stalls). Counts change too: `B2.3`
`sep_optimiser_overstated` 9 against 8, `J3.4` `salt_drift_blind_named_classes_ratio`
1.33 against 0.75. And a spread of zero is weak evidence: 129 metrics are exactly 0 on
all three kernels, but NumPy's own SIMD dispatch (AVX-512 on many CI runners) was never
varied on this machine.

**Prose checks.** 62 pages compare numbers written in the prose with computed ones,
each in its own way (lists of `(name, written, computed, tol)`, bare `assert`s), with
the tolerance set by the last written digit. That fails whenever the written digits go
below the quantity's noise. About 24 `build_page.py` files also read another page's
`agreement.json`.

## 2. Principles

1. **Every published number says what it is:** a value, a bound, or an exact result.
   This decides how it is checked, at run time and in CI.
2. **Classify by role first, by measured spread second, never by size.** A Newton
   residual, an identity or a two-route gap is a bound whatever its size; an error
   that is at round-off (`1e-14`) is a bound too. Only a quantity whose size is set by
   physics, data or discretisation can be a value.
3. **A value must be admissible:** `rtol*|x0| + atol >= 100 * s`, where `s` is the
   measured absolute spread, and `x0 != 0` or `atol > 0`. With `rtol = 0.05` that is a
   relative spread of at most 5e-4.
4. **A conclusion that flips under round-off is a defect of the page, not of the
   check.** A count, branch or classification that changes between kernels gets fixed
   in the numerics or withdrawn. Thresholded results carry their margin (section 3).
5. **One source per number.** The page registers each number once; `agreement.json`,
   the run-time assertion, the prose check and the break-table coverage all read that
   registration. The set of names is fixed: it never depends on values.

## 3. Kinds

| kind | for | run time | CI against committed baseline |
|---|---|---|---|
| `value` | physical, data and discretisation results above round-off | none | `abs(x - x0) <= rtol*abs(x0) + atol`, default `rtol = 0.05`, `atol = 0` |
| `bound` | residuals, closures, identities, two-route gaps, solver-tolerance effects, optimiser stalls, convergence orders (as a range) | `min <= x <= max`, raises | the same bound; the stored `x` is information only |
| `exact` | integers, booleans, labels | margin check when given, raises | equality |

Every bound states its `basis`: `identity` (zero by construction), `tolerance` (set
by a solver or optimiser tolerance), `theory` (a convergence order, a physical limit)
or `chosen` (a judgement, stated in `note`).

**Structural** (`structural: true`) marks a claim no break-table row can move, an
identity named as such (`AGENTS.md`). It is the only exemption from the coverage rule.

**Thresholded results.** A count or label decided by comparing values with a threshold
records its smallest distance to that threshold:

- `claims.count(name, values, threshold, min_margin)` counts `values > threshold`,
  records `margin = min |v - threshold|` (relative to `|threshold|` when that is not
  0), and raises when `margin < min_margin`.
- `claims.exact(name, label, margin=d, min_margin=m)` does the same for a branch or a
  classification computed elsewhere.

`B2.3`, `J3.4` and `D1.1` failed in exactly this way: the result sat inside the noise
of its own threshold, and nothing said so.

**Convergence orders** are a `bound` with a range from theory, plus a `bound` that the
finest-grid error is at least 1e3 times its round-off. Without the second claim, the
order is fitted to noise; that is what `B1.7`'s 2.3 % spread is.

## 4. `agreement.json`, format 2

```json
{
  "page": "B1.6",
  "format": 2,
  "env": {"pymrm": "3.1.0", "numpy": "2.3.4", "scipy": "1.16.2",
          "openblas_coretype": "Haswell", "threads": 1},
  "metrics": {
    "vs_shooting_reference_n800": {"kind": "value", "value": 3.944e-05,
                                   "rtol": 0.05, "unit": "-"},
    "sweep_worst_newton_residual": {"kind": "bound", "value": 4.42e-11, "max": 1e-08,
                                    "basis": "tolerance",
                                    "note": "newton tol 1e-10, scaled by the largest row"},
    "flux_identity": {"kind": "bound", "value": 0.0, "max": 1e-12,
                      "basis": "identity", "structural": true},
    "rows_caught": {"kind": "exact", "value": 9, "margin": 0.31, "min_margin": 0.05}
  },
  "quoted": {"vs_shooting_reference_n800": ["3.9e-5"],
             "sweep_worst_newton_residual": ["below 1e-8"]},
  "breaks": {"vs_shooting_reference_n800": ["closed-form reference off by 1 %"]}
}
```

- `env` lets the check refuse a baseline produced off the pinned environment.
- `quoted` lists the prose literals, so a cut in prose precision is visible in review.
- `breaks` lists, per claim, the break-table rows that fail that claim's own check.
- The measured spread stays out of this file; it changes with every portability run and
  lives in `docs/reproducibility/spread.json`, written by `scripts/measure_spread.py`.

A bare float, as in every file today, reads as format 1: `value` with the old rule
(`rtol = 0.05`, both values below `ABS_FLOOR` skipped). Unmigrated pages keep working.

In format 2, a metric present in the run but absent from the committed file fails, as
does one that disappeared.

**Loosening.** `check_agreement.py --base origin/main` compares the committed file with
the merge-base and fails on a widened bound or tolerance, a changed kind, a removed
claim, a removed break row or a coarser quoted literal, unless `docs/loosened.yaml`
has an entry (page, metric, date, reason). A widened bound would otherwise hide in the
value churn of the diff.

## 5. Page API in `gallery_utils`

```python
from gallery_utils import Claims, load_claims

claims = Claims("B1.6")
claims.value("vs_shooting_reference_n800", err, rtol=0.05)
claims.bound("sweep_worst_newton_residual", res, max=1e-8, basis="tolerance",
             note="newton tol 1e-10, scaled by the largest row")
claims.bound("flux_identity", gap, max=1e-12, basis="identity", structural=True)
claims.count("rows_caught", row_metrics, threshold=1e-6, min_margin=0.05)

claims.breaks("closed-form reference off by 1 %", "vs_shooting_reference_n800", err_b)
claims.says("3.9e-5", "vs_shooting_reference_n800")   # prose check, section 6
claims.report()
```

- Registering a name twice raises.
- `report()` prints the table and writes `agreement.json` (format 2). It raises,
  listing everything at once, on any violated bound or margin, any failing `says`, and
  any non-structural claim with no break row that fails its check.
- `load_claims(page)` returns `{name: float}` from either format. Pages reading another
  page's file switch to it first, before any page writes format 2. A comparison of a
  live number with another page's committed number is a `bound` with a measured basis,
  never "better than 1e-10" (as `F3.2` does with `F3.1` today).
- `report_agreement(page, metrics)` stays, writing format 1.

## 6. Prose numbers

`claims.says(literal, name, scale=1)` checks a literal as it appears in the Markdown.
The prose check is independent of the CI tolerance: the 5 % CI tolerance would let
"0.7186" stand for a computed 0.69.

Grammar, decidable from the literal:
- a number `L`, with `unit(L)` = the place of its last mantissa digit (`0.7186` gives
  1e-4, `3.9e-5` gives 1e-6, `40` gives 1); a trailing `%` sets `scale = 0.01`;
- a bare integer for an `exact` claim is compared exactly;
- a bound phrase from a closed list: `below`, `under`, `less than`, `at most`, `above`,
  `over`, `more than`, `at least`, `between L1 and L2`.

Rules:
- `value` or `bound` with a number: passes iff `|L*scale - x| <= unit(L)*scale`, one
  unit in the last written place, not half (half a unit fails at a rounding boundary).
- `bound` with a phrase: the phrase must be at least as loose as the declared bound
  (prose may say less than the check guarantees, never more).
- `exact`: equality.
- Offline, `scripts/check_claims.py` requires `unit(L)*scale >= 10 * s` for each
  quoted literal with a measured spread `s`.

Phrases such as "an order of magnitude" cannot be decided; they become a bound on the
ratio and are quoted as such.

`README.md`, `meta.yaml` and `models.yaml` quote page numbers too.
`scripts/check_claims.py` checks every number there that is tagged with a metric name
against `agreement.json` under the same rules. It also flags `ABS_FLOOR` in the text of a
migrated page: 57 README/meta files and 13 `models.yaml` passages say metrics are
"outside CI", which stops being true once the page is migrated.

Considered and not chosen: generating the prose numbers from `agreement.json`. The
notebook must render as static Markdown in Jupyter, GitHub and Colab, so a build step
would write the numbers into the Markdown, a second source of truth. It also removes
the check that a person read and agreed with each number.

## 7. Environment and CI

- **The arithmetic is pinned in `scripts/run_pages.py`,** not only in CI:
  `OPENBLAS_CORETYPE=Haswell`, `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`,
  `OPENBLAS_L2_SIZE=512`, and on x86-64
  `NPY_DISABLE_CPU_FEATURES="X86_V4 AVX512_ICL AVX512_SPR"`, unless already set in the
  environment. A local run and CI then compute the same way.
- **`OPENBLAS_L2_SIZE=512` keeps OpenBLAS's blocking fixed** (measured 2026-10-31).
  OpenBLAS 0.3.34, bundled with numpy 2.5.4, raises the GEMM block sizes (DGEMM P and
  Q to 512) on an AMD CPU that exposes AVX-512 and reports a 1024 KB L2, in every core
  table of a `DYNAMIC_ARCH` build, the forced Haswell one included
  (OpenMathLib/OpenBLAS#6013). The Haswell DGEMM kernel copies its B panel, 96 bytes
  per unit of K, into a stack buffer of 28672 bytes plus 0 to 4095 bytes of
  address-dependent slack: K up to 298 is always safe, 299 to 341 overflows in some
  runs, and 342 or more always overwrites the return address. `getrf` blocks at
  min(n/2 rounded to 8, Q), so `H1.7`'s 800 x 800 `numpy.linalg.solve` (K = 400)
  dies whenever the override fires, and `J1.3`'s 600-row dipole ladder (K = 304)
  in about one run in eight. On GitHub's EPYC 9V74 and 9V45 runners with AVX-512
  exposed, the kernels died three times in 2 of the 7 Haswell runs of the last 16 CI
  runs; where `H1.7` survived on such a runner the override cannot have fired (the
  VM presumably reports other cache sizes; the job now records them). The
  Sandybridge and Nehalem kernels have no such buffer and never died, but ran with
  Q = 512 there. The Haswell, Sandybridge and Nehalem tables use the L2 size for
  nothing else, so on a machine that reports 512 KB the pin changes no bit; scipy
  1.18.1's own OpenBLAS (0.3.31) predates the override. Readers are not pinned to a
  core: an unpinned Zen 4 or Zen 5 machine gets the SkylakeX or Cooperlake table,
  for which the larger blocks are meant. The pin becomes redundant once the lock
  carries an OpenBLAS with a78196c5 (no override when the core type is forced).
- **NumPy's SIMD dispatch is part of the arithmetic** (measured 2026-10-30). GitHub's
  runners are a mix of AMD EPYC (AVX2) and Intel Xeon (AVX-512) machines. On an
  AVX-512 runner NumPy computes `exp`, `log`, `pow` and other ufuncs with different
  code. Without the NumPy pin, the locked job on such a runner showed 19 format-1
  metrics drifting beyond 5 % and two format-2 pages failing (`D1.1`, `F3.2`), while
  the same commit passed on an AMD runner. The AVX2 (X86_V3) and baseline (X86_V2)
  paths differ only in `tanh` on the probe used, so this cannot be reproduced on an
  AVX2 machine; the weekly portability run, with the NumPy dispatch left free,
  is the only sample of it.
- **A lock file** `requirements-lock.txt` pins pymrm (released), numpy, scipy and
  the rest. Without it, a new wheel brings a new OpenBLAS and a moved metric no longer
  means changed code. Baselines are produced in a venv built from the lock (the suite's
  shared `.venv` has a development pymrm). `requirements.txt` keeps the floor for
  readers.
- **Jobs**, each changing one thing:
  - `locked`: the lock, pinned kernel. This is the job that defines and checks the
    baselines.
  - `latest`: the newest pymrm, everything else from the lock.
  - `portability`: the lock, `OPENBLAS_CORETYPE=Sandybridge`. It runs on every PR for
    the pages the PR changes, and weekly for all pages. It is the gate, because bounds
    and `says` raise on readers' machines (Colab included). Weekly it leaves NumPy's
    dispatch free, which on an AVX-512 runner gives the fourth spread sample.
- The agreement check runs with `if: always()` (done, `8af9e04`).

## 8. Migration

1. **Infrastructure, no page changes:** `Claims`, `load_claims`, format 2 and the
   loosening check in `check_agreement.py`, `check_claims.py`, `measure_spread.py`, the
   pin in `run_pages.py`, the lock, the CI jobs. Unit tests on fixtures, including
   `says` literals at rounding boundaries.
2. **Readers:** the 24 pages that read another page's `agreement.json` switch to
   `load_claims`, metrics unchanged.
3. **Fourth spread sample:** one CI run on the default kernel; record whether the 129
   zeros stay zero.
4. **The failing pages,** one per commit, each run under all three kernels: `F3.2 J4.1
   J4.6 J4.8 J5.3`.
5. **The other pages with spread above 1e-6** (27 in all, including `A4.2` with 7
   identities and `I1.2` with 3 closures): per claim, a role, then the admissibility
   rule. Where a value is not admissible, it becomes a bound with a basis or the
   numerics are fixed.
6. **`B2.3` and `D1.1`** are numerics problems (conclusions on optimiser stalls; a
   threshold located near multiplicity), handled with `think`, one at a time.
7. **Every other page** moves to `Claims` in lane A of the sweep.
8. **Learnings to the pymrm plugin:** the kinds, role-first classification, margins
   on thresholded results and the `says` rule go into the style guide and the
   verifier's attack list (`refresh-learnings.md`).
