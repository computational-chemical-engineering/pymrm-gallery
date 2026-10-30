# Measured reproducibility of the agreement metrics

`spread.json` holds, per page and metric, how far the metric moved between runs on
different arithmetic. It is written by `scripts/measure_spread.py` and read by
`scripts/check_claims.py`; see [`../claims-design.md`](../claims-design.md).

Samples in the current file:

| label | run |
|---|---|
| `default` | 2026-10-29, Threadripper 2950X (no AVX-512), OpenBLAS kernel chosen by default (Haswell), threads by default, suite `.venv` (pymrm 3.1.1.dev6, numpy 2.4.6, scipy 1.18.0) |
| `sandybridge` | the same with `OPENBLAS_CORETYPE=Sandybridge` |
| `nehalem` | the same with `OPENBLAS_CORETYPE=Nehalem` |

Not yet varied: NumPy's own SIMD dispatch (AVX-512 on many CI runners). Add a CI
runner's run as the next sample before relying on a spread of 0.
