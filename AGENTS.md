# Rules for AI agents that change this repository

## Layout

| Path | Content | May change? |
|---|---|---|
| `SPEC.md` | what the library does, its inputs, outputs and assumptions | yes, together with the code |
| `dfppi/` | the library: `api.py` (public functions), `cli.py`, `crossfit.py`, and the estimators `aipwf.py`, `drf.py`, `rf.py`, `method.py` | yes, under the invariants below |
| `tests/` | `test_paper.py` (the code computes the paper's quantities), `test_api.py` (the library equals the experiment code) | add tests; do not weaken them |
| `experiments/` | the code that produced Sec. 5 and Appendix C, and its figures and tables | only to reproduce or extend the study |
| `manuscript/` | the working manuscript (PDF only) | no |
| `schemas/` | JSON Schemas of the outputs | with `api.py` |

## Commands

```
python3 -m pip install -e .                    # or run from the repository root without installing
python3 -m unittest discover -s tests          # 28 tests, about 3 seconds after the one-time data download
python3 -m dfppi ate --trial t.csv --obs o.csv --trial-propensity 0.5
```

## Invariants (a change that breaks one is a bug)

1. With `trial_propensity`, `covariate_shift` and the rows fixed, `fuse_ate` and `fuse_cate` with `shuffle=False`
   return what `experiments/experiment.replicate` records for the same draw (`tests/test_api.py`).
2. The first rotation (blocks B1, B2, B3 as nuisance, tuning, evaluation) is the single-split estimator of the
   paper's theory.
3. lambda = omega = 0 gives the trial-only AIPW on the evaluation sample; with r_0 = 1 known the calibration
   variance V_cal and B_cal are zero.
4. The trial-only learners never use the OBS: no g in their sieve, no OBS rows in their loss.
5. The CATE of a cross-fitted learner is the mean of the three fitted functions; the coefficients are never
   averaged (each rotation has its own principal components).
6. New estimators in `experiments/experiment.py` are new rows; they must not consume the random numbers of a
   replication, so that the existing rows stay identical.
7. Every random choice is seeded; results do not depend on the number of worker processes.

## Style

Each module starts with a docstring that maps its functions to the paper (definition, equation, algorithm line).
Keep that map current when you change a function. Name variables as the paper does (lambda as `lam`, omega,
`mu_R`, `mu_O`, `g`, `r`). Report numbers only from code that recomputes them from the saved results.
