# DataFusionPPI: prediction-powered data fusion for treatment effect estimation

[![tests](https://github.com/CausalDataScience/DataFusionPPI/actions/workflows/tests.yml/badge.svg)](https://github.com/CausalDataScience/DataFusionPPI/actions/workflows/tests.yml)

Code and working manuscript of "Prediction-Powered Data Fusion for Treatment Effect Estimation" (an earlier
version was submitted to AISTATS 2027). The method estimates the treatment effect in the population of a randomized trial and uses a larger
observational sample, which may be confounded, to reduce the error without changing the target.

* **AI agents:** read [`SPEC.md`](SPEC.md) (how to call the library, assumptions, when not to use it) and
  [`AGENTS.md`](AGENTS.md) (rules for changing the code).
* **Paper:** [`manuscript/`](manuscript/) (the working manuscript, PDF).
* **Reproducing the paper:** [`experiments/README.md`](experiments/README.md).

## Quick start

```
python3 -m pip install -e .
```

```python
import numpy as np, dfppi

rng = np.random.default_rng(0)
x_R, x_O = rng.normal(size=(600, 5)), rng.normal(size=(15000, 5))
a_R, a_O = rng.binomial(1, 0.5, 600), rng.binomial(1, 0.5, 15000)
y_R = x_R[:, 0] + a_R * (1 + x_R[:, 1]) + rng.normal(size=600)
y_O = x_O[:, 0] + a_O * (1 + x_O[:, 1]) + 0.5 * a_O + rng.normal(size=15000)   # a confounded OBS
trial, obs = {"x": x_R, "a": a_R, "y": y_R}, {"x": x_O, "a": a_O, "y": y_O}

out = dfppi.fuse_ate(trial, obs, trial_propensity=0.5, covariate_shift=False)
print(out["estimate"], out["ci"], out["trial_only"]["estimate"], out["warnings"])

fit = dfppi.fuse_cate(trial, obs, trial_propensity=0.5, covariate_shift=False)
print(fit.predict(x_R[:5]), fit.summary()["chosen"])
```

From the shell: `python -m dfppi ate --trial trial.csv --obs obs.csv --trial-propensity 0.5` (JSON on stdout).

## Layout

| Path | Content |
|---|---|
| `SPEC.md`, `AGENTS.md` | specification and change rules, written for AI agents |
| `dfppi/` | the library: `fuse_ate`, `fuse_cate` (`api.py`), the command line (`cli.py`), and the estimators of Algorithms 1 and 2 (`aipwf.py`, `drf.py`, `rf.py`, `method.py`, `crossfit.py`) |
| `schemas/` | JSON Schemas of the outputs |
| `tests/` | unit tests against the paper's definitions and against the experiment code |
| `experiments/` | the simulation study of Sec. 5 and Appendix C: code, figures, tables |
| `manuscript/` | the working manuscript (PDF, arXiv version) |

Requirements: Python 3.9 or later with numpy, scipy and scikit-learn; the experiments also need pandas and
matplotlib (`pip install -e ".[experiments]"`).

## Authors

* Yonghan Jung, University of Illinois Urbana-Champaign, yonghan@illinois.edu
* Shu Yang, North Carolina State University, syang24@ncsu.edu

## Citation

```bibtex
@misc{jung2026dfppi,
  title  = {Prediction-Powered Data Fusion for Treatment Effect Estimation},
  author = {Jung, Yonghan and Yang, Shu},
  year   = {2026},
  note   = {Manuscript; the arXiv identifier will be added once it is posted}
}
```

`CITATION.cff` holds the same information for GitHub's "Cite this repository".

## License

MIT; see `LICENSE`.
