# dfppi specification (for AI agents)

Read this file before calling or changing the code. It is complete for using the library; the paper (`manuscript/`)
has the proofs.

## 1. What it estimates

The treatment effect in the population of a randomized trial (RCT), using a larger observational sample (OBS) to
reduce the error, while the OBS can never move the target away from the trial population.

| Function | Target | Algorithm (paper) |
|---|---|---|
| `dfppi.fuse_ate` | ATE of the trial population, theta_0 = E_R[Y(1) - Y(0)] | Algorithm 1, AIPWF |
| `dfppi.fuse_cate` | CATE tau_0(x) = E_R[Y(1) - Y(0) \| X = x] | Algorithm 2, DRF (`learner="DR"`) or RF (`"R"`) |

## 2. Call

```python
import dfppi
out = dfppi.fuse_ate(trial, obs, trial_propensity=0.5)            # dict, JSON-serializable
fit = dfppi.fuse_cate(trial, obs, trial_propensity=0.5, learner="DR")
fit.predict(x)                # CATE at the rows of x (numpy array, rows by covariates)
fit.predict_trial_only(x)     # the trial-only learner on the same rotations, for comparison
fit.summary()                 # dict, JSON-serializable
```

```
python -m dfppi ate  --trial trial.csv --obs obs.csv --trial-propensity 0.5 [--treatment A] [--outcome Y]
                     [--covariates c1,c2,...] [--no-shift] [--seed 0] [--level 0.95]
python -m dfppi cate --trial trial.csv --obs obs.csv --trial-propensity 0.5 [--learner DR|R]
                     [--predict points.csv] [--out cate.csv]
```

The command line prints JSON: `{"status": "ok", ...}` with the fields of Sec. 4, or `{"status": "error", "error":
"..."}` with exit code 1.

## 3. Inputs

| Argument | Meaning | Required form |
|---|---|---|
| `trial`, `obs` | the two samples | dict with `"x"` (rows by covariates), `"a"` (0 or 1), `"y"` (outcome); the same covariates in the same order; no missing values; at least 30 rows each |
| `trial_propensity` | the known probability of treatment in the trial | a number strictly between 0 and 1 |
| `covariate_shift` | `True` (default): the OBS covariates may follow another law; the density ratio is estimated and calibrated. `False`: the OBS covariates follow the trial law (r_0 = 1) | bool |
| `learner` (CATE) | `"DR"` (DRF) or `"R"` (RF) | string |
| `seed`, `shuffle` | the random split into blocks; `shuffle=False` keeps the row order (used by the tests) | int, bool |

## 4. Outputs

`fuse_ate` returns (schema: `schemas/ate_result.schema.json`):

| Field | Meaning |
|---|---|
| `estimate`, `se`, `ci`, `level` | the cross-fitted AIPWF estimate, its standard error and Wald interval |
| `trial_only` | the cross-fitted trial-only AIPW on the same rotations: `estimate`, `se`, `ci` |
| `variance_ratio_to_trial_only` | (se / trial-only se)^2; below 1 means the OBS helped |
| `coefficients` | lambda and omega of each of the three rotations, and the omega rule |
| `diagnostics` | `effective_obs_fraction` (share of the OBS that carries the trial law), `g_variance_share` |
| `warnings` | list of strings; read them before using the estimate |

`fuse_cate(...).summary()` returns the learner, the sieve, the chosen (lambda, omega) of each rotation, the same
diagnostics and warnings (schema: `schemas/cate_result.schema.json`).

## 5. Assumptions

1. **Trial identification (Assumption 1).** Y = Y(A); the trial treatment is randomized given X; its propensity
   is known (`trial_propensity`) and bounded away from 0 and 1.
2. **Sample separation (Assumption 2).** The rows are independent draws. The code makes the nuisance, tuning and
   evaluation samples by splitting, and cross-fits by rotating them.
3. **Covariate domains (Assumption 3).** Either (i) the OBS covariates follow the trial law, r_0 = 1
   (`covariate_shift=False`), or (ii) the trial covariate law is absolutely continuous with respect to the OBS law
   (every covariate region of the trial occurs in the OBS), and the ratio is estimated (`covariate_shift=True`).

Not assumed: the OBS may be confounded (unmeasured confounding of any form). The target stays the trial
population.

## 6. When it helps, and when not

Evidence: Sec. 5 and Appendix C of the paper, 1,000 replications per cell, synthetic, IHDP and ACIC 2016 designs,
trial sizes n = 300 to 3,000, OBS N = 15,000, cross-fitted.

| Setting | Result relative to the cross-fitted trial-only estimator |
|---|---|
| ATE, n = 300 | MSE 0.47 to 0.60 |
| ATE, n = 3,000 | MSE 0.79 to 1.00 (never above 1 in the 45 cells) |
| CATE, every cell | risk 0.53 to 0.67 (lowest of all compared learners in 45 of 45 cells) |
| Coverage of the 95% interval | 0.936 to 0.960 |

* Helps most with a small trial and a large OBS whose outcome regression resembles the trial's.
* Gives about the trial-only answer with a large trial (ATE), with weak overlap (see `effective_obs_fraction`),
  or with an OBS unrelated to the trial. In two designs with an uninformative OBS, the CATE risk was up to 8%
  above the trial-only learner.
* If the OBS is known to be unconfounded and to share the trial covariate law, pooling all rows is more
  efficient; this method does not use that assumption.

## 7. Errors and warnings

* `ValueError`: malformed input (missing keys, non-binary treatment, different covariates, too few rows, an arm
  with fewer than two units in a block, `trial_propensity` outside (0, 1)).
* `RuntimeError`: the density ratio could not be calibrated (the OBS does not cover the trial covariates). In the
  paper's misspecified-shift design this happened in 2 of 10,000 replications.
* Warnings: treated share of the trial far from `trial_propensity` (Assumption 1 at risk); small trial; OBS
  smaller than the trial; weak overlap; OBS not used (lambda = omega = 0 in every rotation).

## 8. Fixed choices (and where they come from)

| Part | Choice |
|---|---|
| Cross-fitting | 3 trial blocks rotated through (nuisance, tuning, evaluation); OBS split once 60/20/20 |
| Outcome regressions mu_R, mu_O | ridge (penalty 5) on x, x^2, sin x and the pairwise products of the first six covariates, by arm, predictions clipped to the fitted range (`method.fit_outcome_regression`) |
| Density ratio | logistic classifier on five principal components of the trial nuisance block, calibrated to balance g (ATE) or the sieve products (CATE) (`method.balancing_ratio`) |
| ATE coefficients | lambda = -C/A, omega = D/(B + B_cal), projected onto [0, 1] (Eq. (7); B_cal = 0 when r_0 = 1) |
| CATE sieve | (1, five principal components, g), ridge 0.01, grid (lambda, omega) in {0, 0.5, 1}^2 chosen by the validation score |

## 9. Paper to code

| Paper | Code |
|---|---|
| Lemma 1, Eq. (1): AIPW score | `dfppi.method.aipw_score` |
| Def. 1, Eq. (3): AIPWF estimate | `dfppi.aipwf.aipwf` (line 3) |
| Eq. (7): coefficients (lambda, omega) | `dfppi.aipwf.aipwf` (line 2) |
| Thm. 2: standard error (V + V_cal) | `dfppi.aipwf.aipwf` (line 4); cross-fitted: `dfppi.crossfit.crossfit_ate` |
| Def. 2, Eq. (6): calibrated ratio | `dfppi.method.balancing_ratio` |
| Prop. 6: what each ratio balances | `dfppi.aipwf.density_ratio`, `dfppi.drf.density_ratio` |
| Def. 3: DRF and RF pseudo-outcomes | `dfppi.drf.pseudo_outcome`, `dfppi.rf.pseudo_outcome` |
| Def. 5: sieve | `dfppi.drf.sieve` |
| Prop. 4, Eq. (12): normal equation | `dfppi.drf.fit_candidate` |
| Def. 4, Eq. (10): validation score | `dfppi.drf.validation_score` |
