"""Cross-fitted Algorithms 1 (ATE) and 2 (CATE) on user data: the interface for scripts and AI agents.

    import dfppi
    out = dfppi.fuse_ate(trial, obs, trial_propensity=0.5)       # a JSON-serialisable dict
    fit = dfppi.fuse_cate(trial, obs, trial_propensity=0.5)      # learner="DR" (DRF) or "R" (RF)
    tau = fit.predict(x)                                         # the CATE at the rows of x
    fit.summary()                                                # a JSON-serialisable dict

A sample is a dict with "x" (rows by covariates), "a" (treatment, 0 or 1) and "y" (outcome).  The trial is
randomised with the known probability trial_propensity (Assumption 1).  covariate_shift=True estimates the density
ratio of the trial covariates to the OBS covariates and calibrates it (Assumption 3(ii)); covariate_shift=False uses
r_0 = 1, which is right only when the OBS covariates have the trial law (Assumption 3(i)).

Both functions cross-fit as Sec. 5 does: the trial is cut into three blocks that rotate through the roles
(nuisance, tuning, evaluation), the OBS is cut once into 60/20/20 nuisance, tuning and evaluation rows, the ATE is
the mean of the three estimates and the CATE the mean of the three fitted functions.
"""
import numpy as np
from scipy.stats import norm

from . import aipwf as AIPWF
from . import drf as DRF
from . import rf as RF
from .crossfit import ROTATIONS, crossfit_ate
from .method import fit_nuisances, fit_outcome_regression, fused_score, remember, var

__all__ = ["fuse_ate", "fuse_cate", "CATEFit"]
LEARNERS = {"DR": DRF, "R": RF}
CLIP_QUANTILE = 0.005                 # as in fit_nuisances: no regression predicts beyond its own fitted range


# ---------------------------------------------------------------------------- inputs
def _sample(data, name):
    """The arrays of a sample, checked."""
    try:
        x, a, y = (np.asarray(data[key], dtype=float) for key in ("x", "a", "y"))
    except KeyError as missing:
        raise ValueError(f"{name}: a sample needs the keys 'x', 'a' and 'y' (missing {missing})") from None
    x = x[:, None] if x.ndim == 1 else x
    a, y = a.ravel(), y.ravel()
    if x.ndim != 2 or not len(x) == len(a) == len(y):
        raise ValueError(f"{name}: 'x' must have one row per entry of 'a' and 'y'")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError(f"{name}: 'x' and 'y' must be finite (no missing values)")
    if not np.all((a == 0) | (a == 1)):
        raise ValueError(f"{name}: 'a' must be 0 or 1")
    return {"x": x, "a": a, "y": y}


def _split(trial, obs, trial_propensity, seed, shuffle):
    """Three trial blocks and the OBS roles; the checks that every estimator needs."""
    if not 0.0 < trial_propensity < 1.0:
        raise ValueError("trial_propensity must lie strictly between 0 and 1")
    if trial["x"].shape[1] != obs["x"].shape[1]:
        raise ValueError("the trial and the OBS must have the same covariates, in the same order")
    rng = np.random.default_rng(seed)
    n, N = len(trial["y"]), len(obs["y"])
    if n < 30 or N < 30:
        raise ValueError(f"too few rows: {n} trial and {N} OBS rows; each needs at least 30")
    order_R = rng.permutation(n) if shuffle else np.arange(n)
    order_O = rng.permutation(N) if shuffle else np.arange(N)
    blocks = [{k: v[idx] for k, v in trial.items()} for idx in np.array_split(order_R, 3)]
    for b, block in enumerate(blocks):
        if min(block["a"].sum(), (1 - block["a"]).sum()) < 2:
            raise ValueError(f"trial block {b + 1} has fewer than two units in one arm; give more trial rows")
    nuis, tune = round(0.6 * N), round(0.2 * N)
    cuts = {"nuis": order_O[:nuis], "tune": order_O[nuis:nuis + tune], "eval": order_O[nuis + tune:]}
    O = {role: {k: v[idx] for k, v in obs.items()} for role, idx in cuts.items()}
    return blocks, O


def _common_diagnostics(trial, obs, trial_propensity):
    n, N = len(trial["y"]), len(obs["y"])
    share, warnings = float(trial["a"].mean()), []
    z = (share - trial_propensity) / np.sqrt(trial_propensity * (1 - trial_propensity) / n)
    if abs(z) > 4:
        warnings.append(f"the treated share of the trial ({share:.3f}) is {abs(z):.1f} standard errors from "
                        f"trial_propensity ({trial_propensity}); Assumption 1 needs the true randomisation probability")
    if n < 150:
        warnings.append(f"small trial ({n} rows): each block has about {n // 3} rows, so the nuisance fits are rough")
    if N < n:
        warnings.append(f"the OBS ({N} rows) is smaller than the trial ({n} rows); fusion gains need a larger OBS")
    return {"n_trial": n, "n_obs": N, "treated_share_trial": share}, warnings


def _ratio_fraction(r_values):
    """Effective sample fraction (mean r)^2 / mean(r^2) of a ratio on the OBS rows: 1 without a shift, small when
    few OBS rows carry the trial law."""
    return float(np.mean(r_values) ** 2 / np.mean(r_values ** 2))


def _fit_rotation(blocks, roles, O, trial_propensity, covariate_shift, mu_O):
    R_nuis, R_tune, R_eval = (blocks[b] for b in roles)
    try:
        nuis = fit_nuisances(R_nuis, O["nuis"], trial_propensity, same_domain=not covariate_shift, mu_O=mu_O)
        return nuis, R_tune, R_eval
    except RuntimeError as failure:
        raise RuntimeError(f"{failure}.  The density ratio could not be calibrated: the OBS may not cover the trial "
                           "covariates.  Check the overlap, or give more OBS rows.") from None


# ---------------------------------------------------------------------------- ATE
def fuse_ate(trial, obs, trial_propensity, covariate_shift=True, seed=0, shuffle=True, level=0.95):
    """Algorithm 1 (AIPWF), cross-fitted, with its Wald interval, and the cross-fitted trial-only AIPW on the same
    rotations for comparison.  Returns a JSON-serialisable dict (schemas/ate_result.schema.json)."""
    from . import __version__
    trial, obs = _sample(trial, "trial"), _sample(obs, "obs")
    blocks, O = _split(trial, obs, trial_propensity, seed, shuffle)
    sizes, warnings = _common_diagnostics(trial, obs, trial_propensity)
    mu_O = remember(fit_outcome_regression(O["nuis"], clip_quantile=CLIP_QUANTILE))
    fused, trial_only, fractions, shares = [], [], [], []
    for roles in ROTATIONS:
        nuis, R_tune, R_eval = _fit_rotation(blocks, roles, O, trial_propensity, covariate_shift, mu_O)
        r = AIPWF.density_ratio(nuis)
        r_On = r(O["nuis"]["x"])
        fused.append((AIPWF.aipwf(nuis, r, R_tune, O["tune"], R_eval, O["eval"]), r_On))
        trial_only.append((AIPWF.aipwf(nuis, r, R_tune, O["tune"], R_eval, O["eval"], False, False), r_On))
        fractions.append(_ratio_fraction(r_On))
        shares.append(var(nuis.g(R_tune["x"])) / var(fused_score(R_tune, nuis, 0.0)))
    g = lambda x: mu_O(x)[1] - mu_O(x)[0]
    g_blocks, g_On = [g(block["x"]) for block in blocks], g(O["nuis"]["x"])
    estimate, se = crossfit_ate(fused, g_blocks, g_On, calibrated_ratio=covariate_shift)
    estimate_0, se_0 = crossfit_ate(trial_only, g_blocks, g_On, calibrated_ratio=False)
    z = float(norm.ppf(0.5 + level / 2))
    lam, omega = [out["lambda"] for out, _ in fused], [out["omega"] for out, _ in fused]
    if covariate_shift and min(fractions) < 0.1:
        warnings.append(f"weak overlap: the calibrated ratio keeps an effective {min(fractions):.1%} of the OBS rows")
    if max(lam) == 0 and max(omega) == 0:
        warnings.append("the OBS were not used (lambda = omega = 0 in every rotation): the estimate is the "
                        "cross-fitted trial-only AIPW")
    return {
        "estimand": "ATE", "method": "AIPWF (Algorithm 1), 3-fold cross-fitted", "version": __version__,
        "estimate": estimate, "se": se, "ci": [estimate - z * se, estimate + z * se], "level": level,
        "trial_only": {"estimate": estimate_0, "se": se_0, "ci": [estimate_0 - z * se_0, estimate_0 + z * se_0]},
        "variance_ratio_to_trial_only": (se / se_0) ** 2,
        "coefficients": {"lambda": lam, "omega": omega,
                         "omega_rule": "D/(B+B_cal)" if covariate_shift else "D/B (r_0 = 1, so B_cal = 0)"},
        "covariate_shift": bool(covariate_shift),
        "density_ratio": "estimated and calibrated to balance g" if covariate_shift else "known: r_0 = 1",
        "sizes": {**sizes, "trial_blocks": [len(b["y"]) for b in blocks],
                  "obs_nuis_tune_eval": [len(O[k]["y"]) for k in ("nuis", "tune", "eval")]},
        "diagnostics": {"effective_obs_fraction": float(np.mean(fractions)),
                        "g_variance_share": float(np.mean(shares))},
        "warnings": warnings, "seed": seed, "shuffle": bool(shuffle),
    }


# ---------------------------------------------------------------------------- CATE
class CATEFit:
    """The cross-fitted CATE of Algorithm 2: predict(x) is the mean of the three fitted sieve functions.
    predict_trial_only(x) is the trial-only learner of the same type on the same rotations (no OBS input)."""

    def __init__(self, fits, summary):
        self._fits, self._summary = fits, summary

    def predict(self, x):
        x = np.atleast_2d(np.asarray(x, dtype=float))
        return np.mean([DRF.sieve(nuis, x) @ beta for nuis, beta, _ in self._fits], axis=0)

    def predict_trial_only(self, x):
        x = np.atleast_2d(np.asarray(x, dtype=float))
        return np.mean([nuis.components(x) @ beta_0 for nuis, _, beta_0 in self._fits], axis=0)

    def summary(self):
        return dict(self._summary)


def fuse_cate(trial, obs, trial_propensity, learner="DR", covariate_shift=True, seed=0, shuffle=True):
    """Algorithm 2 with the pseudo-outcome of the DR-learner (learner="DR", DRF) or of the R-learner ("R", RF),
    cross-fitted.  Returns a CATEFit."""
    from . import __version__
    if learner not in LEARNERS:
        raise ValueError(f"learner must be one of {sorted(LEARNERS)}")
    estimator = LEARNERS[learner]
    trial, obs = _sample(trial, "trial"), _sample(obs, "obs")
    blocks, O = _split(trial, obs, trial_propensity, seed, shuffle)
    sizes, warnings = _common_diagnostics(trial, obs, trial_propensity)
    mu_O = remember(fit_outcome_regression(O["nuis"], clip_quantile=CLIP_QUANTILE))
    fits, chosen, fractions = [], [], []
    for roles in ROTATIONS:
        nuis, R_tune, R_eval = _fit_rotation(blocks, roles, O, trial_propensity, covariate_shift, mu_O)
        r = DRF.density_ratio(nuis)
        lam, omega, beta = estimator.fit(nuis, r, R_tune, O["tune"], R_eval, O["eval"])
        fits.append((nuis, beta, estimator.rct_only(nuis, R_tune)))
        chosen.append({"lambda": lam, "omega": omega})
        fractions.append(_ratio_fraction(r(O["nuis"]["x"])))
    if covariate_shift and min(fractions) < 0.1:
        warnings.append(f"weak overlap: the calibrated ratio keeps an effective {min(fractions):.1%} of the OBS rows")
    if all(c["lambda"] == 0 and c["omega"] == 0 for c in chosen):
        warnings.append("grid point (0, 0) in every rotation: the OBS enter only through g in the sieve")
    summary = {
        "estimand": "CATE", "method": f"{'DRF' if learner == 'DR' else 'RF'} (Algorithm 2), 3-fold cross-fitted",
        "version": __version__, "learner": learner,
        "sieve": "(1, five principal components of the trial nuisance block, g = mu_O(., 1) - mu_O(., 0)); "
                 "the fitted CATE is the mean of the three rotations' sieve functions",
        "grid": [list(c) for c in DRF.GRID], "ridge": 0.01, "chosen": chosen,
        "covariate_shift": bool(covariate_shift),
        "density_ratio": "estimated and calibrated to balance the sieve products" if covariate_shift
                         else "known: r_0 = 1",
        "sizes": {**sizes, "trial_blocks": [len(b["y"]) for b in blocks],
                  "obs_nuis_tune_eval": [len(O[k]["y"]) for k in ("nuis", "tune", "eval")]},
        "diagnostics": {"effective_obs_fraction": float(np.mean(fractions))},
        "warnings": warnings, "seed": seed, "shuffle": bool(shuffle),
    }
    return CATEFit(fits, summary)
