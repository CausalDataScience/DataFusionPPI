"""Generation adapter for the manuscript-minimal plan.

The structural equations, the effect function and the observational propensity
come from the existing SCM in ``ssem_ate_pilot``.  Only three things change, and
each is required by the plan.

1. The trial propensity is the known constant 0.35.  At one half the two losses
   coincide observation by observation, so 0.35 is the smallest departure that
   still separates the doubly robust and the R-Fusion learners.
2. The trial covariate law is N(mu, Sigma) in both cells.  The observational law
   is N(mu, Sigma) in the shared cell and N(mu + d, Sigma) in the shifted one,
   with d scaled so its Mahalanobis length is one.  The trial population and the
   true effect are therefore identical across cells.
3. One replication draws its latent variable, covariate noise, treatment uniforms
   and outcome noise once, and both cells are built from them.  The two cells are
   paired, and inside a cell every ratio route and both learners see one dataset.

The transport ratio r0 = dP_R^X / dP_O^X for two normals that differ only in
mean is exp{-d' Sigma^-1 (x - mu) + d' Sigma^-1 d / 2}, which the plan states and
``exact_ratio`` implements.  When the shift is zero it returns one exactly.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import ssem_ate_pilot as legacy

TRIAL_PROPENSITY = 0.35
RCT_ROLE_SIZES = {"nuisance": 100, "tuning": 100, "selection": 100}
OBS_ROLE_SIZES = {"nuisance": 3000, "tuning": 1000, "selection": 1000}
TRUTH_GRID = 50_000
SCM_ID = 1


def scm():
    return legacy.generate_scm_parameters(SCM_ID)


def shift_direction(params) -> np.ndarray:
    """d = v / sqrt(v' Sigma^-1 v), so the Mahalanobis length of d is one."""
    v = np.asarray(params.rct_mean_shift, dtype=float)
    if not np.any(v):
        v = np.array([1.0, 0.0, 0.0, 0.0])
    precision = np.linalg.inv(params.covariance)
    return v / np.sqrt(float(v @ precision @ v))


@dataclass(frozen=True)
class Cell:
    name: str                      # "shared" or "shifted"
    obs_offset: np.ndarray

    @property
    def shifted(self) -> bool:
        return bool(np.any(self.obs_offset))


def cells(params) -> tuple[Cell, Cell]:
    d = shift_direction(params)
    return (Cell("shared", np.zeros_like(d)), Cell("shifted", d))


def exact_ratio(x: np.ndarray, params, cell: Cell) -> np.ndarray:
    """r0(x) = dP_R^X / dP_O^X for two normals differing only in mean."""
    if not cell.shifted:
        return np.ones(len(x))
    precision = np.linalg.inv(params.covariance)
    d = cell.obs_offset
    centred = x - params.covariate_mean
    return np.exp(-centred @ precision @ d + 0.5 * float(d @ precision @ d))


# ------------------------------------------------------------------ sampling
@dataclass(frozen=True)
class Draw:
    """The randomness of one source and role, drawn once and reused by both cells."""
    latent: np.ndarray
    covariate_noise: np.ndarray
    treatment_uniform: np.ndarray
    outcome_noise: np.ndarray

    @staticmethod
    def make(n: int, rng: np.random.Generator) -> "Draw":
        return Draw(latent=rng.normal(size=n),
                    covariate_noise=rng.normal(size=(n, 4)),
                    treatment_uniform=rng.uniform(size=n),
                    outcome_noise=rng.normal(size=n))


def realise(draw: Draw, study: str, params, cell: Cell) -> dict[str, np.ndarray]:
    """Build one dataset from shared randomness.

    The covariates move with the cell only for the observational source.  The
    propensity, treatment and outcome are then regenerated at the moved
    covariates, so nothing is carried over from the unshifted draw.
    """
    mean = params.covariate_mean + (cell.obs_offset if study == "OBS" else 0.0)
    x = (mean + np.outer(draw.latent, params.latent_loading)
         + draw.covariate_noise * params.covariate_noise_sd)
    if study == "RCT":
        propensity = np.full(len(x), TRIAL_PROPENSITY)
    else:
        propensity = legacy.obs_propensity(x, draw.latent, params)
    treatment = (draw.treatment_uniform < propensity).astype(float)
    effect = legacy.tau_function(x, params)
    baseline = (params.baseline_intercept
                + legacy.baseline_basis(x) @ params.baseline_coefficients
                + params.outcome_latent_coefficient * draw.latent)
    outcome = baseline + treatment * effect + params.outcome_noise_sd * draw.outcome_noise
    return {"x": x, "a": treatment, "y": outcome, "propensity": propensity,
            "tau": effect}


def truth_grid(params, seed: int = 20260913) -> dict[str, np.ndarray]:
    """Fixed trial-population covariates with the true effect.

    Generated once from its own seed and reused by every cell, learner, route and
    replication, so every risk, bias and variance number is an average over the
    same points.
    """
    rng = np.random.default_rng(seed)
    latent = rng.normal(size=TRUTH_GRID)
    x = (params.covariate_mean + np.outer(latent, params.latent_loading)
         + rng.normal(size=(TRUTH_GRID, 4)) * params.covariate_noise_sd)
    return {"x": x, "tau": legacy.tau_function(x, params)}
