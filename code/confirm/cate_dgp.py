"""Data-generating processes for the confirmatory CATE stages.

Audit plan section 7.1 fixes the bounded core:

    X ~ Unif([-1, 1]^d)
    Y(a) = m(X) + a tau(X) + eta(X, U) + sigma(X) eps
    logit P_O(A = 1 | X, U) = logit e_R(X) + c0 U + c1 h(X) U

The trial propensity is one half throughout, which is the cleanest setting for
the mechanism check because Corollary 7.1's R-Fusion cross term vanishes there.

Two knobs move the two channels, for the reasons M1 established.

``confounding`` (c0, and c1 for the shape) decides how well an observational
outcome regression predicts the trial outcome surface, and so decides the sign
and size of C_p.  ``roughness`` puts part of tau outside the spline span, which
is the only way D_p can be non-zero: when the basis spans tau, zeta_p equals tau
and D_p is exactly zero.

The Gaussian shift family moves the trial covariate mean so the transport ratio
is not one.  The mean separation is set from the target second moment, so
E_O(r0^2) = exp(d2) and the effective sample size fraction is exp(-d2) with
d2 = log 2 for the mild level and log 5 for the strong one.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

SMOOTH_FREQ = 1.0
ROUGH_FREQ = 6.0
SIGMA = 0.5
SHIFT_D2 = {"none": 0.0, "mild": math.log(2.0), "strong": math.log(5.0)}


@dataclass(frozen=True)
class SCMBDesign:
    """Prespecified bounded theorem-check law from handoff Section 6.3."""
    confounding: float = 1.0
    B0 = 9.0
    RBAR = 1.0

    def sample(self, study: str, n: int,
               rng: np.random.Generator) -> dict[str, np.ndarray]:
        x = rng.uniform(-1.0, 1.0, size=(n, 5))
        u = rng.uniform(-1.0, 1.0, size=n)
        epsilon = rng.uniform(-1.0, 1.0, size=n)
        e_r = 0.5 + 0.2 * x[:, 0]
        if study == "RCT":
            e = e_r
        else:
            logit = np.log(e_r / (1.0 - e_r)) + self.confounding * u
            e = 1.0 / (1.0 + np.exp(-logit))
        a = (rng.random(n) < e).astype(float)
        effect = 0.15 + 0.10 * x[:, 0] - 0.05 * x[:, 1] ** 2
        baseline_value = 0.15 * np.sin(np.pi * x[:, 0]) + 0.05 * x[:, 1] * x[:, 2]
        y = baseline_value + 0.10 * u + 0.05 * epsilon + a * effect
        return {"x": x, "a": a, "y": y, "tau": effect,
                "propensity": e if study == "RCT" else None,
                "u": u, "epsilon": epsilon}

    def target_sample(self, n: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
        x = rng.uniform(-1.0, 1.0, size=(n, 5))
        return {"x": x, "tau": 0.15 + 0.10 * x[:, 0] - 0.05 * x[:, 1] ** 2}

    def oracle_ratio(self, x: np.ndarray) -> np.ndarray:
        return np.ones(len(x))


def _direction(dimension: int) -> np.ndarray:
    """A fixed unit direction, so the design carries no hidden randomness."""
    raw = np.array([1.0, -0.6, 0.35, -0.2, 0.15, -0.1])[:dimension]
    return raw / np.linalg.norm(raw)


def index(x: np.ndarray) -> np.ndarray:
    """The single index every structural function is built on."""
    return x @ _direction(x.shape[1])


def tau(x: np.ndarray, roughness: float) -> np.ndarray:
    """A smooth part a cubic spline can follow, plus a rough part it cannot."""
    z = index(x)
    smooth = 0.30 + 0.45 * z - 0.25 * z ** 2 + 0.20 * z ** 3
    return smooth + roughness * np.cos(ROUGH_FREQ * z)


def baseline(x: np.ndarray) -> np.ndarray:
    return 0.40 * index(x) + 0.20 * np.sin(SMOOTH_FREQ * index(x))


def shape(x: np.ndarray) -> np.ndarray:
    """h(X) in the propensity, the term that makes the confounding X-dependent."""
    return np.tanh(1.5 * index(x))


@dataclass(frozen=True)
class Design:
    family: str
    dimension: int = 4
    confounding: float = 0.0
    shape_confounding: float = 0.0
    roughness: float = 0.0
    shift: str = "none"
    obs_baseline_shift: float = 0.0

    # ------------------------------------------------------------- covariates
    @property
    def mean_separation(self) -> float:
        """Trial mean offset along the index, from the target E_O(r0^2)."""
        return math.sqrt(SHIFT_D2[self.shift])

    def draw_covariates(self, study: str, n: int, rng: np.random.Generator) -> np.ndarray:
        if self.family == "gaussian_shift":
            x = rng.normal(size=(n, self.dimension))
            if study == "RCT" and self.mean_separation:
                x = x + self.mean_separation * _direction(self.dimension)
            return x
        return rng.uniform(-1.0, 1.0, size=(n, self.dimension))

    def oracle_ratio(self, x: np.ndarray) -> np.ndarray:
        """r0 = dP_R^X / dP_O^X.  Exactly one when the two laws agree."""
        if self.family != "gaussian_shift" or not self.mean_separation:
            return np.ones(len(x))
        z = index(x)
        d = self.mean_separation
        return np.exp(d * z - 0.5 * d ** 2)

    # ---------------------------------------------------------------- outcome
    def propensity(self, study: str, x: np.ndarray, u: np.ndarray) -> np.ndarray:
        if study == "RCT":
            return np.full(len(x), 0.5)
        logit = self.confounding * u + self.shape_confounding * shape(x) * u
        return 1.0 / (1.0 + np.exp(-logit))

    def sample(self, study: str, n: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
        x = self.draw_covariates(study, n, rng)
        u = rng.normal(size=n)
        e = self.propensity(study, x, u)
        a = (rng.uniform(size=n) < e).astype(float)
        effect = tau(x, self.roughness)
        # the confounder enters the outcome additively; the propensity is what
        # makes its influence depend on X
        # the observational baseline may differ from the trial baseline; the
        # conditional effect is shared, so the estimand is unchanged, but an
        # observational outcome regression then mispredicts the trial surface and
        # the increment it offers stops being worth taking
        offset = self.obs_baseline_shift if study == "OBS" else 0.0
        y = (baseline(x) + offset * np.sin(2.5 * index(x))
             + a * effect + 0.60 * u + SIGMA * rng.normal(size=n))
        out = {"x": x, "a": a, "y": y, "tau": effect}
        # the trial propensity is known by design; the observational one is not,
        # and the engine must estimate it from (X, A) alone
        out["propensity"] = e if study == "RCT" else None
        return out

    def target_sample(self, n: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
        """Trial-population covariates with the true conditional effect."""
        x = self.draw_covariates("RCT", n, rng)
        return {"x": x, "tau": tau(x, self.roughness)}

    def ess_fraction(self) -> float:
        return math.exp(-SHIFT_D2[self.shift])


class StarNotImplemented(NotImplementedError):
    """The STAR source law is a separate blocker (plan P0-4).

    Raised rather than silently falling through to the synthetic sampler, which
    would emit uniform covariates under a STAR label.
    """


def from_cell(cell) -> Design:
    if cell.family == "star":
        raise StarNotImplemented(
            "the frozen STAR source law is not implemented yet; plan P0-4 requires "
            "the exact learner tuple, five-fold deterministic out-of-fold pooling, "
            "exact-cell propensity inside [0.15, 0.85] and a law frozen before any "
            "replication, with source, fold, probability and propensity hashes")
    return Design(family=cell.family, dimension=cell.dimension,
                  confounding=cell.confounding,
                  shape_confounding=cell.shape_confounding,
                  roughness=cell.roughness, shift=cell.shift,
                  obs_baseline_shift=cell.obs_baseline_shift)
