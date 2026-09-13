"""Deterministic mechanism fixtures for the four quadrants of plan section 7.2.

The audit's root cause 5.2 is that the two channels' opportunities were never
moved independently.  These fixtures move them with two knobs that are
algebraically independent, which the module proves rather than asserts.

Construction.  Take a fixed grid of covariate values.  At each grid point place
four rows carrying a fixed two-column pattern with column means zero, unit
column variances and a prescribed correlation kappa.  Scale the columns by
sigma_z and sigma_delta and shift the first by tau(x).  The result is an array
with, exactly and with no sampling error,

    E(Ztilde_0 | x) = tau(x),      E(Delta | x) = 0,
    Cov(Ztilde_0, Delta | x) = sigma_z(x) sigma_delta(x) kappa.

Writing zeta_p for the basis projection of tau, the two moments that decide the
quadrant then separate:

    C_p = tr[Gamma^-1 E{b b' sigma_z sigma_delta kappa}]
    D_p = tr[Gamma^-1 E{b b' (tau - zeta_p)(ghat - zeta_p)}]

C_p depends on kappa and not on how well the basis spans tau.  D_p depends on
the part of tau the basis misses and not on kappa.  So kappa is the lambda knob
and the span of the basis is the omega knob.

A scoped consequence worth stating: for the DRF learner, an exact score, a fixed
basis, exact transport, and population first-order moments, when the basis spans
tau exactly, zeta_p = tau and D_p is exactly zero.  Thus the omega channel has no
first-order opportunity under those conditions.  This fixture does not claim
that conclusion for RF, estimated scores, learned bases, or finite samples.
"""
from __future__ import annotations

import numpy as np

from confirm.moments import Fixture

# A fixed four-row pattern with zero mean and unit variance in each column and
# zero correlation between them.  Correlation is introduced afterwards, so the
# pattern itself carries no design choice.
PATTERN_P = np.array([1.0, 1.0, -1.0, -1.0])
PATTERN_Q = np.array([1.0, -1.0, 1.0, -1.0])
GRID_POINTS = 60
OBS_GRID_POINTS = 90


def _grid(n: int) -> np.ndarray:
    """Midpoints of n equal cells on [-1, 1]; deterministic and symmetric."""
    return (np.arange(n) + 0.5) / n * 2.0 - 1.0


def _basis(x: np.ndarray, degree: int = 3) -> np.ndarray:
    """Polynomial basis with an intercept, columns 1, x, ..., x^degree."""
    return np.vander(x, degree + 1, increasing=True)


def _tau(x: np.ndarray, *, off_span: float) -> np.ndarray:
    """A cubic that the basis spans, plus an optional term that it does not.

    The extra term is a cosine at a frequency the cubic cannot follow, so
    ``off_span = 0`` puts tau inside the span and any positive value puts a
    controlled amount outside it.
    """
    inside = 0.30 + 0.45 * x - 0.25 * x ** 2 + 0.20 * x ** 3
    return inside + off_span * np.cos(6.0 * x)


def _ghat(x: np.ndarray, *, tracks: float) -> np.ndarray:
    """The observational prediction.

    ``tracks`` scales how much of tau the prediction carries.  A prediction that
    equals the projection contributes nothing to D_p, and one that follows the
    off-span part of tau contributes to it.
    """
    return _tau(x, off_span=1.0) * tracks + 0.10 * (1.0 - tracks) * np.sin(2.0 * x)


def build(name: str, *, kappa: float, off_span: float, tracks: float,
          sigma_z: float = 0.8, sigma_delta: float = 0.5,
          n_rct_tune: int = 120, n_obs_tune: int = 900,
          degree: int = 3, ratio_scale: float = 1.0) -> Fixture:
    """One quadrant fixture.  Nothing here is random at call time."""
    grid = _grid(GRID_POINTS)
    b_target = _basis(grid, degree)
    tau_target = _tau(grid, off_span=off_span)

    # four rows per grid point, carrying the prescribed conditional structure
    x_rct = np.repeat(grid, len(PATTERN_P))
    tile_p = np.tile(PATTERN_P, len(grid))
    tile_q = np.tile(PATTERN_Q, len(grid))
    correlated_q = kappa * tile_p + np.sqrt(max(1.0 - kappa ** 2, 0.0)) * tile_q

    b_rct = _basis(x_rct, degree)
    z0_rct = _tau(x_rct, off_span=off_span) + sigma_z * tile_p
    delta_rct = sigma_delta * correlated_q
    ghat_rct = _ghat(x_rct, tracks=tracks)
    weight_rct = np.ones_like(x_rct)

    x_obs = _grid(OBS_GRID_POINTS)
    b_obs = _basis(x_obs, degree)
    ghat_obs = _ghat(x_obs, tracks=tracks)
    ratio_obs = np.full_like(x_obs, ratio_scale)

    return Fixture(name=name, b_target=b_target, tau_target=tau_target,
                   b_rct=b_rct, z0_rct=z0_rct, delta_rct=delta_rct,
                   ghat_rct=ghat_rct, weight_rct=weight_rct,
                   b_obs=b_obs, ghat_obs=ghat_obs, ratio_obs=ratio_obs,
                   n_rct_tune=n_rct_tune, n_obs_tune=n_obs_tune)


# kappa negative gives C_p negative and so lambda* positive, which is the
# direction the manuscript expects when the observational regression helps.
QUADRANTS = {
    "L0O0": dict(kappa=0.0, off_span=0.0, tracks=0.0),
    "L1O0": dict(kappa=-0.6, off_span=0.0, tracks=0.0),
    "L0O1": dict(kappa=0.0, off_span=0.35, tracks=1.0),
    "L1O1": dict(kappa=-0.6, off_span=0.35, tracks=1.0),
}


def quadrant(name: str) -> Fixture:
    return build(name, **QUADRANTS[name])


def rf_cross_term(propensity: float, rows_per_point: int = 10) -> Fixture:
    """An R-Fusion fixture built from the manuscript definitions (4.tex:250).

        w_RF   = chi(A, X) = {A - e(X)}^2 / [e(X){1 - e(X)}],   E_R(chi | X) = 1
        Ztilde_RF,0 = {Y - m_R(X)} / {A - e(X)}
        Delta_RF    = -{m_O(X) - m_R(X)} / {A - e(X)}

    Corollary 7.1 says the lambda-omega cross term vanishes when e is one half.
    The reason is visible once the weight and the increment are multiplied:

        chi * Delta_RF          = -(A - e) delta_m / [e(1-e)]
        chi^2 * Delta_RF        = -(A - e)^3 delta_m / [e(1-e)]^2

    and for a Bernoulli(e) arm, E{(A - e)^3 | X} = e(1-e)(1-2e), which is zero
    exactly at e = 1/2.  The cross term therefore carries a factor (1-2e)/[e(1-e)].
    The fixture places ``rows_per_point`` rows at each grid point with the arm
    split matching e exactly, so the arm law carries no sampling error.
    """
    if not 0.0 < propensity < 1.0:
        raise ValueError("the trial propensity must lie strictly inside (0, 1)")
    treated = propensity * rows_per_point
    if abs(treated - round(treated)) > 1e-12:
        raise ValueError(f"e = {propensity} needs rows_per_point * e to be a whole number")
    treated = int(round(treated))

    grid = _grid(GRID_POINTS)
    b_target = _basis(grid)
    tau_target = _tau(grid, off_span=0.35)

    arm_block = np.r_[np.ones(treated), np.zeros(rows_per_point - treated)]
    x_rct = np.repeat(grid, rows_per_point)
    arm = np.tile(arm_block, len(grid))
    centred_arm = arm - propensity
    chi = centred_arm ** 2 / (propensity * (1.0 - propensity))

    # a mean-zero pattern inside each arm at each grid point, so the conditional
    # mean of the trial pseudo-outcome is exactly tau
    def balanced(count: int) -> np.ndarray:
        if count == 1:
            return np.zeros(1)
        step = np.arange(count) - (count - 1) / 2.0
        return step / np.abs(step).max()

    noise_block = np.r_[balanced(treated), balanced(rows_per_point - treated)]
    noise = np.tile(noise_block, len(grid))

    tau_rct = _tau(x_rct, off_span=0.35)
    delta_m = 0.35 * np.cos(3.0 * x_rct)          # m_O - m_R, a function of X only
    z0_rct = tau_rct + 0.6 * noise / centred_arm
    delta_rct = -delta_m / centred_arm

    x_obs = _grid(OBS_GRID_POINTS)
    return Fixture(name=f"RF_e{propensity:.2f}", b_target=b_target,
                   tau_target=tau_target, b_rct=_basis(x_rct),
                   z0_rct=z0_rct, delta_rct=delta_rct,
                   ghat_rct=_ghat(x_rct, tracks=1.0), weight_rct=chi,
                   b_obs=_basis(x_obs), ghat_obs=_ghat(x_obs, tracks=1.0),
                   ratio_obs=np.ones_like(x_obs),
                   n_rct_tune=120, n_obs_tune=900)
