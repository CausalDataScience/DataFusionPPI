"""Manuscript-faithful DRF and RF fits for the minimal plan.

Both losses are implemented once, from the manuscript, and every rule and route
calls the same functions.  Writing zeta = b' beta and adding rho ||beta||^2:

    L_DRF = P_R (Z_lambda - zeta)^2
            + omega { P_O[r (ghat - zeta)^2] - P_R (ghat - zeta)^2 }

    L_RF  = P_R [Y - m_lambda - (A - e) zeta]^2 / {e (1 - e)}
            + omega { P_O[r (ghat - zeta)^2] - P_R[chi (ghat - zeta)^2] }

Differentiating gives the two normal equations in ``solve``.  With
chi = (A - e)^2 / {e (1 - e)}, the R-Fusion right-hand side first term is
P_R[chi b Ztilde_RF] where Ztilde_RF = (Y - m_lambda) / (A - e), which is the
manuscript's w_j b Ztilde_{j,lambda} at w_RF = chi.  At e = 1/2 both chi and the
two pseudo-outcomes coincide, so the two losses agree observation by observation.

Selection uses the common-reference score for both learners:

    C(lambda, omega) = P_R^sel[(Z_0 - zeta)^2 - omega G] + P_O^sel[omega r G]

with G = (ghat - zeta)^2.  The trial correction carries no chi here; that is a
property of the score, not an oversight, and a mutation test checks it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from fusion_core import fit_outcome, spline_basis_factory

RIDGE = 1e-2
GRID = tuple((lam, om) for lam in (0.0, 0.5, 1.0) for om in (0.0, 0.5, 1.0))
SPLINE_KNOTS = 3


class NumericalFailure(RuntimeError):
    """A solve that is not finite or does not satisfy its own normal equation."""


# ------------------------------------------------------------------ nuisance
@dataclass
class Nuisance:
    mu_r: object
    mu_o: object
    basis: object
    trial_propensity: float


def fit_nuisance(rct_nuis: dict, obs_nuis: dict, trial_propensity: float) -> Nuisance:
    """Identical outcome regressions on both sources; basis from the trial only."""
    return Nuisance(mu_r=fit_outcome(rct_nuis, spec="flexible", alpha=5.0),
                    mu_o=fit_outcome(obs_nuis, spec="flexible", alpha=5.0),
                    basis=spline_basis_factory(rct_nuis["x"], SPLINE_KNOTS),
                    trial_propensity=trial_propensity)


def aipw(data: dict, q0: np.ndarray, q1: np.ndarray, e: np.ndarray) -> np.ndarray:
    a, y = data["a"], data["y"]
    return q1 - q0 + a * (y - q1) / e - (1 - a) * (y - q0) / (1 - e)


@dataclass
class Parts:
    """Everything a candidate needs from one trial sample."""
    b: np.ndarray
    z0: np.ndarray
    z1: np.ndarray
    ghat: np.ndarray
    m_r: np.ndarray
    m_o: np.ndarray
    chi: np.ndarray
    centred_arm: np.ndarray
    y: np.ndarray
    weight_scale: float        # 1 / {e (1 - e)}


def trial_parts(data: dict, nuis: Nuisance) -> Parts:
    e = np.full(len(data["x"]), nuis.trial_propensity)
    r0, r1 = nuis.mu_r.predict(data["x"])
    o0, o1 = nuis.mu_o.predict(data["x"])
    return Parts(b=nuis.basis(data["x"]),
                 z0=aipw(data, r0, r1, e), z1=aipw(data, o0, o1, e),
                 ghat=o1 - o0, m_r=e * r1 + (1 - e) * r0, m_o=e * o1 + (1 - e) * o0,
                 chi=(data["a"] - e) ** 2 / (e * (1 - e)),
                 centred_arm=data["a"] - e, y=data["y"],
                 weight_scale=1.0 / float(nuis.trial_propensity
                                          * (1.0 - nuis.trial_propensity)))


def obs_parts(data: dict, nuis: Nuisance) -> tuple[np.ndarray, np.ndarray]:
    o0, o1 = nuis.mu_o.predict(data["x"])
    return nuis.basis(data["x"]), o1 - o0


# ------------------------------------------------------------- normal equations
def solve(learner: str, tune: Parts, b_obs: np.ndarray, ghat_obs: np.ndarray,
          ratio_obs: np.ndarray, lam: float, omega: float,
          ridge: float = RIDGE) -> np.ndarray:
    """One candidate, from the derivative of the loss it belongs to."""
    gram_obs = (b_obs * ratio_obs[:, None]).T @ b_obs / len(b_obs)
    rhs_obs = (b_obs * (ratio_obs * ghat_obs)[:, None]).mean(0)
    if learner == "DRF":
        weight = np.ones(len(tune.b))
        z_lambda = (1.0 - lam) * tune.z0 + lam * tune.z1
        rhs_trial = (tune.b * z_lambda[:, None]).mean(0)
    elif learner == "RF":
        weight = tune.chi
        m_lambda = (1.0 - lam) * tune.m_r + lam * tune.m_o
        # chi * Ztilde_RF simplifies to (A - e)(Y - m_lambda) / {e (1 - e)}, which
        # never divides by A - e and so never needs a guard
        scaled = tune.weight_scale * tune.centred_arm * (tune.y - m_lambda)
        rhs_trial = (tune.b * scaled[:, None]).mean(0)
    else:
        raise ValueError(learner)

    gram_trial = (tune.b * weight[:, None]).T @ tune.b / len(tune.b)
    trial_ghat = (tune.b * (weight * tune.ghat)[:, None]).mean(0)
    matrix = gram_trial + omega * (gram_obs - gram_trial) + ridge * np.eye(tune.b.shape[1])
    rhs = rhs_trial + omega * (rhs_obs - trial_ghat)
    beta = np.linalg.solve(matrix, rhs)
    residual = matrix @ beta - rhs
    scale = max(float(np.max(np.abs(rhs))), 1e-12)
    if not np.all(np.isfinite(beta)) or float(np.max(np.abs(residual))) / scale > 1e-10:
        raise NumericalFailure(
            f"{learner} candidate ({lam}, {omega}) left a relative normal-equation "
            f"residual of {float(np.max(np.abs(residual))) / scale:.2e}")
    return beta


def selection_score(beta: np.ndarray, sel: Parts, b_obs_sel: np.ndarray,
                    ghat_obs_sel: np.ndarray, ratio_obs_sel: np.ndarray,
                    omega: float) -> float:
    """The common-reference score of Algorithm 1, used by both learners."""
    zeta_trial = sel.b @ beta
    g_trial = (sel.ghat - zeta_trial) ** 2
    zeta_obs = b_obs_sel @ beta
    g_obs = (ghat_obs_sel - zeta_obs) ** 2
    trial = float(np.mean((sel.z0 - zeta_trial) ** 2 - omega * g_trial))
    obs = float(np.mean(omega * ratio_obs_sel * g_obs))
    return trial + obs


# ------------------------------------------------------------------- ratios
def classifier_ratio(rct_x: np.ndarray, obs_x: np.ndarray):
    """Prior-corrected classifier odds, computed on the log-odds scale.

    The fitted probability is never formed, so there is no probability clipping
    and no silent truncation of a large or small ratio.
    """
    x = np.vstack([rct_x, obs_x])
    label = np.r_[np.ones(len(rct_x)), np.zeros(len(obs_x))]
    model = Pipeline([("scale", StandardScaler()),
                      ("logit", LogisticRegression(C=1.0, max_iter=5000))])
    model.fit(x, label)
    prior = len(rct_x) / len(x)
    correction = np.log((1.0 - prior) / prior)

    def ratio(z: np.ndarray) -> np.ndarray:
        log_odds = model.decision_function(z)
        values = np.exp(log_odds + correction)
        if not np.all(np.isfinite(values)):
            raise NumericalFailure("the classifier ratio overflowed on the log-odds scale")
        return values

    converged = bool(np.all(model.named_steps["logit"].n_iter_
                            < model.named_steps["logit"].max_iter))
    return ratio, converged


def ratio_diagnostics(values: np.ndarray) -> dict:
    return {"mean": float(values.mean()),
            "second_moment": float((values ** 2).mean()),
            "ess_fraction": float(values.mean() ** 2 / (values ** 2).mean()),
            "max": float(values.max())}


# --------------------------------------------------------------- rule reading
RULES = {
    "rct_only": lambda lam, om: lam == 0.0 and om == 0.0,
    "lambda_only": lambda lam, om: om == 0.0,
    "omega_only": lambda lam, om: lam == 0.0,
    "joint": lambda lam, om: True,
}


def pick(scores: dict[tuple[float, float], float], rule: str) -> tuple[float, float]:
    """Minimise the selection score inside the rule's subset.

    Ties break on (lambda, omega) in lexicographic order, so the choice does not
    depend on dictionary order or on floating-point noise between equal scores.
    """
    allowed = [key for key in scores if RULES[rule](*key)]
    return min(allowed, key=lambda key: (scores[key], key[0], key[1]))
