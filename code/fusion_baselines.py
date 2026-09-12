"""Baselines of handoff Section 4, written in the manuscript's notation.

Each is a faithful reimplementation of a published idea, not a copy of the
authors' estimator.  All share the roles, the nuisance fits, and the evaluation
samples of the fusion estimators, so the comparison isolates the combination rule.
"""
from __future__ import annotations

import numpy as np

from fusion_core import aipw_score, cov, fit_propensity, project, var


def semi_supervised(z0_tune, gr_tune, gr_obs_tune, r_obs_tune,
                    z0_eval, gr_eval, gr_obs_eval, r_obs_eval, n, n_obs) -> dict:
    """Section 4.1.  Uses the RCT-fitted contrast g_R, so no OBS outcome
    information can enter; this is semi-supervised variance reduction."""
    b = var(gr_tune) + (n / n_obs) * var(r_obs_tune * gr_obs_tune)
    d = cov(z0_tune, gr_tune)
    eta = project(d / b) if b > 1e-12 else 0.0
    estimate = float(z0_eval.mean() + eta * ((r_obs_eval * gr_obs_eval).mean() - gr_eval.mean()))
    psi_r = z0_eval - eta * gr_eval
    psi_o = eta * r_obs_eval * gr_obs_eval
    return {"estimate": estimate, "variance": var(psi_r) / n + var(psi_o) / n_obs, "eta": eta}


def obs_pseudo(obs: dict, mu_o, propensity_model) -> np.ndarray:
    """u_O of handoff Section 4.2: the OBS residual part of the AIPW score."""
    e = np.clip(propensity_model.predict_proba(obs["x"])[:, 1], 0.02, 0.98)
    m0, m1 = mu_o.predict(obs["x"])
    a, y = obs["a"], obs["y"]
    return a * (y - m1) / e - (1 - a) * (y - m0) / (1 - e)


def transported_obs_estimate(ghat_rct: np.ndarray, u_obs: np.ndarray, r_obs: np.ndarray) -> float:
    """theta_O of handoff Section 4.2: the RCT covariate mean of the OBS prediction
    plus the transported OBS residual, so both estimates share the RCT sample."""
    return float(np.mean(ghat_rct) + np.mean(r_obs * u_obs))


def combination_moments(z0_rct: np.ndarray, ghat_rct: np.ndarray, u_obs: np.ndarray,
                        r_obs: np.ndarray, n_rct: int, n_obs: int) -> dict:
    """V_R, V_O and the shared-covariate covariance C_RO of handoff Section 4.2."""
    v_r = var(z0_rct) / n_rct
    v_o = var(ghat_rct) / n_rct + var(r_obs * u_obs) / n_obs
    c_ro = cov(z0_rct, ghat_rct) / n_rct
    return {"V_R": v_r, "V_O": v_o, "C_RO": c_ro,
            "theta_R": float(np.mean(z0_rct)),
            "theta_O": transported_obs_estimate(ghat_rct, u_obs, r_obs)}


def shrinkage_weight(m: dict) -> dict:
    """Handoff Section 4.2, computed entirely on the tuning roles."""
    spread = m["V_R"] + m["V_O"] - 2 * m["C_RO"]
    bias2 = max((m["theta_R"] - m["theta_O"]) ** 2 - spread, 0.0)
    denom = spread + bias2
    gamma = project((m["V_R"] - m["C_RO"]) / denom) if abs(denom) > 1e-12 else 0.0
    return {"gamma": gamma, "bias_squared": bias2, "spread": spread}


def adaptive_weight(m: dict, z: float = 1.959963984540054) -> dict:
    """Handoff Section 4.3: a hard compatibility test, then the correlated-estimator
    minimum-variance weight."""
    spread = m["V_R"] + m["V_O"] - 2 * m["C_RO"]
    detected = abs(m["theta_R"] - m["theta_O"]) > z * np.sqrt(max(spread, 1e-12))
    gamma = 0.0 if detected else (project((m["V_R"] - m["C_RO"]) / spread)
                                  if abs(spread) > 1e-12 else 0.0)
    return {"gamma": float(gamma), "bias_detected": float(detected), "spread": spread}


def combine(theta_r: float, theta_o: float, gamma: float, m_eval: dict) -> dict:
    """The frozen weight applied to the evaluation-role component estimates, with the
    variance of a combination of two correlated estimators."""
    estimate = (1 - gamma) * theta_r + gamma * theta_o
    variance = ((1 - gamma) ** 2 * m_eval["V_R"] + gamma ** 2 * m_eval["V_O"]
                + 2 * gamma * (1 - gamma) * m_eval["C_RO"])
    return {"estimate": float(estimate), "variance": float(max(variance, 0.0)), "gamma": float(gamma)}


def naive_pool_spec(z_rct: np.ndarray, z_obs_transported: np.ndarray,
                    n_rct: int, n_obs: int) -> dict:
    """Handoff Section 4.4: the literal pooled average with denominator n_R + n_O."""
    total = n_rct + n_obs
    estimate = (n_rct * float(np.mean(z_rct)) + n_obs * float(np.mean(z_obs_transported))) / total
    variance = ((n_rct / total) ** 2 * var(z_rct) / n_rct
                + (n_obs / total) ** 2 * var(z_obs_transported) / n_obs)
    return {"estimate": estimate, "variance": variance}


def source_indicator_coefficients(beta: np.ndarray) -> np.ndarray:
    """Handoff Section 4.6: predict with the source indicator set to the trial.

    The indicator contributes its coefficient once.  Adding it to every basis
    coefficient, as an earlier version did, changes the fitted function.
    """
    out = beta[:-1].copy()
    out[0] += beta[-1]            # the first basis column is the intercept
    return out


def fit_pooled_outcome(rct_tune: dict, obs_tune: dict, spec: str = "flexible"):
    """One outcome model on the concatenation, with a source indicator."""
    from fusion_core import fit_outcome
    x = np.vstack([rct_tune["x"], obs_tune["x"]])
    source = np.r_[np.ones(len(rct_tune["x"])), np.zeros(len(obs_tune["x"]))][:, None]
    pooled = {"x": np.column_stack([x, source]),
              "a": np.r_[rct_tune["a"], obs_tune["a"]],
              "y": np.r_[rct_tune["y"], obs_tune["y"]]}
    model = fit_outcome(pooled, spec=spec)

    class Wrapped:
        def predict(self, z):
            return model.predict(np.column_stack([z, np.ones((len(z), 1))]))

    return Wrapped()


def experimental_grounding(basis_fn, rct_tune: dict, s_rct_tune, obs_tune: dict,
                           mu_o, propensity_model, rho: float = 1e-2):
    """Section 4.5: an OBS CATE estimate plus an RCT-fitted correction."""
    e_o = np.clip(propensity_model.predict_proba(obs_tune["x"])[:, 1], 0.02, 0.98)
    m0, m1 = mu_o.predict(obs_tune["x"])
    pseudo_o = aipw_score({"a": obs_tune["a"], "y": obs_tune["y"], "propensity": e_o}, m0, m1)
    b_o = basis_fn(obs_tune["x"])
    beta_o = np.linalg.solve(b_o.T @ b_o / len(b_o) + rho * np.eye(b_o.shape[1]),
                             (b_o * pseudo_o[:, None]).mean(0))
    b_r = basis_fn(rct_tune["x"])
    residual = s_rct_tune.z0 - b_r @ beta_o
    beta_q = np.linalg.solve(b_r.T @ b_r / len(b_r) + rho * np.eye(b_r.shape[1]),
                             (b_r * residual[:, None]).mean(0))
    return beta_o + beta_q


def domain_indicator_pool(basis_fn, rct_tune: dict, s_rct_tune, obs_tune: dict,
                          mu_o, propensity_model, rho: float = 1e-2):
    """Section 4.6: the CATE analogue of naive pooling."""
    e_o = np.clip(propensity_model.predict_proba(obs_tune["x"])[:, 1], 0.02, 0.98)
    m0, m1 = mu_o.predict(obs_tune["x"])
    pseudo_o = aipw_score({"a": obs_tune["a"], "y": obs_tune["y"], "propensity": e_o}, m0, m1)
    b_r, b_o = basis_fn(rct_tune["x"]), basis_fn(obs_tune["x"])
    design = np.vstack([np.column_stack([b_r, np.ones(len(b_r))]),
                        np.column_stack([b_o, np.zeros(len(b_o))])])
    target = np.r_[s_rct_tune.z0, pseudo_o]
    beta = np.linalg.solve(design.T @ design / len(design) + rho * np.eye(design.shape[1]),
                           (design * target[:, None]).mean(0))
    return beta[:-1] + beta[-1] * 0.0, beta
