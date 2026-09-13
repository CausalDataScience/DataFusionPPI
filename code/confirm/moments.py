"""Corollary 7.1 sieve moments, computed by three independent routes.

Notation follows manuscript/main/4.tex.  With b(X) the fixed basis,

    Gamma  = E_R(b b')                                        (4.tex:291)
    beta_p = Gamma^{-1} E_R(b tau),   zeta_p = b' beta_p
    psi_R(lambda, omega) = w b {Ztilde_0 + lambda Delta - omega ghat - (1-omega) zeta_p}
    psi_O(omega)         = omega r0 b (ghat - zeta_p)          (4.tex:296)
    J(lambda, omega)     = tr[Gamma^-1 Var_R psi_R]/n_R + tr[Gamma^-1 Var_O psi_O]/n_O

and Corollary 7.1 (4.tex:334) defines

    A_p = tr[Gamma^-1 Var_R(w b Delta)]
    B_p = tr[Gamma^-1 Var_R{w b (ghat - zeta_p)}]
        + (n_R/n_O) tr[Gamma^-1 Var_O{r0 b (ghat - zeta_p)}]
    C_p = tr[Gamma^-1 Cov_R{w b (Ztilde_0 - zeta_p), w b Delta}]
    D_p = tr[Gamma^-1 Cov_R{w b (Ztilde_0 - zeta_p), w b (ghat - zeta_p)}]

so that n_R {J(l, o) - J(0, 0)} = 2 C_p l - 2 D_p o + A_p l^2 + B_p o^2 - 2 l o E_p,
with the cross term E_p = tr[Gamma^-1 Cov_R(w b Delta, w b (ghat - zeta_p))].  The
corollary states that E_p vanishes for DRF; it does so because the doubly robust
Delta has conditional mean zero given X while (ghat - zeta_p) is a function of X,
and this module measures E_p rather than assuming it.

Three routes are provided on purpose.  ``covariance_route`` forms the p-by-p
covariance matrices and takes a trace, which is what the exploratory
implementation does.  ``whitened_route`` never forms a covariance matrix and
instead averages the scalar v' Gamma^-1 u.  ``direct_risk`` ignores the moments
entirely and evaluates J from its definition.  M1 requires the three to agree on
a fixed array to 1e-10, which is a check on the algebra rather than on sampling.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Fixture:
    """One deterministic mechanism fixture.

    Arrays carry no randomness at read time.  ``weight`` is the learner weight
    w_j, one for the doubly robust loss and the R-learner factor otherwise.
    """
    name: str
    b_target: np.ndarray      # basis on the target law, shape (n_t, p)
    tau_target: np.ndarray    # true CATE on the target law, shape (n_t,)
    b_rct: np.ndarray         # basis on the trial rows, shape (n_r, p)
    z0_rct: np.ndarray        # base pseudo-outcome Ztilde_0
    delta_rct: np.ndarray     # Delta, the increment that lambda multiplies
    ghat_rct: np.ndarray      # observational prediction on the trial rows
    weight_rct: np.ndarray    # learner weight w_j on the trial rows
    b_obs: np.ndarray         # basis on the observational rows, shape (n_o, p)
    ghat_obs: np.ndarray
    ratio_obs: np.ndarray     # r0 on the observational rows
    n_rct_tune: int
    n_obs_tune: int

    ridge: float = 0.0        # P0-2 requires rho = 0 on the Theorem 7 path

    def gram(self) -> np.ndarray:
        """Gamma_R, plus the ridge when one is declared.

        The confirmatory path runs at ridge zero, which is what Corollary 7.1
        assumes.  The exploratory engine adds 1e-10 for conditioning, and the
        parameter exists so the two can be compared on equal terms rather than
        being declared different for an unexamined reason.
        """
        gram = self.b_target.T @ self.b_target / len(self.b_target)
        if self.ridge:
            gram = gram + self.ridge * np.eye(len(gram))
        return gram

    def projection(self) -> np.ndarray:
        """beta_p = Gamma^{-1} E_R(b tau)."""
        moment = (self.b_target * self.tau_target[:, None]).mean(axis=0)
        return np.linalg.solve(self.gram(), moment)


def _contributions(fx: Fixture) -> dict[str, np.ndarray]:
    beta_p = fx.projection()
    zeta_r = fx.b_rct @ beta_p
    zeta_o = fx.b_obs @ beta_p
    w = fx.weight_rct[:, None]
    return {
        "u_base": w * fx.b_rct * (fx.z0_rct - zeta_r)[:, None],
        "u_delta": w * fx.b_rct * fx.delta_rct[:, None],
        "u_gap_rct": w * fx.b_rct * (fx.ghat_rct - zeta_r)[:, None],
        "v_gap_obs": fx.b_obs * (fx.ratio_obs * (fx.ghat_obs - zeta_o))[:, None],
        "zeta_rct": zeta_r,
        "zeta_obs": zeta_o,
        "beta_p": beta_p,
    }


# ------------------------------------------------------------------ route 1
def covariance_route(fx: Fixture) -> dict[str, float]:
    """Form the covariance matrices and take the whitened trace."""
    parts = _contributions(fx)
    inverse = np.linalg.inv(fx.gram())

    def tr_cov(left: np.ndarray, right: np.ndarray) -> float:
        stacked = np.column_stack([left, right])
        block = np.cov(stacked.T, ddof=1)
        p = left.shape[1]
        return float(np.trace(inverse @ block[:p, p:]))

    a_p = tr_cov(parts["u_delta"], parts["u_delta"])
    b_p = (tr_cov(parts["u_gap_rct"], parts["u_gap_rct"])
           + (fx.n_rct_tune / fx.n_obs_tune)
           * tr_cov(parts["v_gap_obs"], parts["v_gap_obs"]))
    c_p = tr_cov(parts["u_base"], parts["u_delta"])
    d_p = tr_cov(parts["u_base"], parts["u_gap_rct"])
    e_p = tr_cov(parts["u_delta"], parts["u_gap_rct"])
    return {"A_p": a_p, "B_p": b_p, "C_p": c_p, "D_p": d_p, "E_p": e_p}


# ------------------------------------------------------------------ route 2
def whitened_route(fx: Fixture) -> dict[str, float]:
    """Average the scalar v' Gamma^{-1} u; no covariance matrix is ever formed.

    Uses tr(Gamma^{-1} E[u v']) = E[v' Gamma^{-1} u], with both arrays centred by
    their own column means and the same 1/(n-1) convention as route 1.  The
    p-by-p covariance block of route 1 is never built.
    """
    parts = _contributions(fx)
    gram = fx.gram()

    def paired(left: np.ndarray, right: np.ndarray) -> float:
        a = left - left.mean(axis=0, keepdims=True)
        c = right - right.mean(axis=0, keepdims=True)
        solved = np.linalg.solve(gram, a.T).T      # rows of Gamma^{-1} u_i
        return float((c * solved).sum() / (len(left) - 1))

    a_p = paired(parts["u_delta"], parts["u_delta"])
    b_p = (paired(parts["u_gap_rct"], parts["u_gap_rct"])
           + (fx.n_rct_tune / fx.n_obs_tune)
           * paired(parts["v_gap_obs"], parts["v_gap_obs"]))
    c_p = paired(parts["u_base"], parts["u_delta"])
    d_p = paired(parts["u_base"], parts["u_gap_rct"])
    e_p = paired(parts["u_delta"], parts["u_gap_rct"])
    return {"A_p": a_p, "B_p": b_p, "C_p": c_p, "D_p": d_p, "E_p": e_p}


# ------------------------------------------------------------------ route 3
def direct_risk(fx: Fixture, lam: float, omega: float) -> float:
    """J(lambda, omega) straight from Equation (4.tex:301), no moments used."""
    parts = _contributions(fx)
    inverse = np.linalg.inv(fx.gram())
    w = fx.weight_rct[:, None]
    psi_r = w * fx.b_rct * (
        fx.z0_rct + lam * fx.delta_rct
        - omega * fx.ghat_rct - (1.0 - omega) * parts["zeta_rct"])[:, None]
    psi_o = omega * parts["v_gap_obs"]
    tr_r = float(np.trace(inverse @ np.cov(psi_r.T, ddof=1)))
    tr_o = float(np.trace(inverse @ np.cov(psi_o.T, ddof=1)))
    return tr_r / fx.n_rct_tune + tr_o / fx.n_obs_tune


def quadratic_risk(moments: dict[str, float], fx: Fixture,
                   lam: float, omega: float) -> float:
    """J(0,0) plus the Corollary 7.1 quadratic, for comparison with route 3."""
    base = direct_risk(fx, 0.0, 0.0)
    increment = (2 * moments["C_p"] * lam - 2 * moments["D_p"] * omega
                 + moments["A_p"] * lam ** 2 + moments["B_p"] * omega ** 2
                 - 2 * lam * omega * moments["E_p"])
    return base + increment / fx.n_rct_tune


def oracle_coordinates(moments: dict[str, float],
                       omega_cap: float = 1.0) -> tuple[float, float]:
    """Coordinatewise projection of the unconstrained minimiser (4.tex:346)."""
    a_p, b_p = moments["A_p"], moments["B_p"]
    lam = -moments["C_p"] / a_p if a_p > 0 else 0.0
    omega = moments["D_p"] / b_p if b_p > 0 else 0.0
    return float(np.clip(lam, 0.0, 1.0)), float(np.clip(omega, 0.0, omega_cap))


def oracle_gain(moments: dict[str, float]) -> float:
    """Equation (4.tex:354): the unconstrained first-order improvement."""
    a_p, b_p = moments["A_p"], moments["B_p"]
    return ((moments["C_p"] ** 2 / a_p if a_p > 0 else 0.0)
            + (moments["D_p"] ** 2 / b_p if b_p > 0 else 0.0))
