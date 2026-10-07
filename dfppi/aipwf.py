"""AIPWF (Algorithm 1): the fused estimator of the average treatment effect theta_0 = E_R tau_0(X).

    theta_hat = P_R{Z^lambda - omega g} + omega P_O(r g)                                  Def. 1, Eq. (3)

Z^lambda is the AIPW score of a trial row with the blended regression (1 - lambda) mu_R + lambda mu_O,
g = mu_O(., 1) - mu_O(., 0) is the OBS effect prediction, and r is the density ratio of the trial covariates
to the OBS covariates.  lambda = omega = 0 is the AIPW estimator of the trial alone.

    nuis = fit_nuisances(R_nuis, O_nuis, e, same_domain)                 line 1, shared with DRF and RF
    r = density_ratio(nuis)                                              line 1, the ratio of the ATE
    out = aipwf(nuis, r, R_tune, O_tune, R_eval, O_eval)                 lines 2-4

Where each object of the paper is computed

    Prop. 6   the ratio balances g exactly                   density_ratio
    Thm. 1    A, B, C, D on the tuning samples               aipwf, line 2
    Eq. (7)   the coefficients (lambda, omega), with B_cal   aipwf, line 2
    Def. 1    the estimate, Eq. (3)                          aipwf, line 3
    Thm. 2    its standard error (V + V_cal)                 aipwf, line 4

The standard error of the cross-fitted average of the three rotations is crossfit.crossfit_ate.

Comparison estimators that pool the trial-only estimate with an estimate from all rows: shrinkage, pretest.
"""
import numpy as np

from . import method as M
from .method import fused_score, var


def density_ratio(nuis):
    """Line 1.  Prop. 6(i): the ATE needs the single function g, balanced exactly (Sec. 3.2, Thm. 2(ii)).
    When r_0 = 1 is known, r = 1."""
    return M.density_ratio(nuis, lambda x: nuis.g(x)[:, None], within_standard_error=False)


def aipwf(nuis, r, R_tune, O_tune, R_eval, O_eval, use_lambda=True, use_omega=True, eps_factor=1e-3,
          omega_cal=True, calibrated=True):
    """Algorithm 1, lines 2-4.  Returns lambda, omega, the estimate, its standard error, and the two score
    vectors psi_R (trial evaluation rows) and psi_O (OBS evaluation rows) whose means make the estimate.

    use_lambda = False or use_omega = False fixes that coefficient at zero.  This gives the comparison
    estimators: RCT only (both False), lambda only, and omega only (PPI++).
    omega_cal = True     (Algorithm 1, Eq. (7)) with an estimated ratio, omega minimises V + V_cal, the
                         variance of Thm. 2(ii), rather than V:
                         omega = D / (B + B_cal), B_cal = n {Var_R(g) / n_nuis + Var_O(r g) / N_nuis}.
                         With r_0 = 1 known, B_cal = 0.
    omega_cal = False    omega = D / B, which ignores V_cal: the rule that PPI++ uses.
    calibrated = False   r is not calibrated to balance g (the bare classifier ratio), so V_cal does not apply.
                         experiment.py calls it with omega_cal = False as well (fusion_cls).
    """
    g = nuis.g
    n, N = len(R_eval["y"]), len(O_eval["y"])

    # Line 2.  A, B, C, D of Thm. 1 on the tuning samples, with r and the evaluation sizes (n, N).
    Z0 = fused_score(R_tune, nuis, 0.0)
    Delta = fused_score(R_tune, nuis, 1.0) - Z0
    g_R, rg_O = g(R_tune["x"]), r(O_tune["x"]) * g(O_tune["x"])
    A, C = var(Delta), np.cov(Z0, Delta)[0, 1]
    B, D = var(g_R) + n / N * var(rg_O), np.cov(Z0, g_R)[0, 1]
    B_cal = 0.0
    if omega_cal and nuis.r_base is not None:
        B_cal = n * (var(g(nuis.x_R)) / len(nuis.x_R) + var(r(nuis.x_O) * g(nuis.x_O)) / len(nuis.x_O))
    # Eq. (7): project the oracle coefficients -C / A and D / B of Thm. 1 onto [0, 1], with thresholds.
    eps_A = eps_B = eps_factor * var(Z0)
    lam = float(np.clip(-C / A, 0.0, 1.0)) if use_lambda and A > eps_A else 0.0
    omega = float(np.clip(D / (B + B_cal), 0.0, 1.0)) if use_omega and B > eps_B else 0.0

    # Line 3.  Def. 1 on the evaluation samples, written as the sum of two sample means, Eq. (3).
    psi_R = fused_score(R_eval, nuis, lam) - omega * g(R_eval["x"])
    psi_O = omega * r(O_eval["x"]) * g(O_eval["x"])
    estimate = psi_R.mean() + psi_O.mean()

    # Line 4.  The variance estimate V.  An estimated r adds the calibration variance of Thm. 2(ii),
    # V_cal = omega^2 {Var_R(g) / n + Var_O(r g) / N} on the nuisance samples, where r balances g.
    V = var(psi_R) / n + var(psi_O) / N
    V_cal = 0.0
    if nuis.r_base is not None and calibrated:
        V_cal = omega ** 2 * (var(g(nuis.x_R)) / len(nuis.x_R) + var(r(nuis.x_O) * g(nuis.x_O)) / len(nuis.x_O))
    return {"lambda": lam, "omega": omega, "estimate": float(estimate), "se": float(np.sqrt(V + V_cal)),
            "psi_R": psi_R, "psi_O": psi_O}


# ---------------------------------------------------------------------------- comparison estimators
def shrinkage(theta_R, se_R, theta_pool):
    """Precision-weighted shrinkage of the trial-only estimate toward a pooled estimate (as in Rosenman et al.):
    the pooled estimate gets the weight se_R^2 / {se_R^2 + (theta_pool - theta_R)^2}, its own variance neglected."""
    weight = se_R ** 2 / (se_R ** 2 + (theta_pool - theta_R) ** 2)
    return theta_R + weight * (theta_pool - theta_R)


def pretest(theta_R, se_R, theta_pool, z=1.959963984540054):
    """Test-then-pool: the pooled estimate if it lies within z standard errors of the trial-only estimate,
    otherwise the trial-only estimate."""
    return theta_pool if abs(theta_pool - theta_R) <= z * se_R else theta_R
