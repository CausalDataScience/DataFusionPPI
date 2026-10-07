"""DRF (Algorithm 2 with the doubly robust pseudo-outcome): the fused sieve learner of the CATE tau_0(x).

A candidate is linear in the sieve, t(x) = zeta(x)'beta with zeta = (1, principal components, g).  For each
(lambda, omega) of a grid it is fitted on the tuning samples, and the candidate with the smallest validation
score on the evaluation samples is returned.  lambda blends the OBS regression into the pseudo-outcome, and
omega moves weight from the trial rows to the OBS rows, reweighted by the density ratio r.

    nuis = fit_nuisances(R_nuis, O_nuis, e, same_domain)                 line 1, shared with AIPWF
    r = density_ratio(nuis)                                              line 1, the ratio of the sieve CATE
    lam, omega, beta = fit(nuis, r, R_tune, O_tune, R_eval, O_eval)      lines 2-4
    tau_hat = sieve(nuis, x) @ beta

Where each object of the paper is computed

    Def. 5    sieve zeta = (1, principal components, g)      sieve
    Prop. 6   the ratio balances the products of the sieve   balance_features, density_ratio
    Def. 3    weight and pseudo-outcome of DRF               pseudo_outcome
    Prop. 4   normal equation of the linear sieve, Eq. (12)  fit_candidate
    Def. 4    validation score, Eq. (10)                     validation_score
    Alg. 2    lines 2-4                                      fit
    comparison learner that uses no OBS input                rct_only
    comparison learner of Kallus, Puli and Shalit (2018)     two_step

rf.py is the same algorithm with another pseudo_outcome.
"""
import numpy as np

from . import method as M
from .method import aipw_score, blended_regression, fused_score

GRID = tuple((lam, omega) for lam in (0.0, 0.5, 1.0) for omega in (0.0, 0.5, 1.0))   # contains (0, 0)


def sieve(nuis, x):
    """Def. 5: zeta(x) = (1, principal components, g)."""
    return np.column_stack([nuis.components(x), nuis.g(x)])


def balance_features(components, g):
    """Prop. 6(ii) for the sieve zeta = (components, g): the products g^2, zeta g and zeta_k zeta_l (k <= l).

    Listing every distinct product once gives g, g^2, and for the principal components Z_k: Z_k, Z_k g and
    Z_k Z_l.  The constant 1 * 1 is left out: a calibrated ratio averages to one, which balances a constant.
    """
    p = components.shape[1]
    products = [components[:, k] * components[:, l] for k in range(p) for l in range(k, p)][1:]
    return np.column_stack([g ** 2, components * g[:, None]] + products)


def density_ratio(nuis):
    """Line 1.  Prop. 6(ii): the sieve CATE needs the products of balance_features, each balanced up to its
    standard error.  When r_0 = 1 is known, r = 1."""
    return M.density_ratio(nuis, lambda x: balance_features(nuis.components(x), nuis.g(x)),
                           within_standard_error=True)


def pseudo_outcome(sample, mu, e):
    """Def. 3, DRF: the weight kappa = 1 and the pseudo-outcome Z = phi(V; mu) of trial rows, for a
    regression mu.  With mu = mu_lambda this is (kappa, Z^lambda)."""
    return np.ones(len(sample["y"])), aipw_score(sample, mu, e)


def fit_candidate(lam, omega, nuis, r, R_tune, O_tune, ridge=0.01, pseudo_outcome=pseudo_outcome):
    """Defs. 3 and 5: the candidate t = zeta'beta that minimises, on the tuning samples, the empirical
    proxy risk Eq. (9) plus ridge ||beta||^2.  By Prop. 4 its coefficients solve the normal equation Eq. (12):

        {(1 - omega) P_R(kappa zeta zeta') + omega P_O(r zeta zeta') + ridge I} beta
            = P_R{kappa zeta (Z^lambda - omega g)} + omega P_O(r zeta g).
    """
    kappa, Z = pseudo_outcome(R_tune, blended_regression(nuis, lam), nuis.e)
    zeta_R, g_R = sieve(nuis, R_tune["x"]), nuis.g(R_tune["x"])
    zeta_O, g_O, r_O = sieve(nuis, O_tune["x"]), nuis.g(O_tune["x"]), r(O_tune["x"])
    H = ((1 - omega) * (zeta_R * kappa[:, None]).T @ zeta_R / len(kappa)
         + omega * (zeta_O * r_O[:, None]).T @ zeta_O / len(r_O))
    rhs = zeta_R.T @ (kappa * (Z - omega * g_R)) / len(kappa) + omega * zeta_O.T @ (r_O * g_O) / len(r_O)
    return np.linalg.solve(H + ridge * np.eye(len(rhs)), rhs)


def validation_score(beta, omega, nuis, r, R_eval, O_eval):
    """Def. 4, Eq. (10), on the evaluation samples, for the candidate t = zeta'beta:

        P_R{(Z^0 - t)^2} + omega P_O{r (g - t)^2} - omega P_R{(g - t)^2}.
    """
    t_R, t_O = sieve(nuis, R_eval["x"]) @ beta, sieve(nuis, O_eval["x"]) @ beta
    Z0 = fused_score(R_eval, nuis, 0.0)
    g_R, g_O, r_O = nuis.g(R_eval["x"]), nuis.g(O_eval["x"]), r(O_eval["x"])
    return float(np.mean((Z0 - t_R) ** 2) + omega * np.mean(r_O * (g_O - t_O) ** 2)
                 - omega * np.mean((g_R - t_R) ** 2))


def fit(nuis, r, R_tune, O_tune, R_eval, O_eval, grid=GRID, ridge=0.01, pseudo_outcome=pseudo_outcome):
    """Algorithm 2, lines 2-4.  Returns the selected (lambda, omega) and the coefficients of its candidate,
    so the fitted CATE is x -> sieve(nuis, x) @ beta.

    A smaller grid restricts the learner: omega = 0 is lambda only, and lambda = 0 is omega only.
    """
    scored = []
    for lam, omega in grid:
        beta = fit_candidate(lam, omega, nuis, r, R_tune, O_tune, ridge, pseudo_outcome)   # line 2: fit on tune
        score = validation_score(beta, omega, nuis, r, R_eval, O_eval)                     # line 3: score on eval
        scored.append((score, lam, omega, beta))
    _, lam, omega, beta = min(scored, key=lambda c: c[:3])    # line 4: smallest score, ties to smaller (lambda, omega)
    return lam, omega, beta


def rct_only(nuis, R_tune, ridge=0.01, pseudo_outcome=pseudo_outcome):
    """The comparison learner that uses no OBS input: the same ridge regression with lambda = omega = 0, the
    trial regression mu_R alone, and the sieve (1, principal components) without g.

    Returns its coefficients, so the fitted CATE is x -> nuis.components(x) @ beta.
    """
    kappa, Z = pseudo_outcome(R_tune, nuis.mu_R, nuis.e)
    basis = nuis.components(R_tune["x"])
    H = (basis * kappa[:, None]).T @ basis / len(kappa)
    return np.linalg.solve(H + ridge * np.eye(basis.shape[1]), basis.T @ (kappa * Z) / len(kappa))


def two_step(nuis, R_tune, ridge=0.01, pseudo_outcome=pseudo_outcome):
    """The comparison learner of Kallus, Puli and Shalit (2018), experimental grounding, with the pseudo-outcome
    and basis of this file: the OBS effect prediction g is corrected by a function fitted on the trial,

        t(x) = g(x) + components(x)'beta,   beta minimises P_R{kappa (Z - g - components'beta)^2} + ridge ||beta||^2,

    where (kappa, Z) are the weight and pseudo-outcome of the trial regression mu_R.  The coefficient of g is
    fixed at one; in the sieve of fit the data choose it.  Returns beta, so the fitted CATE is
    x -> nuis.g(x) + nuis.components(x) @ beta.
    """
    kappa, Z = pseudo_outcome(R_tune, nuis.mu_R, nuis.e)
    basis, g = nuis.components(R_tune["x"]), nuis.g(R_tune["x"])
    H = (basis * kappa[:, None]).T @ basis / len(kappa)
    return np.linalg.solve(H + ridge * np.eye(basis.shape[1]), basis.T @ (kappa * (Z - g)) / len(kappa))
