"""RF (Algorithm 2 with the weight and pseudo-outcome of the R-learner): the fused sieve learner of the CATE.

RF differs from DRF in Def. 3 only.  The sieve, the density ratio, the normal equation, the validation score
and the grid are those of drf.py:

    nuis = fit_nuisances(R_nuis, O_nuis, e, same_domain)                 line 1, shared with AIPWF
    r = density_ratio(nuis)                                              line 1, the ratio of the sieve CATE
    lam, omega, beta = fit(nuis, r, R_tune, O_tune, R_eval, O_eval)      lines 2-4
    tau_hat = sieve(nuis, x) @ beta
"""
import numpy as np

from . import drf as DRF
from .drf import GRID, density_ratio, sieve     # the same for both learners
from .method import principal_components


def pseudo_outcome(sample, mu, e):
    """Def. 3, RF: the weight and the pseudo-outcome of trial rows, for a regression mu,

        kappa = (A - e)^2 / {e (1 - e)},    Z = (Y - m) / (A - e),    m = e mu(., 1) + (1 - e) mu(., 0).

    With mu = mu_lambda this is (kappa, Z^lambda).
    """
    mu0, mu1 = mu(sample["x"])
    a = sample["a"]
    return (a - e) ** 2 / (e * (1 - e)), (sample["y"] - e * mu1 - (1 - e) * mu0) / (a - e)


def fit(nuis, r, R_tune, O_tune, R_eval, O_eval, grid=GRID, ridge=0.01):
    """Algorithm 2, lines 2-4, with the pseudo-outcome of RF.  Returns (lambda, omega, beta) as DRF.fit."""
    return DRF.fit(nuis, r, R_tune, O_tune, R_eval, O_eval, grid, ridge, pseudo_outcome)


def rct_only(nuis, R_tune, ridge=0.01):
    """The comparison learner that uses no OBS input, with the pseudo-outcome of RF."""
    return DRF.rct_only(nuis, R_tune, ridge, pseudo_outcome)


def two_step(nuis, R_tune, ridge=0.01):
    """The comparison learner of Kallus, Puli and Shalit (2018) with the pseudo-outcome of RF; see DRF.two_step."""
    return DRF.two_step(nuis, R_tune, ridge, pseudo_outcome)


# The integrative R-learner of Wu and Yang (2022): the trial and the OBS rows in one R-loss,
#     sum over rows of [Y - m(X, S) - {tau(X) + (1 - S) c(X)} {A - e(X, S)}]^2 + lam_tau |b_tau|^2 + lam_c |b_c|^2,
# with tau in the basis (1, principal components of the trial), the confounding function c in a basis of its
# own, c_basis below: (1, principal components of the OBS, g), so that the OBS prediction g enters c only and c,
# which appears on OBS rows only, takes its directions from the OBS covariates.  m(X, S) = E[Y | X, S] and
# e(X, S) = P(A = 1 | X, S) of each source.  A ridge penalty stands in for their SCAD penalty, and the grid
# is theirs: scale = lam_tau / lam_c in {0, 0.5, 1, 1.5}, lam_c on a log grid, chosen by the R-loss of both
# sources on the evaluation samples (where they use cross-validation).
IR_GRID = tuple((scale, lam_c) for scale in (0.0, 0.5, 1.0, 1.5) for lam_c in (1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0))


def c_basis(nuis, x_O_nuis, n_components=5):
    """The basis of the confounding function in the integrative R-learner: x -> (1, first principal components
    of the OBS nuisance covariates, g)."""
    components = principal_components(x_O_nuis, n_components)
    return lambda x: np.column_stack([components(x), nuis.g(x)])


def integrative(nuis, e_O, R_tune, O_tune, R_eval, O_eval, grid=IR_GRID, c_basis=None):
    """The integrative R-learner (Wu and Yang, 2022) with the regressions of nuis and the OBS propensity e_O.
    tau lies in (1, principal components), c in c_basis (by default the same).  Returns the coefficients of
    tau, so the fitted CATE is x -> nuis.components(x) @ beta."""
    c_basis = c_basis or nuis.components
    def rows(R, O):
        """Design (columns of tau, then of c on OBS rows) and residuals Y - m of the trial and OBS rows."""
        designs, residuals = [], []
        for sample, obs in ((R, False), (O, True)):
            x, a, y = sample["x"], sample["a"], sample["y"]
            mu0, mu1 = (nuis.mu_O if obs else nuis.mu_R)(x)
            e = e_O(x) if obs else np.full(len(y), nuis.e)
            w = (a - e)[:, None]
            c = w * c_basis(x)
            designs.append(np.hstack([w * nuis.components(x), c if obs else np.zeros_like(c)]))
            residuals.append(y - e * mu1 - (1 - e) * mu0)
        return np.vstack(designs), np.concatenate(residuals)

    X, r = rows(R_tune, O_tune)
    X_eval, r_eval = rows(R_eval, O_eval)
    D = nuis.components(R_tune["x"][:1]).shape[1]           # coefficients of tau; the rest belong to c
    gram, moment = X.T @ X / len(r), X.T @ r / len(r)
    best = None
    for scale, lam_c in grid:
        b = np.linalg.solve(gram + np.diag(np.r_[np.full(D, scale * lam_c), np.full(X.shape[1] - D, lam_c)]),
                            moment)
        loss = np.mean((r_eval - X_eval @ b) ** 2)
        if best is None or loss < best[0]:
            best = (loss, b)
    return best[1][:D]
