"""Four-role confirmatory CATE engine.

Plan blocker P0-1 asks for nuisance, tuning, selection and reporting samples that
are independent within each source.  This engine draws each role separately from
the design rather than partitioning one draw, so independence is exact by
construction and a draw fingerprint proves it.

The sieve solve follows the manuscript.  Writing b for the basis, the DRF loss

    P_R(Z_lambda - zeta)^2 + omega [ P_O{r (ghat - zeta)^2} - P_R(ghat - zeta)^2 ]

has normal equations H(omega) beta = v(lambda, omega) with

    H(omega) = P_R(b b') + omega [ P_O(r b b') - P_R(b b') ]
    v        = P_R(b Z_lambda) + omega [ P_O(r b ghat) - P_R(b ghat) ]

and E H(omega) = Gamma_R because E_O(r0 b b') equals E_R(b b'), which is the
statement in Theorem 7 item 2.

Selection and reporting both use the plain trial DR loss at lambda = 0, whose
expectation is the true risk plus a constant that does not depend on the
candidate.  A score that moved with lambda could not rank candidates that differ
in lambda.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from confirm import cate_dgp as dgp
from confirm import cate_protocol as cp
from confirm import moments as mo

PROPENSITY_CLIP = (0.02, 0.98)
TRIAL_TRIM = (0.15, 0.85)
SPLINE_KNOTS = {"spline3": 3, "spline5": 5}


# --------------------------------------------------------------------- basis
def spline_basis(z: np.ndarray, knots: np.ndarray) -> np.ndarray:
    """Truncated cubic power basis: 1, z, z^2, z^3, (z - k)_+^3 for each knot."""
    cols = [np.ones_like(z), z, z ** 2, z ** 3]
    cols += [np.clip(z - k, 0.0, None) ** 3 for k in knots]
    return np.column_stack(cols)


@dataclass(frozen=True)
class Basis:
    knots: np.ndarray
    direction: np.ndarray
    label: str

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return spline_basis(x @ self.direction, self.knots)

    @property
    def dimension(self) -> int:
        return 4 + len(self.knots)


def build_basis(sieve: str, x: np.ndarray) -> Basis:
    """Knots at fixed quantiles of the nuisance-stage index; frozen afterwards."""
    direction = dgp._direction(x.shape[1])
    z = x @ direction
    count = SPLINE_KNOTS[sieve]
    quantiles = np.linspace(0.0, 1.0, count + 2)[1:-1]
    return Basis(knots=np.quantile(z, quantiles), direction=direction, label=sieve)


# ---------------------------------------------------------------- four roles
def _fingerprint(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        digest.update(np.ascontiguousarray(array, dtype=np.float64).tobytes())
    return digest.hexdigest()[:16]


@dataclass
class RoleData:
    data: dict
    seed: int
    fingerprint: str


def draw_roles(design: dgp.Design, seed_key: str, sizes: dict[str, int]) -> dict:
    """Four independent draws per source, each with its own RNG namespace."""
    out: dict[str, dict[str, RoleData]] = {"RCT": {}, "OBS": {}}
    for source in ("RCT", "OBS"):
        for role in cp.ROLES:
            seed = cp.role_seed(seed_key, source, role)
            rng = np.random.default_rng(seed)
            n = sizes[f"{source}/{role}"]
            data = design.sample(source, n, rng)
            out[source][role] = RoleData(data=data, seed=seed,
                                         fingerprint=_fingerprint(data["x"], data["y"]))
    return out


# ------------------------------------------------------------------ nuisance
@dataclass
class Nuisance:
    mu_r: object
    mu_o: object
    propensity_o: object
    basis: Basis
    trial_propensity: float = 0.5


def _outcome_fit(data: dict, degree: int = 3) -> object:
    """One ridge per arm on a polynomial of the index; deterministic."""
    direction = dgp._direction(data["x"].shape[1])
    z = data["x"] @ direction

    def design(zz):
        return np.vander(zz, degree + 1, increasing=True)

    models = {}
    for arm in (0.0, 1.0):
        mask = data["a"] == arm
        if mask.sum() < degree + 2:
            models[arm] = None
            continue
        models[arm] = Ridge(alpha=1e-3).fit(design(z[mask]), data["y"][mask])

    class Fitted:
        def predict(self, x: np.ndarray):
            zz = x @ direction
            grid = design(zz)
            out = []
            for arm in (0.0, 1.0):
                model = models[arm]
                out.append(np.zeros(len(x)) if model is None else model.predict(grid))
            return out[0], out[1]

    return Fitted()


def fit_nuisance(roles: dict, sieve: str) -> Nuisance:
    rct, obs = roles["RCT"]["nuisance"].data, roles["OBS"]["nuisance"].data
    propensity = Pipeline([("scale", StandardScaler()),
                           ("logit", LogisticRegression(C=1.0, max_iter=2000))])
    propensity.fit(obs["x"], obs["a"])
    return Nuisance(mu_r=_outcome_fit(rct), mu_o=_outcome_fit(obs),
                    propensity_o=propensity, basis=build_basis(sieve, rct["x"]))


# -------------------------------------------------------------------- scores
def aipw(y, a, e, q0, q1):
    return q1 - q0 + a / e * (y - q1) - (1 - a) / (1 - e) * (y - q0)


@dataclass
class Scores:
    z0: np.ndarray
    delta: np.ndarray
    ghat: np.ndarray
    keep: np.ndarray


def trial_scores(data: dict, nuis: Nuisance) -> Scores:
    e = np.clip(data["propensity"], *TRIAL_TRIM)
    r0, r1 = nuis.mu_r.predict(data["x"])
    o0, o1 = nuis.mu_o.predict(data["x"])
    z0 = aipw(data["y"], data["a"], e, r0, r1)
    z1 = aipw(data["y"], data["a"], e, o0, o1)
    return Scores(z0=z0, delta=z1 - z0, ghat=o1 - o0,
                  keep=np.ones(len(z0), dtype=bool))


def obs_prediction(x: np.ndarray, nuis: Nuisance) -> np.ndarray:
    o0, o1 = nuis.mu_o.predict(x)
    return o1 - o0


# ------------------------------------------------------------- ratio routes
@dataclass
class Ratio:
    name: str
    values_obs: np.ndarray
    diagnostics: dict
    balance_status: str


def _classifier_ratio(rct_x, obs_x):
    x = np.vstack([rct_x, obs_x])
    label = np.r_[np.ones(len(rct_x)), np.zeros(len(obs_x))]
    model = Pipeline([("scale", StandardScaler()),
                      ("logit", LogisticRegression(C=1.0, max_iter=2000))])
    model.fit(x, label)
    prior = len(rct_x) / len(x)

    def ratio(z):
        p = np.clip(model.predict_proba(z)[:, 1], 1e-6, 1 - 1e-6)
        return (p / (1 - p)) * ((1 - prior) / prior)

    return ratio


def _balanced_ratio(base, features_rct, features_obs, target_weights=None):
    """Exponential tilt of a base ratio onto finite moments (Definition 12).

    Minimises log P_O{base e^{xi' phi}} - xi' P_R{phi}; the gradient is the
    balance residual of Proposition 6, so a converged optimiser leaves the
    observational feature mean equal to the trial feature mean.
    """
    target = features_rct.mean(axis=0)
    w0 = base / base.mean()

    def objective(xi):
        raw = features_obs @ xi
        shift = float(raw.max())
        w = w0 * np.exp(raw - shift)
        s = w.mean()
        value = math.log(s) + shift - float(xi @ target)
        grad = (w[:, None] * features_obs).mean(0) / s - target
        return value, grad

    result = minimize(objective, np.zeros(features_obs.shape[1]), jac=True,
                      method="L-BFGS-B",
                      options={"maxiter": 5000, "ftol": 1e-16, "gtol": 1e-14})
    xi = result.x
    raw = features_obs @ xi
    w = w0 * np.exp(raw - float(raw.max()))
    w = w / w.mean()
    residual = float(np.max(np.abs((w[:, None] * features_obs).mean(0)
                                   - features_rct.mean(axis=0))))
    return w, residual


def build_ratio(route: str, design: dgp.Design, roles: dict, nuis: Nuisance,
                obs_x: np.ndarray) -> Ratio:
    rct_x = roles["RCT"]["nuisance"].data["x"]
    obs_nuis_x = roles["OBS"]["nuisance"].data["x"]
    residual = float("nan")
    status = "N/A"

    if route == "unit":
        values = np.ones(len(obs_x))
    elif route == "oracle":
        values = design.oracle_ratio(obs_x)
    else:
        classifier = _classifier_ratio(rct_x, obs_nuis_x)
        if route == "classifier":
            values = classifier(obs_x)
        else:
            phi_rct = nuis.basis(rct_x)
            phi_obs = nuis.basis(obs_x)
            if route == "balanced_x_ghat":
                phi_rct = np.column_stack([phi_rct, obs_prediction(rct_x, nuis)])
                phi_obs = np.column_stack([phi_obs, obs_prediction(obs_x, nuis)])
            values, residual = _balanced_ratio(classifier(obs_x), phi_rct, phi_obs)
            status = "converged" if residual <= 1e-6 else "not_converged"

    positive = np.clip(values, 1e-12, None)
    diagnostics = {
        "mean": float(values.mean()),
        "second_moment": float((values ** 2).mean()),
        "ess_fraction": float(values.mean() ** 2 / (values ** 2).mean()),
        "normalisation_error": float(values.mean() - 1.0),
        "log_ratio_q05": float(np.quantile(np.log(positive), 0.05)),
        "log_ratio_q95": float(np.quantile(np.log(positive), 0.95)),
        "balance_residual": residual,
    }
    return Ratio(name=route, values_obs=values, diagnostics=diagnostics,
                 balance_status=status)


# -------------------------------------------------------------- sieve solve
def solve_candidate(b_rct, z0, delta, ghat_rct, b_obs, ghat_obs, ratio_obs,
                    lam: float, omega: float, rho: float):
    """One (lambda, omega, rho) candidate from the normal equations above."""
    gram_r = b_rct.T @ b_rct / len(b_rct)
    gram_o = (b_obs * ratio_obs[:, None]).T @ b_obs / len(b_obs)
    hessian = gram_r + omega * (gram_o - gram_r)
    z_lambda = z0 + lam * delta
    rhs = (b_rct * z_lambda[:, None]).mean(0) + omega * (
        (b_obs * (ratio_obs * ghat_obs)[:, None]).mean(0)
        - (b_rct * ghat_rct[:, None]).mean(0))
    matrix = hessian + rho * np.eye(len(hessian))
    condition = float(np.linalg.cond(matrix))
    beta = np.linalg.solve(matrix, rhs)
    return beta, condition


def dr_score(prediction: np.ndarray, z0: np.ndarray) -> float:
    """The plain trial DR loss, whose expectation is the risk plus a constant."""
    return float(np.mean((z0 - prediction) ** 2))


def true_risk(prediction: np.ndarray, tau_values: np.ndarray) -> float:
    return float(np.mean((prediction - tau_values) ** 2))
