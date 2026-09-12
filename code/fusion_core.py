"""Shared building blocks for the DataFusionPPI experiments.

Implements Section 3 of materials/2026-09-09-implementation-handoff-JA.md:
roles, nuisances, pseudo-outcomes, the ATE estimator with Algorithm 1 and its
variance, the three density-ratio routes, the CATE losses with the sieve fit,
and the validation machinery of Algorithm 2.

Every formula carries the manuscript reference it implements.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from scipy.optimize import minimize

BASE_SEED = 190602
PROPENSITY_TRIM = (0.15, 0.85)      # handoff Section 3.1
RECTANGLE = (0.0, 1.0)              # Algorithm 1 coefficient rectangle
THRESHOLD_FACTOR = 1e-3             # t_A = t_B = factor * Var_R(Z0)
CI_Z = 1.959963984540054
ORACLE_DRAWS = 20000                # handoff Section 3.2
RATIO_CLIP = (0.05, 20.0)
DEFAULT_GRID = tuple((lam, om) for lam in (0.0, 0.25, 0.5, 0.75, 1.0)
                     for om in (0.0, 0.25, 0.5, 0.75, 1.0))
RIDGE_AXIS = (1e-2, 1e-1, 1.0)


# ---------------------------------------------------------------- utilities
def rng_for(*parts: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([BASE_SEED, *[int(p) for p in parts]]))


def stable_seed(*parts: object) -> int:
    import hashlib
    payload = "|".join(str(p) for p in parts).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "big")


def take(data: dict[str, np.ndarray], idx: np.ndarray) -> dict[str, np.ndarray]:
    return {k: v[idx] for k, v in data.items()}


def var(x: np.ndarray) -> float:
    return float(np.var(x, ddof=1)) if len(x) > 1 else 0.0


def cov(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.cov(x, y, ddof=1)[0, 1]) if len(x) > 1 else 0.0


def project(u: float, lo: float = RECTANGLE[0], hi: float = RECTANGLE[1]) -> float:
    return float(min(hi, max(lo, u)))


def trim_propensity(data: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], float]:
    """Handoff Section 3.1: keep units with e(X) inside PROPENSITY_TRIM."""
    e = data["propensity"]
    keep = (e >= PROPENSITY_TRIM[0]) & (e <= PROPENSITY_TRIM[1])
    return take(data, np.where(keep)[0]), float(1.0 - keep.mean())


def split_roles(data: dict[str, np.ndarray], fractions: tuple[float, float, float],
                rng: np.random.Generator) -> dict[str, dict[str, np.ndarray]]:
    """Assumption 1: independent nuisance, tuning, and evaluation roles."""
    n = len(data["x"])
    order = rng.permutation(n)
    n_nuis = int(round(fractions[0] * n))
    n_tune = int(round(fractions[1] * n))
    return {
        "nuis": take(data, order[:n_nuis]),
        "tune": take(data, order[n_nuis:n_nuis + n_tune]),
        "eval": take(data, order[n_nuis + n_tune:]),
    }


# ------------------------------------------------------------ outcome models
def flexible_features(x: np.ndarray) -> np.ndarray:
    d = x.shape[1]
    inter = [x[:, i] * x[:, j] for i in range(min(d, 6)) for j in range(i + 1, min(d, 6))]
    pieces = [x, x ** 2, np.sin(x)]
    if inter:
        pieces.append(np.column_stack(inter))
    return np.column_stack(pieces)


def linear_features(x: np.ndarray) -> np.ndarray:
    return x


@dataclass(frozen=True)
class OutcomeModel:
    mean: np.ndarray
    sd: np.ndarray
    model: Ridge
    spec: str
    clip: float | None = None      # a prespecified bound on the prediction

    def predict(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        raw = flexible_features(x) if self.spec == "flexible" else linear_features(x)
        f = (raw - self.mean) / self.sd
        zeros, ones = np.zeros((len(x), 1)), np.ones((len(x), 1))
        mu0 = self.model.predict(np.column_stack([f, zeros, zeros * f]))
        mu1 = self.model.predict(np.column_stack([f, ones, ones * f]))
        if self.clip is not None:
            mu0 = np.clip(mu0, -self.clip, self.clip)
            mu1 = np.clip(mu1, -self.clip, self.clip)
        return mu0, mu1


def fit_outcome(data: dict[str, np.ndarray], spec: str = "flexible", alpha: float = 5.0,
                clip: float | None = None) -> OutcomeModel:
    raw = flexible_features(data["x"]) if spec == "flexible" else linear_features(data["x"])
    mean, sd = raw.mean(0), raw.std(0)
    sd = np.where(sd > 1e-8, sd, 1.0)
    f = (raw - mean) / sd
    a = data["a"][:, None]
    model = Ridge(alpha=alpha, fit_intercept=True, solver="lsqr", tol=1e-10)
    model.fit(np.column_stack([f, data["a"], a * f]), data["y"])
    return OutcomeModel(mean, sd, model, spec, clip)


def fit_propensity(data: dict[str, np.ndarray], seed: int) -> Pipeline:
    model = Pipeline([("scale", StandardScaler()),
                      ("logit", LogisticRegression(C=1.0, max_iter=2000, random_state=seed))])
    model.fit(data["x"], data["a"])
    return model


# -------------------------------------------------------------- pseudo-scores
def aipw_score(data: dict[str, np.ndarray], mu0: np.ndarray, mu1: np.ndarray) -> np.ndarray:
    """Manuscript Equation (5): phi(V; q)."""
    a, y, e = data["a"], data["y"], data["propensity"]
    return mu1 - mu0 + a / e * (y - mu1) - (1.0 - a) / (1.0 - e) * (y - mu0)


@dataclass(frozen=True)
class Scores:
    z0: np.ndarray          # phi(V; muR)
    z1: np.ndarray          # phi(V; muO)
    delta: np.ndarray       # Z1 - Z0
    ghat: np.ndarray        # muO1 - muO0 on these covariates
    m_r: np.ndarray         # e*muR1 + (1-e)*muR0
    m_o: np.ndarray

    @property
    def d_m(self) -> np.ndarray:
        return self.m_o - self.m_r


def build_scores(data: dict[str, np.ndarray], mu_r: OutcomeModel, mu_o: OutcomeModel) -> Scores:
    r0, r1 = mu_r.predict(data["x"])
    o0, o1 = mu_o.predict(data["x"])
    e = data["propensity"]
    z0 = aipw_score(data, r0, r1)
    z1 = aipw_score(data, o0, o1)
    return Scores(z0=z0, z1=z1, delta=z1 - z0, ghat=o1 - o0,
                  m_r=e * r1 + (1 - e) * r0, m_o=e * o1 + (1 - e) * o0)


def ghat_on(x: np.ndarray, mu_o: OutcomeModel) -> np.ndarray:
    o0, o1 = mu_o.predict(x)
    return o1 - o0


# ------------------------------------------------------------ density ratios
@dataclass(frozen=True)
class Ratio:
    name: str
    fn: object
    diagnostics: dict

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return self.fn(x)


def ratio_unit(rct_nuis, obs_nuis, **kw) -> Ratio:
    return Ratio("unit", lambda x: np.ones(len(x)), {})


def ratio_oracle(rct_nuis, obs_nuis, oracle_fn=None, **kw) -> Ratio:
    return Ratio("oracle", oracle_fn, {})


def ratio_classifier(rct_nuis, obs_nuis, seed: int = 0, clip: bool = False, **kw) -> Ratio:
    """Definition 11: prior-corrected classifier odds."""
    x = np.vstack([rct_nuis["x"], obs_nuis["x"]])
    label = np.r_[np.ones(len(rct_nuis["x"])), np.zeros(len(obs_nuis["x"]))]
    model = Pipeline([("scale", StandardScaler()),
                      ("logit", LogisticRegression(C=1.0, max_iter=2000, random_state=seed))])
    model.fit(x, label)
    pi = float(label.mean())

    def fn(z: np.ndarray) -> np.ndarray:
        s = np.clip(model.predict_proba(z)[:, 1], 1e-8, 1 - 1e-8)
        r = (1 - pi) / pi * s / (1 - s)
        return np.clip(r, *RATIO_CLIP) if clip else r

    r_obs = fn(obs_nuis["x"])
    diag = {"normalization_error": float(r_obs.mean() - 1.0),
            "ess": float(r_obs.sum() ** 2 / np.sum(r_obs ** 2)),
            "max": float(r_obs.max()),
            "clip_fraction": float(np.mean((r_obs < RATIO_CLIP[0]) | (r_obs > RATIO_CLIP[1])))}
    return Ratio("classifier" + ("_clipped" if clip else ""), fn, diag)


def ratio_balanced(rct_nuis, obs_nuis, seed: int = 0, phi=None, name: str = "balanced", **kw) -> Ratio:
    """Definition 12: exponential tilt of a base ratio onto finite moments.

    Minimizes log P_O^nuis{r_base e^{xi'phi}} - xi' P_R^nuis{phi}; the gradient
    is the balance residual of Proposition 6.
    """
    base = ratio_classifier(rct_nuis, obs_nuis, seed=seed)
    phi = phi if phi is not None else (lambda z: z)
    phi_o, phi_r = phi(obs_nuis["x"]), phi(rct_nuis["x"])
    mu, sd = phi_o.mean(0), phi_o.std(0)
    sd = np.where(sd > 1e-8, sd, 1.0)
    phi_o, phi_r = (phi_o - mu) / sd, (phi_r - mu) / sd
    target = phi_r.mean(0)
    w0 = base(obs_nuis["x"])

    def objective(xi: np.ndarray) -> tuple[float, np.ndarray]:
        raw = phi_o @ xi
        shift = float(raw.max())          # keep the shift; it belongs in the value
        w = w0 * np.exp(raw - shift)
        s = w.mean()
        value = math.log(s) + shift - float(xi @ target)
        grad = (w[:, None] * phi_o).mean(0) / s - target
        return value, grad

    res = minimize(objective, np.zeros(phi_o.shape[1]), jac=True, method="L-BFGS-B",
                   options={"maxiter": 500, "ftol": 1e-14, "gtol": 1e-12})
    xi = res.x
    t_o = phi_o @ xi
    shift = float(t_o.max())
    log_norm = math.log(float((w0 * np.exp(t_o - shift)).mean())) + shift

    def fn(z: np.ndarray) -> np.ndarray:
        return base(z) * np.exp(((phi(z) - mu) / sd) @ xi - log_norm)

    r_obs = fn(obs_nuis["x"])
    residual = ((r_obs[:, None] * phi_o).mean(0) - target)
    diag = {"normalization_error": float(r_obs.mean() - 1.0),
            "balance_residual": float(np.max(np.abs(residual))),
            "ess": float(r_obs.sum() ** 2 / np.sum(r_obs ** 2)),
            "max": float(r_obs.max()),
            "converged": float(res.success)}
    return Ratio(name, fn, diag)


# ------------------------------------------------------- ATE: Algorithm 1
def ate_estimate(z_lambda: np.ndarray, ghat_r: np.ndarray, r_ghat_o: np.ndarray, omega: float) -> float:
    """Manuscript Equation (15)."""
    return float(z_lambda.mean() + omega * (r_ghat_o.mean() - ghat_r.mean()))


def tuning_moments(s_tune_r: Scores, ghat_o: np.ndarray, r_o: np.ndarray, n: int, n_obs: int) -> dict[str, float]:
    """Definition 5 evaluated on the tuning samples, with the evaluation-size ratio."""
    return {"A": var(s_tune_r.delta),
            "B": var(s_tune_r.ghat) + (n / n_obs) * var(r_o * ghat_o),
            "C": cov(s_tune_r.z0, s_tune_r.delta),
            "D": cov(s_tune_r.z0, s_tune_r.ghat),
            "var_z0": var(s_tune_r.z0)}


def algorithm1(m: dict[str, float]) -> dict[str, float]:
    """Equation (ate-separated-plugin) with thresholds and interval projection."""
    t = THRESHOLD_FACTOR * m["var_z0"]
    lam_raw = -m["C"] / m["A"] if m["A"] > t else 0.0
    om_raw = m["D"] / m["B"] if m["B"] > t else 0.0
    lam, om = project(lam_raw), project(om_raw)
    return {"lambda": lam, "omega": om, "lambda_raw": lam_raw, "omega_raw": om_raw,
            "lambda_fallback": float(m["A"] <= t), "omega_fallback": float(m["B"] <= t),
            "lambda_boundary": float(lam in RECTANGLE), "omega_boundary": float(om in RECTANGLE)}


def evaluation_variance(z_lambda: np.ndarray, ghat_r: np.ndarray, r_ghat_o: np.ndarray,
                        omega: float) -> float:
    """Equation (ate-evaluation-variance)."""
    psi_r = z_lambda - omega * ghat_r
    psi_o = omega * r_ghat_o
    return var(psi_r) / len(psi_r) + var(psi_o) / len(psi_o)


def theorem1_variance(z0: np.ndarray, delta: np.ndarray, ghat_r: np.ndarray,
                      r_ghat_o: np.ndarray, lam: float, om: float, n: int, n_obs: int) -> float:
    """Manuscript Equation (13), the exact variance at fixed coefficients."""
    return var(z0 + lam * delta - om * ghat_r) / n + (om ** 2) * var(r_ghat_o) / n_obs


# ----------------------------------------------------------- CATE machinery
def chi_weight(data: dict[str, np.ndarray]) -> np.ndarray:
    e = data["propensity"]
    return (data["a"] - e) ** 2 / (e * (1 - e))


def learner_parts(data: dict[str, np.ndarray], s: Scores, learner: str
                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Handoff Section 3.4: (w_j, Ztilde_{j,0}, Deltatilde_j)."""
    if learner == "DRF":
        return np.ones(len(s.z0)), s.z0, s.delta
    e = data["propensity"]
    denom = data["a"] - e
    return chi_weight(data), (data["y"] - s.m_r) / denom, -s.d_m / denom


def sieve_solver(b_r: np.ndarray, w: np.ndarray, zt0: np.ndarray, dzt: np.ndarray,
                 ghat_r: np.ndarray, b_o: np.ndarray, r_o: np.ndarray, ghat_o: np.ndarray,
                 omega: float, rho: float):
    """Proposition 10 items 2 and 3: one factorization per omega, affine in lambda."""
    n_r, n_o = len(b_r), len(b_o)
    h = (1 - omega) * (b_r * w[:, None]).T @ b_r / n_r + omega * (b_o * r_o[:, None]).T @ b_o / n_o
    h = h + rho * np.eye(b_r.shape[1])
    rhs0 = (b_r * (w * (zt0 - omega * ghat_r))[:, None]).mean(0) + omega * (b_o * (r_o * ghat_o)[:, None]).mean(0)
    rhs1 = (b_r * (w * dzt)[:, None]).mean(0)
    solve_fn = _symmetric_solver(h)

    def beta(lam: float) -> np.ndarray:
        return solve_fn(rhs0 + lam * rhs1)

    return beta, h


def _symmetric_solver(h: np.ndarray):
    """Factorize once and reuse.  At omega = 1 the RCT block drops out and the
    remaining Gram matrix can be singular on one-hot covariates, so add a
    trace-scaled jitter and fall back to a pseudo-inverse."""
    p = len(h)
    scale = max(float(np.trace(h)) / p, 1e-12)
    for power in range(-10, 1):
        candidate = h + (10.0 ** power) * scale * np.eye(p)
        try:
            factor = np.linalg.cholesky(candidate)
        except np.linalg.LinAlgError:
            continue
        return lambda rhs: np.linalg.solve(factor.T, np.linalg.solve(factor, rhs))
    inverse = np.linalg.pinv(h + 1e-6 * scale * np.eye(p), hermitian=True)
    return lambda rhs: inverse @ rhs


def validation_score(z0_eval: np.ndarray, zeta_r: np.ndarray, ghat_r: np.ndarray,
                     zeta_o: np.ndarray, ghat_o: np.ndarray, r_o: np.ndarray,
                     omega: float) -> float:
    """Definition 18 / Equation (63)."""
    g_r = (ghat_r - zeta_r) ** 2
    g_o = (ghat_o - zeta_o) ** 2
    return float(np.mean((z0_eval - zeta_r) ** 2 - omega * g_r) + np.mean(omega * r_o * g_o))


def eval_radius(bound_b: float, bar_r: float, n: int, n_obs: int, m_candidates: int,
                alpha: float = 0.05) -> float:
    """Theorem 6 Equation (64)."""
    return 2 * bound_b ** 2 * math.sqrt(2 * math.log(4 * m_candidates / alpha)) * (
        2 / math.sqrt(n) + bar_r / math.sqrt(n_obs))


def cate_risk(zeta: np.ndarray, tau: np.ndarray) -> float:
    return float(np.mean((zeta - tau) ** 2))


# ------------------------------------------------------------- sieve bases
def spline_basis_factory(x_ref: np.ndarray, n_knots: int, categorical: np.ndarray | None = None):
    """Additive cubic B-splines at empirical quantiles of the reference sample."""
    from scipy.interpolate import BSpline
    d = x_ref.shape[1]
    designs = []
    for j in range(d):
        col = x_ref[:, j]
        if categorical is not None and categorical[j]:
            levels = np.unique(col)[1:]
            designs.append(("cat", j, levels))
            continue
        qs = np.quantile(col, np.linspace(0, 1, n_knots + 2))
        lo, hi = col.min() - 1e-6, col.max() + 1e-6
        knots = np.r_[[lo] * 4, qs[1:-1], [hi] * 4]
        designs.append(("spline", j, knots))

    def basis(x: np.ndarray) -> np.ndarray:
        cols = [np.ones((len(x), 1))]
        for kind, j, payload in designs:
            if kind == "cat":
                cols.append(np.column_stack([(x[:, j] == lv).astype(float) for lv in payload])
                            if len(payload) else np.zeros((len(x), 0)))
            else:
                knots = payload
                n_basis = len(knots) - 4
                mat = np.zeros((len(x), n_basis))
                xv = np.clip(x[:, j], knots[0], knots[-1])
                for k in range(n_basis):
                    c = np.zeros(n_basis)
                    c[k] = 1.0
                    mat[:, k] = BSpline(knots, c, 3, extrapolate=False)(xv)
                mat = np.nan_to_num(mat)
                cols.append(mat[:, 1:])
        return np.column_stack(cols)

    return basis


def neural_basis_factory(obs_nuis: dict[str, np.ndarray], seed: int, width: int = 64,
                         epochs: int = 400, holdout: float = 0.1):
    """Frozen neural trunk trained on the OBS nuisance sample (handoff Section 3.4).

    Returns (basis_fn, ghat_fn).  The trunk never sees the tuning or evaluation
    samples, so Assumption 1 is preserved.
    """
    import torch

    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)
    x = torch.tensor(obs_nuis["x"], dtype=torch.float32)
    a = torch.tensor(obs_nuis["a"], dtype=torch.float32)[:, None]
    y = torch.tensor(obs_nuis["y"], dtype=torch.float32)[:, None]
    mx, sx = x.mean(0, keepdim=True), x.std(0, keepdim=True).clamp_min(1e-8)
    my, sy = y.mean(), y.std().clamp_min(1e-8)
    xn, yn = (x - mx) / sx, (y - my) / sy
    n = len(x)
    g = torch.Generator().manual_seed(seed)
    perm = torch.randperm(n, generator=g)
    n_hold = max(8, int(holdout * n))
    hold, train = perm[:n_hold], perm[n_hold:]

    trunk = torch.nn.Sequential(torch.nn.Linear(x.shape[1], width), torch.nn.ReLU(),
                                torch.nn.Linear(width, width), torch.nn.ReLU())
    head = torch.nn.Linear(width, 2)
    params = list(trunk.parameters()) + list(head.parameters())
    opt = torch.optim.Adam(params, lr=3e-3)

    def forward(idx):
        h = trunk(xn[idx])
        out = head(h)
        return out.gather(1, a[idx].long())

    best, best_state, patience = float("inf"), None, 0
    for _ in range(epochs):
        opt.zero_grad()
        loss = torch.nn.functional.mse_loss(forward(train), yn[train])
        loss.backward()
        opt.step()
        with torch.no_grad():
            v = float(torch.nn.functional.mse_loss(forward(hold), yn[hold]))
        if v < best - 1e-5:
            best, patience = v, 0
            best_state = [p.detach().clone() for p in params]
        else:
            patience += 1
            if patience >= 40:
                break
    if best_state is not None:
        with torch.no_grad():
            for p, q in zip(params, best_state):
                p.copy_(q)

    @torch.no_grad()
    def basis(z: np.ndarray) -> np.ndarray:
        h = trunk((torch.tensor(z, dtype=torch.float32) - mx) / sx).numpy()
        return np.column_stack([np.ones(len(z)), h])

    @torch.no_grad()
    def ghat(z: np.ndarray) -> np.ndarray:
        h = trunk((torch.tensor(z, dtype=torch.float32) - mx) / sx)
        out = head(h).numpy()
        return (out[:, 1] - out[:, 0]) * float(sy)

    return basis, ghat


def standardize_factory(basis, x_ref: np.ndarray):
    """Standardize a basis on the RCT tuning sample, keeping the intercept."""
    b = basis(x_ref)
    mu, sd = b.mean(0), b.std(0)
    keep = sd > 1e-8
    mu = np.where(keep, mu, 0.0)
    sd = np.where(keep, sd, 1.0)
    mu[0], sd[0] = 0.0, 1.0

    def fn(z: np.ndarray) -> np.ndarray:
        return (basis(z) - mu) / sd

    return fn

def prediction_bound(z0_tune: np.ndarray, factor: float = 3.0) -> float:
    """A prespecified range bound for the candidate functions and the scores.

    Theorem 6 assumes the validation scores lie in intervals of a fixed length.
    That fails if a candidate or a fitted prediction is unbounded: the score
    subtracts a term in the squared gap, so one exploding quantity drags it to
    minus infinity and captures the selector.  Clipping every candidate, every
    fitted prediction, and the pseudo-outcome to this bound is a deterministic,
    prespecified transformation, so the theorem applies to the clipped objects
    and the range condition holds by construction.

    The bound is read off the tuning-sample pseudo-outcome only.  Deriving it
    from the fitted prediction as well would let a diverged nuisance fit inflate
    the bound and defeat its purpose, which was observed with a neural trunk
    whose training diverged on about one replication in a hundred.
    """
    return factor * max(float(np.max(np.abs(z0_tune))), 1e-6)


def clip_prediction(values: np.ndarray, bound: float) -> np.ndarray:
    return np.clip(values, -bound, bound)

def variance_standard_error(values: np.ndarray, bootstrap: int = 2000, seed: int = 7) -> float:
    """Bootstrap standard error of a sample variance.

    The normal-theory expression var * sqrt(2/(K-1)) assumes a normal sampling
    distribution.  AIPW-based estimates are heavy tailed at small trial sizes,
    where that expression is optimistic, so the bootstrap is used instead.
    """
    values = np.asarray(values, dtype=float)
    k = len(values)
    if k < 3:
        return float("nan")
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, k, size=(bootstrap, k))
    return float(np.std(np.var(values[idx], axis=1, ddof=1), ddof=1))
