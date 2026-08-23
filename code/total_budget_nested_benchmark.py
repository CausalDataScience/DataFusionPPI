#!/usr/bin/env python3
"""Total-budget nested benchmark for RCT/OBS ATE fusion."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import ssem_ate_pilot as legacy

SEED = 190602
N_RCT_VALUES = (20, 30)
N_OBS_FIXED = 5000
OBS_MULTIPLIER: int | None = None
REGIMES = ("shared", "shifted")
OUTER_FOLDS = 5
INNER_FOLDS = 2
DEFAULT_REPLICATIONS = 20
BOX = (-2.0, 2.0)
RATIO_CLIP = (0.05, 20.0)
STAR_URL = "https://vincentarelbundock.github.io/Rdatasets/csv/AER/STAR.csv"
STAR_SHA256 = "0e8b179ea3d883730b25008ca1293c4c8f886aa8d5d299c03ab6ff05d0123ae6"
STAR_COLUMNS = ("gender", "ethnicity", "birth", "lunchk", "schoolk")
PROJECT_DIR = Path(__file__).resolve().parents[1]
MATERIALS_DIR = PROJECT_DIR / "materials"
METHODS = (
    "rct_aipw",
    "haipw_only",
    "ppi_power_tuned_only",
    "full_datafusionppi",
)


def obs_size_for(n_rct: int) -> int:
    return N_OBS_FIXED if OBS_MULTIPLIER is None else OBS_MULTIPLIER * n_rct


def parse_n_rct_values(raw: str) -> tuple[int, ...]:
    values = tuple(int(value.strip()) for value in raw.split(",") if value.strip())
    if not values or any(value < OUTER_FOLDS for value in values) or len(set(values)) != len(values):
        raise argparse.ArgumentTypeError(
            f"n-RCT values must be distinct integers at least {OUTER_FOLDS}"
        )
    return values


def output_paths(stem: str) -> tuple[Path, Path, Path]:
    if not stem or Path(stem).name != stem or any(character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for character in stem):
        raise ValueError("output stem must be a simple alphanumeric, underscore, or hyphen name")
    return (
        MATERIALS_DIR / f"{stem}_replications.csv",
        MATERIALS_DIR / f"{stem}_summary.csv",
        MATERIALS_DIR / f"{stem}.png",
    )


@dataclass(frozen=True)
class OutcomeModel:
    mean: np.ndarray
    sd: np.ndarray
    model: Ridge

    def predict(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        f = (feature_map(x) - self.mean) / self.sd
        d0 = np.column_stack([f, np.zeros(len(x)), np.zeros_like(f)])
        d1 = np.column_stack([f, np.ones(len(x)), f])
        return self.model.predict(d0), self.model.predict(d1)


@dataclass(frozen=True)
class RatioModel:
    model: Pipeline
    prior: float

    def raw(self, x: np.ndarray) -> np.ndarray:
        p = np.clip(self.model.predict_proba(x)[:, 1], 1e-10, 1 - 1e-10)
        return p / (1 - p) * (1 - self.prior) / self.prior


@dataclass(frozen=True)
class StarSupport:
    x: np.ndarray
    u: np.ndarray
    y0: np.ndarray
    y1: np.ndarray
    tau: np.ndarray
    e_rct: np.ndarray
    e_obs: np.ndarray
    direction: np.ndarray

    def probabilities(self, regime: str) -> tuple[np.ndarray, np.ndarray]:
        size = len(self.x)
        p_obs = np.full(size, 1.0 / size)
        if regime == "shared":
            return p_obs.copy(), p_obs
        score = self.x @ self.direction / math.sqrt(self.x.shape[1])
        log_weight = 0.35 * score
        weight = np.exp(log_weight - log_weight.max())
        return weight / weight.sum(), p_obs


def rng_for(*parts: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([SEED, *parts]))


def stable_seed(*parts: object) -> int:
    return int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:4], "little")


def take(data: dict[str, np.ndarray], idx: np.ndarray) -> dict[str, np.ndarray]:
    return {key: value[idx] for key, value in data.items()}


def feature_map(x: np.ndarray) -> np.ndarray:
    d = min(x.shape[1], 4)
    interactions = [x[:, i] * x[:, j] for i in range(d) for j in range(i + 1, d)]
    pieces = [x, x**2, np.sin(x)]
    if interactions:
        pieces.append(np.column_stack(interactions))
    return np.column_stack(pieces)


def fit_outcome(data: dict[str, np.ndarray]) -> OutcomeModel:
    raw = feature_map(data["x"])
    mean, sd = raw.mean(0), raw.std(0)
    sd = np.where(sd > 1e-8, sd, 1.0)
    f = (raw - mean) / sd
    a = data["a"][:, None]
    design = np.column_stack([f, data["a"], a * f])
    model = Ridge(alpha=5.0, fit_intercept=True, solver="lsqr", tol=1e-8)
    model.fit(design, data["y"])
    return OutcomeModel(mean, sd, model)


def fit_ratio(rct: dict[str, np.ndarray], obs: dict[str, np.ndarray], seed: int) -> RatioModel:
    x = np.vstack([rct["x"], obs["x"]])
    label = np.r_[np.ones(len(rct["x"])), np.zeros(len(obs["x"]))]
    model = Pipeline([("scale", StandardScaler()), ("logit", LogisticRegression(C=1.0, max_iter=1000, solver="lbfgs", random_state=seed))])
    model.fit(x, label)
    return RatioModel(model, float(label.mean()))


def pseudo(data: dict[str, np.ndarray], mu0: np.ndarray, mu1: np.ndarray) -> np.ndarray:
    a, y, e = data["a"], data["y"], data["propensity"]
    return mu1 - mu0 + a / e * (y - mu1) - (1 - a) / (1 - e) * (y - mu0)


def var(x: np.ndarray) -> float:
    return float(np.var(x, ddof=1)) if len(x) > 1 else 0.0


def cov(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.cov(x, y, ddof=1)[0, 1]) if len(x) > 1 else 0.0


def tuning_moments(
    z0: np.ndarray,
    delta: np.ndarray,
    gr: np.ndarray,
    rgo: np.ndarray,
    rho: float,
) -> dict[str, float]:
    return {
        "a": var(delta),
        "b": var(gr) + rho * var(rgo),
        "c": cov(z0, delta),
        "d": cov(z0, gr),
        "h": cov(delta, gr),
    }


def coefficient_record(
    raw_lambda: float,
    raw_omega: float,
    lambda_used: float,
    omega_used: float,
    restriction: str,
    uses_ratio: float,
    fallback: float,
) -> dict[str, float | str]:
    raw = np.array([raw_lambda, raw_omega], dtype=float)
    raw = np.where(np.isfinite(raw), raw, 0.0)
    used = np.array([lambda_used, omega_used], dtype=float)
    return {
        "raw_lambda_used": float(raw[0]),
        "raw_omega_used": float(raw[1]),
        "lambda_used": float(used[0]),
        "omega_used": float(used[1]),
        "lambda_constrained_from_raw": float(not np.isclose(raw[0], used[0])),
        "omega_constrained_from_raw": float(not np.isclose(raw[1], used[1])),
        "lambda_at_boundary": float(np.isclose(abs(used[0]), BOX[1])),
        "omega_at_boundary": float(np.isclose(abs(used[1]), BOX[1])),
        "coefficient_fallback": fallback,
        "tuning_restriction": restriction,
        "uses_ratio_channel": uses_ratio,
    }


def tuning_objective(lam: float, omega: float, moments: dict[str, float]) -> float:
    return float(
        2.0 * lam * moments["c"]
        + lam**2 * moments["a"]
        - 2.0 * omega * moments["d"]
        - 2.0 * lam * omega * moments["h"]
        + omega**2 * moments["b"]
    )


def exact_box_minimizer(
    moments: dict[str, float],
    axis_candidates: list[tuple[float, float]],
) -> tuple[np.ndarray, np.ndarray, float]:
    lower, upper = BOX
    a_value, b_value = moments["a"], moments["b"]
    c_value, d_value, h_value = moments["c"], moments["d"], moments["h"]
    matrix = np.array([[a_value, -h_value], [-h_value, b_value]], dtype=float)
    target = np.array([-c_value, d_value], dtype=float)
    eigenvalues = np.linalg.eigvalsh(matrix)
    if eigenvalues.min() < -1e-8:
        raise AssertionError("Empirical tuning quadratic is not convex")
    fallback = float(eigenvalues.min() <= 1e-10 or np.linalg.cond(matrix) > 1e10)
    raw = np.linalg.pinv(matrix, rcond=1e-12) @ target
    raw = np.where(np.isfinite(raw), raw, 0.0)

    candidates: list[tuple[float, float]] = []
    if lower <= raw[0] <= upper and lower <= raw[1] <= upper:
        candidates.append((float(raw[0]), float(raw[1])))
    for fixed_lambda in BOX:
        omega = (d_value + h_value * fixed_lambda) / b_value if b_value > 1e-10 else 0.0
        candidates.append((fixed_lambda, float(np.clip(omega, lower, upper))))
    for fixed_omega in BOX:
        lam = (h_value * fixed_omega - c_value) / a_value if a_value > 1e-10 else 0.0
        candidates.append((float(np.clip(lam, lower, upper)), fixed_omega))
    candidates.extend((lam, omega) for lam in BOX for omega in BOX)
    candidates.extend(axis_candidates)
    candidates.append((0.0, 0.0))
    objective_values = [tuning_objective(lam, omega, moments) for lam, omega in candidates]
    best_index = min(range(len(candidates)), key=lambda index: (objective_values[index], index))
    return raw, np.asarray(candidates[best_index], dtype=float), fallback


def tune_all(moments: dict[str, float]) -> dict[str, dict[str, float | str]]:
    a_value, b_value = moments["a"], moments["b"]
    c_value, d_value, h_value = moments["c"], moments["d"], moments["h"]
    haipw_fallback = float(a_value <= 1e-10)
    ppi_fallback = float(b_value <= 1e-10)
    raw_haipw_lambda = -c_value / a_value if not haipw_fallback else 0.0
    raw_ppi_omega = d_value / b_value if not ppi_fallback else 0.0
    haipw_lambda = float(np.clip(raw_haipw_lambda, *BOX))
    ppi_omega = float(np.clip(raw_ppi_omega, *BOX))
    raw_full, used_full, fallback = exact_box_minimizer(
        moments,
        axis_candidates=[(haipw_lambda, 0.0), (0.0, ppi_omega)],
    )
    coefficients = {
        "rct_aipw": coefficient_record(
            0.0, 0.0, 0.0, 0.0, "fixed_lambda_0_omega_0", 0.0, 0.0
        ),
        "haipw_only": coefficient_record(
            raw_haipw_lambda,
            0.0,
            haipw_lambda,
            0.0,
            "restricted_omega_0",
            0.0,
            haipw_fallback,
        ),
        "ppi_power_tuned_only": coefficient_record(
            0.0,
            raw_ppi_omega,
            0.0,
            ppi_omega,
            "restricted_lambda_0",
            1.0,
            ppi_fallback,
        ),
        "full_datafusionppi": coefficient_record(
            float(raw_full[0]),
            float(raw_full[1]),
            float(used_full[0]),
            float(used_full[1]),
            "joint_box_lambda_omega",
            1.0,
            fallback,
        ),
    }
    if coefficients["rct_aipw"]["lambda_used"] != 0.0 or coefficients["rct_aipw"]["omega_used"] != 0.0:
        raise AssertionError("AIPW coefficient restriction failed")
    if coefficients["haipw_only"]["omega_used"] != 0.0:
        raise AssertionError("lambda-only outcome-model fusion omega restriction failed")
    if coefficients["ppi_power_tuned_only"]["lambda_used"] != 0.0:
        raise AssertionError("omega-only PPI correction lambda restriction failed")
    if not haipw_fallback and not np.isclose(coefficients["haipw_only"]["raw_lambda_used"], -c_value / a_value):
        raise AssertionError("lambda-only outcome-model fusion did not solve its restricted problem")
    if not ppi_fallback and not np.isclose(coefficients["ppi_power_tuned_only"]["raw_omega_used"], d_value / b_value):
        raise AssertionError("omega-only PPI correction did not solve its restricted problem")
    objectives = {
        method: tuning_objective(
            float(coefficient["lambda_used"]),
            float(coefficient["omega_used"]),
            moments,
        )
        for method, coefficient in coefficients.items()
    }
    full_objective = objectives["full_datafusionppi"]
    for reference in ("rct_aipw", "haipw_only", "ppi_power_tuned_only"):
        if full_objective > objectives[reference] + 1e-8:
            raise AssertionError(f"Full box objective exceeds {reference}")
    return coefficients


def folds(size: int, count: int, rng: np.random.Generator) -> list[np.ndarray]:
    return [np.asarray(x, int) for x in np.array_split(rng.permutation(size), count)]


def synthetic_sample(scm: legacy.SCMParameters, regime: str, n: int, rep: int) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], float]:
    rc = REGIMES.index(regime)
    rct = legacy.sample_dataset(rng_for(100, scm.scm_id, rc, n, rep), n, "RCT", regime, scm)
    obs = legacy.sample_dataset(
        rng_for(101, scm.scm_id, rc, n, rep), obs_size_for(n), "OBS", regime, scm
    )
    rct["true_ratio"] = legacy.oracle_ratio(rct["x"], regime, scm)
    obs["true_ratio"] = legacy.oracle_ratio(obs["x"], regime, scm)
    return rct, obs, legacy.true_theta(scm, regime)


def load_star_support() -> StarSupport:
    download = subprocess.run(
        ["curl", "-fsSL", STAR_URL],
        check=True,
        capture_output=True,
        timeout=30,
    )
    payload = download.stdout
    digest = hashlib.sha256(payload).hexdigest()
    if digest != STAR_SHA256:
        raise RuntimeError(f"STAR source hash mismatch: {digest}")
    frame = pd.read_csv(io.BytesIO(payload))
    if len(frame) != 11598 or not set(STAR_COLUMNS).issubset(frame.columns):
        raise RuntimeError("STAR source shape/schema mismatch")
    numeric = pd.DataFrame(index=frame.index)
    birth = pd.to_numeric(frame["birth"], errors="coerce")
    numeric["birth"] = birth.fillna(birth.median())
    numeric["birth_missing"] = birth.isna().astype(float)
    categorical = frame[["gender", "ethnicity", "lunchk", "schoolk"]].astype("string").fillna("Missing")
    encoded = pd.concat([numeric, pd.get_dummies(categorical, prefix=categorical.columns, dtype=float)], axis=1)
    raw = encoded.to_numpy(dtype=float)
    sd = raw.std(axis=0)
    raw = raw[:, sd > 1e-10]
    x = (raw - raw.mean(axis=0)) / raw.std(axis=0)
    rng = rng_for(6000)
    phi = feature_map(x)
    phi = (phi - phi.mean(axis=0)) / np.where(phi.std(axis=0) > 1e-8, phi.std(axis=0), 1.0)
    beta_m = rng.normal(0.0, 0.55 / math.sqrt(phi.shape[1]), phi.shape[1])
    beta_tau = rng.normal(0.0, 0.22 / math.sqrt(phi.shape[1]), phi.shape[1])
    u, noise = rng.normal(size=len(x)), rng.normal(size=len(x))
    baseline = phi @ beta_m + 0.7 * u + noise
    tau = 1.0 + phi @ beta_tau
    beta_r = rng.normal(0.0, 0.25 / math.sqrt(x.shape[1]), x.shape[1])
    beta_o = rng.normal(0.0, 0.35 / math.sqrt(x.shape[1]), x.shape[1])
    e_rct = np.clip(1.0 / (1.0 + np.exp(np.clip(-(0.1 + x @ beta_r), -35, 35))), 0.2, 0.8)
    e_obs = np.clip(1.0 / (1.0 + np.exp(np.clip(-(-0.2 + x @ beta_o + 0.9 * u), -35, 35))), 0.05, 0.95)
    direction = rng.normal(size=x.shape[1])
    support = StarSupport(x=x, u=u, y0=baseline, y1=baseline + tau, tau=tau, e_rct=e_rct, e_obs=e_obs, direction=direction)
    for regime in REGIMES:
        p_r, p_o = support.probabilities(regime)
        ratio = p_r / p_o
        if abs(float(np.sum(p_o * ratio)) - 1.0) > 1e-12:
            raise AssertionError("STAR exact ratio does not normalize")
        if not np.isfinite(float(np.sum(p_r * tau))):
            raise AssertionError("STAR exact target is not finite")
    return support


def star_sample(support: StarSupport, regime: str, n: int, rep: int) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], float]:
    p_r, p_o = support.probabilities(regime)
    regime_code = REGIMES.index(regime)
    r_rng, o_rng = rng_for(6100, regime_code, n, rep), rng_for(6101, regime_code, n, rep)
    r_id = r_rng.choice(len(support.x), size=n, replace=True, p=p_r)
    o_id = o_rng.choice(len(support.x), size=obs_size_for(n), replace=True, p=p_o)

    def observed(ids: np.ndarray, study: str, rng: np.random.Generator) -> dict[str, np.ndarray]:
        propensity = support.e_rct[ids] if study == "RCT" else support.e_obs[ids]
        treatment = rng.binomial(1, propensity).astype(float)
        outcome = np.where(treatment == 1, support.y1[ids], support.y0[ids])
        true_ratio = (p_r / p_o)[ids]
        return {"x": support.x[ids], "a": treatment, "y": outcome, "propensity": propensity, "tau": support.tau[ids], "true_ratio": true_ratio}

    rct = observed(r_id, "RCT", r_rng)
    obs = observed(o_id, "OBS", o_rng)
    theta = float(np.sum(p_r * support.tau))
    return rct, obs, theta


def inner_scores(rct: dict[str, np.ndarray], obs: dict[str, np.ndarray], seed: int) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    rf, of = folds(len(rct["x"]), 2, rng_for(seed, 1)), folds(len(obs["x"]), 2, rng_for(seed, 2))
    parts = {key: [] for key in ("z0", "delta", "gr", "rgo")}
    min_r, min_o = len(rct["x"]), len(obs["x"])
    for k in range(2):
        ri, oi = rf[k], of[k]
        rt, ot = np.setdiff1d(np.arange(len(rct["x"])), ri), np.setdiff1d(np.arange(len(obs["x"])), oi)
        tr, to, vr, vo = take(rct, rt), take(obs, ot), take(rct, ri), take(obs, oi)
        min_r = min(min_r, int(tr["a"].sum()), int(len(tr["a"]) - tr["a"].sum()))
        min_o = min(min_o, int(to["a"].sum()), int(len(to["a"]) - to["a"].sum()))
        mr, mo, ratio = fit_outcome(tr), fit_outcome(to), fit_ratio(tr, to, stable_seed(seed, k))
        r0, r1 = mr.predict(vr["x"]); o0r, o1r = mo.predict(vr["x"]); o0o, o1o = mo.predict(vo["x"])
        z0, z1 = pseudo(vr, r0, r1), pseudo(vr, o0r, o1r)
        parts["z0"].append(z0); parts["delta"].append(z1 - z0); parts["gr"].append(o1r - o0r)
        parts["rgo"].append(np.clip(ratio.raw(vo["x"]), *RATIO_CLIP) * (o1o - o0o))
    return {key: np.concatenate(value) for key, value in parts.items()}, {"min_rct_train_arm": float(min_r), "min_obs_train_arm": float(min_o)}


def evaluate(
    benchmark: str,
    scenario: int,
    regime: str,
    n: int,
    rep: int,
    rct: dict[str, np.ndarray],
    obs: dict[str, np.ndarray],
    theta: float,
    source_sha256: str = "not_applicable",
    support_rows: int = 0,
) -> list[dict[str, object]]:
    seed = stable_seed(benchmark, scenario, regime, n, rep)
    paired_design_id = f"{seed:08x}"
    n_obs = len(obs["x"])
    rf, of = folds(n, 5, rng_for(seed, 10)), folds(n_obs, 5, rng_for(seed, 11))
    fold_estimates = {method: [] for method in METHODS}
    variance_proxies = {method: 0.0 for method in METHODS}
    coefficient_keys = (
        "raw_lambda_used",
        "raw_omega_used",
        "lambda_used",
        "omega_used",
        "lambda_constrained_from_raw",
        "omega_constrained_from_raw",
        "lambda_at_boundary",
        "omega_at_boundary",
        "coefficient_fallback",
    )
    coefficient_diagnostics = {
        method: {key: [] for key in coefficient_keys} for method in METHODS
    }
    restrictions: dict[str, str] = {}
    ratio_usage: dict[str, float] = {}
    common_keys = (
        "tuning_a",
        "tuning_b",
        "tuning_c",
        "tuning_d",
        "tuning_h",
        "ratio_normalization_error_raw",
        "ratio_normalization_error_clipped",
        "ratio_ess",
        "ratio_max_raw",
        "ratio_clip_fraction",
        "true_r_product_bias_raw",
        "true_r_product_bias_clipped",
        "min_rct_train_arm",
        "min_obs_train_arm",
    )
    common_diagnostics = {key: [] for key in common_keys}
    for k in range(5):
        ri, oi = rf[k], of[k]
        rd = np.setdiff1d(np.arange(n), ri)
        od = np.setdiff1d(np.arange(n_obs), oi)
        dr, do, er, eo = take(rct, rd), take(obs, od), take(rct, ri), take(obs, oi)
        scores, arms = inner_scores(dr, do, stable_seed(seed, k, 20))
        moments = tuning_moments(
            scores["z0"], scores["delta"], scores["gr"], scores["rgo"], len(rd) / len(od)
        )
        coefficients = tune_all(moments)
        mr, mo = fit_outcome(dr), fit_outcome(do)
        ratio = fit_ratio(dr, do, stable_seed(seed, k, 21))
        r0, r1 = mr.predict(er["x"])
        o0r, o1r = mo.predict(er["x"])
        o0o, o1o = mo.predict(eo["x"])
        z0, z1 = pseudo(er, r0, r1), pseudo(er, o0r, o1r)
        delta, gr, go = z1 - z0, o1r - o0r, o1o - o0o
        raw = ratio.raw(eo["x"])
        clipped = np.clip(raw, *RATIO_CLIP)
        rg = clipped * go
        correction = float(rg.mean() - gr.mean())

        for method, coefficient in coefficients.items():
            lam = float(coefficient["lambda_used"])
            omega = float(coefficient["omega_used"])
            fold_estimates[method].append(
                float(np.mean(z0 + lam * delta) + omega * correction)
            )
            variance_proxies[method] += (
                var(z0 + lam * delta - omega * gr) / len(z0)
                + omega**2 * var(rg) / len(rg)
            ) / 25
            for key in coefficient_keys:
                coefficient_diagnostics[method][key].append(float(coefficient[key]))
            restrictions[method] = str(coefficient["tuning_restriction"])
            ratio_usage[method] = float(coefficient["uses_ratio_channel"])

        for key, value in moments.items():
            common_diagnostics[f"tuning_{key}"].append(float(value))
        for key, value in arms.items():
            common_diagnostics[key].append(float(value))
        common_diagnostics["ratio_normalization_error_raw"].append(float(raw.mean() - 1))
        common_diagnostics["ratio_normalization_error_clipped"].append(float(clipped.mean() - 1))
        common_diagnostics["ratio_ess"].append(float(clipped.sum() ** 2 / np.sum(clipped**2)))
        common_diagnostics["ratio_max_raw"].append(float(raw.max()))
        common_diagnostics["ratio_clip_fraction"].append(
            float(np.mean((raw < RATIO_CLIP[0]) | (raw > RATIO_CLIP[1])))
        )
        common_diagnostics["true_r_product_bias_raw"].append(
            float(np.mean((raw - eo["true_ratio"]) * go))
        )
        common_diagnostics["true_r_product_bias_clipped"].append(
            float(np.mean((clipped - eo["true_ratio"]) * go))
        )

    rows = []
    for method in METHODS:
        estimate = float(np.mean(fold_estimates[method]))
        row: dict[str, object] = {
            "benchmark": benchmark,
            "scenario_id": scenario,
            "regime": regime,
            "n_rct_total": n,
            "n_obs_total": n_obs,
            "outer_replication": rep,
            "paired_design_id": paired_design_id,
            "estimator": method,
            "tuning_restriction": restrictions[method],
            "uses_ratio_channel": ratio_usage[method],
            "estimate": estimate,
            "true_theta": theta,
            "error": estimate - theta,
            "score_variance_proxy_nonexact": variance_proxies[method],
            "unique_rct_generated": n,
            "unique_obs_generated": n_obs,
            "source_sha256": source_sha256,
            "support_rows": support_rows,
        }
        for key, values in coefficient_diagnostics[method].items():
            row[f"mean_{key}"] = float(np.mean(values))
        for key, values in common_diagnostics.items():
            row[f"common_min_{key[4:]}" if key.startswith("min_") else f"common_mean_{key}"] = float(
                np.min(values) if key.startswith("min_") else np.mean(values)
            )
        rows.append(row)
    return rows


def run_synthetic(b: int) -> list[dict[str, object]]:
    rows = []
    for scm in legacy.SCMS:
        for regime in REGIMES:
            for n in N_RCT_VALUES:
                for rep in range(b):
                    rct, obs, theta = synthetic_sample(scm, regime, n, rep)
                    rows += evaluate("arbitrary_scm", scm.scm_id, regime, n, rep, rct, obs, theta)
    return rows


def run_star(b: int, support: StarSupport) -> list[dict[str, object]]:
    rows = []
    for regime in REGIMES:
        for n in N_RCT_VALUES:
            for rep in range(b):
                rct, obs, theta = star_sample(support, regime, n, rep)
                rows += evaluate("star_real_x", 0, regime, n, rep, rct, obs, theta, STAR_SHA256, len(support.x))
    return rows


def summarize(
    frame: pd.DataFrame,
    interpretation: str = "total_budget_nested_debugging_benchmark",
) -> pd.DataFrame:
    group_keys = ["benchmark", "scenario_id", "regime", "n_rct_total"]
    diagnostic_columns = [
        "mean_raw_lambda_used", "mean_raw_omega_used", "mean_lambda_used",
        "mean_omega_used",
        "mean_lambda_constrained_from_raw", "mean_omega_constrained_from_raw",
        "mean_lambda_at_boundary", "mean_omega_at_boundary",
        "mean_coefficient_fallback",
        "common_mean_tuning_a", "common_mean_tuning_b", "common_mean_tuning_c",
        "common_mean_tuning_d", "common_mean_tuning_h",
        "common_mean_ratio_normalization_error_raw",
        "common_mean_ratio_normalization_error_clipped", "common_mean_ratio_ess",
        "common_mean_ratio_max_raw", "common_mean_ratio_clip_fraction",
        "common_mean_true_r_product_bias_raw",
        "common_mean_true_r_product_bias_clipped", "common_min_rct_train_arm",
        "common_min_obs_train_arm",
        "score_variance_proxy_nonexact",
    ]
    rows: list[dict[str, object]] = []
    for cell_key, cell in frame.groupby(group_keys, sort=True):
        estimates = {
            method: group.sort_values("outer_replication")["estimate"].to_numpy()
            for method, group in cell.groupby("estimator")
        }
        base = estimates["rct_aipw"]
        full = estimates["full_datafusionppi"]
        base_variance = var(base)
        full_variance = var(full)
        for method, group in cell.groupby("estimator", sort=True):
            group = group.sort_values("outer_replication")
            estimate = group["estimate"].to_numpy()
            theta = float(group["true_theta"].iloc[0])
            empirical_variance = var(estimate)
            bias = float(estimate.mean() - theta)
            error = estimate - theta
            base_error = base - theta
            full_error = full - theta
            paired_base = estimate - base
            paired_full = estimate - full
            squared_error_base = error**2 - base_error**2
            squared_error_full = error**2 - full_error**2
            b = len(estimate)
            paired_variance_contribution = b / (b - 1) * (
                (full - full.mean()) ** 2 - (estimate - estimate.mean()) ** 2
            )
            rmse = float(math.sqrt(np.mean(error**2)))
            base_rmse = float(math.sqrt(np.mean(base_error**2)))
            full_rmse = float(math.sqrt(np.mean(full_error**2)))
            row: dict[str, object] = dict(zip(group_keys, cell_key))
            row.update({
                "estimator": method,
                "tuning_restriction": str(group["tuning_restriction"].iloc[0]),
                "uses_ratio_channel": float(group["uses_ratio_channel"].iloc[0]),
                "replications": b,
                "true_theta": theta,
                "mean_estimate": float(estimate.mean()),
                "bias": bias,
                "bias_mcse": float(np.std(estimate, ddof=1) / math.sqrt(b)),
                "empirical_variance": empirical_variance,
                "empirical_variance_mcse_normal": float(empirical_variance * math.sqrt(2.0 / (b - 1))),
                "rmse": rmse,
                "variance_ratio_to_aipw": float(empirical_variance / base_variance) if base_variance > 1e-14 else 0.0,
                "variance_ratio_to_full": float(empirical_variance / full_variance) if full_variance > 1e-14 else 0.0,
                "rmse_ratio_to_aipw": float(rmse / base_rmse) if base_rmse > 1e-14 else 0.0,
                "rmse_ratio_to_full": float(rmse / full_rmse) if full_rmse > 1e-14 else 0.0,
                "paired_mean_difference_vs_aipw": float(paired_base.mean()),
                "paired_mean_difference_mcse_vs_aipw": float(np.std(paired_base, ddof=1) / math.sqrt(b)),
                "paired_mean_difference_vs_full": float(paired_full.mean()),
                "paired_mean_difference_mcse_vs_full": float(np.std(paired_full, ddof=1) / math.sqrt(b)),
                "paired_squared_error_difference_vs_aipw": float(squared_error_base.mean()),
                "paired_squared_error_difference_mcse_vs_aipw": float(np.std(squared_error_base, ddof=1) / math.sqrt(b)),
                "paired_squared_error_difference_vs_full": float(squared_error_full.mean()),
                "paired_squared_error_difference_mcse_vs_full": float(np.std(squared_error_full, ddof=1) / math.sqrt(b)),
                "paired_variance_difference_full_minus_method": float(paired_variance_contribution.mean()),
                "paired_variance_difference_mcse_full_minus_method": float(np.std(paired_variance_contribution, ddof=1) / math.sqrt(b)),
                "source_sha256": str(group["source_sha256"].iloc[0]),
                "support_rows": int(group["support_rows"].iloc[0]),
                "interpretation": interpretation,
            })
            for column in diagnostic_columns:
                row[f"average_{column}"] = float(group[column].mean())
            rows.append(row)
    return pd.DataFrame(rows)


def validate_replications(frame: pd.DataFrame, b: int, include_star: bool = True) -> None:
    scenario_count = 4 if include_star else 3
    expected_rows = scenario_count * len(REGIMES) * len(N_RCT_VALUES) * b * len(METHODS)
    if len(frame) != expected_rows:
        raise AssertionError(f"replication shape mismatch: {frame.shape}")
    if set(frame["estimator"]) != set(METHODS):
        raise AssertionError("primary method set mismatch")
    counts = frame.groupby(["benchmark", "scenario_id", "regime", "n_rct_total", "estimator"]).size()
    if not (counts == b).all():
        raise AssertionError("dropped or duplicated replication")
    design_counts = frame.groupby(
        ["benchmark", "scenario_id", "regime", "n_rct_total", "outer_replication"]
    )["paired_design_id"].nunique()
    if not (design_counts == 1).all():
        raise AssertionError("methods do not share the same paired design")
    restrictions = {
        "rct_aipw": "fixed_lambda_0_omega_0",
        "haipw_only": "restricted_omega_0",
        "ppi_power_tuned_only": "restricted_lambda_0",
        "full_datafusionppi": "joint_box_lambda_omega",
    }
    for method, restriction in restrictions.items():
        if set(frame.loc[frame.estimator == method, "tuning_restriction"]) != {restriction}:
            raise AssertionError(f"{method} tuning restriction mismatch")
    if not (frame.loc[frame.estimator == "rct_aipw", ["mean_lambda_used", "mean_omega_used"]] == 0).all().all():
        raise AssertionError("AIPW coefficient restriction failed")
    if not (frame.loc[frame.estimator == "haipw_only", "mean_omega_used"] == 0).all():
        raise AssertionError("lambda-only outcome-model fusion omega restriction failed")
    if not (frame.loc[frame.estimator == "ppi_power_tuned_only", "mean_lambda_used"] == 0).all():
        raise AssertionError("omega-only PPI correction lambda restriction failed")
    expected_obs = frame["n_rct_total"].map(obs_size_for)
    if not (frame["n_obs_total"] == expected_obs).all() or not (frame["unique_obs_generated"] == expected_obs).all():
        raise AssertionError("OBS total-budget accounting failed")
    if not (frame["n_rct_total"] == frame["unique_rct_generated"]).all():
        raise AssertionError("RCT total-budget accounting failed")
    if frame.isna().any().any() or not np.isfinite(frame.select_dtypes(include=[np.number]).to_numpy()).all():
        raise AssertionError("missing or non-finite replication output")


def validate(frame: pd.DataFrame, summary: pd.DataFrame, b: int, include_star: bool = True) -> None:
    validate_replications(frame, b, include_star)
    scenario_count = 4 if include_star else 3
    expected_summary = scenario_count * len(REGIMES) * len(N_RCT_VALUES) * len(METHODS)
    if len(summary) != expected_summary:
        raise AssertionError(f"summary shape mismatch: {summary.shape}")
    if summary.isna().any().any() or not np.isfinite(summary.select_dtypes(include=[np.number]).to_numpy()).all():
        raise AssertionError("missing or non-finite summary output")
    aipw = summary[summary.estimator == "rct_aipw"]
    full = summary[summary.estimator == "full_datafusionppi"]
    if not np.allclose(aipw["variance_ratio_to_aipw"], 1.0) or not np.allclose(aipw["rmse_ratio_to_aipw"], 1.0):
        raise AssertionError("AIPW self-comparison failed")
    if not np.allclose(full["variance_ratio_to_full"], 1.0) or not np.allclose(full["rmse_ratio_to_full"], 1.0):
        raise AssertionError("Full self-comparison failed")
    if not np.allclose(full["paired_mean_difference_vs_full"], 0.0) or not np.allclose(full["paired_variance_difference_full_minus_method"], 0.0):
        raise AssertionError("Full paired self-comparison failed")
    full_variance_by_cell = full.set_index(
        ["benchmark", "scenario_id", "regime", "n_rct_total"]
    )["empirical_variance"]
    expected_difference = np.array(
        [
            full_variance_by_cell.loc[(row.benchmark, row.scenario_id, row.regime, row.n_rct_total)]
            - row.empirical_variance
            for row in summary.itertuples()
        ]
    )
    if not np.allclose(
        summary["paired_variance_difference_full_minus_method"], expected_difference
    ):
        raise AssertionError("Paired variance-difference diagnostic is inconsistent")


def make_figure(summary: pd.DataFrame, figure_path: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(22, 9), constrained_layout=True)
    panels = [
        ("arbitrary_scm", "Synthetic realisations"),
        ("star_real_x", "STAR real-X semi-synthetic"),
    ]
    colors = {
        "rct_aipw": "#4C78A8",
        "haipw_only": "#F58518",
        "ppi_power_tuned_only": "#54A24B",
        "full_datafusionppi": "#E45756",
    }
    labels = {
        "rct_aipw": "RCT-AIPW",
        "haipw_only": r"$\lambda$-only outcome-model fusion",
        "ppi_power_tuned_only": r"$\omega$-only PPI correction",
        "full_datafusionppi": "Full DataFusionPPI",
    }
    for row_index, (benchmark, title) in enumerate(panels):
        data = summary[summary["benchmark"] == benchmark].copy()
        data["regime_order"] = data["regime"].map({"shared": 0, "shifted": 1})
        cells = data[["scenario_id", "regime", "regime_order", "n_rct_total"]].drop_duplicates().sort_values(
            ["scenario_id", "regime_order", "n_rct_total"]
        )
        cells["cell_position"] = np.arange(len(cells))
        cells["cell_label"] = cells.apply(
            lambda row: (
                f"S{int(row.scenario_id)}-{'same' if row.regime == 'shared' else 'shift'}-{int(row.n_rct_total)}"
                if benchmark == "arbitrary_scm"
                else f"{'same' if row.regime == 'shared' else 'shift'}-{int(row.n_rct_total)}"
            ),
            axis=1,
        )
        data = data.merge(cells, on=["scenario_id", "regime", "regime_order", "n_rct_total"])
        bias_axis, variance_axis, rmse_axis = axes[row_index]
        for method in METHODS:
            group = data[data.estimator == method].sort_values("cell_position")
            x = group["cell_position"].to_numpy(dtype=float)
            bias_axis.errorbar(
                x,
                group["bias"],
                yerr=2.0 * group["bias_mcse"],
                marker="o",
                color=colors[method],
                linestyle="-",
                capsize=2,
                alpha=0.85,
                label=labels[method],
            )
            variance_axis.plot(
                x, group["variance_ratio_to_aipw"], marker="o", color=colors[method], label=labels[method]
            )
            rmse_axis.plot(
                x, group["rmse_ratio_to_aipw"], marker="o", color=colors[method], label=labels[method]
            )
        bias_axis.axhline(0, color="black", linewidth=0.8)
        variance_axis.axhline(1, color="black", linewidth=0.8)
        rmse_axis.axhline(1, color="black", linewidth=0.8)
        bias_axis.set_title(f"{title}: bias ± 2 MCSE")
        variance_axis.set_title(f"{title}: variance / AIPW variance")
        rmse_axis.set_title(f"{title}: RMSE / AIPW RMSE")
        for axis in (bias_axis, variance_axis, rmse_axis):
            axis.set_xticks(cells["cell_position"], cells["cell_label"], rotation=45, ha="right")
            axis.grid(alpha=0.2)
        bias_axis.set_ylabel("bias")
        variance_axis.set_ylabel("variance ratio")
        rmse_axis.set_ylabel("RMSE ratio")
    handles, legend_labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="outside upper center", ncol=4, fontsize=9)
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)


def write_outputs(
    frame: pd.DataFrame,
    summary: pd.DataFrame,
    replication_path: Path,
    summary_path: Path,
    figure_path: Path,
    overwrite: bool,
) -> None:
    existing = [path for path in (replication_path, summary_path, figure_path) if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(f"refusing to overwrite existing outputs: {existing}")
    frame.to_csv(replication_path, index=False, float_format="%.12g")
    summary.to_csv(summary_path, index=False, float_format="%.12g")
    make_figure(summary, figure_path)


def whi_check(data_dir: Path | None, crosswalk: Path | None) -> int:
    if data_dir is None or crosswalk is None or not data_dir.is_dir() or not crosswalk.is_file():
        print("WHI_BLOCKED: authorized BioLINCC directory and semantic crosswalk required")
        return 3
    mapping = json.loads(crosswalk.read_text())
    required = {"study_membership", "uterus_status", "ep_trial_arm", "index_date", "chd_event", "chd_event_date", "last_followup_date"}
    files = sorted(data_dir.rglob("*.csv"))
    inventory: dict[str, list[str]] = {}
    for path in files:
        relative = str(path.relative_to(data_dir))
        try:
            inventory[relative] = list(pd.read_csv(path, nrows=0).columns)
        except Exception as error:
            inventory[relative] = [f"READ_ERROR:{type(error).__name__}"]
    invalid = []
    for semantic, specification in mapping.items():
        if not isinstance(specification, dict):
            invalid.append(semantic)
            continue
        file_name, column = specification.get("file"), specification.get("column")
        if file_name not in inventory or column not in inventory.get(file_name, []):
            invalid.append(semantic)
    missing = sorted(required - set(mapping))
    status = "blocked" if missing or invalid else "schema_review_required"
    print(json.dumps({"status": status, "missing_semantics": missing, "invalid_mappings": sorted(invalid), "csv_file_count": len(files), "inventory": inventory}, sort_keys=True))
    return 3 if status == "blocked" else 2


def main() -> int:
    global N_RCT_VALUES, OBS_MULTIPLIER
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replications", type=int, default=20)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--n-rct-values", default=",".join(map(str, N_RCT_VALUES)))
    parser.add_argument("--obs-multiplier", type=int)
    parser.add_argument("--output-stem", default="total_budget_nested_benchmark")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--whi-check-dir", type=Path); parser.add_argument("--whi-crosswalk", type=Path)
    args = parser.parse_args()
    if args.whi_check_dir is not None or args.whi_crosswalk is not None: return whi_check(args.whi_check_dir, args.whi_crosswalk)
    N_RCT_VALUES = parse_n_rct_values(args.n_rct_values)
    if args.obs_multiplier is not None and args.obs_multiplier <= 0:
        parser.error("--obs-multiplier must be positive")
    OBS_MULTIPLIER = args.obs_multiplier
    replication_path, summary_path, figure_path = output_paths(args.output_stem)
    b = 1 if args.smoke else args.replications
    if b <= 0:
        parser.error("--replications must be positive")
    started = time.perf_counter()
    rows = run_synthetic(b)
    star_status = "complete"
    try:
        support = load_star_support()
        rows += run_star(b, support)
    except Exception as error:
        star_status = f"blocked:{type(error).__name__}:{error}"
        raise RuntimeError(star_status) from error
    frame = pd.DataFrame(rows)
    interpretation = (
        "proportional_budget_nested_debugging_benchmark"
        if OBS_MULTIPLIER is not None
        else "total_budget_nested_debugging_benchmark"
    )
    summary = summarize(frame, interpretation) if b > 1 else pd.DataFrame()
    if b == 1:
        validate_replications(frame, b, include_star=True)
    else:
        validate(frame, summary, b, include_star=star_status == "complete")
        write_outputs(
            frame,
            summary,
            replication_path,
            summary_path,
            figure_path,
            args.overwrite,
        )
    print(json.dumps({
        "phase": "combined",
        "rows": len(frame),
        "summary_rows": len(summary),
        "replications": b,
        "n_rct_values": N_RCT_VALUES,
        "obs_multiplier": OBS_MULTIPLIER,
        "output_stem": args.output_stem,
        "star_status": star_status,
        "runtime_seconds": round(time.perf_counter() - started, 3),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
