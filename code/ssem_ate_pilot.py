#!/usr/bin/env python3
"""Debugging pilot for feasible prediction-powered RCT/OBS ATE fusion.

Three illustrative nonlinear SCMs are drawn once from a fixed random generator.
For every SCM, regime, RCT size, and outer replication, independent training,
calibration, and evaluation samples are regenerated.  RCT sizes 20 and 30 and
OBS size 5000 are per-stage sizes: one outer replication therefore generates
3*n RCT and 15000 OBS observations.  This is not a total-budget comparison,
not a random-SCM population study, and not theorem verification.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


SEED = 190602
SCM_COUNT = 3
N_RCT_VALUES = (20, 30)
N_OBS = 5000
REGIMES = ("shared", "shifted")
N_OUTER_REPLICATIONS = 20
OUTCOME_RIDGE_ALPHA = 5.0
RATIO_LOGISTIC_C = 1.0
COEFFICIENT_BOX = (-2.0, 2.0)
RATIO_CLIP = (0.05, 20.0)
DENOMINATOR_TOLERANCE = 1e-10
CI_Z = 1.959963984540054

PROJECT_DIR = Path(__file__).resolve().parents[1]
MATERIALS_DIR = PROJECT_DIR / "materials"
REPLICATION_PATH = MATERIALS_DIR / "ssem_ate_pilot_replications.csv"
SUMMARY_PATH = MATERIALS_DIR / "ssem_ate_pilot_summary.csv"
FIGURE_PATH = MATERIALS_DIR / "ssem_ate_pilot.png"


@dataclass(frozen=True)
class SCMParameters:
    scm_id: int
    covariate_mean: np.ndarray
    latent_loading: np.ndarray
    covariate_noise_sd: np.ndarray
    rct_mean_shift: np.ndarray
    baseline_intercept: float
    baseline_coefficients: np.ndarray
    outcome_latent_coefficient: float
    treatment_intercept: float
    treatment_coefficients: np.ndarray
    outcome_noise_sd: float
    rct_propensity_intercept: float
    rct_propensity_coefficients: np.ndarray
    obs_propensity_intercept: float
    obs_propensity_coefficients: np.ndarray
    obs_latent_coefficient: float

    @property
    def covariance(self) -> np.ndarray:
        return np.outer(self.latent_loading, self.latent_loading) + np.diag(
            self.covariate_noise_sd**2
        )

    def rct_mean(self, regime: str) -> np.ndarray:
        if regime == "shared":
            return self.covariate_mean
        if regime == "shifted":
            return self.covariate_mean + self.rct_mean_shift
        raise ValueError(f"Unknown regime: {regime}")


@dataclass(frozen=True)
class JointOutcomeModel:
    feature_mean: np.ndarray
    feature_sd: np.ndarray
    baseline_coefficients: np.ndarray
    treatment_coefficients: np.ndarray

    def design(self, x: np.ndarray) -> np.ndarray:
        raw = generic_outcome_features(x)
        standardized = (raw - self.feature_mean) / self.feature_sd
        return np.column_stack([np.ones(x.shape[0]), standardized])

    def predict(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        design = self.design(x)
        mu0 = design @ self.baseline_coefficients
        mu1 = mu0 + design @ self.treatment_coefficients
        return mu0, mu1


@dataclass(frozen=True)
class EstimatedRatioModel:
    classifier: Pipeline
    rct_prior: float

    def predict(self, x: np.ndarray) -> np.ndarray:
        probability = self.classifier.predict_proba(x)[:, 1]
        probability = np.clip(probability, 1e-8, 1.0 - 1e-8)
        odds = probability / (1.0 - probability)
        prior_correction = (1.0 - self.rct_prior) / self.rct_prior
        return np.clip(odds * prior_correction, *RATIO_CLIP)


def expit(value: np.ndarray) -> np.ndarray:
    value = np.clip(value, -35.0, 35.0)
    return 1.0 / (1.0 + np.exp(-value))


def random_sign(rng: np.random.Generator) -> float:
    return float(rng.choice(np.array([-1.0, 1.0])))


def generate_scm_parameters(scm_id: int) -> SCMParameters:
    rng = np.random.default_rng(np.random.SeedSequence([SEED, 7000, scm_id]))
    covariate_mean = rng.normal(0.0, 0.15, 4)
    latent_loading = rng.normal(0.0, 0.30, 4)
    covariate_noise_sd = rng.uniform(0.75, 1.05, 4)
    rct_mean_shift = rng.normal(0.0, 0.18, 4)
    baseline_intercept = float(rng.normal(0.0, 0.2))
    baseline_coefficients = rng.normal(0.0, 0.30, 9)
    outcome_latent_coefficient = random_sign(rng) * float(rng.uniform(0.4, 0.8))
    treatment_intercept = float(rng.uniform(0.8, 1.2))
    treatment_coefficients = rng.normal(0.0, 0.22, 5)
    outcome_noise_sd = float(rng.uniform(0.8, 1.2))
    rct_propensity_intercept = float(rng.normal(0.0, 0.1))
    rct_propensity_coefficients = rng.normal(0.0, 0.18, 5)
    obs_propensity_intercept = float(rng.normal(-0.2, 0.1))
    obs_propensity_coefficients = rng.normal(0.0, 0.30, 5)
    obs_latent_coefficient = random_sign(rng) * float(rng.uniform(0.7, 1.1))
    return SCMParameters(
        scm_id=scm_id,
        covariate_mean=covariate_mean,
        latent_loading=latent_loading,
        covariate_noise_sd=covariate_noise_sd,
        rct_mean_shift=rct_mean_shift,
        baseline_intercept=baseline_intercept,
        baseline_coefficients=baseline_coefficients,
        outcome_latent_coefficient=outcome_latent_coefficient,
        treatment_intercept=treatment_intercept,
        treatment_coefficients=treatment_coefficients,
        outcome_noise_sd=outcome_noise_sd,
        rct_propensity_intercept=rct_propensity_intercept,
        rct_propensity_coefficients=rct_propensity_coefficients,
        obs_propensity_intercept=obs_propensity_intercept,
        obs_propensity_coefficients=obs_propensity_coefficients,
        obs_latent_coefficient=obs_latent_coefficient,
    )


SCMS = tuple(generate_scm_parameters(scm_id) for scm_id in range(SCM_COUNT))


def baseline_basis(x: np.ndarray) -> np.ndarray:
    return np.column_stack(
        [
            x[:, 0],
            x[:, 1],
            x[:, 2],
            x[:, 3],
            np.sin(x[:, 0]),
            np.cos(x[:, 1]),
            x[:, 0] * x[:, 1],
            x[:, 2] * x[:, 3],
            x[:, 0] ** 2,
        ]
    )


def treatment_effect_basis(x: np.ndarray) -> np.ndarray:
    return np.column_stack(
        [
            x[:, 0],
            x[:, 1],
            np.sin(x[:, 2]),
            x[:, 0] * x[:, 3],
            x[:, 1] ** 2,
        ]
    )


def treatment_assignment_basis(x: np.ndarray) -> np.ndarray:
    return np.column_stack(
        [
            x[:, 0],
            x[:, 1],
            np.tanh(x[:, 2]),
            np.sin(x[:, 3]),
            x[:, 0] * x[:, 1],
        ]
    )


def generic_outcome_features(x: np.ndarray) -> np.ndarray:
    """Generic 22-feature map, deliberately not the DGP basis itself."""
    linear = [x[:, column] for column in range(4)]
    squares = [x[:, column] ** 2 for column in range(4)]
    interactions = [
        x[:, first] * x[:, second]
        for first in range(4)
        for second in range(first + 1, 4)
    ]
    sine = [np.sin(x[:, column]) for column in range(4)]
    cosine = [np.cos(x[:, column]) for column in range(4)]
    return np.column_stack(linear + squares + interactions + sine + cosine)


def tau_function(x: np.ndarray, scm: SCMParameters) -> np.ndarray:
    return scm.treatment_intercept + treatment_effect_basis(
        x
    ) @ scm.treatment_coefficients


def true_theta(scm: SCMParameters, regime: str) -> float:
    mean = scm.rct_mean(regime)
    covariance = scm.covariance
    coefficient = scm.treatment_coefficients
    expectation = np.array(
        [
            mean[0],
            mean[1],
            math.exp(-covariance[2, 2] / 2.0) * math.sin(mean[2]),
            covariance[0, 3] + mean[0] * mean[3],
            covariance[1, 1] + mean[1] ** 2,
        ]
    )
    return float(scm.treatment_intercept + coefficient @ expectation)


def rct_propensity(x: np.ndarray, scm: SCMParameters) -> np.ndarray:
    linear_predictor = (
        scm.rct_propensity_intercept
        + treatment_assignment_basis(x) @ scm.rct_propensity_coefficients
    )
    return np.clip(expit(linear_predictor), 0.20, 0.80)


def obs_propensity(
    x: np.ndarray, u: np.ndarray, scm: SCMParameters
) -> np.ndarray:
    linear_predictor = (
        scm.obs_propensity_intercept
        + treatment_assignment_basis(x) @ scm.obs_propensity_coefficients
        + scm.obs_latent_coefficient * u
    )
    return np.clip(expit(linear_predictor), 0.05, 0.95)


def seed_rng(*components: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([SEED, *components]))


def regime_index(regime: str) -> int:
    return REGIMES.index(regime)


def sample_rng(
    scm_id: int,
    regime: str,
    n_rct: int,
    outer_replication: int,
    stage_code: int,
    study_code: int,
) -> np.random.Generator:
    return seed_rng(
        scm_id,
        regime_index(regime),
        n_rct,
        outer_replication,
        stage_code,
        study_code,
    )


def sample_dataset(
    rng: np.random.Generator,
    sample_size: int,
    study: str,
    regime: str,
    scm: SCMParameters,
) -> dict[str, np.ndarray]:
    u = rng.normal(size=sample_size)
    mean = scm.covariate_mean if study == "OBS" else scm.rct_mean(regime)
    x = (
        mean
        + np.outer(u, scm.latent_loading)
        + rng.normal(size=(sample_size, 4)) * scm.covariate_noise_sd
    )
    if study == "RCT":
        propensity = rct_propensity(x, scm)
    elif study == "OBS":
        propensity = obs_propensity(x, u, scm)
    else:
        raise ValueError(f"Unknown study: {study}")
    treatment = rng.binomial(1, propensity, size=sample_size).astype(float)
    tau = tau_function(x, scm)
    baseline = (
        scm.baseline_intercept
        + baseline_basis(x) @ scm.baseline_coefficients
        + scm.outcome_latent_coefficient * u
    )
    outcome = (
        baseline
        + treatment * tau
        + scm.outcome_noise_sd * rng.normal(size=sample_size)
    )
    return {
        "x": x,
        "a": treatment,
        "y": outcome,
        "propensity": propensity,
        "tau": tau,
    }


def fit_joint_outcome_model(
    data: dict[str, np.ndarray], minimum_arm_count: int
) -> JointOutcomeModel:
    treatment = data["a"]
    arm1 = int(treatment.sum())
    arm0 = int(len(treatment) - arm1)
    if min(arm0, arm1) < minimum_arm_count:
        raise RuntimeError(
            f"Outcome-model arm gate failed: arm0={arm0}, arm1={arm1}, "
            f"minimum={minimum_arm_count}"
        )
    raw = generic_outcome_features(data["x"])
    feature_mean = raw.mean(axis=0)
    feature_sd = raw.std(axis=0, ddof=0)
    feature_sd = np.where(feature_sd > 1e-8, feature_sd, 1.0)
    standardized = (raw - feature_mean) / feature_sd
    base = np.column_stack([np.ones(len(treatment)), standardized])
    joint = np.column_stack([base, treatment[:, None] * base])
    penalty = np.eye(joint.shape[1]) * OUTCOME_RIDGE_ALPHA
    penalty[0, 0] = 0.0
    system = joint.T @ joint + penalty
    target = joint.T @ data["y"]
    try:
        coefficients = np.linalg.solve(system, target)
    except np.linalg.LinAlgError:
        coefficients = np.linalg.lstsq(system, target, rcond=None)[0]
    feature_count = base.shape[1]
    return JointOutcomeModel(
        feature_mean=feature_mean,
        feature_sd=feature_sd,
        baseline_coefficients=coefficients[:feature_count],
        treatment_coefficients=coefficients[feature_count:],
    )


def fit_estimated_ratio(
    rct_data: dict[str, np.ndarray],
    obs_data: dict[str, np.ndarray],
    random_state: int,
) -> EstimatedRatioModel:
    x = np.vstack([rct_data["x"], obs_data["x"]])
    label = np.concatenate(
        [np.ones(len(rct_data["x"])), np.zeros(len(obs_data["x"]))]
    )
    prior = float(label.mean())
    classifier = Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "logit",
                LogisticRegression(
                    C=RATIO_LOGISTIC_C,
                    max_iter=2000,
                    solver="lbfgs",
                    random_state=random_state,
                ),
            ),
        ]
    )
    classifier.fit(x, label)
    logit = classifier.named_steps["logit"]
    if logit.coef_.shape != (1, 4):
        raise AssertionError("Ratio model is not the prespecified 4-variable logit")
    if int(logit.n_iter_[0]) >= int(logit.max_iter):
        raise RuntimeError("Density-ratio logistic regression did not converge")
    return EstimatedRatioModel(classifier=classifier, rct_prior=prior)


def gaussian_log_density(
    x: np.ndarray, mean: np.ndarray, covariance: np.ndarray
) -> np.ndarray:
    sign, log_determinant = np.linalg.slogdet(covariance)
    if sign <= 0:
        raise ValueError("Covariance matrix is not positive definite")
    difference = x - mean
    inverse = np.linalg.inv(covariance)
    quadratic = np.einsum("ni,ij,nj->n", difference, inverse, difference)
    return -0.5 * (
        x.shape[1] * math.log(2.0 * math.pi) + log_determinant + quadratic
    )


def oracle_ratio(x: np.ndarray, regime: str, scm: SCMParameters) -> np.ndarray:
    if regime == "shared":
        return np.ones(x.shape[0])
    covariance = scm.covariance
    log_ratio = gaussian_log_density(
        x, scm.rct_mean(regime), covariance
    ) - gaussian_log_density(x, scm.covariate_mean, covariance)
    return np.exp(log_ratio)


def pseudo_outcome(
    rct_data: dict[str, np.ndarray], mu0: np.ndarray, mu1: np.ndarray
) -> np.ndarray:
    treatment = rct_data["a"]
    outcome = rct_data["y"]
    propensity = rct_data["propensity"]
    return (
        mu1
        - mu0
        + treatment / propensity * (outcome - mu1)
        - (1.0 - treatment) / (1.0 - propensity) * (outcome - mu0)
    )


def sample_variance(value: np.ndarray) -> float:
    return float(np.var(value, ddof=1))


def sample_covariance(first: np.ndarray, second: np.ndarray) -> float:
    return float(np.cov(first, second, ddof=1)[0, 1])


def safe_correlation(first: np.ndarray, second: np.ndarray) -> float:
    denominator = math.sqrt(sample_variance(first) * sample_variance(second))
    if denominator <= DENOMINATOR_TOLERANCE:
        return 0.0
    return sample_covariance(first, second) / denominator


def clipped_ratio(numerator: float, denominator: float) -> float:
    if denominator <= DENOMINATOR_TOLERANCE:
        return 0.0
    return float(np.clip(numerator / denominator, *COEFFICIENT_BOX))


def weight_ess(weight: np.ndarray) -> float:
    return float(weight.sum() ** 2 / np.sum(weight**2))


def prepare_cell(
    scm: SCMParameters,
    regime: str,
    n_rct: int,
    outer_replication: int,
) -> dict[str, object]:
    train_rct = sample_dataset(
        sample_rng(scm.scm_id, regime, n_rct, outer_replication, 10, 1),
        n_rct,
        "RCT",
        regime,
        scm,
    )
    train_obs = sample_dataset(
        sample_rng(scm.scm_id, regime, n_rct, outer_replication, 10, 0),
        N_OBS,
        "OBS",
        regime,
        scm,
    )
    calibration_rct = sample_dataset(
        sample_rng(scm.scm_id, regime, n_rct, outer_replication, 20, 1),
        n_rct,
        "RCT",
        regime,
        scm,
    )
    calibration_obs = sample_dataset(
        sample_rng(scm.scm_id, regime, n_rct, outer_replication, 20, 0),
        N_OBS,
        "OBS",
        regime,
        scm,
    )

    rct_model = fit_joint_outcome_model(train_rct, minimum_arm_count=3)
    obs_model = fit_joint_outcome_model(train_obs, minimum_arm_count=30)
    ratio_model = fit_estimated_ratio(
        train_rct,
        train_obs,
        random_state=(
            SEED
            + 100_000 * scm.scm_id
            + 10_000 * regime_index(regime)
            + 100 * outer_replication
            + n_rct
        ),
    )

    rct_mu0_r, rct_mu1_r = rct_model.predict(calibration_rct["x"])
    obs_mu0_r, obs_mu1_r = obs_model.predict(calibration_rct["x"])
    obs_mu0_o, obs_mu1_o = obs_model.predict(calibration_obs["x"])
    z0 = pseudo_outcome(calibration_rct, rct_mu0_r, rct_mu1_r)
    z1 = pseudo_outcome(calibration_rct, obs_mu0_r, obs_mu1_r)
    delta = z1 - z0
    g_rct = obs_mu1_r - obs_mu0_r
    g_obs = obs_mu1_o - obs_mu0_o

    a_hat = sample_variance(delta)
    c_hat = sample_covariance(z0, delta)
    d_hat = sample_covariance(z0, g_rct)
    lambda_hat = clipped_ratio(-c_hat, a_hat)
    estimated_weight = ratio_model.predict(calibration_obs["x"])
    b_hat = sample_variance(g_rct) + (n_rct / N_OBS) * sample_variance(
        estimated_weight * g_obs
    )
    omega_hat = clipped_ratio(d_hat, b_hat)
    oracle_weight = oracle_ratio(calibration_obs["x"], regime, scm)
    ratio_bias_functional = float(
        np.mean((estimated_weight - oracle_weight) * g_obs)
    )

    nominal_lambda_gain = (
        c_hat**2 / (n_rct * a_hat)
        if a_hat > DENOMINATOR_TOLERANCE
        else 0.0
    )
    nominal_omega_gain = (
        d_hat**2 / (n_rct * b_hat)
        if b_hat > DENOMINATOR_TOLERANCE
        else 0.0
    )
    clipped_lambda_gain = max(
        -(2.0 * lambda_hat * c_hat + lambda_hat**2 * a_hat) / n_rct,
        0.0,
    )
    clipped_omega_gain = max(
        (2.0 * omega_hat * d_hat - omega_hat**2 * b_hat) / n_rct,
        0.0,
    )
    train_arm1 = int(train_rct["a"].sum())
    train_arm0 = int(n_rct - train_arm1)
    return {
        "rct_model": rct_model,
        "obs_model": obs_model,
        "ratio_model": ratio_model,
        "lambda": lambda_hat,
        "omega": omega_hat,
        "a_hat": a_hat,
        "b_hat": b_hat,
        "c_hat": c_hat,
        "d_hat": d_hat,
        "cross_covariance": sample_covariance(delta, g_rct),
        "calibration_ratio_bias_functional": ratio_bias_functional,
        "nominal_unconstrained_lambda_gain": nominal_lambda_gain,
        "nominal_unconstrained_omega_gain": nominal_omega_gain,
        "predicted_clipped_lambda_gain": clipped_lambda_gain,
        "predicted_clipped_omega_gain": clipped_omega_gain,
        "g_tau_correlation": safe_correlation(g_rct, calibration_rct["tau"]),
        "between_x_share": sample_variance(calibration_rct["tau"])
        / max(sample_variance(z0), DENOMINATOR_TOLERANCE),
        "train_rct_arm0": train_arm0,
        "train_rct_arm1": train_arm1,
    }


def method_specifications(cell: dict[str, object]) -> list[dict[str, object]]:
    return [
        {
            "estimator": "rct_aipw",
            "lambda": 0.0,
            "omega": 0.0,
            "ratio": "neutral",
            "variance_scope": "evaluation_variance_given_rct_nuisance_fit",
        },
        {
            "estimator": "estimated_r_datafusion",
            "lambda": float(cell["lambda"]),
            "omega": float(cell["omega"]),
            "ratio": "estimated",
            "variance_scope": "evaluation_variance_excludes_fit_tuning_uncertainty",
        },
    ]


def evaluate_cell(
    scm: SCMParameters,
    regime: str,
    n_rct: int,
    n_outer_replications: int,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for replication in range(n_outer_replications):
        cell = prepare_cell(scm, regime, n_rct, replication)
        rct_model: JointOutcomeModel = cell["rct_model"]
        obs_model: JointOutcomeModel = cell["obs_model"]
        ratio_model: EstimatedRatioModel = cell["ratio_model"]
        evaluation_rct = sample_dataset(
            sample_rng(scm.scm_id, regime, n_rct, replication, 100, 1),
            n_rct,
            "RCT",
            regime,
            scm,
        )
        evaluation_obs = sample_dataset(
            sample_rng(scm.scm_id, regime, n_rct, replication, 100, 0),
            N_OBS,
            "OBS",
            regime,
            scm,
        )

        rct_mu0_r, rct_mu1_r = rct_model.predict(evaluation_rct["x"])
        obs_mu0_r, obs_mu1_r = obs_model.predict(evaluation_rct["x"])
        obs_mu0_o, obs_mu1_o = obs_model.predict(evaluation_obs["x"])
        z0 = pseudo_outcome(evaluation_rct, rct_mu0_r, rct_mu1_r)
        z1 = pseudo_outcome(evaluation_rct, obs_mu0_r, obs_mu1_r)
        delta = z1 - z0
        g_rct = obs_mu1_r - obs_mu0_r
        g_obs = obs_mu1_o - obs_mu0_o

        estimated_weight = ratio_model.predict(evaluation_obs["x"])
        oracle_weight = oracle_ratio(evaluation_obs["x"], regime, scm)
        ratio_difference = estimated_weight - oracle_weight
        fitted_ratio_diagnostics = {
            "fitted_ratio_mean_minus_one": float(estimated_weight.mean() - 1.0),
            "fitted_ratio_ess": weight_ess(estimated_weight),
            "fitted_ratio_max": float(estimated_weight.max()),
            "fitted_ratio_clip_fraction": float(
                np.mean(
                    np.isclose(estimated_weight, RATIO_CLIP[0])
                    | np.isclose(estimated_weight, RATIO_CLIP[1])
                )
            ),
            "fitted_ratio_l2_oracle": float(
                math.sqrt(np.mean(ratio_difference**2))
            ),
            "fitted_ratio_bias_functional": float(
                np.mean(ratio_difference * g_obs)
            ),
        }
        evaluation_ratio = {
            "neutral": np.ones(N_OBS),
            "estimated": estimated_weight,
        }

        for method in method_specifications(cell):
            ratio_name = str(method["ratio"])
            weight = evaluation_ratio[ratio_name]
            lambda_value = float(method["lambda"])
            omega_value = float(method["omega"])
            z_lambda = z0 + lambda_value * delta
            estimate = float(
                np.mean(z_lambda)
                + omega_value * (np.mean(weight * g_obs) - np.mean(g_rct))
            )
            estimated_variance = (
                sample_variance(z_lambda - omega_value * g_rct) / n_rct
                + omega_value**2 * sample_variance(weight * g_obs) / N_OBS
            )
            if estimated_variance < 0.0:
                raise AssertionError("Estimated variance is negative")
            standard_error = math.sqrt(estimated_variance)
            target = true_theta(scm, regime)
            ci_lower = estimate - CI_Z * standard_error
            ci_upper = estimate + CI_Z * standard_error
            row: dict[str, object] = {
                "scm_id": scm.scm_id,
                "regime": regime,
                "n": n_rct,
                "N": N_OBS,
                "replication": replication,
                "estimator": method["estimator"],
                "true_theta": target,
                "estimate": estimate,
                "error": estimate - target,
                "estimated_variance": estimated_variance,
                "standard_error": standard_error,
                "ci_lower": ci_lower,
                "ci_upper": ci_upper,
                "coverage": int(ci_lower <= target <= ci_upper),
                "lambda": lambda_value,
                "omega": omega_value,
                "calibration_A": float(cell["a_hat"]),
                "calibration_B": float(cell["b_hat"]),
                "calibration_C": float(cell["c_hat"]),
                "calibration_D": float(cell["d_hat"]),
                "calibration_cov_delta_g": float(cell["cross_covariance"]),
                "calibration_ratio_bias_functional": float(
                    cell["calibration_ratio_bias_functional"]
                ),
                "nominal_unconstrained_lambda_gain": float(
                    cell["nominal_unconstrained_lambda_gain"]
                ),
                "nominal_unconstrained_omega_gain": float(
                    cell["nominal_unconstrained_omega_gain"]
                ),
                "predicted_clipped_lambda_gain": float(
                    cell["predicted_clipped_lambda_gain"]
                ),
                "predicted_clipped_omega_gain": float(
                    cell["predicted_clipped_omega_gain"]
                ),
                "calibration_g_tau_correlation": float(
                    cell["g_tau_correlation"]
                ),
                "calibration_between_x_share": float(cell["between_x_share"]),
                "train_rct_arm0": int(cell["train_rct_arm0"]),
                "train_rct_arm1": int(cell["train_rct_arm1"]),
                "ratio_mode": ratio_name,
                "variance_scope": method["variance_scope"],
                "feasible_minus_oracle_ratio_channel": float(
                    omega_value * np.mean(ratio_difference * g_obs)
                ),
            }
            row.update(fitted_ratio_diagnostics)
            rows.append(row)
    return rows


def clipping_count(value: pd.Series) -> int:
    return int(
        np.isclose(value, COEFFICIENT_BOX[0]).sum()
        + np.isclose(value, COEFFICIENT_BOX[1]).sum()
    )


def build_summary(replications: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    keys = ["scm_id", "regime", "n", "estimator"]
    for group_key, group in replications.groupby(keys, sort=True):
        scm_id, regime, n_rct, estimator = group_key
        count = len(group)
        error = group["error"].to_numpy()
        estimate = group["estimate"].to_numpy()
        coverage = float(group["coverage"].mean())
        empirical_variance = float(np.var(estimate, ddof=1))
        lambda_clip_count = clipping_count(group["lambda"])
        omega_clip_count = clipping_count(group["omega"])
        records.append(
            {
                "scm_id": int(scm_id),
                "regime": regime,
                "n": int(n_rct),
                "N": int(group["N"].iloc[0]),
                "estimator": estimator,
                "replications": count,
                "true_theta": float(group["true_theta"].iloc[0]),
                "bias": float(np.mean(error)),
                "bias_mcse": float(np.std(error, ddof=1) / math.sqrt(count)),
                "empirical_variance": empirical_variance,
                "empirical_variance_mcse": empirical_variance
                * math.sqrt(2.0 / (count - 1)),
                "mean_estimated_variance": float(group["estimated_variance"].mean()),
                "rmse": float(np.sqrt(np.mean(error**2))),
                "mean_standard_error": float(group["standard_error"].mean()),
                "coverage_95": coverage,
                "coverage_mcse": math.sqrt(coverage * (1.0 - coverage) / count),
                "mean_lambda": float(group["lambda"].mean()),
                "sd_lambda": float(group["lambda"].std(ddof=1)),
                "lambda_clipping_count": lambda_clip_count,
                "lambda_clipping_rate": lambda_clip_count / count,
                "mean_omega": float(group["omega"].mean()),
                "sd_omega": float(group["omega"].std(ddof=1)),
                "omega_clipping_count": omega_clip_count,
                "omega_clipping_rate": omega_clip_count / count,
                "mean_nominal_unconstrained_lambda_gain": float(
                    group["nominal_unconstrained_lambda_gain"].mean()
                ),
                "mean_nominal_unconstrained_omega_gain": float(
                    group["nominal_unconstrained_omega_gain"].mean()
                ),
                "mean_predicted_clipped_lambda_gain": float(
                    group["predicted_clipped_lambda_gain"].mean()
                ),
                "mean_predicted_clipped_omega_gain": float(
                    group["predicted_clipped_omega_gain"].mean()
                ),
                "mean_calibration_g_tau_correlation": float(
                    group["calibration_g_tau_correlation"].mean()
                ),
                "mean_calibration_between_x_share": float(
                    group["calibration_between_x_share"].mean()
                ),
                "mean_calibration_ratio_bias_functional": float(
                    group["calibration_ratio_bias_functional"].mean()
                ),
                "mean_fitted_ratio_mean_minus_one": float(
                    group["fitted_ratio_mean_minus_one"].mean()
                ),
                "mean_fitted_ratio_ess": float(group["fitted_ratio_ess"].mean()),
                "mean_fitted_ratio_max": float(group["fitted_ratio_max"].mean()),
                "mean_fitted_ratio_clip_fraction": float(
                    group["fitted_ratio_clip_fraction"].mean()
                ),
                "mean_fitted_ratio_l2_oracle": float(
                    group["fitted_ratio_l2_oracle"].mean()
                ),
                "mean_fitted_ratio_bias_functional": float(
                    group["fitted_ratio_bias_functional"].mean()
                ),
                "mean_feasible_minus_oracle_ratio_channel": float(
                    group["feasible_minus_oracle_ratio_channel"].mean()
                ),
                "minimum_train_rct_arm_count": int(
                    min(group["train_rct_arm0"].min(), group["train_rct_arm1"].min())
                ),
                "interpretation": "illustrative_random_scm_debugging_pilot",
            }
        )
    summary = pd.DataFrame.from_records(records)
    baseline = summary.loc[
        summary["estimator"] == "rct_aipw",
        ["scm_id", "regime", "n", "empirical_variance"],
    ].rename(columns={"empirical_variance": "baseline_empirical_variance"})
    summary = summary.merge(baseline, on=["scm_id", "regime", "n"], how="left")
    summary["variance_ratio_vs_aipw"] = (
        summary["empirical_variance"] / summary["baseline_empirical_variance"]
    )
    summary = summary.drop(columns=["baseline_empirical_variance"])

    paired = replications.pivot(
        index=["scm_id", "regime", "n", "replication"],
        columns="estimator",
        values="estimate",
    ).reset_index()
    paired["paired_difference"] = (
        paired["estimated_r_datafusion"] - paired["rct_aipw"]
    )
    paired_records: list[dict[str, object]] = []
    for key, group in paired.groupby(["scm_id", "regime", "n"], sort=True):
        scm_id, regime, n_rct = key
        difference = group["paired_difference"].to_numpy()
        paired_records.append(
            {
                "scm_id": int(scm_id),
                "regime": regime,
                "n": int(n_rct),
                "paired_mean_fusion_minus_aipw": float(np.mean(difference)),
                "paired_difference_mcse": float(
                    np.std(difference, ddof=1) / math.sqrt(len(difference))
                ),
            }
        )
    summary = summary.merge(
        pd.DataFrame.from_records(paired_records),
        on=["scm_id", "regime", "n"],
        how="left",
    )
    return summary.sort_values(keys).reset_index(drop=True)


def plot_summary(summary: pd.DataFrame) -> None:
    method_order = ["rct_aipw", "estimated_r_datafusion"]
    method_labels = {"rct_aipw": "RCT AIPW", "estimated_r_datafusion": "Fusion"}
    colors = {"rct_aipw": plt.cm.tab10.colors[0], "estimated_r_datafusion": plt.cm.tab10.colors[1]}
    figure, axes = plt.subplots(2, 4, figsize=(16, 8), sharex=True)
    column_spec = [
        (20, "bias", "bias_mcse", "Bias"),
        (20, "variance_ratio_vs_aipw", None, "Variance ratio"),
        (30, "bias", "bias_mcse", "Bias"),
        (30, "variance_ratio_vs_aipw", None, "Variance ratio"),
    ]
    for row_index, regime in enumerate(REGIMES):
        for column_index, (n_rct, metric, error_metric, label) in enumerate(column_spec):
            axis = axes[row_index, column_index]
            cell = summary[(summary["regime"] == regime) & (summary["n"] == n_rct)]
            for method in method_order:
                method_data = cell[cell["estimator"] == method].sort_values("scm_id")
                error = method_data[error_metric] if error_metric is not None else None
                axis.errorbar(
                    method_data["scm_id"],
                    method_data[metric],
                    yerr=error,
                    marker="o",
                    capsize=3,
                    linewidth=1.5,
                    color=colors[method],
                    label=method_labels[method],
                )
            reference = 0.0 if metric == "bias" else 1.0
            axis.axhline(reference, color="black", linewidth=0.8, linestyle="--")
            axis.set_title(f"{regime}, per-stage n={n_rct}: {label}")
            axis.set_xticks(range(SCM_COUNT), [f"SCM {k}" for k in range(SCM_COUNT)])
            axis.grid(alpha=0.25)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
    figure.suptitle(
        "Illustrative random-SCM debugging pilot (20 outer refits per cell)\n"
        "per stage: RCT n / OBS 5000; total generated per replicate: 3n / 15000",
        fontsize=14,
    )
    figure.tight_layout(rect=(0, 0.06, 1, 0.92))
    figure.savefig(FIGURE_PATH, dpi=180, bbox_inches="tight")
    plt.close(figure)


def scm_payload() -> list[dict[str, object]]:
    payload: list[dict[str, object]] = []
    for scm in SCMS:
        payload.append(
            {
                "scm_id": scm.scm_id,
                "covariate_mean": scm.covariate_mean.tolist(),
                "latent_loading": scm.latent_loading.tolist(),
                "covariate_noise_sd": scm.covariate_noise_sd.tolist(),
                "rct_mean_shift": scm.rct_mean_shift.tolist(),
                "baseline_intercept": scm.baseline_intercept,
                "baseline_coefficients": scm.baseline_coefficients.tolist(),
                "outcome_latent_coefficient": scm.outcome_latent_coefficient,
                "treatment_intercept": scm.treatment_intercept,
                "treatment_coefficients": scm.treatment_coefficients.tolist(),
                "outcome_noise_sd": scm.outcome_noise_sd,
                "rct_propensity_intercept": scm.rct_propensity_intercept,
                "rct_propensity_coefficients": scm.rct_propensity_coefficients.tolist(),
                "obs_propensity_intercept": scm.obs_propensity_intercept,
                "obs_propensity_coefficients": scm.obs_propensity_coefficients.tolist(),
                "obs_latent_coefficient": scm.obs_latent_coefficient,
                "theta_shared": true_theta(scm, "shared"),
                "theta_shifted": true_theta(scm, "shifted"),
            }
        )
    return payload


def target_smoke() -> dict[str, object]:
    results: dict[str, object] = {}
    sample_size = 250_000
    for scm in SCMS:
        for regime in REGIMES:
            rng = seed_rng(9900, scm.scm_id, regime_index(regime))
            x = rng.multivariate_normal(
                scm.rct_mean(regime), scm.covariance, size=sample_size
            )
            tau = tau_function(x, scm)
            estimate = float(tau.mean())
            mcse = float(tau.std(ddof=1) / math.sqrt(sample_size))
            analytic = true_theta(scm, regime)
            difference = abs(estimate - analytic)
            tolerance = max(0.01, 5.0 * mcse)
            if difference > tolerance:
                raise AssertionError(
                    f"Analytic target smoke failed for SCM {scm.scm_id}, {regime}: "
                    f"difference={difference}, tolerance={tolerance}"
                )
            results[f"scm{scm.scm_id}_{regime}"] = {
                "analytic": analytic,
                "monte_carlo": estimate,
                "difference": difference,
                "tolerance": tolerance,
            }
    return results


def oracle_ratio_smoke() -> dict[str, object]:
    results: dict[str, object] = {}
    sample_size = 150_000
    for scm in SCMS:
        eigenvalues = np.linalg.eigvalsh(scm.covariance)
        if eigenvalues.min() <= 0.0:
            raise AssertionError(f"SCM {scm.scm_id} covariance is not positive definite")
        if np.linalg.norm(scm.rct_mean_shift) <= 1e-8:
            raise AssertionError(f"SCM {scm.scm_id} shifted regime is not shifted")
        shared_x = seed_rng(9800, scm.scm_id, 0).multivariate_normal(
            scm.covariate_mean, scm.covariance, size=1000
        )
        if not np.array_equal(
            oracle_ratio(shared_x, "shared", scm), np.ones(1000)
        ):
            raise AssertionError("Shared oracle ratio is not exactly one")
        rng = seed_rng(9800, scm.scm_id, 1)
        x = rng.multivariate_normal(
            scm.covariate_mean, scm.covariance, size=sample_size
        )
        weight = oracle_ratio(x, "shifted", scm)
        weighted_mean = np.mean(weight[:, None] * x, axis=0)
        weighted_second = np.mean(weight[:, None] * x**2, axis=0)
        target_mean = scm.rct_mean("shifted")
        target_second = np.diag(scm.covariance) + target_mean**2
        errors = np.concatenate(
            [weighted_mean - target_mean, weighted_second - target_second]
        )
        maximum_error = float(np.max(np.abs(errors)))
        if abs(float(weight.mean()) - 1.0) > 0.03 or maximum_error > 0.08:
            raise AssertionError(
                f"Oracle ratio smoke failed for SCM {scm.scm_id}: "
                f"mean={weight.mean()}, max_error={maximum_error}"
            )
        results[f"scm{scm.scm_id}"] = {
            "minimum_covariance_eigenvalue": float(eigenvalues.min()),
            "shift_norm": float(np.linalg.norm(scm.rct_mean_shift)),
            "oracle_ratio_mean": float(weight.mean()),
            "weighted_moment_max_error": maximum_error,
            "oracle_ratio_second_moment": float(np.mean(weight**2)),
        }
    return results


def arm_gate_smoke() -> dict[str, int]:
    minimum = 10**9
    failures = 0
    for scm in SCMS:
        for regime in REGIMES:
            for n_rct in N_RCT_VALUES:
                for replication in range(N_OUTER_REPLICATIONS):
                    data = sample_dataset(
                        sample_rng(
                            scm.scm_id, regime, n_rct, replication, 10, 1
                        ),
                        n_rct,
                        "RCT",
                        regime,
                        scm,
                    )
                    arm1 = int(data["a"].sum())
                    arm0 = n_rct - arm1
                    minimum = min(minimum, arm0, arm1)
                    failures += int(min(arm0, arm1) < 3)
    if failures:
        raise AssertionError(
            f"RCT training arm preflight failed in {failures} cells; no redraw allowed"
        )
    return {"minimum_rct_training_arm_count": minimum, "arm_gate_failures": failures}


def variance_identity_smoke() -> float:
    rng = seed_rng(9700)
    n_rct = 300
    n_obs = 3000
    z0 = rng.normal(size=n_rct)
    delta = 0.4 * z0 + rng.normal(scale=0.8, size=n_rct)
    g_rct = rng.normal(size=n_rct)
    weighted_g_obs = rng.normal(size=n_obs)
    lambda_value = 0.37
    omega_value = -0.28
    compact = (
        sample_variance(z0 + lambda_value * delta - omega_value * g_rct) / n_rct
        + omega_value**2 * sample_variance(weighted_g_obs) / n_obs
    )
    expanded = (
        sample_variance(z0) / n_rct
        + 2.0 * lambda_value * sample_covariance(z0, delta) / n_rct
        + lambda_value**2 * sample_variance(delta) / n_rct
        - 2.0 * omega_value * sample_covariance(z0, g_rct) / n_rct
        - 2.0
        * lambda_value
        * omega_value
        * sample_covariance(delta, g_rct)
        / n_rct
        + omega_value**2 * sample_variance(g_rct) / n_rct
        + omega_value**2 * sample_variance(weighted_g_obs) / n_obs
    )
    difference = abs(compact - expanded)
    if difference > 1e-12:
        raise AssertionError(f"Variance identity mismatch: {difference}")
    return difference


def deterministic_smoke() -> None:
    first = pd.DataFrame(evaluate_cell(SCMS[0], "shared", 20, 1))
    second = pd.DataFrame(evaluate_cell(SCMS[0], "shared", 20, 1))
    pd.testing.assert_frame_equal(first, second, check_exact=True)


def run_smoke_tests() -> dict[str, object]:
    if generic_outcome_features(np.zeros((2, 4))).shape[1] != 22:
        raise AssertionError("Generic outcome feature dimension is not 22")
    target_results = target_smoke()
    ratio_results = oracle_ratio_smoke()
    arm_results = arm_gate_smoke()
    variance_difference = variance_identity_smoke()
    deterministic_smoke()
    return {
        "target_checks": target_results,
        "oracle_ratio_checks": ratio_results,
        **arm_results,
        "generic_outcome_feature_count": 22,
        "joint_outcome_parameter_count": 46,
        "ratio_parameter_count": 5,
        "compact_expanded_variance_difference": variance_difference,
        "deterministic_miniature_run": "pass",
    }


def validate_outputs(replications: pd.DataFrame, summary: pd.DataFrame) -> None:
    if len(replications) != 480:
        raise AssertionError(f"Expected 480 replication rows, got {len(replications)}")
    if len(summary) != 24:
        raise AssertionError(f"Expected 24 summary rows, got {len(summary)}")
    expected_estimators = {"rct_aipw", "estimated_r_datafusion"}
    if set(replications["estimator"]) != expected_estimators:
        raise AssertionError("Primary output does not contain exactly two estimators")
    group_sizes = replications.groupby(
        ["scm_id", "regime", "n", "estimator"]
    ).size()
    if not (group_sizes == N_OUTER_REPLICATIONS).all():
        raise AssertionError("A primary group does not have 20 replications")
    if min(
        replications["train_rct_arm0"].min(),
        replications["train_rct_arm1"].min(),
    ) < 3:
        raise AssertionError("An RCT training arm failed the count gate")
    numeric_replications = replications.select_dtypes(include=[np.number])
    numeric_summary = summary.select_dtypes(include=[np.number])
    if not np.isfinite(numeric_replications.to_numpy()).all():
        raise AssertionError("Replication output contains NaN or infinity")
    if not np.isfinite(numeric_summary.to_numpy()).all():
        raise AssertionError("Summary output contains NaN or infinity")
    if replications.isna().any().any() or summary.isna().any().any():
        raise AssertionError("Output contains missing values")
    if not (
        replications["fitted_ratio_clip_fraction"].between(0.0, 1.0).all()
    ):
        raise AssertionError("Ratio clipping fraction is outside [0,1]")
    if not FIGURE_PATH.exists() or FIGURE_PATH.stat().st_size <= 0:
        raise AssertionError("Figure output is missing or empty")


def run_full_simulation() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    for scm in SCMS:
        for regime in REGIMES:
            for n_rct in N_RCT_VALUES:
                cell_rows = evaluate_cell(
                    scm, regime, n_rct, N_OUTER_REPLICATIONS
                )
                rows.extend(cell_rows)
                print(
                    f"completed scm={scm.scm_id} regime={regime} "
                    f"per_stage_n={n_rct}: {len(cell_rows)} rows",
                    flush=True,
                )
    replications = pd.DataFrame.from_records(rows).sort_values(
        ["scm_id", "regime", "n", "replication", "estimator"]
    )
    replications = replications.reset_index(drop=True)
    summary = build_summary(replications)
    replications.to_csv(REPLICATION_PATH, index=False, float_format="%.10g")
    summary.to_csv(SUMMARY_PATH, index=False, float_format="%.10g")
    plot_summary(summary)
    validate_outputs(replications, summary)
    return replications, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--smoke-only",
        action="store_true",
        help="Run mathematical, design, and reproducibility checks only.",
    )
    arguments = parser.parse_args()
    start = time.perf_counter()
    print("SCM_PARAMETERS=" + json.dumps(scm_payload(), sort_keys=True), flush=True)
    smoke_results = run_smoke_tests()
    print("SMOKE_RESULTS=" + json.dumps(smoke_results, sort_keys=True), flush=True)
    if arguments.smoke_only:
        print(f"SMOKE_RUNTIME_SECONDS={time.perf_counter() - start:.3f}")
        return
    replications, summary = run_full_simulation()
    payload = {
        "replication_rows": len(replications),
        "summary_rows": len(summary),
        "replication_path": str(REPLICATION_PATH),
        "summary_path": str(SUMMARY_PATH),
        "figure_path": str(FIGURE_PATH),
        "runtime_seconds": time.perf_counter() - start,
    }
    print("RUN_RESULTS=" + json.dumps(payload, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
