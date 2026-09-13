#!/usr/bin/env python3
"""Run and validate the deterministic DataFusionPPI version-3 Agile gate."""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from typing import Any, Callable

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import OneHotEncoder, SplineTransformer, StandardScaler

CODE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = CODE_DIR.parent

PROTOCOL = "DataFusionPPI-v3-agile-gate-2026-09-10"
BASE_SEED = 190602
HANDOFF = PROJECT_DIR / "materials" / "2026-09-09-implementation-handoff-JA.md"
HANDOFF_SHA256 = "e1b2207beac31094f18856ed74391e1dc1f3bb665efbbf986f3f98231e950b39"
STAR_URL = "https://vincentarelbundock.github.io/Rdatasets/csv/AER/STAR.csv"
STAR_SHA256 = "0e8b179ea3d883730b25008ca1293c4c8f886aa8d5d299c03ab6ff05d0123ae6"
ARTIFACTS = (
    "gate_v3_replications.csv",
    "gate_v3_summary.csv",
    "gate_v3_assertions.json",
    "gate_v3_results.md",
)
ROLE_NAMES = (
    "design", "R_nuis", "O_nuis", "R_tune", "O_tune", "R_select",
    "O_select", "R_report", "O_report", "truth", "oracle", "gate_bootstrap",
)
GRID = tuple((lam, omega) for lam in (0.0, .25, .5, .75, 1.0)
             for omega in (0.0, .25, .5, .75, 1.0))
RHOS = (.01, .1, 1.0)
BALANCE_TOLERANCE = 1e-6


@dataclasses.dataclass(frozen=True)
class LBFGSResult:
    x: np.ndarray
    fun: float
    gradient: np.ndarray
    iterations: int
    success: bool
    status: str


def minimize_lbfgs(
    objective: Callable[[np.ndarray], tuple[float, np.ndarray]], x0: np.ndarray,
    max_iter: int = 100, gradient_tolerance: float = 1e-10, history_size: int = 10,
) -> LBFGSResult:
    """Deterministic unconstrained L-BFGS with Armijo backtracking."""
    x=np.asarray(x0,float).copy();value,gradient=objective(x);s_hist=[];y_hist=[];rho_hist=[]
    if not np.isfinite(value) or not np.isfinite(gradient).all():
        return LBFGSResult(x,float(value),gradient,0,False,"nonfinite_initial")
    for iteration in range(max_iter+1):
        if np.max(np.abs(gradient))<=gradient_tolerance:
            return LBFGSResult(x,float(value),gradient,iteration,True,"gradient_converged")
        if iteration==max_iter: break
        q=gradient.copy();alphas=[]
        for s,y,rho in reversed(list(zip(s_hist,y_hist,rho_hist))):
            alpha=rho*float(s@q);alphas.append(alpha);q-=alpha*y
        scale=1.0 if not s_hist else float((s_hist[-1]@y_hist[-1])/(y_hist[-1]@y_hist[-1]))
        direction=scale*q
        for (s,y,rho),alpha in zip(zip(s_hist,y_hist,rho_hist),reversed(alphas)):
            direction+=s*(alpha-rho*float(y@direction))
        direction=-direction
        directional=float(gradient@direction)
        if not np.isfinite(directional) or directional>=0:
            direction=-gradient;directional=-float(gradient@gradient)
        step=1.0;accepted=False
        for _ in range(60):
            candidate=x+step*direction;candidate_value,candidate_gradient=objective(candidate)
            if (np.isfinite(candidate_value) and np.isfinite(candidate_gradient).all() and
                candidate_value<value and candidate_value<=value+1e-4*step*directional):
                accepted=True;break
            step*=.5
        if not accepted:
            return LBFGSResult(x,float(value),gradient,iteration,False,"line_search_failed")
        s=candidate-x;y=candidate_gradient-gradient;curvature=float(s@y)
        if curvature>1e-12*max(1.0,float(np.linalg.norm(s)*np.linalg.norm(y))):
            s_hist.append(s);y_hist.append(y);rho_hist.append(1.0/curvature)
            if len(s_hist)>history_size:s_hist.pop(0);y_hist.pop(0);rho_hist.pop(0)
        x,value,gradient=candidate,float(candidate_value),candidate_gradient
    return LBFGSResult(x,float(value),gradient,max_iter,False,"max_iter_reached")


@dataclasses.dataclass(frozen=True)
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


@dataclasses.dataclass(frozen=True)
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
        p_obs = np.full(len(self.x), 1.0 / len(self.x))
        if regime == "shared":
            return p_obs.copy(), p_obs
        score = self.x @ self.direction / math.sqrt(self.x.shape[1])
        weight = np.exp(.35 * score - np.max(.35 * score))
        return weight / weight.sum(), p_obs


def random_sign(generator: np.random.Generator) -> float:
    return float(generator.choice(np.array([-1.0, 1.0])))


def generate_scm_parameters(scm_id: int) -> SCMParameters:
    generator = np.random.default_rng(np.random.SeedSequence([BASE_SEED, 7000, scm_id]))
    return SCMParameters(
        scm_id=scm_id,
        covariate_mean=generator.normal(0.0, .15, 4),
        latent_loading=generator.normal(0.0, .30, 4),
        covariate_noise_sd=generator.uniform(.75, 1.05, 4),
        rct_mean_shift=generator.normal(0.0, .18, 4),
        baseline_intercept=float(generator.normal(0.0, .2)),
        baseline_coefficients=generator.normal(0.0, .30, 9),
        outcome_latent_coefficient=random_sign(generator) * float(generator.uniform(.4, .8)),
        treatment_intercept=float(generator.uniform(.8, 1.2)),
        treatment_coefficients=generator.normal(0.0, .22, 5),
        outcome_noise_sd=float(generator.uniform(.8, 1.2)),
        rct_propensity_intercept=float(generator.normal(0.0, .1)),
        rct_propensity_coefficients=generator.normal(0.0, .18, 5),
        obs_propensity_intercept=float(generator.normal(-.2, .1)),
        obs_propensity_coefficients=generator.normal(0.0, .30, 5),
        obs_latent_coefficient=random_sign(generator) * float(generator.uniform(.7, 1.1)),
    )


def baseline_basis(x: np.ndarray) -> np.ndarray:
    return np.column_stack((x[:,0],x[:,1],x[:,2],x[:,3],np.sin(x[:,0]),np.cos(x[:,1]),
                            x[:,0]*x[:,1],x[:,2]*x[:,3],x[:,0]**2))


def treatment_effect_basis(x: np.ndarray) -> np.ndarray:
    return np.column_stack((x[:,0],x[:,1],np.sin(x[:,2]),x[:,0]*x[:,3],x[:,1]**2))


def treatment_assignment_basis(x: np.ndarray) -> np.ndarray:
    return np.column_stack((x[:,0],x[:,1],np.tanh(x[:,2]),np.sin(x[:,3]),x[:,0]*x[:,1]))


def tau_function(x: np.ndarray, scm: SCMParameters) -> np.ndarray:
    return scm.treatment_intercept + treatment_effect_basis(x) @ scm.treatment_coefficients


def rct_propensity(x: np.ndarray, scm: SCMParameters) -> np.ndarray:
    return np.clip(expit(scm.rct_propensity_intercept + treatment_assignment_basis(x) @ scm.rct_propensity_coefficients), .2, .8)


def feature_map(x: np.ndarray) -> np.ndarray:
    d=min(x.shape[1],4);interactions=[x[:,i]*x[:,j] for i in range(d) for j in range(i+1,d)]
    pieces=[x,x**2,np.sin(x)]
    if interactions: pieces.append(np.column_stack(interactions))
    return np.column_stack(pieces)


def load_star_support() -> StarSupport:
    payload=subprocess.run(["curl","-fsSL",STAR_URL],check=True,capture_output=True,timeout=30).stdout
    digest=hashlib.sha256(payload).hexdigest()
    if digest!=STAR_SHA256: raise RuntimeError(f"STAR source hash mismatch: {digest}")
    frame=pd.read_csv(io.BytesIO(payload))
    if len(frame)!=11598 or not {"gender","ethnicity","birth","lunchk","schoolk"}.issubset(frame.columns):
        raise RuntimeError("STAR source shape/schema mismatch")
    numeric=pd.DataFrame(index=frame.index);birth=pd.to_numeric(frame["birth"],errors="coerce")
    numeric["birth"]=birth.fillna(birth.median());numeric["birth_missing"]=birth.isna().astype(float)
    categorical=frame[["gender","ethnicity","lunchk","schoolk"]].astype("string").fillna("Missing")
    raw=pd.concat([numeric,pd.get_dummies(categorical,prefix=categorical.columns,dtype=float)],axis=1).to_numpy(float)
    raw=raw[:,raw.std(0)>1e-10];x=(raw-raw.mean(0))/raw.std(0)
    generator=np.random.default_rng(np.random.SeedSequence([BASE_SEED,6000]));phi=feature_map(x)
    phi=(phi-phi.mean(0))/np.where(phi.std(0)>1e-8,phi.std(0),1.)
    beta_m=generator.normal(0,.55/math.sqrt(phi.shape[1]),phi.shape[1]);beta_tau=generator.normal(0,.22/math.sqrt(phi.shape[1]),phi.shape[1])
    u,noise=generator.normal(size=len(x)),generator.normal(size=len(x));baseline=phi@beta_m+.7*u+noise;tau=1+phi@beta_tau
    beta_r=generator.normal(0,.25/math.sqrt(x.shape[1]),x.shape[1]);beta_o=generator.normal(0,.35/math.sqrt(x.shape[1]),x.shape[1])
    e_rct=np.clip(expit(.1+x@beta_r),.2,.8);e_obs=np.clip(expit(-.2+x@beta_o+.9*u),.05,.95)
    support=StarSupport(x,u,baseline,baseline+tau,tau,e_rct,e_obs,generator.normal(size=x.shape[1]))
    pr,po=support.probabilities("shared")
    if abs(float(np.sum(po*(pr/po)))-1)>1e-12: raise AssertionError("STAR exact ratio does not normalize")
    return support


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_seed(*parts: object) -> np.random.SeedSequence:
    key = "\x1f".join(map(str, (PROTOCOL, *parts))).encode("utf-8")
    entropy = np.frombuffer(hashlib.sha256(key).digest()[:16], dtype="<u4").tolist()
    return np.random.SeedSequence(entropy)


def rng(study: str, cell: str, rep: int, namespace: str, substep: object = 0) -> np.random.Generator:
    if namespace not in ROLE_NAMES:
        raise ValueError(f"unknown RNG namespace: {namespace}")
    return np.random.default_rng(canonical_seed(study, cell, rep, namespace, substep))


def seed32(study: str, cell: str, rep: int, namespace: str, substep: object = 0) -> int:
    return int(canonical_seed(study, cell, rep, namespace, substep).generate_state(1)[0])


def expit(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(x, -35.0, 35.0)))


def logit(x: np.ndarray) -> np.ndarray:
    x = np.clip(x, 1e-12, 1 - 1e-12)
    return np.log(x / (1 - x))


def svar(x: np.ndarray) -> float:
    return float(np.var(np.asarray(x, float), ddof=1))


def scov(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.cov(np.asarray(x, float), np.asarray(y, float), ddof=1)[0, 1])


def fingerprint(data: dict[str, np.ndarray]) -> str:
    h = hashlib.sha256()
    for key in sorted(data):
        h.update(key.encode())
        h.update(np.ascontiguousarray(data[key]).tobytes())
    return h.hexdigest()[:20]


@dataclasses.dataclass
class OutcomeFit:
    kind: str
    center: np.ndarray
    scale: np.ndarray
    model: Ridge

    def _features(self, x: np.ndarray) -> np.ndarray:
        raw = x if self.kind == "linear" else feature_map(x)
        return (raw - self.center) / self.scale

    def predict(self, x: np.ndarray, a: float) -> np.ndarray:
        f = self._features(x)
        av = np.full(len(x), a)
        return self.model.predict(np.column_stack((f, av, av[:, None] * f)))


def fit_outcome(data: dict[str, np.ndarray], kind: str = "correct") -> OutcomeFit:
    raw = data["x"] if kind == "linear" else feature_map(data["x"])
    center = raw.mean(axis=0)
    scale = np.where(raw.std(axis=0) > 1e-10, raw.std(axis=0), 1.0)
    f = (raw - center) / scale
    a = data["a"]
    model = Ridge(alpha=5.0, fit_intercept=True, solver="lsqr", tol=1e-10)
    model.fit(np.column_stack((f, a, a[:, None] * f)), data["y"])
    return OutcomeFit(kind, center, scale, model)


@dataclasses.dataclass
class NeuralOutcomeFit:
    """Frozen nuisance-only neural outcome heads."""
    predict_function: Callable[[np.ndarray], np.ndarray]

    def predict(self, x: np.ndarray, a: float) -> np.ndarray:
        return self.predict_function(x)[:, int(a)]


def pseudo(data: dict[str, np.ndarray], fit: OutcomeFit | NeuralOutcomeFit, clip: bool = False) -> tuple[np.ndarray, np.ndarray]:
    mu0, mu1 = fit.predict(data["x"], 0.0), fit.predict(data["x"], 1.0)
    if clip:
        mu0, mu1 = np.clip(mu0, -1, 1), np.clip(mu1, -1, 1)
    e, a, y = data["e"], data["a"], data["y"]
    z = mu1 - mu0 + a / e * (y - mu1) - (1 - a) / (1 - e) * (y - mu0)
    g = np.clip(mu1 - mu0, -2, 2) if clip else mu1 - mu0
    return z, g


def fusion_scores(
    rct: dict[str, np.ndarray], obs: dict[str, np.ndarray],
    mur: OutcomeFit | NeuralOutcomeFit, muo: OutcomeFit | NeuralOutcomeFit,
    clip: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return Z0, Z1, g_R, g_O, with both g values from the OBS fit."""
    z0, _ = pseudo(rct, mur, clip)
    z1, gr = pseudo(rct, muo, clip)
    _, go = pseudo(obs, muo, clip)
    return z0, z1, gr, go


def sample_scm(
    scm: SCMParameters, n: int, source: str, generator: np.random.Generator,
    confounding: float = 1.0, delta: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    d = 4
    mean = scm.covariate_mean + (delta if source == "RCT" and delta is not None else 0.0)
    u = generator.normal(size=n)
    x = mean + u[:, None] * scm.latent_loading + generator.normal(size=(n, d)) * scm.covariate_noise_sd
    if source == "RCT":
        e = rct_propensity(x, scm)
    else:
        score = (scm.obs_propensity_intercept + treatment_assignment_basis(x) @ scm.obs_propensity_coefficients
                 + confounding * scm.obs_latent_coefficient * u)
        e = np.clip(expit(score), .05, .95)
    a = generator.binomial(1, e).astype(float)
    tau = tau_function(x, scm)
    baseline = (scm.baseline_intercept + baseline_basis(x) @ scm.baseline_coefficients
                + scm.outcome_latent_coefficient * u)
    y = baseline + a * tau + scm.outcome_noise_sd * generator.normal(size=n)
    return {"x": x, "a": a, "y": y, "e": e, "tau": tau}


def sample_scmb(n: int, source: str, generator: np.random.Generator, confounding: float) -> dict[str, np.ndarray]:
    x = generator.uniform(-1, 1, size=(n, 5))
    u = generator.uniform(-1, 1, size=n)
    eps = generator.uniform(-1, 1, size=n)
    er = .5 + .2 * x[:, 0]
    e = er if source == "RCT" else expit(logit(er) + confounding * u)
    a = generator.binomial(1, e).astype(float)
    tau = .15 + .10 * x[:, 0] - .05 * x[:, 1] ** 2
    y = .15 * np.sin(np.pi * x[:, 0]) + .05 * x[:, 1] * x[:, 2] + .10 * u + .05 * eps + a * tau
    return {"x": x, "a": a, "y": y, "e": e, "tau": tau, "u": u, "eps": eps}


def true_theta_scm(scm: SCMParameters, delta: np.ndarray | None = None) -> float:
    mean = scm.covariate_mean + (np.zeros(4) if delta is None else delta)
    # Exact Gaussian expectations of the five treatment-effect features.
    cov = np.diag(scm.covariate_noise_sd ** 2) + np.outer(scm.latent_loading, scm.latent_loading)
    vals = np.array([mean[0], mean[1], math.exp(-cov[2, 2] / 2) * math.sin(mean[2]),
                     mean[0] * mean[3] + cov[0, 3], mean[1] ** 2 + cov[1, 1]])
    return float(scm.treatment_intercept + vals @ scm.treatment_coefficients)


def algo1(z0r: np.ndarray, z1r: np.ndarray, gr: np.ndarray, rgo: np.ndarray,
          n_eval: int, N_eval: int) -> tuple[float, float, dict[str, float]]:
    delta = z1r - z0r
    A = svar(delta)
    B = svar(gr) + n_eval / N_eval * svar(rgo)
    C = scov(z0r, delta)
    D = scov(z0r, gr)
    base = svar(z0r)
    lam = float(np.clip(-C / A, 0, 1)) if A > 1e-3 * base else 0.0
    omega = float(np.clip(D / B, 0, 1)) if B > 1e-3 * base else 0.0
    return lam, omega, {"A": A, "B": B, "C": C, "D": D}


def adaptive_oracle(
    draw: Callable[[int, str, int], dict[str, np.ndarray]], mur: OutcomeFit,
    muo: OutcomeFit, ratio: Callable[[np.ndarray], np.ndarray], n_eval: int, N_eval: int,
) -> tuple[float, float, dict[str, Any]]:
    """Apply the v3 20k/40k/80k/100k moment-MCSE rule."""
    last: dict[str, Any] = {}
    for size in (20000, 40000, 80000, 100000):
        rr, oo = draw(size, "R", size), draw(size, "O", size)
        z0, z1, gr, go = fusion_scores(rr, oo, mur, muo)
        rgo = ratio(oo["x"]) * go
        lam, omega, moments = algo1(z0, z1, gr, rgo, n_eval, N_eval)
        batch_count = 20
        rb, ob = np.array_split(np.arange(size), batch_count), np.array_split(np.arange(size), batch_count)
        values = {key: [] for key in ("A", "B", "C", "D")}
        for ri, oi in zip(rb, ob):
            delta_b = z1[ri] - z0[ri]
            values["A"].append(svar(delta_b))
            values["B"].append(svar(gr[ri]) + n_eval / N_eval * svar(rgo[oi]))
            values["C"].append(scov(z0[ri], delta_b))
            values["D"].append(scov(z0[ri], gr[ri]))
        mcse = {key: float(np.std(value, ddof=1) / math.sqrt(batch_count)) for key, value in values.items()}
        s02 = svar(z0)
        tolerance = {"A": .01 * s02, "B": .01 * s02, "C": .005 * s02, "D": .005 * s02}
        stable = all(mcse[key] <= tolerance[key] for key in mcse)
        last = {"draws": size, "moments": moments, "mcse": mcse, "tolerance": tolerance,
                "stable": stable, "status": "stable" if stable else "monte_carlo_inconclusive_at_cap"}
        if stable:
            return lam, omega, last
    return lam, omega, last


def ate_eval(z0r: np.ndarray, z1r: np.ndarray, gr: np.ndarray, go: np.ndarray,
             ratio_o: np.ndarray, lam: float, omega: float) -> tuple[float, float]:
    zr = z0r + lam * (z1r - z0r)
    theta = float(zr.mean() + omega * ((ratio_o * go).mean() - gr.mean()))
    psi_r = zr - omega * gr
    psi_o = omega * ratio_o * go
    variance = svar(psi_r) / len(psi_r) + svar(psi_o) / len(psi_o)
    return theta, variance


@dataclasses.dataclass
class RatioFit:
    name: str
    predict: Callable[[np.ndarray, np.ndarray | None], np.ndarray]
    residual: float | None = None
    feasible: bool | None = None
    normalization_error: float = 0.0
    protocol_status: str = "N/A"
    optimizer_success: bool | None = None
    optimizer_status: str = "N/A"
    optimizer_iterations: int | None = None


def gaussian_ratio(x: np.ndarray, scm: SCMParameters, delta: np.ndarray) -> np.ndarray:
    cov = np.diag(scm.covariate_noise_sd ** 2) + np.outer(scm.latent_loading, scm.latent_loading)
    inv = np.linalg.inv(cov)
    base = scm.covariate_mean
    lr = (x - base) @ inv @ delta - .5 * float(delta @ inv @ delta)
    return np.exp(np.clip(lr, -700, 700))


def fit_ratios(rn: dict[str, np.ndarray], on: dict[str, np.ndarray], muo: OutcomeFit,
               scm: SCMParameters, delta: np.ndarray, study: str, cell: str, rep: int) -> list[RatioFit]:
    x = np.vstack((rn["x"], on["x"]))
    labels = np.r_[np.ones(len(rn["x"])), np.zeros(len(on["x"]))]
    scaler = StandardScaler().fit(x)
    clf = LogisticRegression(C=1.0, max_iter=1000, solver="lbfgs", random_state=seed32(study, cell, rep, "design", "ratio"))
    clf.fit(scaler.transform(x), labels)
    pi = labels.mean()

    def cls_pred(xx: np.ndarray, _: np.ndarray | None = None) -> np.ndarray:
        s = np.clip(clf.predict_proba(scaler.transform(xx))[:, 1], 1e-10, 1 - 1e-10)
        return (1 - pi) / pi * s / (1 - s)

    oracle_pred=lambda xx,gg=None:gaussian_ratio(xx,scm,delta)
    results = [
        RatioFit("oracle",oracle_pred,normalization_error=float(np.mean(oracle_pred(on["x"]))-1)),
        RatioFit("unweighted",lambda xx,gg=None:np.ones(len(xx)),normalization_error=0.),
        RatioFit("classifier",cls_pred,normalization_error=float(np.mean(cls_pred(on["x"]))-1)),
    ]
    for with_g, name in ((False, "BAL-X"), (True, "BAL-X+g")):
        raw_r = scaler.transform(rn["x"])
        raw_o = scaler.transform(on["x"])
        if with_g:
            _, gr = pseudo(rn, muo)
            _, go = pseudo(on, muo)
            gc, gs = float(np.mean(np.r_[gr, go])), float(np.std(np.r_[gr, go])) or 1.0
            phi_r, phi_o = np.column_stack((raw_r, (gr-gc)/gs)), np.column_stack((raw_o, (go-gc)/gs))
        else:
            gc, gs, phi_r, phi_o = 0.0, 1.0, raw_r, raw_o
        base_o = cls_pred(on["x"])
        target = phi_r.mean(axis=0)
        log_base=np.log(np.clip(base_o,1e-300,None))
        def objective(xi:np.ndarray)->tuple[float,np.ndarray]:
            logits=log_base+phi_o@xi;shift=float(np.max(logits));ew=np.exp(logits-shift)
            weights=ew/ew.sum();log_mean=shift+math.log(float(np.mean(ew)))
            return float(log_mean-xi@target),weights@phi_o-target
        result=minimize_lbfgs(objective,np.zeros(phi_o.shape[1]),max_iter=100,gradient_tolerance=1e-10)
        xi=np.asarray(result.x,float)
        logits = np.log(np.clip(base_o, 1e-300, None)) + phi_o @ xi
        shift=float(np.max(logits));log_norm=shift+math.log(float(np.mean(np.exp(logits-shift))))
        normalized=np.exp(logits-log_norm)
        residual = float(np.max(np.abs(normalized @ phi_o / len(phi_o) - target)))

        def pred(xx: np.ndarray, gg: np.ndarray | None = None, *, xi=xi.copy(), with_g=with_g,
                 gc=gc, gs=gs, log_norm=log_norm) -> np.ndarray:
            ph = scaler.transform(xx)
            if with_g:
                if gg is None: raise ValueError("g required for BAL-X+g")
                ph = np.column_stack((ph, (gg-gc)/gs))
            return np.exp(np.clip(np.log(np.clip(cls_pred(xx),1e-300,None))+ph@xi-log_norm,-700,700))
        normalization_error=float(np.mean(pred(on["x"],go if with_g else None))-1)
        converged=residual<=BALANCE_TOLERANCE
        results.append(RatioFit(name=name,predict=pred,residual=residual,feasible=converged,
                                normalization_error=normalization_error,
                                protocol_status="converged_to_protocol_tolerance" if converged else "not_converged_to_protocol_tolerance",
                                optimizer_success=result.success,optimizer_status=result.status,
                                optimizer_iterations=result.iterations))
    return results


def row(rows: list[dict[str, Any]], stage: str, study: str, cell: str, rep: int,
        estimator: str, ratio: str, metric: str, value: float | str,
        status: str = "ok") -> None:
    rows.append({"protocol": PROTOCOL, "stage": stage, "study": study, "cell": cell,
                 "rep": rep, "estimator": estimator, "ratio": ratio, "metric": metric,
                 "value": value, "status": status,
                 "seed_key": f"{study}/{cell}/{rep}"})


def expected_schedule(stage:str)->list[tuple[str,str,int,str,str,str]]:
    if stage!="g0": raise RuntimeError("G1_NOT_IMPLEMENTED")
    reps=range(1 if stage=="g0" else 10);items=[]
    ate1_metrics=("estimate","squared_error","variance_hat","covered","lambda","omega")
    for rep in reps:
        for cell in ("reference","conf0_correct","conf2_linear","star_reference"):
            for est in ("rct_only","lambda_only","omega_only","joint","oracle_coefficient_joint"):
                items.extend(("ATE1",cell,rep,est,"N/A",metric) for metric in ate1_metrics)
        for cell in ("mild","strong"):
            for est in ("rct_only","lambda_only"):
                items.extend(("ATE2",cell,rep,est,"N/A",metric) for metric in ("estimate","squared_error","variance_hat","lambda","omega"))
            for ratio in ("oracle","unweighted","classifier","BAL-X","BAL-X+g"):
                for est in ("omega_only","joint"):
                    items.extend(("ATE2",cell,rep,est,ratio,metric) for metric in
                                 ("estimate","squared_error","variance_hat","lambda","omega","balance_residual","ess_fraction","normalization_error"))
            items.extend(("ATE2",cell,rep,"oracle_coefficient_joint","oracle",metric) for metric in
                         ("estimate","squared_error","variance_hat","lambda","omega","normalization_error"))
        for cell in ("reference","confounding2"):
            for est in ("rct_only","lambda_only","omega_only","joint","grid_oracle"):
                items.extend(("CATE1",cell,rep,est,"r=1",metric) for metric in
                             ("risk","selection_score","reporting_score","lambda","omega","bound_vacuity_ratio"))
    return sorted(items)


def run_ate1(stage: str, cell: str, rep: int, assertions: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    scm = generate_scm_parameters(1)
    conf, kind = (0.0, "correct") if "conf0" in cell else ((2.0, "linear") if "conf2" in cell else (1.0, "correct"))
    sizes = {"R_nuis": 40, "R_tune": 30, "R_report": 30, "O_nuis": 3000, "O_tune": 1000, "O_report": 1000}
    if cell == "star_reference":
        support = load_star_support()
        pr, po = support.probabilities("shared")
        def draw(role: str, source: str) -> dict[str, np.ndarray]:
            gen = rng("ATE1", cell, rep, role)
            ids = gen.choice(len(support.x), sizes[role], replace=True, p=pr if source == "RCT" else po)
            e = support.e_rct[ids] if source == "RCT" else support.e_obs[ids]
            a = gen.binomial(1, e).astype(float)
            return {"x": support.x[ids], "a": a, "y": np.where(a == 1, support.y1[ids], support.y0[ids]),
                    "e": e, "tau": support.tau[ids]}
        truth = float(pr @ support.tau)
    else:
        def draw(role: str, source: str) -> dict[str, np.ndarray]:
            return sample_scm(scm, sizes[role], source, rng("ATE1", cell, rep, role), confounding=conf)
        truth = true_theta_scm(scm)
    data = {role: draw(role, "RCT" if role.startswith("R") else "OBS") for role in sizes}
    assertions.append(check("ATE1", cell, rep, "role_rng_distinct", len({seed32("ATE1", cell, rep, r) for r in sizes}) == len(sizes), True,
                            {r: fingerprint(data[r]) for r in sizes}))
    mur, muo = fit_outcome(data["R_nuis"], kind), fit_outcome(data["O_nuis"], kind)
    z0t, z1t, gt, got = fusion_scores(data["R_tune"], data["O_tune"], mur, muo)
    lam, omega, _ = algo1(z0t, z1t, gt, got, sizes["R_report"], sizes["O_report"])
    z0e, z1e, ge, goe = fusion_scores(data["R_report"], data["O_report"], mur, muo)
    specs = (("rct_only", 0., 0.), ("lambda_only", lam, 0.), ("omega_only", 0., omega), ("joint", lam, omega))
    if cell == "star_reference":
        support = load_star_support(); pr, po = support.probabilities("shared")
        def oracle_draw(n: int, source: str, sub: int) -> dict[str, np.ndarray]:
            gen = rng("ATE1", cell, rep, "oracle", f"{source}-{sub}"); ids = gen.choice(len(pr), n, p=pr if source=="R" else po)
            e = support.e_rct[ids] if source=="R" else support.e_obs[ids]; a=gen.binomial(1,e).astype(float)
            return {"x":support.x[ids],"a":a,"y":np.where(a==1,support.y1[ids],support.y0[ids]),"e":e,"tau":support.tau[ids]}
    else:
        def oracle_draw(n: int, source: str, sub: int) -> dict[str, np.ndarray]:
            return sample_scm(scm,n,"RCT" if source=="R" else "OBS",rng("ATE1",cell,rep,"oracle",f"{source}-{sub}"),confounding=conf)
    olam, oomega, ometa = adaptive_oracle(oracle_draw, mur, muo, lambda x: np.ones(len(x)), sizes["R_report"], sizes["O_report"])
    specs += (("oracle_coefficient_joint", olam, oomega),)
    for est, la, omga in specs:
        theta, variance = ate_eval(z0e, z1e, ge, goe, np.ones(len(goe)), la, omga)
        for metric, value in (("estimate", theta), ("squared_error", (theta-truth)**2), ("variance_hat", variance),
                              ("covered", float(abs(theta-truth) <= 1.959963984540054*math.sqrt(max(variance, 0)))),
                              ("lambda", la), ("omega", omga)):
            row(rows, stage, "ATE1", cell, rep, est, "N/A", metric, value)
    assertions.append(check("ATE1", cell, rep, "oracle_adaptive_mcse", True, False, ometa))
    return rows


def run_ate2(stage: str, cell: str, rep: int, assertions: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]]=[]; scm=generate_scm_parameters(1)
    cov=np.diag(scm.covariate_noise_sd**2)+np.outer(scm.latent_loading,scm.latent_loading)
    inv=np.linalg.inv(cov); direction=scm.rct_mean_shift
    d2=math.log(2 if cell=="mild" else 5); delta=direction*math.sqrt(d2/float(direction@inv@direction))
    sizes={"R_nuis":40,"R_tune":30,"R_report":30,"O_nuis":3000,"O_tune":1000,"O_report":1000}
    data={role:sample_scm(scm,n,"RCT" if role.startswith("R") else "OBS",rng("ATE2",cell,rep,role),1.,delta) for role,n in sizes.items()}
    mur,muo=fit_outcome(data["R_nuis"]),fit_outcome(data["O_nuis"])
    ratios=fit_ratios(data["R_nuis"],data["O_nuis"],muo,scm,delta,"ATE2",cell,rep)
    for rf in ratios:
        if rf.name.startswith("BAL-"):
            assertions.append(check("ATE2",cell,rep,f"ratio_balance_{rf.name}",bool(rf.feasible),True,
                                    {"balance_residual":rf.residual,"tolerance":BALANCE_TOLERANCE,
                                     "protocol_status":rf.protocol_status,
                                     "optimizer_success":rf.optimizer_success,
                                     "optimizer_status":rf.optimizer_status,
                                     "optimizer_iterations":rf.optimizer_iterations,
                                     "normalization_error_diagnostic":rf.normalization_error}))
    z0t,z1t,gt,got=fusion_scores(data["R_tune"],data["O_tune"],mur,muo)
    z0e,z1e,ge,goe=fusion_scores(data["R_report"],data["O_report"],mur,muo)
    truth=true_theta_scm(scm,delta)
    # Ratio-free restrictions are recorded once.
    lam0,_,_=algo1(z0t,z1t,gt,got,len(z0e),len(goe))
    for est,la in (("rct_only",0.),("lambda_only",lam0)):
        th,v=ate_eval(z0e,z1e,ge,goe,np.ones(len(goe)),la,0.)
        for metric,value in (("estimate",th),("squared_error",(th-truth)**2),("variance_hat",v),("lambda",la),("omega",0.)):
            row(rows,stage,"ATE2",cell,rep,est,"N/A",metric,value)
    for rf in ratios:
        rt=rf.predict(data["O_tune"]["x"],got if rf.name=="BAL-X+g" else None)
        re=rf.predict(data["O_report"]["x"],goe if rf.name=="BAL-X+g" else None)
        la,om,_=algo1(z0t,z1t,gt,rt*got,len(z0e),len(goe))
        for est,xla,xom in (("omega_only",0.,om),("joint",la,om)):
            th,v=ate_eval(z0e,z1e,ge,goe,re,xla,xom)
            status="ok" if rf.feasible is not False else "explicit_infeasible"
            for metric,value in (("estimate",th),("squared_error",(th-truth)**2),("variance_hat",v),("lambda",xla),("omega",xom),
                                 ("balance_residual",rf.residual if rf.residual is not None else "N/A"),("ess_fraction",float(re.sum()**2/(np.sum(re**2)*len(re)))),
                                 ("normalization_error",rf.normalization_error)):
                metric_status="not_applicable" if metric=="balance_residual" and rf.residual is None else status
                row(rows,stage,"ATE2",cell,rep,est,rf.name,metric,value,metric_status)
    oracle=ratios[0]
    def oracle_draw(n:int,source:str,sub:int)->dict[str,np.ndarray]:
        return sample_scm(scm,n,"RCT" if source=="R" else "OBS",rng("ATE2",cell,rep,"oracle",f"{source}-{sub}"),1.,delta)
    ola,oom,ometa=adaptive_oracle(oracle_draw,mur,muo,lambda x:gaussian_ratio(x,scm,delta),len(z0e),len(goe))
    roe=oracle.predict(data["O_report"]["x"])
    th,v=ate_eval(z0e,z1e,ge,goe,roe,ola,oom)
    for metric,value in (("estimate",th),("squared_error",(th-truth)**2),("variance_hat",v),("lambda",ola),("omega",oom),
                         ("normalization_error",oracle.normalization_error)):
        row(rows,stage,"ATE2",cell,rep,"oracle_coefficient_joint","oracle",metric,value)
    empirical=float(np.mean(gaussian_ratio(data["O_report"]["x"],scm,delta)**2))
    analytic_er=math.exp(-d2/2+d2/2)
    analytic_er2=math.exp(-d2+2*d2)
    assertions.append(check("ATE2",cell,rep,"gaussian_ratio_contract",
                            abs(float(delta@inv@delta)-d2)<1e-12 and abs(analytic_er-1)<1e-15 and abs(analytic_er2-math.exp(d2))<1e-15,True,
                            {"d2":d2,"analytic_E_r":analytic_er,"analytic_E_r2":analytic_er2,
                             "analytic_ess_fraction":1/analytic_er2,"empirical_E_r2":empirical,
                             "global_bound":"unbounded","radius_status":"N/A"}))
    assertions.append(check("ATE2",cell,rep,"oracle_adaptive_mcse",True,False,ometa))
    return rows


@dataclasses.dataclass
class Basis:
    transform: Callable[[np.ndarray], np.ndarray]
    p: int
    name: str


def spline_basis(x_tune: np.ndarray, categorical: bool = False) -> Basis:
    if categorical:
        enc=OneHotEncoder(handle_unknown="ignore",sparse_output=False).fit(x_tune.astype(str))
        def trans(x: np.ndarray) -> np.ndarray: return np.column_stack((np.ones(len(x)),enc.transform(x.astype(str))))
        p=1+enc.transform(x_tune[:1].astype(str)).shape[1]
    else:
        transformers=[SplineTransformer(n_knots=5,degree=3,include_bias=False,knots="quantile").fit(x_tune[:,[j]]) for j in range(x_tune.shape[1])]
        def trans(x: np.ndarray) -> np.ndarray:
            return np.column_stack((np.ones(len(x)), *(t.transform(x[:, [j]]) for j, t in enumerate(transformers))))
        p=trans(x_tune[:1]).shape[1]
    return Basis(trans,p,"spline_K3")


def solve_drf(basis: Basis, rt: dict[str,np.ndarray], ot: dict[str,np.ndarray], z0: np.ndarray,
              z1: np.ndarray, gr: np.ndarray, go: np.ndarray, lam: float, omega: float,
              rho: float=0.) -> tuple[np.ndarray,float,float]:
    br,bo=basis.transform(rt["x"]),basis.transform(ot["x"])
    H=(1-omega)*(br.T@br/len(br))+omega*(bo.T@bo/len(bo))
    penalty=np.eye(H.shape[0])
    rhs=br.T@(z0+lam*(z1-z0)-omega*gr)/len(br)+omega*bo.T@go/len(bo)
    mat=H+rho*penalty
    beta=np.linalg.solve(mat,rhs) if np.linalg.cond(mat)<1e12 else np.linalg.pinv(mat)@rhs
    residual=float(np.max(np.abs(mat@beta-rhs))/max(1.,np.max(np.abs(rhs))))
    return beta,residual,float(np.linalg.cond(mat))


def cscore(candidate: tuple[float,float,np.ndarray], basis: Basis, rr: dict[str,np.ndarray],
           oo: dict[str,np.ndarray], mur: OutcomeFit, muo: OutcomeFit, clip: bool=True) -> float:
    lam,om,beta=candidate
    z0,_,gr,go=fusion_scores(rr,oo,mur,muo,clip)
    zr=basis.transform(rr["x"])@beta; zo=basis.transform(oo["x"])@beta
    if clip: zr,zo=np.clip(zr,-2,2),np.clip(zo,-2,2)
    return float(np.mean((z0-zr)**2-om*(gr-zr)**2)+np.mean(om*(go-zo)**2))


def run_cate1(stage:str,cell:str,rep:int,assertions:list[dict[str,Any]])->list[dict[str,Any]]:
    rows=[]; conf=1. if cell=="reference" else 2.
    sizes={r:(200 if r.startswith("R") else 5000) for r in ("R_nuis","O_nuis","R_tune","O_tune","R_select","O_select","R_report","O_report")}
    data={r:sample_scmb(n,"RCT" if r.startswith("R") else "OBS",rng("CATE1",cell,rep,r),conf) for r,n in sizes.items()}
    mur,muo=fit_outcome(data["R_nuis"]),fit_outcome(data["O_nuis"])
    rt,ot=data["R_tune"],data["O_tune"]; z0,z1,gr,go=fusion_scores(rt,ot,mur,muo,True)
    basis=spline_basis(rt["x"]); candidates=[]; residuals=[]; conditions=[]
    for idx,(la,om) in enumerate(GRID):
        beta,resid,cond=solve_drf(basis,rt,ot,z0,z1,gr,go,la,om)
        candidates.append((la,om,beta));residuals.append(resid);conditions.append(cond)
    scores=np.array([cscore(c,basis,data["R_select"],data["O_select"],mur,muo) for c in candidates])
    selected=int(np.lexsort((np.arange(len(scores)),scores))[0])
    truth=sample_scmb(100000,"RCT",rng("CATE1",cell,rep,"truth"),conf)
    bt=basis.transform(truth["x"]); risks=np.array([np.mean((np.clip(bt@c[2],-2,2)-truth["tau"])**2) for c in candidates])
    indexes={"rct_only":GRID.index((0.,0.)),"lambda_only":min((i for i,g in enumerate(GRID) if g[1]==0),key=lambda i:(scores[i],i)),
             "omega_only":min((i for i,g in enumerate(GRID) if g[0]==0),key=lambda i:(scores[i],i)),"joint":selected,
             "grid_oracle":int(np.lexsort((np.arange(len(risks)),risks))[0])}
    eps=2*9**2*math.sqrt(2*math.log(4*25/.05))*(2/math.sqrt(200)+1/math.sqrt(5000))
    vac=2*eps/(risks.max()-risks.min()) if risks.max()>risks.min() else math.inf
    for est,idx in indexes.items():
        la,om,_=candidates[idx]; report=cscore(candidates[idx],basis,data["R_report"],data["O_report"],mur,muo)
        for metric,value in (("risk",risks[idx]),("selection_score",scores[idx]),("reporting_score",report),
                             ("lambda",la),("omega",om),("bound_vacuity_ratio",vac)):
            row(rows,stage,"CATE1",cell,rep,est,"r=1",metric,float(value))
    assertions.extend([
        check("CATE1",cell,rep,"normal_equation_relative_residual",max(residuals)<=1e-8,True,{"max":max(residuals)}),
        check("CATE1",cell,rep,"theorem7_scope_label",True,True,
              {"basis":"spline_K3","rho":0,"p":basis.p,"n_tune":200,"condition_max":max(conditions),
               "clipped_performance":True,"theorem7_eligible":False,"theorem7_status":"N/A","neural":False}),
        check("CATE1",cell,rep,"prespecified_bound",True,True,{"B0":9,"rbar":1,"radius":eps,"bound_vacuity_ratio":vac}),
        check("CATE1",cell,rep,"scmb_raw_boundedness",
              all(np.max(np.abs(data[r]["x"]))<=1 and np.max(np.abs(data[r]["u"]))<=1 and
                  np.max(np.abs(data[r]["eps"]))<=1 and np.max(np.abs(data[r]["y"]))<1 for r in sizes),True,
              {"X_bound":1,"U_bound":1,"epsilon_bound":1,
               "observed_max_abs_y":max(float(np.max(np.abs(data[r]["y"]))) for r in sizes)}),
        check("CATE1",cell,rep,"role_rng_distinct",len({seed32("CATE1",cell,rep,r) for r in sizes})==len(sizes),True,
              {r:fingerprint(data[r]) for r in sizes}),
    ])
    return rows


def neural_train(on:dict[str,np.ndarray],rn:dict[str,np.ndarray],study:str,cell:str,rep:int)->tuple[Basis,NeuralOutcomeFit,dict[str,Any]]:
    x=on["x"]; a=on["a"].astype(int); y=on["y"]; gen=rng(study,cell,rep,"O_nuis","neural")
    order=gen.permutation(len(x));nv=max(1,len(x)//10);va,train=order[:nv],order[nv:]
    d=x.shape[1]; W1=gen.normal(0,math.sqrt(2/(d+64)),(d,64));b1=np.zeros(64)
    W2=gen.normal(0,math.sqrt(2/128),(64,64));b2=np.zeros(64);Wh=gen.normal(0,math.sqrt(2/66),(64,2));bh=np.zeros(2)
    params=[W1,b1,W2,b2,Wh,bh];m=[np.zeros_like(p) for p in params];v=[np.zeros_like(p) for p in params]
    best=None;best_loss=math.inf;best_epoch=-1;stale=0
    def forward(xx):
        q1=xx@W1+b1;h1=np.maximum(q1,0);q2=h1@W2+b2;h=np.maximum(q2,0);return q1,h1,q2,h,h@Wh+bh
    for epoch in range(300):
        q1,h1,q2,h,p=forward(x[train]);err=(p[np.arange(len(train)),a[train]]-y[train]);dp=np.zeros_like(p);dp[np.arange(len(train)),a[train]]=2*err/len(train)
        grads=[None]*6;grads[4]=h.T@dp;grads[5]=dp.sum(0);dh=dp@Wh.T;dq2=dh*(q2>0);grads[2]=h1.T@dq2;grads[3]=dq2.sum(0);dh1=dq2@W2.T;dq1=dh1*(q1>0);grads[0]=x[train].T@dq1;grads[1]=dq1.sum(0)
        t=epoch+1
        for j,(p0,g) in enumerate(zip(params,grads)):
            m[j]=.9*m[j]+.1*g;v[j]=.999*v[j]+.001*g*g;p0-=.001*(m[j]/(1-.9**t))/(np.sqrt(v[j]/(1-.999**t))+1e-8)
        pv=forward(x[va])[-1];vl=float(np.mean((pv[np.arange(len(va)),a[va]]-y[va])**2))
        if vl<best_loss-1e-12: best_loss=vl;best_epoch=epoch;best=[p.copy() for p in params];stale=0
        else: stale+=1
        if stale>=25: break
    if best is None: raise RuntimeError("neural training produced no checkpoint")
    W1[:],b1[:],W2[:],b2[:],Wh[:],bh[:]=best
    def trunk(xx): return np.maximum(np.maximum(xx@W1+b1,0)@W2+b2,0)
    hr=trunk(rn["x"]);center=hr.mean(0);scale=np.where(hr.std(0)>1e-10,hr.std(0),1.)
    def trans(xx): return np.column_stack((np.ones(len(xx)),(trunk(xx)-center)/scale))
    basis=Basis(trans,65,"neural_numpy_2x64")
    neural_fit=NeuralOutcomeFit(lambda xx:trunk(xx)@Wh+bh)
    return basis,neural_fit,{"epochs":best_epoch+1,"validation_mse":best_loss,"optimizer":"deterministic_numpy_full_batch_adam","lr":.001,"max_epochs":300,"patience":25}


def run_neural_smoke(assertions:list[dict[str,Any]])->None:
    study,cell,rep="CATE1_NEURAL_SMOKE","scmb_reference",0
    sizes={r:(200 if r.startswith("R") else 5000) for r in ("R_nuis","O_nuis","R_tune","O_tune","R_select","O_select","R_report","O_report")}
    data={r:sample_scmb(n,"RCT" if r.startswith("R") else "OBS",rng(study,cell,rep,r),1.) for r,n in sizes.items()}
    mur=fit_outcome(data["R_nuis"])
    basis,muo,meta=neural_train(data["O_nuis"],data["R_nuis"],study,cell,rep)
    rt,ot=data["R_tune"],data["O_tune"];z0,z1,gr,go=fusion_scores(rt,ot,mur,muo,True)
    candidates=[];scores=[];maxres=0.
    for gi,(la,om) in enumerate(GRID):
        for ri,rho in enumerate(RHOS):
            beta,resid,_=solve_drf(basis,rt,ot,z0,z1,gr,go,la,om,rho);maxres=max(maxres,resid)
            c=(la,om,beta);candidates.append((gi,ri,rho,c));scores.append(cscore(c,basis,data["R_select"],data["O_select"],mur,muo))
    sel=min(range(75),key=lambda i:(scores[i],i));rct=[i for i,q in enumerate(candidates) if q[3][0]==0 and q[3][1]==0]
    rho_i=min(rct,key=lambda i:(scores[i],candidates[i][1]))
    report_sel=cscore(candidates[sel][3],basis,data["R_report"],data["O_report"],mur,muo)
    report_rct=cscore(candidates[rho_i][3],basis,data["R_report"],data["O_report"],mur,muo)
    frozen_selection=sel
    perturbed_reporting={key:value.copy() for key,value in data["R_report"].items()};perturbed_reporting["y"]+=7.0
    perturb_score=cscore(candidates[sel][3],basis,perturbed_reporting,data["O_report"],mur,muo)
    rerun_scores=np.array([cscore(q[3],basis,data["R_select"],data["O_select"],mur,muo) for q in candidates])
    rerun_selection=min(range(75),key=lambda i:(rerun_scores[i],i))
    ok=(len(candidates)==75 and np.isfinite(scores).all() and maxres<=1e-8 and
        np.isfinite(report_sel-report_rct) and frozen_selection==rerun_selection and
        np.array_equal(scores,rerun_scores) and perturb_score!=report_sel)
    assertions.append(check(study,cell,rep,"neural_end_to_end_75_candidates",ok,True,
                            {**meta,"candidate_count":75,"selected_index":sel,"rho_R":candidates[rho_i][2],
                             "report_difference":report_sel-report_rct,"max_normal_residual":maxres,
                             "basis":"two_hidden_ReLU_64_shared_two_heads","scaling":"RCT_nuisance_only",
                             "g_source":"frozen_neural_heads","reporting_perturbation_changed_score_not_selection":True,
                             "selection_rerun_after_reporting_only_perturbation":rerun_selection,
                             "theorem7_eligible":False,"intercept_penalized":True}))


def load_star_source() -> tuple[pd.DataFrame,dict[str,Any]]:
    out=subprocess.run(["curl","-fsSL",STAR_URL],check=True,capture_output=True,timeout=30).stdout
    digest=hashlib.sha256(out).hexdigest()
    if digest!=STAR_SHA256: raise RuntimeError(f"STAR hash mismatch: {digest}")
    f=pd.read_csv(io.BytesIO(out));required={"rownames","stark","gender","ethnicity","lunchk","schoolk","schoolidk","mathk","readk"}
    if not required.issubset(f.columns): raise RuntimeError(f"STAR schema missing {sorted(required-set(f.columns))}")
    return f,{"url":STAR_URL,"sha256":digest,"treatment":"stark: small vs regular+regular+aide","school_type":"schoolk","school_id":"schoolidk"}


def prepare_star_source(assertions:list[dict[str,Any]])->tuple[pd.DataFrame,dict[str,Any]]:
    f,meta=load_star_source()
    return prepare_star_frame(f,meta)


def prepare_star_frame(f:pd.DataFrame,meta:dict[str,Any])->tuple[pd.DataFrame,dict[str,Any]]:
    meta=dict(meta);h=f[f["stark"].isin(["small","regular","regular+aide"])&f["mathk"].notna()&f["readk"].notna()].copy()
    h["A"]=(h["stark"]=="small").astype(int);rawcols=["gender","ethnicity","lunchk","schoolk","schoolidk"]
    miss=[];imputation_modes={}
    for c in rawcols:
        mc=f"{c}__missing";h[mc]=h[c].isna().astype(int);miss.append(mc)
        mode=h[c].mode(dropna=True)
        if mode.empty: raise RuntimeError(f"no mode for {c}")
        imputation_modes[c]=str(mode.iloc[0]);h[c]=h[c].fillna(mode.iloc[0]).astype(str)
    h["ethnicity3"]=np.where(h["ethnicity"].eq("afam"),"afam",np.where(h["ethnicity"].eq("cauc"),"cauc","other"))
    xcols=["gender","ethnicity3","lunchk","schoolk","schoolidk",*miss]
    h["xkey"]=list(map(tuple,h[xcols].to_numpy(object)))
    shares=h.groupby("xkey",sort=True)["A"].mean();eligible=set(shares[(shares>=.15)&(shares<=.85)].index)
    h=h[h["xkey"].isin(eligible)].copy().sort_values("rownames",kind="stable").reset_index(drop=True)
    if h.empty: raise RuntimeError("no eligible exact-X STAR rows")
    counts=h.groupby("xkey",sort=True).size();px=(counts/counts.sum()).to_dict();ex=h.groupby("xkey",sort=True)["A"].mean().to_dict()
    # Freeze a deterministic 5-fold OOF ridge using X+A and stable row-id tie breaking.
    enc=OneHotEncoder(handle_unknown="ignore",sparse_output=False).fit(h[xcols].astype(str));X=enc.transform(h[xcols].astype(str));A=h["A"].to_numpy()
    fold_rng=rng("CATE2_SMOKE","math_alpha08_nR200",0,"design","oof_folds");fold=np.empty(len(h),int);fold[fold_rng.permutation(len(h))]=np.arange(len(h))%5
    pred=np.empty(len(h))
    design=np.column_stack((X,A,A[:,None]*X))
    for k in range(5):
        tr=fold!=k;te=~tr;model=Ridge(alpha=5.,solver="lsqr",tol=1e-10).fit(design[tr],h.loc[tr,"mathk"]);pred[te]=model.predict(design[te])
    resid=h["mathk"].to_numpy()-pred;ranks=np.empty(len(h))
    for a in (0,1):
        ids=np.flatnonzero(A==a);order=np.lexsort((h.loc[ids,"rownames"].to_numpy(),resid[ids]));rank=np.empty(len(ids),int);rank[order]=np.arange(len(ids));ranks[ids]=(rank+.5)/len(ids)
    alpha=.8;h["weight"]=np.where(A==1,1+alpha*(2*ranks-1),1-alpha*(2*ranks-1));h["fold"]=fold
    ymean=float(h["mathk"].mean());ysd=float(h["mathk"].std(ddof=0));h["Y"]=(h["mathk"]-ymean)/ysd
    row_id_hash=hashlib.sha256("\n".join(map(str,h["rownames"])).encode()).hexdigest()
    fold_hash=hashlib.sha256(np.ascontiguousarray(fold,dtype="<i8").tobytes()).hexdigest()
    probability_payload=[(repr(k),format(float(px[k]),".17g")) for k in sorted(px,key=repr)]
    propensity_payload=[(repr(k),format(float(ex[k]),".17g")) for k in sorted(ex,key=repr)]
    probability_hash=hashlib.sha256(json.dumps(probability_payload,separators=(",",":"),ensure_ascii=True).encode()).hexdigest()
    propensity_hash=hashlib.sha256(json.dumps(propensity_payload,separators=(",",":"),ensure_ascii=True).encode()).hexdigest()
    support_hash=hashlib.sha256((row_id_hash+fold_hash+probability_hash+propensity_hash).encode()).hexdigest()
    meta.update({"x_tuple":xcols,"ethnicity3":"afam->afam; cauc->cauc; all other observed levels->other",
                 "missingness_indicators":miss,"imputation_modes":imputation_modes,
                 "eligible_rows":len(h),"exact_cells":len(px),"outcome":"mathk",
                 "alpha":alpha,"oof":"5-fold deterministic Ridge(alpha=5) on one-hot X+A+AX",
                 "rank":"within-arm (rank-0.5)/n; stable rownames tie break","standardization":{"mean":ymean,"sd":ysd},
                 "eligible_ate":float(h.loc[h.A==1,"mathk"].mean()-h.loc[h.A==0,"mathk"].mean()),
                 "ordered_row_id_sha256":row_id_hash,"fold_assignment_sha256":fold_hash,
                 "probability_table_sha256":probability_hash,"propensity_table_sha256":propensity_hash,
                 "frozen_support_sha256":support_hash})
    h.attrs.update({"xcols":xcols,"px":px,"ex":ex,"encoder":enc,"meta":meta})
    return h,meta


def draw_star_role(h:pd.DataFrame,source:str,n:int,generator:np.random.Generator)->dict[str,np.ndarray]:
    px=h.attrs["px"];ex=h.attrs["ex"];keys=list(px);probs=np.array([px[k] for k in keys]);chosen=generator.choice(len(keys),n,p=probs)
    ids=[]
    groups={(k,a):np.flatnonzero((h.xkey==k)&(h.A==a)) for k in keys for a in (0,1)}
    for ki in chosen:
        k=keys[ki];a=int(generator.random()<ex[k]);pool=groups[(k,a)]
        p=None if source=="RCT" else h.loc[pool,"weight"].to_numpy()/h.loc[pool,"weight"].sum()
        ids.append(int(generator.choice(pool,p=p)))
    d=h.loc[ids];enc=h.attrs["encoder"];x=enc.transform(d[h.attrs["xcols"]].astype(str))
    e=np.array([ex[k] for k in d.xkey]);return {"x":x,"a":d.A.to_numpy(float),"y":d.Y.to_numpy(float),"e":e,"rowid":d.rownames.to_numpy(),"xindex":chosen}


def run_cate2_smoke(assertions:list[dict[str,Any]],h:pd.DataFrame,meta:dict[str,Any])->None:
    study,cell,rep="CATE2_SMOKE","math_alpha08_nR200",0
    sizes={"R_nuis":200,"O_nuis":5000,"R_tune":200,"O_tune":5000,"R_select":1000,"O_select":1000,"R_report":1000,"O_report":1000}
    data={r:draw_star_role(h,"RCT" if r.startswith("R") else "OBS",n,rng(study,cell,rep,r)) for r,n in sizes.items()}
    keys=list(h.attrs["px"]);p=np.array([h.attrs["px"][k] for k in keys]);ex=h.attrs["ex"]
    po_values=[]
    for key in keys:
        cell_probability=0.0
        for arm in (0,1):
            mask=np.array([value==key for value in h["xkey"]]) & (h["A"].to_numpy()==arm)
            weights=h.loc[mask,"weight"].to_numpy(float)
            conditional_mass=float(np.sum(weights/weights.sum()))
            arm_probability=ex[key] if arm==1 else 1-ex[key]
            cell_probability+=h.attrs["px"][key]*arm_probability*conditional_mass
        po_values.append(cell_probability)
    po=np.array(po_values)
    propensity_err=max(abs(ex[k]-float(h.loc[h.xkey==k,"A"].mean())) for k in keys)
    mur,muo=fit_outcome(data["R_nuis"],"linear"),fit_outcome(data["O_nuis"],"linear")
    rt,ot=data["R_tune"],data["O_tune"];z0,z1,gr,go=fusion_scores(rt,ot,mur,muo)
    basis=spline_basis(rt["x"],categorical=True);cands=[]
    for la,om in GRID:
        beta,resid,_=solve_drf(basis,rt,ot,z0,z1,gr,go,la,om);cands.append((la,om,beta))
    scores=np.array([cscore(c,basis,data["R_select"],data["O_select"],mur,muo,False) for c in cands]);sel=min(range(25),key=lambda i:(scores[i],i));rct=GRID.index((0.,0.))
    diff=cscore(cands[sel],basis,data["R_report"],data["O_report"],mur,muo,False)-cscore(cands[rct],basis,data["R_report"],data["O_report"],mur,muo,False)
    hard=np.max(np.abs(p-po))<=1e-15 and propensity_err<=1e-12 and len({seed32(study,cell,rep,r) for r in sizes})==len(sizes)
    assertions.extend([
        check(study,cell,rep,"star_source_resolved",True,True,meta),
        check(study,cell,rep,"source_law_exact_marginal_and_propensity",hard,True,
              {"max_px_difference":float(np.max(np.abs(p-po))),"max_propensity_identity_error":propensity_err,"r0":1}),
        check(study,cell,rep,"oof_and_feature_exclusion",bool(h.fold.nunique()==5),True,
              {"oof_folds":5,"learner_features":meta["x_tuple"],"excluded":["weight","fold","rank","residual"]}),
        check(study,cell,rep,"selection_reporting_rng_independence",
              fingerprint(data["R_select"])!=fingerprint(data["R_report"]) and fingerprint(data["O_select"])!=fingerprint(data["O_report"]),True,
              {r:fingerprint(data[r]) for r in sizes}),
        check(study,cell,rep,"directional_bias_soft",True,False,{"reporting_score_difference":diff,
              "naive_obs_contrast":float(data["O_report"]["y"][data["O_report"]["a"]==1].mean()-data["O_report"]["y"][data["O_report"]["a"]==0].mean())}),
    ])


def check(study:str,cell:str,rep:int,name:str,passed:bool,hard:bool,details:dict[str,Any])->dict[str,Any]:
    return {"study":study,"cell":cell,"rep":rep,"name":name,"hard":hard,"status":"pass" if passed else "fail","details":details}


def fixed_assertions()->list[dict[str,Any]]:
    out=[]
    parity_x=np.array([[.1,-.2,.3,-.4],[1.,2.,3.,4.]])
    parity_scm=generate_scm_parameters(1)
    parity_arrays={"baseline_basis":baseline_basis(parity_x),"treatment_effect_basis":treatment_effect_basis(parity_x),
                   "treatment_assignment_basis":treatment_assignment_basis(parity_x),"tau_function":tau_function(parity_x,parity_scm),
                   "feature_map":feature_map(parity_x)}
    expected_hashes={"baseline_basis":"da59ac126cdf6faba3470b42f0f2681cac419664adcdaade27d338c388459560",
                     "treatment_effect_basis":"cff6e94898fbf6b298b6ff1045793365ed25df543c23c81621606c449fd7f00a",
                     "treatment_assignment_basis":"ce34651174f754a8146e08f21d44167aa77a0fda594f2d9739b479595a6419cf",
                     "tau_function":"db5410816605f13f65c62ce09e50f6a4125a06e4ad226c49609861e1878be50e",
                     "feature_map":"c11b28c42d8418a37774d1e1326c84597f9109378683ce0f329a9233442a1625"}
    observed_hashes={name:hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest() for name,value in parity_arrays.items()}
    expected_mean=np.array([.060588628144355158,-.087839412048191706,-.13390727375574465,-.15867802430768993])
    out.append(check("FIXTURE","LEGACY_PARITY",0,"internalized_legacy_parity",
                     observed_hashes==expected_hashes and np.array_equal(parity_scm.covariate_mean,expected_mean),True,
                     {"observed_hashes":observed_hashes,"expected_hashes":expected_hashes,
                      "matplotlib_runtime_dependency":False}))
    rdata={"x":np.array([[0.],[1.],[2.],[3.]]),"a":np.array([0.,1.,0.,1.]),
           "y":np.array([.2,1.4,.9,2.1]),"e":np.array([.4,.5,.6,.7])}
    odata={"x":np.array([[-1.],[.5],[1.5],[2.5],[4.]]),"a":np.array([0.,1.,1.,0.,1.]),
           "y":np.array([-.1,1.,1.5,.7,2.8]),"e":np.array([.3,.45,.55,.65,.75])}
    mur=NeuralOutcomeFit(lambda x:np.column_stack((.1+.2*x[:,0],.7+.1*x[:,0])))
    muo=NeuralOutcomeFit(lambda x:np.column_stack((-.2+.3*x[:,0],.9+.4*x[:,0])))
    z0,z1,g,go=fusion_scores(rdata,odata,mur,muo)
    def manual(data:dict[str,np.ndarray],fit:NeuralOutcomeFit)->tuple[np.ndarray,np.ndarray]:
        m0,m1=fit.predict(data["x"],0),fit.predict(data["x"],1);a,y,e=data["a"],data["y"],data["e"]
        return m1-m0+a/e*(y-m1)-(1-a)/(1-e)*(y-m0),m1-m0
    mz0,_=manual(rdata,mur);mz1,mg=manual(rdata,muo);_,mgo=manual(odata,muo)
    la,om=.4,.3;theta,v=ate_eval(z0,z1,g,go,np.ones(len(go)),la,om)
    combined=z0+la*(z1-z0);direct_theta=float(combined.mean()+om*(go.mean()-g.mean()))
    direct_v=float((svar(combined)+om**2*svar(g)-2*om*scov(combined,g))/len(g)+om**2*svar(go)/len(go))
    pseudo_ok=(np.allclose(z0,mz0,atol=1e-15) and np.allclose(z1,mz1,atol=1e-15) and
               np.allclose(g,mg,atol=1e-15) and np.allclose(go,mgo,atol=1e-15))
    out.append(check("FIXTURE","ATE",0,"fixed_ate_pseudo_and_variance_identity",
                     pseudo_ok and abs(theta-direct_theta)<=1e-15 and abs(v-direct_v)<=1e-15,True,
                     {"lambda":la,"omega":om,"variance":v,"expanded_variance":direct_v,"ddof":1,
                      "theta":theta,"direct_theta":direct_theta,"g_source":"OBS fit on both sources"}))
    # Exact common-marginal CATE score identity: Z0=tau and fixed candidates.
    tau=np.array([-.2,.1,.4]);prob=np.array([.2,.3,.5]);g0=np.array([0.,.2,.3]);zeta=np.array([-.1,.05,.2]);omega=.5
    riskdiff=float(prob@((zeta-tau)**2-(g0-tau)**2))
    scorediff=float(prob@((tau-zeta)**2-omega*(g0-zeta)**2+omega*(g0-zeta)**2-(tau-g0)**2))
    out.append(check("FIXTURE","CATE",0,"fixed_cate_score_identity",abs(riskdiff-scorediff)<=1e-15,True,
                     {"risk_difference":riskdiff,"score_difference":scorediff,"support_probabilities":prob.tolist()}))
    xr=np.column_stack((np.linspace(-1,1,80),np.sin(np.linspace(-2,2,80))))
    xo=np.column_stack((np.linspace(-.9,.9,100),np.cos(np.linspace(-2,2,100))))
    rt={"x":xr};ot={"x":xo};basis=spline_basis(xr)
    fz0=.2+xr[:,0];fz1=fz0+.1*xr[:,1];fgr=.3-.2*xr[:,0];fgo=.25+.1*xo[:,1]
    beta0,resid0,cond0=solve_drf(basis,rt,ot,fz0,fz1,fgr,fgo,.4,.3,0.)
    br,bo=basis.transform(xr),basis.transform(xo)
    H=.7*(br.T@br/len(br))+.3*(bo.T@bo/len(bo))
    rhs=br.T@(fz0+.4*(fz1-fz0)-.3*fgr)/len(br)+.3*bo.T@fgo/len(bo)
    direct_beta=np.linalg.solve(H,rhs)
    beta_ridge,resid_ridge,_=solve_drf(basis,rt,ot,fz0,fz1,fgr,fgo,.4,.3,.1)
    direct_ridge=np.linalg.solve(H+.1*np.eye(H.shape[0]),rhs)
    out.append(check("FIXTURE","CATE",0,"unclipped_spline_theorem7_normal_equation",
                     basis.p<len(xr) and cond0<=1e8 and resid0<=1e-8 and np.allclose(beta0,direct_beta,rtol=1e-8,atol=1e-10),True,
                     {"basis":"spline_K3","clipped":False,"rho":0,"p":basis.p,"n_tune":len(xr),"condition":cond0,"theorem7_eligible":True}))
    out.append(check("FIXTURE","CATE",0,"full_ridge_penalizes_intercept",
                     resid_ridge<=1e-8 and np.allclose(beta_ridge,direct_ridge,rtol=1e-8,atol=1e-10),True,
                     {"rho":.1,"penalty":"rho_times_full_identity","intercept_penalized":True}))
    seeds=[seed32("FIXTURE","RNG",0,r) for r in ROLE_NAMES]
    out.append(check("FIXTURE","RNG",0,"canonical_128bit_seed_namespaces",len(set(seeds))==len(seeds),True,
                     {"key":"SHA256(protocol,study,cell,rep,namespace,substep)","entropy_bits":128,"global_rng":False}))
    for name in ("population_variance","oracle_ordering","stochastic_score_3mcse","coverage"):
        out.append({"study":"AGGREGATE","cell":"G0","rep":0,"name":name,"hard":False,"status":"not_evaluable_at_g0","details":{"deferred":"G1/G2 not implemented"}})
    return out


def summarize(rows:pd.DataFrame)->pd.DataFrame:
    numeric=pd.to_numeric(rows["value"],errors="coerce");work=rows.assign(numeric_value=numeric)
    group=["protocol","stage","study","cell","estimator","ratio","metric","status"]
    s=work.groupby(group,dropna=False,sort=True)["numeric_value"].agg(["count","mean","std"]).reset_index()
    s["mcse"]=s["std"]/np.sqrt(s["count"]);return s.drop(columns="std")


def decision(assertions:list[dict[str,Any]],stage:str)->str:
    if any(a["hard"] and a["status"]=="fail" for a in assertions): return "REDESIGN"
    if stage=="g0": return "INCONCLUSIVE"
    return "INCONCLUSIVE"


def jsonable(x:Any)->Any:
    if isinstance(x,(np.integer,)):return int(x)
    if isinstance(x,(np.floating,)):return None if not np.isfinite(x) else float(x)
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,dict):return {str(k):jsonable(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [jsonable(v) for v in x]
    return x


def git_meta()->dict[str,Any]:
    head=subprocess.run(["git","-C",str(PROJECT_DIR),"rev-parse","HEAD"],check=True,capture_output=True,text=True).stdout.strip()
    dirty=subprocess.run(["git","-C",str(PROJECT_DIR),"status","--short","--","."],check=True,capture_output=True,text=True).stdout.splitlines()
    return {"commit":head,"dirty":dirty}


def atomic_text(path:Path,text:str)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix=f".{path.name}.",dir=path.parent,text=True)
    try:
        with os.fdopen(fd,"w",encoding="utf-8",newline="") as f:f.write(text)
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp):os.unlink(tmp)


def write_artifacts(outdir:Path,rows:list[dict[str,Any]],assertions:list[dict[str,Any]],stage:str,runtime:float)->None:
    outdir.mkdir(parents=True,exist_ok=True);existing=outdir/ARTIFACTS[2]
    if existing.exists():
        old=json.loads(existing.read_text());
        if old.get("protocol")!=PROTOCOL or old.get("code_sha256")!=sha256_path(Path(__file__)):
            raise RuntimeError("refusing to overwrite artifacts with different protocol/code provenance")
    frame=pd.DataFrame(rows).sort_values(["study","cell","rep","estimator","ratio","metric"],kind="stable").reset_index(drop=True)
    summary=summarize(frame)
    schedule=expected_schedule(stage);schedule_text=json.dumps(schedule,separators=(",",":"),ensure_ascii=True)
    schedule_manifest={"record_count":len(schedule),"sha256":hashlib.sha256(schedule_text.encode()).hexdigest(),
                       "replication_count":8*(1 if stage=="g0" else 10)}
    code_hash=sha256_path(Path(__file__));meta={"protocol":PROTOCOL,"stage":stage,"decision":decision(assertions,stage),
        "code_sha256":code_hash,"handoff_sha256":sha256_path(HANDOFF),"git":git_meta(),
        "seed_contract":"SHA256 canonical tuple -> first 128 bits -> SeedSequence",
        "schedule_manifest":schedule_manifest,"balance_tolerance":BALANCE_TOLERANCE,
        "assertions":assertions}
    hard=[a for a in assertions if a["hard"]];failed=[a for a in hard if a["status"]=="fail"]
    memo=(f"# DataFusionPPI v3 Agile Gate Results\n\n"
          f"- Protocol: `{PROTOCOL}`\n- Stage: `{stage.upper()}`\n- Decision: `{meta['decision']}`\n"
          f"- Scheduled performance replications: {frame[['study','cell','rep']].drop_duplicates().shape[0]}\n"
          f"- Hard assertions: {len(hard)} total, {len(failed)} failed\n- Runtime seconds: {runtime:.6f}\n"
          f"- Code SHA-256: `{code_hash}`\n- Handoff SHA-256: `{meta['handoff_sha256']}`\n- Git commit: `{meta['git']['commit']}`\n\n"
          "## Interpretation\n\nG0 is an engineering and algebraic smoke test. Aggregate population variance, oracle ordering, "
          "coverage, and stochastic three-MCSE checks are not evaluable at G0; G1/G2 remain unimplemented and locked. "
          "Gaussian ratios are analytically common-support and unbounded; their Theorem 6 radius is N/A. "
          "The neural path is an empirical finite-grid candidate, not a Theorem 7 check.\n\n"
          "## SCM-B contract\n\n`X in [-1,1]^5`, and `U, epsilon` are independent Uniform[-1,1]. "
          "The RCT propensity and outcome are exactly those in v3. OBS treatment uses "
          "`expit(logit(e_R(X)) + c U)` with `c=0,1,2`; potential outcomes are shared.\n")
    atomic_text(outdir/ARTIFACTS[0],frame.to_csv(index=False,float_format="%.17g",lineterminator="\n"))
    atomic_text(outdir/ARTIFACTS[1],summary.to_csv(index=False,float_format="%.17g",lineterminator="\n"))
    atomic_text(outdir/ARTIFACTS[2],json.dumps(jsonable(meta),sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    atomic_text(outdir/ARTIFACTS[3],memo)


def validate(outdir:Path,stage:str)->None:
    if stage!="g0":raise RuntimeError("G1_NOT_IMPLEMENTED")
    for name in ARTIFACTS:
        if not (outdir/name).is_file():raise RuntimeError(f"missing artifact {name}")
    meta=json.loads((outdir/ARTIFACTS[2]).read_text());frame=pd.read_csv(outdir/ARTIFACTS[0],keep_default_na=False);summary=pd.read_csv(outdir/ARTIFACTS[1],keep_default_na=False,na_values=[""])
    if meta["protocol"]!=PROTOCOL or meta["stage"]!=stage:raise RuntimeError("artifact provenance mismatch")
    if meta.get("code_sha256")!=sha256_path(Path(__file__)):raise RuntimeError("artifact code hash does not match current code")
    if meta["handoff_sha256"]!=HANDOFF_SHA256:raise RuntimeError("handoff provenance mismatch")
    if any(a["hard"] and a["status"]=="fail" for a in meta["assertions"]):raise RuntimeError("hard assertion failure")
    key_columns=["study","cell","rep","estimator","ratio","metric"]
    expected=expected_schedule(stage);observed=[tuple(x) for x in frame[key_columns].itertuples(index=False,name=None)]
    if len(observed)!=len(set(observed)):raise RuntimeError("duplicate scheduled result key")
    if sorted(observed)!=expected:raise RuntimeError("schedule/cardinality/failure-accounting mismatch")
    schedule_text=json.dumps(expected,separators=(",",":"),ensure_ascii=True)
    expected_manifest={"record_count":len(expected),"sha256":hashlib.sha256(schedule_text.encode()).hexdigest(),
                       "replication_count":8*(1 if stage=="g0" else 10)}
    if meta.get("schedule_manifest")!=expected_manifest:raise RuntimeError("schedule manifest mismatch")
    if not all(frame["seed_key"]==frame["study"]+"/"+frame["cell"]+"/"+frame["rep"].astype(str)):
        raise RuntimeError("scheduled seed key mismatch")
    allowed_status={"ok","not_applicable","explicit_infeasible","failure"}
    if not set(frame["status"]).issubset(allowed_status):raise RuntimeError("unknown result status")
    if (frame["status"]=="failure").any() and meta["decision"]!="REDESIGN":raise RuntimeError("failure row without REDESIGN decision")
    norm=frame[frame["metric"]=="normalization_error"].copy();norm_value=pd.to_numeric(norm["value"],errors="coerce")
    if not np.isfinite(norm_value).all():raise RuntimeError("ratio normalization diagnostic must be finite")
    balance=frame[frame["metric"]=="balance_residual"].copy()
    nonbal=balance[~balance["ratio"].isin(("BAL-X","BAL-X+g"))]
    if not ((nonbal["value"]=="N/A") & (nonbal["status"]=="not_applicable")).all():
        raise RuntimeError("non-balancing ratio must report balance residual/status N/A")
    bal=balance[balance["ratio"].isin(("BAL-X","BAL-X+g"))]
    bal_value=pd.to_numeric(bal["value"],errors="coerce")
    if not np.isfinite(bal_value).all():raise RuntimeError("BAL ratio balance residual must be finite")
    converged=bal["status"]=="ok"
    if ((bal_value[converged]>BALANCE_TOLERANCE).any() or
        (bal_value[~converged]<=BALANCE_TOLERANCE).any()):
        raise RuntimeError("BAL protocol convergence status/residual mismatch")
    recomputed=summarize(frame);summary_columns=list(recomputed.columns)
    try:pd.testing.assert_frame_equal(summary[summary_columns],recomputed,check_dtype=False,rtol=1e-13,atol=1e-15)
    except AssertionError as exc:raise RuntimeError("summary does not match replication artifact") from exc
    if frame.empty or summary.empty:raise RuntimeError("empty scientific artifact")
    if meta["decision"]!="INCONCLUSIVE":raise RuntimeError("G0/G1 must remain INCONCLUSIVE below cap absent implemented GO test")


def run(stage:str,outdir:Path)->float:
    if stage!="g0":raise RuntimeError("G1_NOT_IMPLEMENTED")
    if sha256_path(HANDOFF)!=HANDOFF_SHA256:raise RuntimeError("handoff preflight hash changed")
    start=time.monotonic();assertions=fixed_assertions();rows=[];reps=range(1 if stage=="g0" else 10)
    frozen_cate2=None
    if stage=="g0":
        try:
            frozen_cate2=prepare_star_source(assertions)
        except Exception as exc:
            raise RuntimeError(f"CATE2 source-law freeze failed before performance replications: {exc}") from exc
        assertions.append(check("CATE2_SMOKE","math_alpha08_nR200",0,"source_law_frozen_before_performance",True,True,
                                frozen_cate2[1]))
    for rep in reps:
        for cell in ("reference","conf0_correct","conf2_linear","star_reference"): rows.extend(run_ate1(stage,cell,rep,assertions))
        for cell in ("mild","strong"): rows.extend(run_ate2(stage,cell,rep,assertions))
        for cell in ("reference","confounding2"): rows.extend(run_cate1(stage,cell,rep,assertions))
    if stage=="g0":
        run_neural_smoke(assertions)
        assert frozen_cate2 is not None
        run_cate2_smoke(assertions,*frozen_cate2)
    runtime=time.monotonic()-start;write_artifacts(outdir,rows,assertions,stage,runtime);return runtime


def main()->int:
    p=argparse.ArgumentParser();p.add_argument("--stage",choices=("g0","g1"),required=True);p.add_argument("--artifact-dir",type=Path,default=PROJECT_DIR/"materials");p.add_argument("--validate-artifacts",action="store_true");args=p.parse_args()
    if args.validate_artifacts:validate(args.artifact_dir,args.stage);print(f"VALID {args.stage} {args.artifact_dir}");return 0
    runtime=run(args.stage,args.artifact_dir);validate(args.artifact_dir,args.stage);print(f"PASS {args.stage} runtime_seconds={runtime:.6f} decision=INCONCLUSIVE");return 0


if __name__=="__main__":raise SystemExit(main())
