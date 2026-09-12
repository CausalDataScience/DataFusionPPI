"""Conditional-variance mini-experiment.

The exact variance of Theorem 1 is conditional on the fitted objects and on the
coefficients.  A study that refits everything in every replication therefore
cannot check that formula directly.  Here the nuisances, the transport weight
and the coefficients are fitted once per training draw and then frozen, and only
the evaluation data is regenerated, so the replicated spread is exactly the
quantity the theorem describes.

Every coefficient pair sees the same evaluation draw inside a replication, and
the omega comparisons are paired on lambda by construction.
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd

import fusion_data as fd
from fusion_core import (CI_Z, algorithm1, ate_estimate, build_scores, evaluation_variance,
                         fit_outcome, ghat_on, ratio_classifier, ratio_unit, rng_for,
                         stable_seed, theorem1_variance, tuning_moments, trim_propensity,
                         var, variance_standard_error)
from fusion_io import summarize, write_outputs

REFERENCE = {"n_rct": 100, "n_obs": 5000, "confounding": 1.0, "spec": "flexible"}
N_EVAL_RCT = 30       # the evaluation share of the reference trial under the honest split
N_EVAL_OBS = 1000


def frozen_fit(scm, regime: str, cfg: dict, seed: int, ratio_builder, oracle_fn):
    """Fit the nuisances, the weight and the coefficients once, then freeze them."""
    nuis_r = fd.sample_scm(scm, "RCT", int(0.4 * cfg["n_rct"]), regime, rng_for(seed, 1))
    tune_r = fd.sample_scm(scm, "RCT", int(0.3 * cfg["n_rct"]), regime, rng_for(seed, 2))
    nuis_o = fd.sample_scm(scm, "OBS", int(0.6 * cfg["n_obs"]), regime, rng_for(seed, 3))
    tune_o = fd.sample_scm(scm, "OBS", int(0.2 * cfg["n_obs"]), regime, rng_for(seed, 4))
    nuis_r, _ = trim_propensity(nuis_r)
    tune_r, _ = trim_propensity(tune_r)

    mu_r = fit_outcome(nuis_r, spec="flexible")
    mu_o = fit_outcome(nuis_o, spec=cfg["spec"])
    ratio = ratio_builder(nuis_r, nuis_o, seed=seed, oracle_fn=oracle_fn)

    s_tune = build_scores(tune_r, mu_r, mu_o)
    g_tune_o = ghat_on(tune_o["x"], mu_o)
    moments = tuning_moments(s_tune, g_tune_o, ratio(tune_o["x"]), N_EVAL_RCT, N_EVAL_OBS)
    coef = algorithm1(moments)

    big_r = fd.sample_scm(scm, "RCT", 200000, regime, rng_for(seed, 5))
    big_o = fd.sample_scm(scm, "OBS", 200000, regime, rng_for(seed, 6))
    s_big = build_scores(big_r, mu_r, mu_o)
    g_big = ghat_on(big_o["x"], mu_o)
    r_big = ratio(big_o["x"])
    a = var(s_big.delta)
    b = var(s_big.ghat) + (N_EVAL_RCT / N_EVAL_OBS) * var(r_big * g_big)
    c = float(np.cov(s_big.z0, s_big.delta, ddof=1)[0, 1])
    d = float(np.cov(s_big.z0, s_big.ghat, ddof=1)[0, 1])
    oracle = {"lambda": float(np.clip(-c / a, 0, 1)) if a > 1e-12 else 0.0,
              "omega": float(np.clip(d / b, 0, 1)) if b > 1e-12 else 0.0}

    def theorem1(lam: float, om: float) -> float:
        return theorem1_variance(s_big.z0, s_big.delta, s_big.ghat, r_big * g_big,
                                 lam, om, N_EVAL_RCT, N_EVAL_OBS)

    return {"mu_r": mu_r, "mu_o": mu_o, "ratio": ratio, "coef": coef, "oracle": oracle,
            "theorem1": theorem1}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--training-draws", type=int, default=5)
    parser.add_argument("--evaluations", type=int, default=100)
    parser.add_argument("--stem", type=str, default="conditional_variance")
    parser.add_argument("--regime", type=str, default="shared")
    args = parser.parse_args()

    rows, checks, started = [], [], time.time()
    builder = ratio_unit if args.regime == "shared" else ratio_classifier
    for scm_id in (1, 2, 3):
        scm = fd.scm_family(scm_id, confounding=REFERENCE["confounding"])
        theta = fd.scm_true_theta(scm, args.regime)
        oracle_fn = fd.scm_oracle_ratio(scm, args.regime)
        for draw in range(args.training_draws):
            seed = stable_seed("cv", scm_id, draw)
            frozen = frozen_fit(scm, args.regime, REFERENCE, seed, builder, oracle_fn)
            lam, om = frozen["coef"]["lambda"], frozen["coef"]["omega"]
            olam, oom = frozen["oracle"]["lambda"], frozen["oracle"]["omega"]
            pairs = {"rct_only": (0.0, 0.0), "omega_only": (0.0, om),
                     "lambda_only": (lam, 0.0), "joint": (lam, om),
                     "oracle_omega_only": (0.0, oom), "oracle_lambda_only": (olam, 0.0),
                     "oracle_joint": (olam, oom)}
            for rep in range(args.evaluations):
                s = stable_seed("cv-eval", scm_id, draw, rep)
                eval_r = fd.sample_scm(scm, "RCT", N_EVAL_RCT, args.regime, rng_for(s, 1))
                eval_o = fd.sample_scm(scm, "OBS", N_EVAL_OBS, args.regime, rng_for(s, 2))
                eval_r, _ = trim_propensity(eval_r)
                sc = build_scores(eval_r, frozen["mu_r"], frozen["mu_o"])
                g_o = ghat_on(eval_o["x"], frozen["mu_o"])
                rg = frozen["ratio"](eval_o["x"]) * g_o
                for name, (l, w) in pairs.items():
                    z = sc.z0 + l * sc.delta
                    est = ate_estimate(z, sc.ghat, rg, w)
                    v = evaluation_variance(z, sc.ghat, rg, w)
                    half = CI_Z * np.sqrt(max(v, 0.0))
                    rows.append({"family": f"scm{scm_id}", "training_draw": draw,
                                 "replication": rep, "estimator": name,
                                 "estimate": est, "error": est - theta,
                                 "variance_hat": v, "ci_length": 2 * half,
                                 "covered": float(abs(est - theta) <= half),
                                 "lambda_used": l, "omega_used": w,
                                 "theorem1_variance": frozen["theorem1"](l, w),
                                 "failed": float(not np.isfinite(est)),
                                 "true_theta": theta, "seed": s})
            checks.append({"family": f"scm{scm_id}", "training_draw": draw,
                           "lambda_hat": lam, "omega_hat": om,
                           "oracle_lambda": olam, "oracle_omega": oom})
            print(f"  scm{scm_id} training draw {draw}: lambda {lam:.3f} omega {om:.3f}, "
                  f"oracle {olam:.3f} {oom:.3f}  [{time.time() - started:.0f}s]", flush=True)

    frame = pd.DataFrame(rows)
    # acceptance: the frozen coefficients never move across evaluation draws
    moved = frame.groupby(["family", "training_draw", "estimator"])[["lambda_used", "omega_used"]].nunique()
    assert int(moved.max().max()) == 1, "a frozen coefficient changed across evaluation draws"

    group = ["family", "training_draw", "estimator"]
    agg = []
    for keys, block in frame.groupby(group):
        est = block["estimate"].to_numpy(dtype=float)
        k = len(est)
        empirical = float(np.var(est, ddof=1))
        theory = float(block["theorem1_variance"].iloc[0])
        se = variance_standard_error(est)
        agg.append({**dict(zip(group, keys)),
                    "empirical_variance": empirical, "theorem1_variance": theory,
                    "gap": empirical - theory, "monte_carlo_se": se,
                    "z_gap": (empirical - theory) / se if se > 0 else np.nan,
                    "mean_estimated_variance": float(block["variance_hat"].mean()),
                    "coverage": float(block["covered"].mean()),
                    "bias": float(block["error"].mean()),
                    "replications": k})
    summary = pd.DataFrame(agg)
    paths = write_outputs(args.stem, frame, summary)
    pd.DataFrame(checks).to_csv(paths["summary"].parent / f"{args.stem}_frozen_coefficients.csv",
                                index=False)
    within = int((summary["z_gap"].abs() <= 3).sum())
    print(f"wrote {paths['replications'].name}; {len(frame)} rows in {time.time() - started:.0f}s; "
          f"exact variance matched in {within}/{len(summary)} frozen fits")
    return 0


if __name__ == "__main__":
    sys.exit(main())
