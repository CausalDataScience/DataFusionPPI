"""ATE-1: when fusion helps, and whether Algorithm 1 finds it (handoff 5.1).

Checks Theorems 1, 2, 2.1, 3, and 4.  Three axes around a reference
configuration give 11 configurations per data family and 44 cells.
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd

import fusion_data as fd
from fusion_ate import ate_replication, exact_variance_check, RCT_FRACTIONS, OBS_FRACTIONS
from fusion_core import (algorithm1, ate_estimate, build_scores, evaluation_variance,
                         fit_outcome, ghat_on, ratio_unit, rng_for, split_roles,
                         stable_seed, take, tuning_moments)
from fusion_io import summarize, write_outputs

FAMILIES = ("scm1", "scm2", "scm3", "star")
REFERENCE = {"n_rct": 100, "n_obs": 5000, "confounding": 1.0, "spec": "flexible"}


def configurations() -> list[dict]:
    out = []
    for n_rct in (50, 100, 200, 400):
        out.append({**REFERENCE, "n_rct": n_rct, "axis": "rct_size"})
    for n_obs in (1000, 20000):
        out.append({**REFERENCE, "n_obs": n_obs, "axis": "obs_size"})
    for conf, spec in ((0.0, "flexible"), (2.0, "flexible"),
                       (0.0, "linear"), (1.0, "linear"), (2.0, "linear")):
        out.append({**REFERENCE, "confounding": conf, "spec": spec, "axis": "quality"})
    return out


def build_cell(family: str, cfg: dict, star_support=None):
    """Return (draw, sampler, theta, test_sampler) for one design cell."""
    if family.startswith("scm"):
        scm = fd.scm_family(int(family[-1]), confounding=cfg["confounding"])
        theta = fd.scm_true_theta(scm, "shared")

        def draw(rep: int):
            s = stable_seed(family, tuple(sorted(cfg.items())), rep)
            return (fd.sample_scm(scm, "RCT", cfg["n_rct"], "shared", rng_for(s, 1)),
                    fd.sample_scm(scm, "OBS", cfg["n_obs"], "shared", rng_for(s, 2)))

        def sampler(n: int, seed: int):
            return (fd.sample_scm(scm, "RCT", n, "shared", rng_for(seed, 11)),
                    fd.sample_scm(scm, "OBS", n, "shared", rng_for(seed, 12)))

        return draw, sampler, theta

    support = star_support

    def draw(rep: int):
        rct, obs, _ = fd.star_realx_sample(support, "shared", cfg["n_rct"], cfg["n_obs"], rep)
        return rct, obs

    def sampler(n: int, seed: int):
        rct, obs, _ = fd.star_realx_sample(support, "shared", n, n, int(seed) % 100000)
        return rct, obs

    _, _, theta = fd.star_realx_sample(support, "shared", 50, 100, 0)
    return draw, sampler, theta


def cross_fitted_joint(rct, obs, theta, spec, seed, folds: int = 5):
    """The scheme used by the existing benchmark, kept as a comparison panel."""
    n, n_obs = len(rct["x"]), len(obs["x"])
    r_idx = rng_for(seed, 31).permutation(n)
    o_idx = rng_for(seed, 32).permutation(n_obs)
    r_folds = np.array_split(r_idx, folds)
    o_folds = np.array_split(o_idx, folds)
    estimates, variances = [], []
    for k in range(folds):
        r_eval, o_eval = take(rct, r_folds[k]), take(obs, o_folds[k])
        r_rest = take(rct, np.setdiff1d(r_idx, r_folds[k]))
        o_rest = take(obs, np.setdiff1d(o_idx, o_folds[k]))
        half_r = len(r_rest["x"]) // 2
        half_o = len(o_rest["x"]) // 2
        mu_r = fit_outcome(take(r_rest, np.arange(half_r)), spec="flexible")
        mu_o = fit_outcome(take(o_rest, np.arange(half_o)), spec=spec)
        s_tune = build_scores(take(r_rest, np.arange(half_r, len(r_rest["x"]))), mu_r, mu_o)
        g_tune = ghat_on(o_rest["x"][half_o:], mu_o)
        m = tuning_moments(s_tune, g_tune, np.ones(len(g_tune)), len(r_eval["x"]), len(o_eval["x"]))
        coef = algorithm1(m)
        s_eval = build_scores(r_eval, mu_r, mu_o)
        g_eval = ghat_on(o_eval["x"], mu_o)
        z = s_eval.z0 + coef["lambda"] * s_eval.delta
        estimates.append(ate_estimate(z, s_eval.ghat, g_eval, coef["omega"]))
        variances.append(evaluation_variance(z, s_eval.ghat, g_eval, coef["omega"]))
    est = float(np.mean(estimates))
    v = float(np.mean(variances)) / folds
    half = 1.959963984540054 * np.sqrt(max(v, 0.0))
    return {"estimator": "joint_cross_fitted", "estimate": est, "error": est - theta,
            "variance_hat": v, "ci_length": 2 * half,
            "covered": float(abs(est - theta) <= half),
            "lambda_used": np.nan, "omega_used": np.nan, "split": "cross_fitted"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replications", type=int, default=100)
    parser.add_argument("--families", type=str, default=",".join(FAMILIES))
    parser.add_argument("--stem", type=str, default="ate1_when_fusion_helps")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()

    families = tuple(args.families.split(","))
    configs = configurations()
    if args.smoke:
        configs = configs[:2]
    star_support = fd.star_realx_support() if "star" in families else None

    rows, checks = [], []
    started = time.time()
    for family in families:
        for cfg in configs:
            draw, sampler, theta = build_cell(family, cfg, star_support)
            label = {"family": family, "axis": cfg["axis"], "n_rct": cfg["n_rct"],
                     "n_obs": cfg["n_obs"], "confounding": cfg["confounding"],
                     "spec": cfg["spec"], "true_theta": theta}
            for rep in range(args.replications):
                rct, obs = draw(rep)
                seed = stable_seed(family, tuple(sorted(cfg.items())), rep)
                got = ate_replication(rct, obs, theta, ratio_unit, "unit", cfg["spec"],
                                      seed, oracle_sampler=sampler)
                if cfg["axis"] == "rct_size" and cfg["n_rct"] == REFERENCE["n_rct"]:
                    got.append(cross_fitted_joint(rct, obs, theta, cfg["spec"], seed))
                for row in got:
                    row.update(label)
                    row["replication"] = rep
                    row["seed"] = seed
                rows.extend(got)
            if family.startswith("scm"):
                check = exact_variance_check(sampler, theta, cfg["spec"],
                                             stable_seed(family, cfg["axis"], "ev"),
                                             cfg["n_rct"], cfg["n_obs"],
                                             repeats=200 if not args.smoke else 40)
                checks.append({**label, **check})
            print(f"  {family} {cfg['axis']} n_R={cfg['n_rct']} N_O={cfg['n_obs']} "
                  f"conf={cfg['confounding']} spec={cfg['spec']}  "
                  f"[{time.time() - started:.0f}s]", flush=True)

    frame = pd.DataFrame(rows)
    group = ["family", "axis", "n_rct", "n_obs", "confounding", "spec", "estimator"]
    summary = summarize(frame, group,
                        {"rmse": "error", "bias": "error", "coverage": "covered",
                         "ci_length": "ci_length", "variance_hat": "variance_hat",
                         "lambda_used": "lambda_used", "omega_used": "omega_used"})
    ratios = []
    for keys, block in summary[summary["metric"] == "rmse"].groupby(group[:-1]):
        base = block[block["estimator"] == "rct_only"]["value"]
        if len(base) == 0 or float(base.iloc[0]) == 0:
            continue
        for _, r in block.iterrows():
            ratios.append({**r.to_dict(), "metric": "rmse_ratio_to_rct",
                           "value": float(r["value"]) / float(base.iloc[0])})
    summary = pd.concat([summary, pd.DataFrame(ratios)], ignore_index=True)
    paths = write_outputs(args.stem, frame, summary)
    if checks:
        pd.DataFrame(checks).to_csv(paths["summary"].parent / f"{args.stem}_variance_check.csv", index=False)
    print(f"wrote {paths['replications'].name}, {paths['summary'].name}; "
          f"{len(frame)} rows in {time.time() - started:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
