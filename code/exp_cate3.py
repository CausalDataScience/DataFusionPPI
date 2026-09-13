"""CATE-3: the conditional effect when the two covariate laws differ.

CATE-1 draws both sources in the ``shared`` regime, where the transport ratio
r0(x) = p_R(x) / p_O(x) is exactly one, and hands the pipeline the constant one.
This study is the conditional counterpart of the average-effect shift study in
``exp_ate2.py``: the trial covariate mean is shifted away from the observational
mean, so r0 is not one, and the same five ratio routes compete.

The target is the conditional effect on the TRIAL population, so the test sample
and the anchor are drawn with the shifted trial mean.

Run: python3 exp_cate3.py --replications 100 --stem cate3_shift
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd

import fusion_cate as fc
import fusion_data as fd
from exp_ate2 import SHIFT_LEVELS, ratio_specs
from fusion_core import rng_for, stable_seed
from fusion_io import summarize, write_outputs

REGIME = "shifted"
FAMILIES = ("scm1", "scm2", "scm3")
N_OBS = 5000
TEST_DRAWS = 50000
ANCHOR_SIZE = 500


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replications", type=int, default=100)
    parser.add_argument("--families", type=str, default=",".join(FAMILIES))
    parser.add_argument("--sieves", type=str, default="spline3")
    parser.add_argument("--sizes", type=str, default="200")
    parser.add_argument("--levels", type=str, default=",".join(SHIFT_LEVELS))
    parser.add_argument("--stem", type=str, default="cate3_shift")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()

    families = tuple(args.families.split(","))
    sieves = tuple(args.sieves.split(","))
    sizes = tuple(int(v) for v in args.sizes.split(","))
    levels = tuple(args.levels.split(","))
    routes = ratio_specs()
    if args.smoke:
        families, sieves, sizes, levels = families[:1], sieves[:1], sizes[:1], levels[:1]
        routes = routes[:2]
        args.replications = min(args.replications, 3)

    rows, decomposition, started = [], [], time.time()
    for family in families:
        scm_id = int(family[-1])
        for level in levels:
            scale = SHIFT_LEVELS[level]
            scm = fd.scm_family(scm_id, shift=scale)
            oracle_fn = fd.scm_oracle_ratio(scm, REGIME)
            overlap = fd.overlap_check(scm, REGIME)
            # the estimand lives on the trial population, so the test covariates
            # carry the shifted trial mean
            test = fd.scm_test_sample(scm, REGIME, TEST_DRAWS, rng_for(9100, scm_id))
            anchor_x = test["x"][:ANCHOR_SIZE]
            anchor_tau = test["tau"][:ANCHOR_SIZE]

            for n_rct in sizes:
                def sampler(n, seed, _scm=scm):
                    return (fd.sample_scm(_scm, "RCT", n, REGIME, rng_for(seed, 11)),
                            fd.sample_scm(_scm, "OBS", n, REGIME, rng_for(seed, 12)))

                label = {"family": family, "shift": level, "shift_scale": scale,
                         "n_rct": n_rct, "n_obs": N_OBS,
                         "oracle_r_min": overlap["min"], "oracle_r_max": overlap["max"]}
                for sieve in sieves:
                    kind = "neural" if sieve == "neural" else "spline"
                    knots = 5 if sieve == "spline5" else 3
                    for route_name, builder in routes:
                        predictions: dict = {}
                        for rep in range(args.replications):
                            # one draw per replication, shared by every route, so
                            # the routes are compared on identical data
                            s = stable_seed("cate3", family, level, n_rct, rep)
                            rct = fd.sample_scm(scm, "RCT", n_rct, REGIME, rng_for(s, 1))
                            obs = fd.sample_scm(scm, "OBS", N_OBS, REGIME, rng_for(s, 2))
                            seed = stable_seed("cate3", family, level, n_rct, sieve, rep)
                            got = fc.cate_replication(
                                rct, obs, test, kind, seed, builder, route_name,
                                spline_knots=knots,
                                oracle_sampler=sampler if kind == "spline" else None,
                                oracle_fn=oracle_fn,
                                anchor=anchor_x, predictions=predictions)
                            for row in got:
                                row.update(label)
                                row["sieve_label"] = sieve
                                row["route"] = route_name
                                row["replication"] = rep
                                row["seed"] = seed
                            rows.extend(got)
                        for (learner, tag), stack in predictions.items():
                            matrix = np.vstack(stack)
                            k = len(matrix)
                            pred_var = float(np.mean(np.var(matrix, axis=0, ddof=1)))
                            sq_bias = float(np.mean((matrix.mean(axis=0) - anchor_tau) ** 2))
                            mean_risk = float(np.mean(np.mean((matrix - anchor_tau) ** 2, axis=1)))
                            decomposition.append({
                                **label, "sieve_label": sieve, "route": route_name,
                                "learner": learner, "rule": tag, "replications": k,
                                "prediction_variance": pred_var, "squared_bias": sq_bias,
                                "mean_risk_on_anchor": mean_risk,
                                "decomposition_gap":
                                    mean_risk - ((k - 1) / k * pred_var + sq_bias),
                                "anchor_points": matrix.shape[1]})
                        print(f"  {family} {level} n_R={n_rct} {sieve} {route_name}"
                              f"  [{time.time() - started:.0f}s]", flush=True)

    frame = pd.DataFrame(rows)
    group = ["family", "shift", "n_rct", "sieve_label", "route", "learner"]
    summary = summarize(frame, group,
                        {"risk_selected": "risk_selected", "risk_rct_only": "risk_rct_only",
                         "risk_grid_oracle": "risk_grid_oracle",
                         "risk_lambda_only": "risk_lambda_only",
                         "risk_omega_only": "risk_omega_only",
                         "selected_lambda": "selected_lambda",
                         "selected_omega": "selected_omega",
                         "risk_selected_minus_rct": "risk_selected_minus_rct",
                         "score_selected_minus_rct": "score_selected_minus_rct",
                         "bar_r": "bar_r"})
    paths = write_outputs(args.stem, frame, summary)
    decomp = pd.DataFrame(decomposition)
    decomp.to_csv(paths["replications"].parent / f"{args.stem}_decomposition.csv", index=False)
    worst = float(np.nanmax(np.abs(decomp["decomposition_gap"]))) if len(decomp) else float("nan")
    print(f"wrote {paths['replications'].name}; {len(frame)} rows in "
          f"{time.time() - started:.0f}s; largest risk decomposition gap {worst:.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
