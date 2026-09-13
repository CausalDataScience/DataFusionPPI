"""CATE-1: honest validation on the sieve and the oracle calculus (handoff 5.3).

Checks Theorem 6, Theorem 7, and Corollary 7.1.  Three axes around a reference
configuration give 6 configurations per data family and 24 cells.
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd

import fusion_data as fd
import fusion_cate as fc
from fusion_core import ratio_unit, rng_for, stable_seed
from fusion_io import summarize, write_outputs

FAMILIES = ("scm1", "scm2", "scm3", "star", "bounded")
REFERENCE = {"n_rct": 200, "n_obs": 5000, "confounding": 1.0, "heterogeneity": 1.0}
TEST_DRAWS = 50000


def configurations() -> list[dict]:
    out = [{**REFERENCE, "n_rct": n, "axis": "rct_size"} for n in (100, 200, 400)]
    out += [{**REFERENCE, "confounding": c, "axis": "quality"} for c in (0.0, 2.0)]
    out += [{**REFERENCE, "heterogeneity": 2.0, "axis": "heterogeneity"}]
    return out


def build_cell(family: str, cfg: dict, star_support=None):
    """Returns (draw, sampler, test, score_range).

    ``score_range`` is the (bound, ratio cap) pair that Theorem 6 needs, and is
    supplied only by the bounded design, whose range holds before any data is
    drawn.  The Gaussian and STAR families return None, so their radius is
    reported as not applicable rather than read off the sample.
    """
    if family == "bounded":
        # the confounding axis scales the latent loading of the observational
        # assignment and the heterogeneity axis scales the effect; nothing else moves
        fam = fd.bounded_family(1, confounding=cfg["confounding"],
                                heterogeneity=cfg["heterogeneity"])
        test = fam.test_sample(TEST_DRAWS, rng_for(8900, 1))

        def draw(rep):
            s = stable_seed(family, tuple(sorted(cfg.items())), rep)
            return (fam.sample("RCT", cfg["n_rct"], rng_for(s, 1)),
                    fam.sample("OBS", cfg["n_obs"], rng_for(s, 2)))

        def sampler(n, seed):
            return (fam.sample("RCT", n, rng_for(seed, 11)), fam.sample("OBS", n, rng_for(seed, 12)))

        return draw, sampler, test, fam.score_range

    if family.startswith("scm"):
        scm = fd.scm_family(int(family[-1]), confounding=cfg["confounding"],
                            heterogeneity=cfg["heterogeneity"])
        test = fd.scm_test_sample(scm, "shared", TEST_DRAWS, rng_for(8800, int(family[-1])))

        def draw(rep):
            s = stable_seed(family, tuple(sorted(cfg.items())), rep)
            return (fd.sample_scm(scm, "RCT", cfg["n_rct"], "shared", rng_for(s, 1)),
                    fd.sample_scm(scm, "OBS", cfg["n_obs"], "shared", rng_for(s, 2)))

        def sampler(n, seed):
            return (fd.sample_scm(scm, "RCT", n, "shared", rng_for(seed, 11)),
                    fd.sample_scm(scm, "OBS", n, "shared", rng_for(seed, 12)))

        return draw, sampler, test, None

    support = star_support
    test = fd.star_realx_test(support, "shared", TEST_DRAWS, 0)

    def draw(rep):
        rct, obs, _ = fd.star_realx_sample(support, "shared", cfg["n_rct"], cfg["n_obs"], rep)
        return rct, obs

    def sampler(n, seed):
        rct, obs, _ = fd.star_realx_sample(support, "shared", n, n, int(seed) % 100000)
        return rct, obs

    return draw, sampler, test, None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replications", type=int, default=100)
    parser.add_argument("--families", type=str, default=",".join(FAMILIES))
    parser.add_argument("--sieves", type=str, default="spline3,spline5,neural")
    parser.add_argument("--stem", type=str, default="cate1_sieve_validation")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()

    families = tuple(args.families.split(","))
    sieves = tuple(args.sieves.split(","))
    configs = configurations()[:1] if args.smoke else configurations()
    star_support = fd.star_realx_support() if "star" in families else None

    rows, decomposition, started = [], [], time.time()
    anchor_size = 500
    for family in families:
        for cfg in configs:
            draw, sampler, test, score_range = build_cell(family, cfg, star_support)
            bounded = family == "bounded"
            label = {"family": family, "axis": cfg["axis"], "n_rct": cfg["n_rct"],
                     "n_obs": cfg["n_obs"], "confounding": cfg["confounding"],
                     "heterogeneity": cfg["heterogeneity"]}
            anchor_x = test["x"][:anchor_size]
            anchor_tau = test["tau"][:anchor_size]
            for sieve in sieves:
                kind = "neural" if sieve == "neural" else "spline"
                knots = 5 if sieve == "spline5" else 3
                predictions: dict = {}
                for rep in range(args.replications):
                    rct, obs = draw(rep)
                    seed = stable_seed(family, tuple(sorted(cfg.items())), sieve, rep)
                    got = fc.cate_replication(rct, obs, test, kind, seed, ratio_unit, "unit",
                                              spline_knots=knots,
                                              oracle_sampler=sampler if kind == "spline" else None,
                                              anchor=anchor_x, predictions=predictions,
                                              score_range=score_range,
                                              candidate_bound=fd.BOUNDED_CANDIDATE if bounded else None,
                                              outcome_clip=fd.BOUNDED_PREDICTION if bounded else None)
                    for row in got:
                        row.update(label)
                        row["sieve_label"] = sieve
                        row["replication"] = rep
                        row["seed"] = seed
                    rows.extend(got)
                for (learner, tag), stack in predictions.items():
                    matrix = np.vstack(stack)
                    k = len(matrix)
                    pred_var = float(np.mean(np.var(matrix, axis=0, ddof=1)))
                    mean_pred = matrix.mean(axis=0)
                    sq_bias = float(np.mean((mean_pred - anchor_tau) ** 2))
                    mean_risk = float(np.mean(np.mean((matrix - anchor_tau) ** 2, axis=1)))
                    decomposition.append({**label, "sieve_label": sieve, "learner": learner,
                                          "rule": tag, "replications": k,
                                          "prediction_variance": pred_var,
                                          "squared_bias": sq_bias,
                                          "mean_risk_on_anchor": mean_risk,
                                          "decomposition_gap":
                                              mean_risk - ((k - 1) / k * pred_var + sq_bias),
                                          "anchor_points": matrix.shape[1]})
                print(f"  {family} {cfg['axis']} n_R={cfg['n_rct']} conf={cfg['confounding']} "
                      f"het={cfg['heterogeneity']} {sieve}  [{time.time() - started:.0f}s]", flush=True)

    frame = pd.DataFrame(rows)
    group = ["family", "axis", "n_rct", "confounding", "heterogeneity", "sieve_label", "learner"]
    summary = summarize(frame, group,
                        {"risk_selected": "risk_selected", "risk_rct_only": "risk_rct_only",
                         "risk_grid_oracle": "risk_grid_oracle",
                         "risk_lambda_only": "risk_lambda_only",
                         "risk_omega_only": "risk_omega_only",
                         "regret": "regret", "radius": "radius",
                         "regret_within_radius": "regret_within_radius",
                         "regret_over_radius": "regret_over_radius",
                         "selected_lambda": "selected_lambda", "selected_omega": "selected_omega",
                         "oracle_sieve_lambda": "oracle_sieve_lambda",
                         "oracle_sieve_omega": "oracle_sieve_omega",
                         "risk_selected_minus_rct": "risk_selected_minus_rct",
                         "score_selected_minus_rct": "score_selected_minus_rct",
                         "a_p2": "a_p2",
                         "risk_at_0.0_0.0": "risk_at_0.0_0.0", "risk_at_0.5_0.0": "risk_at_0.5_0.0",
                         "risk_at_0.0_0.5": "risk_at_0.0_0.5", "risk_at_0.5_0.5": "risk_at_0.5_0.5"})
    paths = write_outputs(args.stem, frame, summary)
    decomp = pd.DataFrame(decomposition)
    decomp.to_csv(paths["replications"].parent / f"{args.stem}_decomposition.csv", index=False)
    worst = float(np.nanmax(np.abs(decomp["decomposition_gap"]))) if len(decomp) else float("nan")
    print(f"wrote {paths['replications'].name}; {len(frame)} rows in {time.time() - started:.0f}s; "
          f"largest risk decomposition gap {worst:.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
