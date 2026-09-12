"""ATE-only STAR real-outcome stress test.

This is the ATE projection of ``exp_cate2.py``.  It preserves the same cohort,
source-law sampler, seed, estimand, estimator call, labels, and output schema,
without fitting or writing any CATE object.
"""
from __future__ import annotations

import argparse
import sys
import time

import pandas as pd

import fusion_data as fd
from fusion_ate import ate_replication
from fusion_core import ratio_unit, stable_seed
from fusion_io import MATERIALS, summarize


OUTCOMES = ("mathk", "readk")
SIZES = (200, 400)
ALPHAS = (0.4, 0.8)


def one_replication(outcome: str, n_rct: int, alpha: float, rep: int,
                    cohort=None) -> list[dict]:
    """Return the exact ATE rows used by the mixed legacy driver."""
    cohort = cohort or fd.star_cohort(outcome, "pooled")
    _, obs, extra = fd.star_cate2_replication(
        cohort, n_rct, alpha, rep, n_rct_ate=n_rct + 2000)
    seed = stable_seed("cate2", outcome, n_rct, alpha, rep)
    rows = ate_replication(extra["ate_rct"], obs, cohort.ate_reference,
                           ratio_unit, "unit", "linear", seed,
                           with_baselines=True)
    label = {"outcome": outcome, "n_rct": n_rct, "alpha": alpha,
             "replication": rep, **extra["diagnostics"]}
    for row in rows:
        row.update(label)
        row["n_trial_for_ate"] = len(extra["ate_rct"]["x"])
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replications", type=int, default=100)
    parser.add_argument("--stem", default="cate2_star_real_ate")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    outcomes = OUTCOMES[:1] if args.smoke else OUTCOMES
    sizes = SIZES[:1] if args.smoke else SIZES
    alphas = ALPHAS[:1] if args.smoke else ALPHAS
    reps = min(args.replications, 2) if args.smoke else args.replications
    rows, started = [], time.time()
    for outcome in outcomes:
        cohort = fd.star_cohort(outcome, "pooled")
        for n_rct in sizes:
            for alpha in alphas:
                for rep in range(reps):
                    rows.extend(one_replication(outcome, n_rct, alpha, rep, cohort))
    frame = pd.DataFrame(rows)
    summary = summarize(frame, ["outcome", "n_rct", "alpha", "estimator"],
                        {"rmse": "error", "bias": "error", "coverage": "covered",
                         "estimate": "estimate", "ci_length": "ci_length"})
    rep_path = MATERIALS / f"{args.stem}_replications.csv"
    sum_path = MATERIALS / f"{args.stem}_summary.csv"
    frame.to_csv(rep_path, index=False)
    summary.to_csv(sum_path, index=False)
    print(f"wrote {rep_path.name}; {len(frame)} rows in {time.time()-started:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
