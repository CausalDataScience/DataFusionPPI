"""ATE-2: covariate shift and the density ratio (handoff 5.2).

Checks the drift identity of Theorem 1, Definition 12 with Proposition 6, and
Theorem 5 with Corollary 5.1.  Two panels: shifted synthetic families with a
known ratio, and the NSW trial against the CPS and PSID comparison groups.

Shift levels were fixed by the overlap requirement of the handoff: scale 0.5
gives r0 in [0.31, 3.55] and scale 1.0 gives r0 in [0.105, 13.5]; every larger
scale exceeds the cap of 20 and is therefore inadmissible.
"""
from __future__ import annotations

import argparse
import functools
import sys
import time

import numpy as np
import pandas as pd

import fusion_data as fd
from fusion_ate import ate_replication
from fusion_core import (fit_outcome, ratio_balanced, ratio_classifier, ratio_oracle,
                         ratio_unit, ghat_on, rng_for, stable_seed)
from fusion_io import summarize, write_outputs

SHIFT_LEVELS = {"mild": 0.5, "strong": 1.0}


def _balanced_with_prediction(rct_nuis, obs_nuis, seed: int = 0, **kw):
    """The loss-targeted balancing map: the covariates together with the fitted
    observational prediction, which is what the manuscript asks for and what the
    earlier version silently dropped."""
    mu_o = fit_outcome(obs_nuis, spec="flexible")

    def phi(x):
        return np.column_stack([x, ghat_on(x, mu_o)])

    return ratio_balanced(rct_nuis, obs_nuis, seed=seed, phi=phi, name="balanced_x_ghat")


def ratio_specs():
    """Five routes.  The loss-targeted map appends the OBS prediction to phi."""
    return [("oracle", ratio_oracle),
            ("unit", ratio_unit),
            ("classifier", ratio_classifier),
            ("balanced_x", functools.partial(ratio_balanced, name="balanced_x")),
            ("balanced_x_ghat", _balanced_with_prediction)]


def drift_reference(scm, regime, ratio, mu_o, draws: int = 100000) -> float:
    """Theorem 1: omega * E_O[(r - r0) ghat], the exact mean drift per unit omega."""
    obs = fd.sample_scm(scm, "OBS", draws, regime, rng_for(4242, 7))
    r = ratio(obs["x"])
    return float(np.mean((r - obs["true_ratio"]) * ghat_on(obs["x"], mu_o)))


def synthetic_panel(replications: int, smoke: bool) -> pd.DataFrame:
    rows = []
    levels = list(SHIFT_LEVELS.items())[:1] if smoke else list(SHIFT_LEVELS.items())
    families = (1,) if smoke else (1, 2, 3)
    sizes = (100,) if smoke else (100, 400)
    for scm_id in families:
        for level_name, scale in levels:
            scm = fd.scm_family(scm_id, shift=scale)
            theta = fd.scm_true_theta(scm, "shifted")
            oracle_fn = fd.scm_oracle_ratio(scm, "shifted")
            overlap = fd.overlap_check(scm, "shifted")
            for n_rct in sizes:
                for name, builder in ratio_specs():
                    for rep in range(replications):
                        seed = stable_seed("ate2", scm_id, level_name, n_rct, rep)
                        rct = fd.sample_scm(scm, "RCT", n_rct, "shifted", rng_for(seed, 1))
                        obs = fd.sample_scm(scm, "OBS", 5000, "shifted", rng_for(seed, 2))
                        got = ate_replication(rct, obs, theta, builder, name, "flexible",
                                              seed, oracle_fn=oracle_fn, with_baselines=False)
                        for row in got:
                            row.update({"panel": "synthetic", "family": f"scm{scm_id}",
                                        "shift": level_name, "shift_scale": scale,
                                        "n_rct": n_rct, "n_obs": 5000, "replication": rep,
                                        "true_theta": theta, "seed": seed,
                                        "oracle_r_min": overlap["min"], "oracle_r_max": overlap["max"]})
                        rows.extend(got)
                print(f"  synthetic scm{scm_id} {level_name} done", flush=True)
    return pd.DataFrame(rows)


def nsw_panel(replications: int, smoke: bool) -> pd.DataFrame:
    rows = []
    comparisons = ("cps",) if smoke else ("cps", "psid")
    sizes = (200,) if smoke else (100, 200, None)
    for comparison in comparisons:
        design = fd.nsw_design(comparison)
        for n_rct in sizes:
            for name, builder in ratio_specs():
                if name == "oracle":
                    continue
                for rep in range(replications):
                    seed = stable_seed("nsw", comparison, n_rct or 0, rep)
                    rct, obs, reference = fd.nsw_replication(design, n_rct, rep)
                    got = ate_replication(rct, obs, reference, builder, name, "flexible",
                                          seed, with_baselines=False)
                    for row in got:
                        row.update({"panel": "nsw", "family": comparison, "shift": "real",
                                    "shift_scale": np.nan, "n_rct": n_rct or len(rct["x"]),
                                    "n_obs": len(obs["x"]), "replication": rep,
                                    "true_theta": reference, "seed": seed})
                    rows.extend(got)
            print(f"  nsw {comparison} n_R={n_rct} done", flush=True)
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replications", type=int, default=100)
    parser.add_argument("--stem", type=str, default="ate2_shift")
    parser.add_argument("--panel", type=str, default="both")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    started = time.time()
    frames = []
    if args.panel in ("both", "synthetic"):
        frames.append(synthetic_panel(args.replications, args.smoke))
    if args.panel in ("both", "nsw"):
        frames.append(nsw_panel(args.replications, args.smoke))
    frame = pd.concat(frames, ignore_index=True)
    group = ["panel", "family", "shift", "n_rct", "ratio", "estimator"]
    summary = summarize(frame, group,
                        {"rmse": "error", "bias": "error", "coverage": "covered",
                         "ci_length": "ci_length", "omega_used": "omega_used",
                         "ratio_balance_residual": "ratio_balance_residual",
                         "ratio_normalization_error": "ratio_normalization_error",
                         "ratio_ess": "ratio_ess"})
    paths = write_outputs(args.stem, frame, summary)
    print(f"wrote {paths['replications'].name}; {len(frame)} rows in {time.time() - started:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
