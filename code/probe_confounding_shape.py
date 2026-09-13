"""Diagnostic: is the observational confounding bias a constant shift in x?

The source-indicator pooled learner absorbs any part of the observational bias
that does not vary with the covariates, because the indicator is a constant.
This script measures how much of the bias varies.  It writes one row per family,
confounding level and covariate quintile.

Run: python3 probe_confounding_shape.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

import fusion_data as fd
from fusion_io import MATERIALS

N_DRAW = 200_000
QUINTILES = 5


def one_family(fid: int, confounding: float, seed: int) -> list[dict]:
    scm = fd.scm_family(fid, confounding=confounding)
    obs = fd.sample_scm(scm, "OBS", N_DRAW, "shared", np.random.default_rng(seed))
    tau = fd.scm_tau(obs["x"], scm)

    # the propensity the pipeline can actually fit sees x only, never the confounder
    model = LogisticRegression(max_iter=2000).fit(obs["x"], obs["a"])
    e = np.clip(model.predict_proba(obs["x"])[:, 1], 0.02, 0.98)
    contrast = obs["a"] * obs["y"] / e - (1 - obs["a"]) * obs["y"] / (1 - e)
    bias = contrast - tau

    z = obs["x"][:, 0]
    edges = np.quantile(z, np.linspace(0.0, 1.0, QUINTILES + 1))
    rows = []
    for i in range(QUINTILES):
        lo, hi = edges[i], edges[i + 1]
        keep = (z >= lo) & (z <= hi if i == QUINTILES - 1 else z < hi)
        rows.append({"family": f"scm{fid}", "confounding": confounding,
                     "quintile": i + 1,
                     "bias_in_quintile": float(bias[keep].mean()),
                     "mean_bias": float(bias.mean()),
                     "tau_sd": float(tau.std()),
                     "draws": int(keep.sum())})
    return rows


def main() -> None:
    rows = []
    for fid in (1, 2, 3):
        for conf in (0.0, 1.0, 2.0):
            rows.extend(one_family(fid, conf, seed=1 + fid))
    frame = pd.DataFrame(rows)
    grouped = frame.groupby(["family", "confounding"])["bias_in_quintile"]
    frame = frame.merge(
        (grouped.max() - grouped.min()).rename("bias_spread").reset_index(),
        on=["family", "confounding"])
    frame["spread_over_tau_sd"] = frame["bias_spread"] / frame["tau_sd"]
    path = MATERIALS / "confounding_shape.csv"
    frame.to_csv(path, index=False)
    summary = frame.drop_duplicates(["family", "confounding"])
    print(f"wrote {path.name}; {len(frame)} rows from {N_DRAW} draws per cell")
    for _, r in summary.iterrows():
        print(f"  {r['family']} confounding={r['confounding']:.0f}: "
              f"mean bias {r['mean_bias']:+.3f}, spread across quintiles "
              f"{r['bias_spread']:.3f}, tau sd {r['tau_sd']:.3f}, "
              f"ratio {r['spread_over_tau_sd']:.3f}")


if __name__ == "__main__":
    main()
