"""Variance-first analysis of the ATE studies.

The omega channel is defined to reduce variance, so it is judged on variance
directly.  Every comparison pairs two estimators that share the same lambda
inside the same replication, and every interval comes from a paired bootstrap
that resamples replication indices once and reuses them for both arms.
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

from fusion_core import variance_standard_error
from fusion_io import MATERIALS

BOOTSTRAP = 2000

# (name, omega-off arm, omega-on arm, what the pair isolates)
PAIRS = [
    ("omega_alone", "rct_only", "omega_only", "omega with lambda held at zero"),
    ("omega_given_lambda", "lambda_only", "joint", "omega on top of the selected lambda"),
    ("whole_method", "rct_only", "joint", "both channels against the trial-only estimator"),
    ("omega_alone_oracle", "rct_only", "oracle_omega_only", "oracle omega with lambda at zero"),
    ("omega_given_lambda_oracle", "oracle_lambda_only", "oracle_joint",
     "oracle omega on top of the oracle lambda"),
]


def _variance(values: np.ndarray) -> float:
    return float(np.var(values, ddof=1)) if len(values) > 1 else np.nan


def paired_variance_comparison(off: np.ndarray, on: np.ndarray, rng: np.random.Generator,
                               bootstrap: int = BOOTSTRAP) -> dict:
    """Variance reduction and its rate, with a paired bootstrap interval."""
    k = len(off)
    v_off, v_on = _variance(off), _variance(on)
    reduction = v_off - v_on
    rate = 1.0 - v_on / v_off if v_off > 0 else np.nan
    idx = rng.integers(0, k, size=(bootstrap, k))
    boot_off = np.var(off[idx], axis=1, ddof=1)
    boot_on = np.var(on[idx], axis=1, ddof=1)
    boot_red = boot_off - boot_on
    with np.errstate(divide="ignore", invalid="ignore"):
        boot_rate = np.where(boot_off > 0, 1.0 - boot_on / boot_off, np.nan)
    finite = boot_rate[np.isfinite(boot_rate)]
    return {"variance_off": v_off, "variance_on": v_on,
            "variance_reduction": reduction,
            "reduction_low": float(np.percentile(boot_red, 2.5)),
            "reduction_high": float(np.percentile(boot_red, 97.5)),
            "reduction_rate": rate,
            "rate_low": float(np.percentile(finite, 2.5)) if len(finite) else np.nan,
            "rate_high": float(np.percentile(finite, 97.5)) if len(finite) else np.nan,
            "rate_excludes_zero": float(len(finite) > 0 and np.percentile(finite, 2.5) > 0),
            "replications": k}


def check_uniqueness(frame: pd.DataFrame, group: list[str]) -> None:
    """One row per cell, replication and estimator.  Duplicates would be hidden by
    the mean inside a pivot, so they are an error rather than something to average."""
    keys = group + ["replication", "estimator"]
    counts = frame.groupby(keys, dropna=False).size()
    duplicated = int((counts > 1).sum())
    if duplicated:
        raise AssertionError(f"{duplicated} duplicated (cell, replication, estimator) rows")



def pooled_variance_comparison(cells: list[tuple[np.ndarray, np.ndarray]], weights: np.ndarray,
                               rng: np.random.Generator, bootstrap: int = BOOTSTRAP) -> dict:
    """A fixed-weight average over cells, with the whole aggregate recomputed in
    every resample.

    Cells are a fixed design, so they are not resampled.  Replications are
    resampled inside each cell with the two arms kept paired, and the weighted
    average rate is recomputed from the resampled cells.  Averaging the endpoints
    of per-cell intervals would not give an interval for the average.
    """
    weights = np.asarray(weights, dtype=float)
    weights = weights / weights.sum()
    rates = np.array([1.0 - _variance(on) / _variance(off) if _variance(off) > 0 else np.nan
                      for off, on in cells])
    point = float(np.nansum(weights * rates))
    draws = np.empty(bootstrap)
    for b in range(bootstrap):
        vals = np.empty(len(cells))
        for i, (off, on) in enumerate(cells):
            idx = rng.integers(0, len(off), size=len(off))
            v_off, v_on = _variance(off[idx]), _variance(on[idx])
            vals[i] = 1.0 - v_on / v_off if v_off > 0 else np.nan
        draws[b] = np.nansum(weights * vals)
    return {"reduction_rate": point,
            "rate_low": float(np.percentile(draws, 2.5)),
            "rate_high": float(np.percentile(draws, 97.5)),
            "rate_excludes_zero": float(np.percentile(draws, 2.5) > 0),
            "cells": len(cells)}



def pooled_mean(frame: pd.DataFrame, group: list[str], column: str, estimator: str,
                rng: np.random.Generator, bootstrap: int = BOOTSTRAP) -> dict:
    """Fixed-weight average of a per-replication quantity over design cells.

    Cells carry equal declared weight and are not resampled; replications are
    resampled inside each cell and the weighted average is recomputed each time.
    """
    cells = []
    for _, block in frame[frame.estimator == estimator].groupby(group, dropna=False):
        values = block[column].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if len(values) >= 2:
            cells.append(values)
    if not cells:
        return {"value": np.nan, "low": np.nan, "high": np.nan, "cells": 0}
    weights = np.full(len(cells), 1.0 / len(cells))
    point = float(np.sum(weights * np.array([v.mean() for v in cells])))
    draws = np.empty(bootstrap)
    for b in range(bootstrap):
        draws[b] = float(np.sum(weights * np.array(
            [v[rng.integers(0, len(v), size=len(v))].mean() for v in cells])))
    return {"value": point, "low": float(np.percentile(draws, 2.5)),
            "high": float(np.percentile(draws, 97.5)), "cells": len(cells)}


def pooled_reporting(frame: pd.DataFrame, group: list[str], columns: list[str],
                     seed: int = 20260912) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for estimator in sorted(frame["estimator"].dropna().unique()):
        for column in columns:
            if column not in frame:
                continue
            rows.append({"estimator": estimator, "quantity": column,
                         **pooled_mean(frame, group, column, estimator, rng)})
    return pd.DataFrame(rows)


def reporting_table(frame: pd.DataFrame, group: list[str]) -> pd.DataFrame:
    """The seven required quantities per cell and estimator, with the identity check.

    One mask decides which replications count, and every column is read from the
    same masked rows, so variance, bias and error are aggregated over exactly the
    same replications."""
    rows = []
    for keys, block in frame.groupby(group + ["estimator"], dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        base = dict(zip(group + ["estimator"], keys))
        failed = block["failed"].fillna(0).to_numpy(dtype=float) if "failed" in block else np.zeros(len(block))
        valid = (failed == 0) & np.isfinite(block["error"].to_numpy(dtype=float)) \
            & np.isfinite(block["estimate"].to_numpy(dtype=float))
        attempted = len(block)
        ok = block[valid]
        err = ok["error"].to_numpy(dtype=float)
        k = len(err)
        if k < 2:
            continue
        est = ok["estimate"].to_numpy(dtype=float)
        s2 = _variance(est)
        bias = float(np.mean(err))
        mse = float(np.mean(err ** 2))
        vhat = ok["variance_hat"].to_numpy(dtype=float)
        vhat = vhat[np.isfinite(vhat)]
        cov = ok["covered"].to_numpy(dtype=float)
        cov = cov[np.isfinite(cov)]
        length = ok["ci_length"].to_numpy(dtype=float)
        length = length[np.isfinite(length)]
        identity = (k - 1) / k * s2 + bias ** 2
        rows.append({**base,
                     "mean_estimate": float(np.mean(est)),
                     "bias": bias, "bias_se": float(np.std(err, ddof=1) / np.sqrt(k)),
                     "empirical_variance": s2,
                     "empirical_variance_se": variance_standard_error(est),
                     "empirical_sd": float(np.sqrt(s2)),
                     "mean_estimated_variance": float(np.mean(vhat)) if len(vhat) else np.nan,
                     "variance_ratio_hat_to_empirical": float(np.mean(vhat)) / s2 if s2 > 0 and len(vhat) else np.nan,
                     "mse": mse, "mse_se": float(np.std(err ** 2, ddof=1) / np.sqrt(k)),
                     "rmse": float(np.sqrt(mse)),
                     "coverage": float(np.mean(cov)) if len(cov) else np.nan,
                     "coverage_se": float(np.sqrt(np.mean(cov) * (1 - np.mean(cov)) / len(cov))) if len(cov) else np.nan,
                     "mean_ci_length": float(np.mean(length)) if len(length) else np.nan,
                     "effective_replications": k,
                     "attempted_replications": attempted,
                     "failure_rate": float(1.0 - k / attempted),
                     "mse_identity_gap": mse - identity})
    return pd.DataFrame(rows)


def variance_table(frame: pd.DataFrame, group: list[str], seed: int = 20260911) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for keys, block in frame.groupby(group, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        base = dict(zip(group, keys))
        wide = block.pivot_table(index="replication", columns="estimator", values="estimate")
        for name, off_name, on_name, purpose in PAIRS:
            if off_name not in wide or on_name not in wide:
                continue
            sub = wide[[off_name, on_name]].dropna()
            if len(sub) < 10:
                continue
            result = paired_variance_comparison(sub[off_name].to_numpy(),
                                                sub[on_name].to_numpy(), rng)
            rows.append({**base, "comparison": name, "off": off_name, "on": on_name,
                         "purpose": purpose, **result})
    table = pd.DataFrame(rows)
    if table.empty:
        return table
    table["scope"] = "cell"

    pooled_rows = []
    # a pooled row per family as well as overall, so a figure can show a real
    # interval for the average instead of an envelope of per-cell intervals
    partitions = [("pooled", frame)]
    if "family" in frame.columns:
        partitions += [(f"pooled:{fam}", block) for fam, block in frame.groupby("family")]
    for scope_name, part in partitions:
        pooled_rows.extend(_pooled_for(part, group, rng, scope_name))
    table = pd.concat([table, pd.DataFrame(pooled_rows)], ignore_index=True)
    return table


def _pooled_for(frame: pd.DataFrame, group: list[str], rng, scope_name: str) -> list[dict]:
    rows = []
    for name, off_name, on_name, purpose in PAIRS:
        cells, weights = [], []
        for _, block in frame.groupby(group, dropna=False):
            wide = block.pivot_table(index="replication", columns="estimator", values="estimate")
            if off_name not in wide or on_name not in wide:
                continue
            sub = wide[[off_name, on_name]].dropna()
            if len(sub) < 10:
                continue
            cells.append((sub[off_name].to_numpy(), sub[on_name].to_numpy()))
            weights.append(1.0)          # equal weight per design cell, declared here
        if cells:
            rows.append({"scope": scope_name, "comparison": name, "off": off_name,
                         "on": on_name, "purpose": purpose,
                         **pooled_variance_comparison(cells, np.array(weights), rng)})
    return rows
    # recovery of the oracle reduction, reported only where the oracle gain is real
    keyed = table.set_index(group + ["comparison"])
    recovery = []
    for keys, _ in table.groupby(group, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        base = dict(zip(group, keys))
        for est_name, orc_name in (("omega_alone", "omega_alone_oracle"),
                                   ("omega_given_lambda", "omega_given_lambda_oracle")):
            try:
                est = keyed.loc[tuple(keys) + (est_name,)]
                orc = keyed.loc[tuple(keys) + (orc_name,)]
            except KeyError:
                continue
            orc_rate = float(orc["reduction_rate"])
            usable = float(orc["rate_low"]) > 0.005
            recovery.append({**base, "comparison": f"{est_name}_recovery",
                             "reduction_rate": float(est["reduction_rate"]) / orc_rate if usable else np.nan,
                             "oracle_reduction_rate": orc_rate,
                             "recovery_reportable": float(usable)})
    return pd.concat([table, pd.DataFrame(recovery)], ignore_index=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stem", type=str, default="ate1_when_fusion_helps")
    parser.add_argument("--group", type=str,
                        default="family,axis,n_rct,n_obs,confounding,spec")
    args = parser.parse_args()
    group = [g for g in args.group.split(",") if g]
    frame = pd.read_csv(MATERIALS / f"{args.stem}_replications.csv")
    group = [g for g in group if g in frame.columns]
    check_uniqueness(frame, group)
    report = reporting_table(frame, group)
    variance = variance_table(frame, group)
    report.to_csv(MATERIALS / f"{args.stem}_reporting.csv", index=False)
    pooled_reporting(frame, group, ["covered", "ci_length", "variance_hat", "error"]).to_csv(
        MATERIALS / f"{args.stem}_pooled.csv", index=False)
    variance.to_csv(MATERIALS / f"{args.stem}_variance.csv", index=False)
    worst = float(np.nanmax(np.abs(report["mse_identity_gap"]))) if len(report) else np.nan
    print(f"{args.stem}: {len(report)} reporting rows, {len(variance)} variance rows; "
          f"largest MSE identity gap {worst:.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
