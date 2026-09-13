"""Aggregation, paired bootstrap, two figures and one table.

Everything is computed from the sufficient statistics the run stored, which the
plan defines as

    G[k, l] = P_T[f_k f_l],   v[k] = P_T[f_k tau],   t = P_T[tau^2].

For resampling weights a summing to one,

    mean risk = sum_k a_k (G[k,k] - 2 v[k] + t)
    squared bias = a' G a - 2 a' v + t
    variance component = sum_k a_k G[k,k] - a' G a

and the third is the first minus the second, so the plan's identity

    mean risk = squared bias + (K - 1) / K * variance

holds at every resample rather than only at the point estimate.  No model is
refitted inside the bootstrap and the fifty thousand truth points are never
revisited.

The bootstrap resamples whole seed bundles.  One draw of replication indices is
applied to every cell, learner, route and rule at once, which is what makes a
difference between two of them paired.

Run: python3 -m confirm.mm_analyze --stage <path>
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from confirm import mm_run as mr
from confirm import provenance as pv

BOOTSTRAP = 2000
BOOTSTRAP_SEED = 91307
PANELS = (("shared", "unit_exact", "shared, exact r0 = 1"),
          ("shifted", "exact", "shifted, true ratio"),
          ("shifted", "classifier", "shifted, estimated ratio"),
          ("shifted", "wrong_unit", "shifted, ratio ignored"))
RULE_LABEL = {"rct_only": "RCT only", "lambda_only": r"$\lambda$ only",
              "omega_only": r"$\omega$ only", "joint": "joint",
              "grid_oracle": "grid oracle"}


def pieces(block: dict, t: float, weights: np.ndarray) -> dict:
    gram = np.asarray(block["G"])
    v = np.asarray(block["v"])
    k = block["K"]
    diagonal = np.diag(gram)
    mean_risk = float(weights @ (diagonal - 2 * v + t))
    bias2 = float(weights @ gram @ weights - 2 * weights @ v + t)
    component = mean_risk - bias2
    return {"mean_risk": mean_risk, "squared_bias": bias2,
            "variance_component": component,
            "variance": component * k / max(k - 1, 1), "K": k}


def load(stage: Path) -> tuple[dict, dict]:
    stats = json.loads((stage / "sufficient_statistics.json").read_text())
    config = json.loads((stage / "config.json").read_text())
    return stats, config


def analyse(stats: dict, rules: tuple[str, ...]) -> dict:
    t = stats["t"]
    groups = stats["groups"]
    any_block = next(iter(groups.values()))
    k = any_block["K"]
    uniform = np.full(k, 1.0 / k)

    point = {key: pieces(block, t, uniform) for key, block in groups.items()}

    rng = np.random.default_rng(BOOTSTRAP_SEED)
    draws = {key: [] for key in groups}
    ratios = {key: [] for key in groups}
    differences = {key: [] for key in groups}
    for _ in range(BOOTSTRAP):
        index = rng.integers(0, k, size=k)
        weights = np.bincount(index, minlength=k).astype(float) / k
        resampled = {key: pieces(block, t, weights) for key, block in groups.items()}
        for key, value in resampled.items():
            cell, learner, route, rule = key.split("|")
            base = resampled[mr.group_key(cell, learner, route, "rct_only")]
            draws[key].append(value["mean_risk"])
            ratios[key].append(value["mean_risk"] / base["mean_risk"]
                               if base["mean_risk"] > 0 else np.nan)
            differences[key].append(value["mean_risk"] - base["mean_risk"])

    out = {}
    for key, value in point.items():
        cell, learner, route, rule = key.split("|")
        base = point[mr.group_key(cell, learner, route, "rct_only")]
        out[key] = {
            **value, "cell": cell, "learner": learner, "route": route, "rule": rule,
            "risk_ratio": value["mean_risk"] / base["mean_risk"],
            "risk_difference": value["mean_risk"] - base["mean_risk"],
            "variance_ratio": (value["variance"] / base["variance"]
                               if base["variance"] > 0 else float("nan")),
            "ratio_low": float(np.quantile(ratios[key], 0.025)),
            "ratio_high": float(np.quantile(ratios[key], 0.975)),
            "difference_low": float(np.quantile(differences[key], 0.025)),
            "difference_high": float(np.quantile(differences[key], 0.975)),
            "difference_mcse": float(np.std(differences[key], ddof=1)),
        }

    # the omega increment, joint against lambda only, on the same resamples
    for cell, route, _ in PANELS:
        for learner in mr.LEARNERS:
            joint = mr.group_key(cell, learner, route, "joint")
            lam = mr.group_key(cell, learner, route, "lambda_only")
            if joint not in out:
                continue
            gaps = np.array(draws[joint]) - np.array(draws[lam])
            out[joint]["omega_increment"] = (out[joint]["mean_risk"]
                                             - out[lam]["mean_risk"])
            out[joint]["omega_increment_low"] = float(np.quantile(gaps, 0.025))
            out[joint]["omega_increment_high"] = float(np.quantile(gaps, 0.975))
    return out


def write_table(result: dict, stage: Path, rules: tuple[str, ...]) -> None:
    import csv
    fields = ["cell", "learner", "route", "rule", "K", "mean_risk", "risk_ratio",
              "ratio_low", "ratio_high", "risk_difference", "difference_low",
              "difference_high", "difference_mcse", "variance", "variance_ratio",
              "squared_bias"]
    path = stage / "table1_results.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for cell, route, _ in PANELS:
            for learner in mr.LEARNERS:
                for rule in (*rules, "grid_oracle"):
                    key = mr.group_key(cell, learner, route, rule)
                    if key in result:
                        writer.writerow(result[key])
    print(f"  wrote {path.name}")


def write_figures(result: dict, stage: Path, rules: tuple[str, ...]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"rct_only": "#11233F", "lambda_only": "#0E7490",
              "omega_only": "#53657A", "joint": "#173B7A"}
    order = [r for r in rules]

    fig, axes = plt.subplots(2, 4, figsize=(11.0, 5.0), sharey="row")
    for row, learner in enumerate(mr.LEARNERS):
        for col, (cell, route, title) in enumerate(PANELS):
            ax = axes[row][col]
            xs = np.arange(len(order))
            for i, rule in enumerate(order):
                key = mr.group_key(cell, learner, route, rule)
                entry = result[key]
                ax.bar(i, entry["risk_ratio"], 0.62, color=colors[rule])
                ax.plot([i, i], [entry["ratio_low"], entry["ratio_high"]],
                        color="#31415A", lw=1.2)
            ax.axhline(1.0, color="#999", lw=0.8, ls=":")
            ax.set_xticks(xs)
            ax.set_xticklabels([RULE_LABEL[r] for r in order], fontsize=6.5,
                               rotation=25, ha="right")
            if row == 0:
                ax.set_title(title, fontsize=8.5)
            if col == 0:
                ax.set_ylabel(f"{learner}\nrisk ratio to RCT only", fontsize=8)
            ax.tick_params(labelsize=7)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(stage / "figure1_risk_comparison.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(2, 4, figsize=(11.0, 5.0), sharey="row")
    for row, learner in enumerate(mr.LEARNERS):
        for col, (cell, route, title) in enumerate(PANELS):
            ax = axes[row][col]
            base = result[mr.group_key(cell, learner, route, "rct_only")]["mean_risk"]
            for i, rule in enumerate(order):
                entry = result[mr.group_key(cell, learner, route, rule)]
                bias = entry["squared_bias"] / base
                var = entry["variance_component"] / base
                ax.bar(i, bias, 0.62, color="#9A3412", label="squared bias" if i == 0 else "")
                ax.bar(i, var, 0.62, bottom=bias, color="#173B7A",
                       label="variance" if i == 0 else "")
            ax.axhline(1.0, color="#999", lw=0.8, ls=":")
            ax.set_xticks(np.arange(len(order)))
            ax.set_xticklabels([RULE_LABEL[r] for r in order], fontsize=6.5,
                               rotation=25, ha="right")
            if row == 0:
                ax.set_title(title, fontsize=8.5)
            if col == 0:
                ax.set_ylabel(f"{learner}\nrisk share of RCT only", fontsize=8)
            if row == 0 and col == 3:
                ax.legend(fontsize=6.5, frameon=False, loc="upper right")
            ax.tick_params(labelsize=7)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(stage / "figure2_risk_decomposition.png", dpi=200)
    plt.close(fig)
    print("  wrote figure1_risk_comparison.png, figure2_risk_decomposition.png")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True)
    args = parser.parse_args(argv)
    stage = Path(args.stage)
    if not stage.is_absolute():
        stage = pv.CONFIRM / stage
    stats, config = load(stage)
    rules = tuple(config["rules"])
    result = analyse(stats, rules)
    pv.write_json(stage / "analysis.json", result)
    write_table(result, stage, rules)
    write_figures(result, stage, rules)

    print(f"\n{'cell':8s}{'learner':8s}{'route':12s}{'rule':12s}"
          f"{'risk':>10s}{'ratio':>8s}{'95% interval':>18s}{'var ratio':>11s}")
    for cell, route, _ in PANELS:
        for learner in mr.LEARNERS:
            for rule in rules:
                e = result[mr.group_key(cell, learner, route, rule)]
                print(f"{cell:8s}{learner:8s}{route:12s}{rule:12s}"
                      f"{e['mean_risk']:10.4f}{e['risk_ratio']:8.3f}"
                      f"   [{e['ratio_low']:.3f}, {e['ratio_high']:.3f}]"
                      f"{e['variance_ratio']:11.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
