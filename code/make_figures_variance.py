"""The four variance-first figures.

Figure 1 judges omega on variance reduction directly, figure 2 splits accuracy
into variance and squared bias, figure 3 checks that the variance estimate and
the interval keep up, and figure 4 does the same split for the conditional
effect.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MATERIALS = Path(__file__).resolve().parents[1] / "materials"
NAVY, TEAL, RUST, MUTED = "#173B7A", "#0E7490", "#9A3412", "#53657A"


def _clean(ax):
    ax.tick_params(labelsize=7.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def figure_variance_reduction() -> None:
    """Figure 1: the reduction rate against the omega-off variance, with paired
    bootstrap intervals, estimated against oracle, conditional against full."""
    full = pd.read_csv(MATERIALS / "ate1_when_fusion_helps_variance.csv")
    cond = pd.read_csv(MATERIALS / "conditional_variance_400_variance.csv")
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.0), sharex=True)

    panels = [
        (axes[0], "conditional on the fitted objects", cond, "family"),
        (axes[1], "whole procedure refitted", full, "family"),
    ]
    for ax, title, frame, key in panels:
        rows = []
        for comparison, color, label in (("omega_given_lambda", NAVY, "estimated omega"),
                                         ("omega_given_lambda_oracle", RUST, "oracle omega")):
            block = frame[frame.comparison == comparison]
            if block.empty:
                continue
            # one pooled row per family, so each bar is a real interval for that
            # family's average rather than an envelope of its cells
            per_family = block[block["scope"].astype(str).str.startswith("pooled:")].copy()
            per_family["family"] = per_family["scope"].str.split(":").str[1]
            agg = per_family.set_index("family")[["reduction_rate", "rate_low", "rate_high"]]
            agg.columns = ["rate", "low", "high"]
            pooled = block[block["scope"] == "pooled"]
            rows.append((agg.sort_index(), color, label, pooled))
        if not rows:
            continue
        labels = list(rows[0][0].index) + ["pooled"]
        offsets = np.linspace(-0.16, 0.16, len(rows))
        for (agg, color, label, pooled), off in zip(rows, offsets):
            rate = list(agg["rate"].to_numpy())
            low = list(agg["low"].to_numpy())
            high = list(agg["high"].to_numpy())
            if len(pooled):
                rate.append(float(pooled["reduction_rate"].iloc[0]))
                low.append(float(pooled["rate_low"].iloc[0]))
                high.append(float(pooled["rate_high"].iloc[0]))
            else:
                rate.append(np.nan); low.append(np.nan); high.append(np.nan)
            rate, low, high = np.array(rate), np.array(low), np.array(high)
            y = np.arange(len(labels)) + off
            ax.errorbar(100 * rate, y,
                        xerr=[100 * (rate - low), 100 * (high - rate)],
                        fmt="o", ms=3.5, lw=1.1, capsize=2, color=color, label=label)
        ax.axvline(0.0, color=MUTED, lw=0.8, ls=":")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=7.5)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("variance reduction against omega off, percent\n(each bar is a bootstrap interval for that average)", fontsize=7.5)
        _clean(ax)
    handles, labels_ = axes[0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.legend(handles, labels_, fontsize=7.5, frameon=False, ncol=2,
               loc="upper center", bbox_to_anchor=(0.5, 0.995))
    fig.savefig(MATERIALS / "figure1_omega_variance_reduction.png", dpi=200)
    plt.close(fig)


def figure_decomposition() -> None:
    """Figure 2: empirical variance, squared bias, and RMSE for the same cells."""
    rep = pd.read_csv(MATERIALS / "ate1_when_fusion_helps_reporting.csv")
    block = rep[(rep.axis == "rct_size")]
    order = ["rct_only", "omega_only", "lambda_only", "joint", "shrinkage", "obs_transported"]
    names = ["trial only", "omega only", "lambda only", "joint", "shrinkage", "OBS AIPW"]
    sizes = sorted(block.n_rct.unique())
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.5))
    width = 0.8 / len(sizes)
    for ax, column, title in zip(axes,
                                 ["empirical_variance", "squared_bias", "rmse"],
                                 ["empirical variance", "squared bias", "RMSE"]):
        work = block.copy()
        work["squared_bias"] = work["bias"] ** 2
        for j, n in enumerate(sizes):
            sub = work[work.n_rct == n].groupby("estimator")[column].mean()
            vals = [sub.get(e, np.nan) for e in order]
            ax.bar(np.arange(len(order)) + j * width - 0.4 + width / 2, vals, width,
                   color=plt.cm.Blues(0.35 + 0.18 * j), label=f"$n_R$={n}" if column == "empirical_variance" else None)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(names, rotation=35, ha="right", fontsize=7)
        ax.set_title(title, fontsize=9)
        ax.set_yscale("log")
        _clean(ax)
    axes[0].legend(fontsize=6.5, frameon=False)
    fig.tight_layout()
    fig.savefig(MATERIALS / "figure2_variance_bias_rmse.png", dpi=200)
    plt.close(fig)


def figure_inference() -> None:
    """Figure 3: does the variance estimate track the spread, and does the
    interval keep its coverage?"""
    rep = pd.read_csv(MATERIALS / "ate1_when_fusion_helps_reporting.csv")
    block = rep[rep.axis == "rct_size"]
    order = ["rct_only", "omega_only", "lambda_only", "joint", "shrinkage", "adaptive",
             "naive_pool", "obs_transported"]
    names = ["trial only", "omega only", "lambda only", "joint", "shrinkage", "adaptive",
             "pooling", "OBS AIPW"]
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.6))
    ratio = block.groupby("estimator")["variance_ratio_hat_to_empirical"].mean()
    axes[0].barh(range(len(order)), [ratio.get(e, np.nan) for e in order], color=NAVY, alpha=0.85)
    axes[0].axvline(1.0, color=RUST, lw=1.0)
    axes[0].set_yticks(range(len(order)))
    axes[0].set_yticklabels(names, fontsize=7.5)
    axes[0].set_xlabel("mean estimated variance over empirical variance", fontsize=8)
    pooled = pd.read_csv(MATERIALS / "ate1_when_fusion_helps_pooled.csv")
    cov_p = pooled[pooled.quantity == "covered"].set_index("estimator")
    len_p = pooled[pooled.quantity == "ci_length"].set_index("estimator")
    centre = np.array([cov_p["value"].get(e, np.nan) for e in order])
    lo = np.array([cov_p["low"].get(e, np.nan) for e in order])
    hi = np.array([cov_p["high"].get(e, np.nan) for e in order])
    axes[1].barh(range(len(order)), centre,
                 xerr=[np.nan_to_num(centre - lo), np.nan_to_num(hi - centre)],
                 color=TEAL, alpha=0.85, error_kw={"lw": 0.8})
    axes[1].axvline(0.95, color=RUST, lw=1.0)
    axes[1].set_yticks(range(len(order)))
    axes[1].set_yticklabels([])
    axes[1].set_xlabel("coverage of the 95 percent interval", fontsize=8)
    for i, e in enumerate(order):
        length = len_p["value"].get(e, np.nan)
        if np.isfinite(length):
            axes[1].text(0.02, i, f"length {length:.2f}", fontsize=6, va="center", color=MUTED)
    for ax in axes:
        _clean(ax)
    fig.tight_layout()
    fig.savefig(MATERIALS / "figure3_variance_and_coverage.png", dpi=200)
    plt.close(fig)


def figure_cate_decomposition() -> None:
    """Figure 4: the conditional risk split into prediction variance and squared bias."""
    frames = []
    for stem in ("cate1_spline", "cate1_neural"):
        path = MATERIALS / f"{stem}_decomposition.csv"
        if path.exists():
            frames.append(pd.read_csv(path))
    if not frames:
        return
    dec = pd.concat(frames, ignore_index=True)
    dec = dec[(dec.axis == "rct_size") & (dec.learner == "DRF")]
    sieves = [s for s in ("spline3", "spline5", "neural") if s in set(dec.sieve_label)]
    rules = [("rct_only", "trial only"), ("selected", "selected"), ("grid_oracle", "grid oracle")]
    fig, axes = plt.subplots(1, len(sieves), figsize=(2.3 * len(sieves) + 0.8, 2.9), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, sieve in zip(axes, sieves):
        sub = dec[dec.sieve_label == sieve]
        sizes = sorted(sub.n_rct.unique())
        width = 0.8 / len(rules)
        for j, (rule, label) in enumerate(rules):
            work = sub[sub.rule == rule].copy()
            # the empirical risk carries (K-1)/K times the variance, so the bar uses
            # the same correction the aggregate file applies
            work["corrected_variance"] = work["prediction_variance"] * (work["replications"] - 1) / work["replications"]
            block = work.groupby("n_rct")[["corrected_variance", "squared_bias"]].mean()
            xs = np.arange(len(sizes)) + j * width - 0.4 + width / 2
            var = [block["corrected_variance"].get(n, np.nan) for n in sizes]
            bias = [block["squared_bias"].get(n, np.nan) for n in sizes]
            ax.bar(xs, var, width, color=NAVY, alpha=0.85,
                   label="prediction variance" if j == 0 and sieve == sieves[0] else None)
            ax.bar(xs, bias, width, bottom=var, color=RUST, alpha=0.85,
                   label="squared bias" if j == 0 and sieve == sieves[0] else None)
            for x in xs:
                ax.text(x, -0.02, label, rotation=90, fontsize=5.5, ha="center",
                        va="top", color=MUTED, transform=ax.get_xaxis_transform(),
                        clip_on=False)
        ax.set_xticks(range(len(sizes)))
        ax.set_xticklabels(sizes, fontsize=7.5)
        ax.tick_params(axis="x", pad=26)
        ax.set_xlabel("trial size $n_R$", fontsize=8, labelpad=2)
        ax.set_title(sieve, fontsize=9)
        ax.set_yscale("log")
        _clean(ax)
    axes[0].set_ylabel("conditional risk on the anchor", fontsize=8)
    axes[0].legend(fontsize=6.5, frameon=False)
    fig.tight_layout()
    fig.savefig(MATERIALS / "figure4_cate_variance_bias.png", dpi=200)
    plt.close(fig)


def judgment_table() -> None:
    """The three judgment axes of the revised plan, one row per comparison."""
    rows = []
    for stem, scope in (("conditional_variance_400", "conditional"),
                        ("ate1_when_fusion_helps", "whole procedure")):
        v = pd.read_csv(MATERIALS / f"{stem}_variance.csv")
        r = pd.read_csv(MATERIALS / f"{stem}_reporting.csv")
        for comparison in ("omega_alone", "omega_given_lambda", "whole_method"):
            block = v[v.comparison == comparison]
            if block.empty:
                continue
            orc = v[v.comparison == f"{comparison}_oracle"]
            on = block["on"].iloc[0]
            acc = r[r.estimator == on]
            base = r[r.estimator == block["off"].iloc[0]]
            rows.append({
                "scope": scope, "comparison": comparison,
                "mean_reduction_rate": round(float(block["reduction_rate"].mean()), 4),
                "cells_with_interval_above_zero": round(float(block["rate_excludes_zero"].mean()), 3),
                "oracle_mean_reduction_rate": round(float(orc["reduction_rate"].mean()), 4) if len(orc) else np.nan,
                "bias_change": round(float(acc["bias"].abs().mean() - base["bias"].abs().mean()), 5),
                "rmse_change": round(float(acc["rmse"].mean() - base["rmse"].mean()), 5),
                "coverage_change": round(float(acc["coverage"].mean() - base["coverage"].mean()), 4),
                "ci_length_change": round(float(acc["mean_ci_length"].mean() - base["mean_ci_length"].mean()), 5)})
    pd.DataFrame(rows).to_csv(MATERIALS / "judgment_axes.csv", index=False)


def main() -> int:
    figure_variance_reduction()
    figure_decomposition()
    figure_inference()
    figure_cate_decomposition()
    judgment_table()
    made = sorted(p.name for p in MATERIALS.glob("figure[1-4]_*.png"))
    print("wrote", ", ".join(made))
    return 0


if __name__ == "__main__":
    sys.exit(main())
