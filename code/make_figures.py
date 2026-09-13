"""Figures and tables for Section 5 of the paper (handoff Section 8.3)."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MATERIALS = Path(__file__).resolve().parents[1] / "materials"
PALETTE = {"rct_only": "#11233F", "joint": "#173B7A", "lambda_only": "#0E7490",
           "omega_only": "#53657A", "oracle_joint": "#9A3412",
           "shrinkage": "#7A7A7A", "adaptive": "#B0B0B0", "obs_transported": "#C0C0C0"}
LABELS = {"rct_only": "trial only", "joint": "joint", "lambda_only": r"$\lambda$ only",
          "omega_only": r"$\omega$ only", "oracle_joint": "oracle coefficients",
          "shrinkage": "shrinkage", "adaptive": "adaptive", "obs_transported": "observational AIPW"}


def figure1() -> None:
    s = pd.read_csv(MATERIALS / "ate1_when_fusion_helps_summary.csv")
    r = s[(s.metric == "rmse_ratio_to_rct") & (s.axis.isin(["rct_size", "quality"]))]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.9), sharey=True)
    size = r[r.axis == "rct_size"]
    for est in ("rct_only", "joint", "lambda_only", "omega_only", "oracle_joint",
                "shrinkage", "obs_transported"):
        block = size[size.estimator == est].groupby("n_rct")["value"].mean()
        if block.empty:
            continue
        axes[0].plot(block.index, block.values, marker="o", ms=3, lw=1.3,
                     color=PALETTE.get(est, "#888"), label=LABELS.get(est, est))
    axes[0].set_xscale("log")
    axes[0].set_xticks([50, 100, 200, 400])
    axes[0].set_xticklabels([50, 100, 200, 400])
    axes[0].minorticks_off()
    axes[0].set_xlabel("trial size $n_R$")
    axes[0].set_ylabel("RMSE ratio to trial only")
    axes[0].set_title("correctly specified OBS model", fontsize=9)
    axes[0].axhline(1.0, color="#999", lw=0.7, ls=":")

    qual = r[r.axis == "quality"]
    order = [(0.0, "flexible"), (1.0, "linear"), (2.0, "flexible")]
    labels = ["no confounding", "misspecified", "strong confounding"]
    for est in ("joint", "lambda_only", "oracle_joint", "shrinkage", "obs_transported"):
        vals = []
        for conf, spec in order:
            b = qual[(qual.estimator == est) & (qual.confounding == conf) & (qual.spec == spec)]
            vals.append(float(b["value"].mean()) if len(b) else np.nan)
        axes[1].plot(range(len(order)), vals, marker="s", ms=3, lw=1.3,
                     color=PALETTE.get(est, "#888"), label=LABELS.get(est, est))
    axes[1].set_xticks(range(len(order)))
    axes[1].set_xticklabels(labels, fontsize=7)
    axes[1].set_title(r"prediction quality at $n_R=100$", fontsize=9)
    axes[1].axhline(1.0, color="#999", lw=0.7, ls=":")
    handles, labels_ = axes[0].get_legend_handles_labels()
    _fig1_legend = (handles, labels_)
    for ax in axes:
        ax.tick_params(labelsize=7.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    fig.legend(*_fig1_legend, fontsize=7, frameon=False, ncol=4,
               loc="upper center", bbox_to_anchor=(0.5, 0.995))
    fig.savefig(MATERIALS / "figure1_ate_rmse.png", dpi=200)
    plt.close(fig)


def figure_shift_routes() -> None:
    """Figure 6: what the transport weight costs when the two covariate laws
    differ, for the average effect and for the conditional effect."""
    routes = ["oracle", "classifier", "balanced_x", "unit"]
    names = ["oracle r", "classifier", "balanced", "ignore (r=1)"]
    colors = ["#9A3412", "#173B7A", "#0E7490", "#7A7A7A"]

    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.0))
    xs = np.arange(len(routes))

    ax = axes[0]
    path = MATERIALS / "ate2_shift_summary.csv"
    if path.exists():
        d = pd.read_csv(path)
        d = d[(d.panel == "synthetic") & (d.estimator == "joint") & (d.metric == "bias")]
        w = 0.38
        for off, level, alpha in ((-w / 2, "mild", 0.45), (w / 2, "strong", 1.0)):
            vals = [abs(float(d[(d["shift"] == level) & (d.ratio == r)]["value"].mean()))
                    for r in routes]
            ax.bar(xs + off, vals, w, color=colors, alpha=alpha, edgecolor="none")
        ax.set_ylabel("absolute bias", fontsize=8)
        ax.set_title("ATE under covariate shift", fontsize=9)
        ax.text(0.02, 0.93, "pale = mild shift,  solid = strong shift",
                transform=ax.transAxes, fontsize=6.8, color="#31415A")

    for ax, stem, title, level in ((axes[1], "cate3_shift_spline",
                                    "CATE, spline basis", "strong"),
                                   (axes[2], "cate3_shift_neural",
                                    "CATE, neural basis", "strong")):
        path = MATERIALS / f"{stem}_replications.csv"
        if not path.exists():
            continue
        r = pd.read_csv(path)
        r = r[r.learner.isin(["DRF", "RF"]) & (r["shift"] == level) & (r.n_rct == 200)]
        vals, chan = [], []
        for route in routes:
            b = r[r.route == route]
            vals.append(float(np.median(b["risk_selected"] / b["risk_rct_only"])))
            chan.append(float(np.median(b["risk_omega_only"] / b["risk_rct_only"])))
        ax.bar(xs, vals, 0.56, color=colors, edgecolor="none", label="both channels")
        ax.plot(xs, chan, marker="o", ms=4, lw=1.2, color="#31415A", ls="--",
                label=r"$\omega$ channel alone")
        ax.set_ylabel("median risk ratio to RCT only", fontsize=8)
        ax.set_title(f"{title}, strong shift", fontsize=9)
        ax.set_ylim(0, max(max(vals), max(chan)) * 1.35)
        ax.legend(fontsize=6.5, frameon=False, loc="upper left")

    for ax in axes:
        ax.set_xticks(xs)
        ax.set_xticklabels(names, fontsize=7, rotation=28, ha="right")
        ax.tick_params(labelsize=7.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(MATERIALS / "figure6_shift_routes.png", dpi=200)
    plt.close(fig)


def figure_five_rules() -> None:
    """Figure 0: the five coefficient rules and their two oracle versions, side by
    side, for the average effect and for the conditional effect.  Everything is a
    ratio to the trial-only estimator in the same cell, so 1.0 is the baseline."""
    order = ["rct_only", "lambda_only", "omega_only", "joint",
             "oracle_lambda_only", "oracle_omega_only", "oracle_joint"]
    names = ["RCT only", r"$\lambda$ est", r"$\omega$ est", r"$\lambda+\omega$ est",
             r"$\lambda$ oracle", r"$\omega$ oracle", r"$\lambda+\omega$ oracle"]
    fill = ["#11233F", "#0E7490", "#53657A", "#173B7A", "#C2765A", "#C0A08C", "#9A3412"]

    rep = pd.read_csv(MATERIALS / "ate1_when_fusion_helps_reporting.csv")
    cell = ["family", "axis", "n_rct", "n_obs", "confounding", "spec"]
    base_r = rep[rep.estimator == "rct_only"].set_index(cell)["rmse"]
    base_v = rep[rep.estimator == "rct_only"].set_index(cell)["empirical_variance"]
    idx = rep.set_index(cell)
    rmse, var = [], []
    for est in order:
        b = idx[idx.estimator == est]
        rmse.append(float((b["rmse"] / base_r.reindex(b.index)).mean()))
        var.append(float((b["empirical_variance"] / base_v.reindex(b.index)).mean()))

    cols = {"rct_only": "risk_rct_only", "lambda_only": "risk_lambda_only",
            "omega_only": "risk_omega_only", "joint": "risk_selected",
            "oracle_lambda_only": "risk_lambda_only_oracle",
            "oracle_omega_only": "risk_omega_only_oracle",
            "oracle_joint": "risk_grid_oracle"}
    frames = [pd.read_csv(MATERIALS / f"{s}_replications.csv")
              for s in ("cate1_spline", "cate1_neural")
              if (MATERIALS / f"{s}_replications.csv").exists()]
    cate = None
    if frames:
        c = pd.concat(frames, ignore_index=True)
        c = c[c.learner.isin(["DRF", "RF"])]
        cate = [float(np.median(c[cols[e]] / c["risk_rct_only"])) for e in order]

    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.2))
    xs = np.arange(len(order))

    ax = axes[0]
    w = 0.38
    ax.bar(xs - w / 2, rmse, w, color=fill, edgecolor="none")
    ax.bar(xs + w / 2, var, w, color=fill, edgecolor="none", alpha=0.45, hatch="///")
    ax.set_title("ATE: error and variance", fontsize=9)
    ax.set_ylabel("ratio to RCT only", fontsize=8)
    ax.text(0.02, 0.93, "solid bar = RMSE,  hatched bar = variance",
            transform=ax.transAxes, fontsize=6.8, color="#31415A")

    ax = axes[1]
    if cate is not None:
        ax.bar(xs, cate, 0.62, color=fill, edgecolor="none")
    ax.set_title("CATE: risk on the anchor", fontsize=9)
    ax.set_ylabel("median ratio to RCT only", fontsize=8)

    for ax in axes[:2]:
        ax.axhline(1.0, color="#999", lw=0.8, ls=":")
        ax.set_xticks(xs)
        ax.set_xticklabels(names, fontsize=7, rotation=32, ha="right")
        ax.set_ylim(0, 1.2)

    ax = axes[2]
    dec_paths = [MATERIALS / f"{s}_decomposition.csv" for s in ("cate1_spline", "cate1_neural")]
    dec = [pd.read_csv(p) for p in dec_paths if p.exists()]
    if dec:
        d = pd.concat(dec, ignore_index=True)
        piv = d.pivot_table(index="rule", values=["prediction_variance", "squared_bias"],
                            aggfunc="mean")
        piv = piv.reindex([r for r in ("rct_only", "selected", "grid_oracle") if r in piv.index])
        z = np.arange(len(piv))
        ax.bar(z, piv["prediction_variance"], 0.55, color="#173B7A", label="prediction variance")
        ax.bar(z, piv["squared_bias"], 0.55, bottom=piv["prediction_variance"],
               color="#9A3412", label="squared bias")
        ax.set_xticks(z)
        ax.set_xticklabels(["RCT only", r"$\lambda+\omega$ est", r"$\lambda+\omega$ oracle"],
                           fontsize=7, rotation=32, ha="right")
        for zi, (v, b) in enumerate(zip(piv["prediction_variance"], piv["squared_bias"])):
            ax.text(zi, v + b, f"{v + b:.2f}", ha="center", va="bottom", fontsize=6.8)
        ax.set_ylabel("CATE risk on the anchor", fontsize=8)
        ax.set_title("CATE: what the risk is made of", fontsize=9)
        ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    for a in axes:
        a.tick_params(labelsize=7.5)
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(MATERIALS / "figure0_five_rules.png", dpi=200)
    plt.close(fig)


def figure_pooled_comparator() -> None:
    """Figure 5: the source-indicator pooled learner against the fusion rule, and
    the reason it wins.  Its indicator is a constant, so it absorbs whatever part
    of the observational bias does not vary with the covariates."""
    frames = []
    for stem in ("cate1_spline", "cate1_neural"):
        path = MATERIALS / f"{stem}_replications.csv"
        if path.exists():
            frames.append(pd.read_csv(path))
    if not frames:
        return
    rep = pd.concat(frames, ignore_index=True)
    keys = ["family", "axis", "n_rct", "confounding", "heterogeneity",
            "sieve_label", "replication"]
    fusion = rep[rep.learner == "DRF"].set_index(keys)
    pooled = rep[rep.learner == "domain_indicator_pool"].set_index(keys)
    joined = fusion[["risk_selected", "risk_rct_only"]].join(
        pooled[["risk_selected"]], rsuffix="_pooled", how="inner")
    joined.columns = ["fusion", "trial_only", "pooled"]
    strong = joined.reset_index()
    strong = strong[strong.confounding == strong.confounding.max()]
    order = ["bounded", "scm1", "scm2", "scm3", "star"]
    med = strong.groupby("family")[["trial_only", "fusion", "pooled"]].median()
    med = med.reindex([f for f in order if f in med.index])

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9))
    ax = axes[0]
    width, xs = 0.26, np.arange(len(med))
    for off, col, color, label in ((-width, "trial_only", "#11233F", "trial only"),
                                   (0.0, "fusion", "#173B7A", "fusion, selected"),
                                   (width, "pooled", "#9A3412", "source-indicator pooled")):
        ax.bar(xs + off, med[col].to_numpy(), width, color=color, label=label)
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(med.index, fontsize=7.5)
    ax.set_ylabel("median risk on the anchor", fontsize=8)
    ax.set_title("strongest confounding in each design", fontsize=9)
    ax.legend(fontsize=6.5, frameon=False, loc="lower left", ncol=1)
    ax.set_ylim(top=med.to_numpy().max() * 3.0)

    ax = axes[1]
    shape_path = MATERIALS / "confounding_shape.csv"
    if shape_path.exists():
        shape = pd.read_csv(shape_path)
        shape = shape[shape.confounding == shape.confounding.max()]
        for fam, g in shape.groupby("family"):
            left = g["bias_in_quintile"] - g["mean_bias"]
            ax.plot(g["quintile"], left, marker="o", ms=3, lw=1.2, label=fam)
        span = float(shape["tau_sd"].mean())
        ax.axhspan(-span, span, color="#C9D3E0", alpha=0.55, zorder=0)
        ax.text(3.0, span * 0.72, "one standard deviation of the true effect",
                ha="center", fontsize=6.5, color="#31415A")
        ax.axhline(0.0, color="#999", lw=0.7, ls=":")
        ax.set_xticks([1, 2, 3, 4, 5])
        ax.set_xlabel("quintile of the leading covariate", fontsize=8)
        ax.set_ylabel("observational bias left\nafter removing its average", fontsize=8)
        ax.set_title("what an indicator cannot absorb", fontsize=9)
        ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    for a in axes:
        a.tick_params(labelsize=7.5)
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(MATERIALS / "figure5_pooled_comparator.png", dpi=200)
    plt.close(fig)


def _cate_frames() -> pd.DataFrame:
    frames = []
    for stem in ("cate1_spline", "cate1_neural"):
        path = MATERIALS / f"{stem}_replications.csv"
        if path.exists():
            frames.append(pd.read_csv(path))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def figure2() -> None:
    frame = _cate_frames()
    if frame.empty:
        return
    core = frame[frame.learner.isin(["DRF", "RF"])]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.9), sharey=True)
    for ax, learner in zip(axes, ("DRF", "RF")):
        block = core[(core.learner == learner) & (core.axis == "rct_size")]
        for sieve, color in (("spline3", "#173B7A"), ("spline5", "#0E7490"), ("neural", "#9A3412")):
            sub = block[block.sieve_label == sieve]
            if sub.empty:
                continue
            sel = sub.groupby("n_rct")["risk_selected"].median()
            base = sub.groupby("n_rct")["risk_rct_only"].median()
            ax.plot(sel.index, sel.values, marker="o", ms=3, lw=1.3, color=color,
                    label=f"{sieve}, selected")
            ax.plot(base.index, base.values, marker="o", ms=3, lw=1.0, ls="--",
                    color=color, alpha=0.6, label=f"{sieve}, trial only")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks([100, 200, 400])
        ax.set_xticklabels([100, 200, 400])
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.set_xlabel("trial size $n_R$")
        ax.set_title(f"{learner} learner", fontsize=9)
        ax.tick_params(labelsize=7.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("median CATE risk")
    axes[0].legend(fontsize=6, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(MATERIALS / "figure2_cate_risk.png", dpi=200)
    plt.close(fig)


def figure_regret() -> None:
    frame = _cate_frames()
    if frame.empty:
        return
    core = frame[frame.learner.isin(["DRF", "RF"])].dropna(subset=["regret_over_radius"])
    if core.empty:
        return
    fig, ax = plt.subplots(figsize=(3.3, 2.3))
    values = core["regret_over_radius"].to_numpy()
    values = values[values > 0]
    ax.hist(np.log10(values), bins=40, color="#173B7A", alpha=0.85)
    ax.axvline(0.0, color="#9A3412", lw=1.2)
    ax.set_xlabel(r"$\log_{10}$ of regret divided by $2\varepsilon_{\mathrm{eval}}(0.05)$", fontsize=8)
    ax.set_ylabel("replications", fontsize=8)
    ax.tick_params(labelsize=7.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(MATERIALS / "figureA_regret_vs_radius.png", dpi=200)
    plt.close(fig)


def figure_coefficients() -> None:
    path = MATERIALS / "ate1_when_fusion_helps_replications.csv"
    if not path.exists():
        return
    frame = pd.read_csv(path)
    block = frame[(frame.estimator == "joint") & frame.oracle_lambda.notna()]
    if block.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(5.4, 2.4))
    for ax, (est, orc, name) in zip(axes, [("lambda_used", "oracle_lambda", r"$\lambda$"),
                                           ("omega_used", "oracle_omega", r"$\omega$")]):
        ax.scatter(block[orc], block[est], s=2, alpha=0.12, color="#173B7A", edgecolors="none")
        ax.plot([0, 1], [0, 1], color="#9A3412", lw=1.0)
        ax.set_xlabel(f"oracle {name}", fontsize=8)
        ax.set_ylabel(f"selected {name}", fontsize=8)
        ax.tick_params(labelsize=7.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(MATERIALS / "figureA_coefficient_recovery.png", dpi=200)
    plt.close(fig)


def tables() -> None:
    rows = []
    s = pd.read_csv(MATERIALS / "ate1_when_fusion_helps_summary.csv")
    cov = s[(s.metric == "coverage") & (s.axis == "rct_size")]
    for est in ("rct_only", "joint", "shrinkage", "adaptive", "naive_pool", "obs_transported"):
        block = cov[cov.estimator == est].groupby("n_rct")["value"].mean()
        for n, v in block.items():
            rows.append({"table": "1_coverage", "estimator": est, "n_rct": n, "value": round(float(v), 3)})
    frame = _cate_frames()
    if not frame.empty:
        core = frame[frame.learner.isin(["DRF", "RF"])]
        for sieve in core.sieve_label.dropna().unique():
            sub = core[core.sieve_label == sieve]
            rows.append({"table": "1_regret", "estimator": sieve, "n_rct": np.nan,
                         "value": round(float(sub["regret_within_radius"].mean()), 3)})
            rows.append({"table": "1_regret_ratio_median", "estimator": sieve, "n_rct": np.nan,
                         "value": round(float(sub["regret_over_radius"].median()), 6)})
    pd.DataFrame(rows).to_csv(MATERIALS / "paper_tables.csv", index=False)


def main() -> int:
    figure1()
    figure2()
    figure_five_rules()
    figure_shift_routes()
    figure_pooled_comparator()
    figure_regret()
    figure_coefficients()
    tables()
    made = sorted(p.name for p in MATERIALS.glob("figure*.png"))
    print("wrote", ", ".join(made) if made else "no figures")
    return 0


if __name__ == "__main__":
    sys.exit(main())
