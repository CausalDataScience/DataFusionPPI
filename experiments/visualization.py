"""Draws the figures of Sec. 5 and Appendix C and the ATE table of Appendix C (Tables 13 and 14) from the results
of simulation.py, and prints the numbers they show.

    python3 visualization.py                       reads results/, writes figures/
    python3 visualization.py --zoom --out figures/zoom
    python3 visualization.py --law weak_shift
    python3 visualization.py --baseline full       trial-only references with nuisance on two thirds of the trial

Style: small panels with large text, thick lines, a 95% interval around
every curve (shaded band and whiskers), and the legends as their own files.  Each panel is also saved on its
own (ate_<column>, cate_<column>, sweep_n<size>, obs_<column>) so that LaTeX can place the panels with \\subfloat.

ATE panel      top: mean squared error; bottom: its ratio to RCT only.  Columns MAIN: SCM 1 with r_0 = 1, SCM 1
               with r_0 != 1, IHDP and ACIC 2016, the last three under the OBS covariate law LAW.
CATE panel     DRF and RF with their comparisons, every risk divided by the mean risk of the RCT-only DR-learner.
Sweep panel    SCM 1 under LAW as its confounding bias grows (horizontal axis: root mean square of delta in units
               of sd tau_0); top: ATE, bottom: CATE, both relative to RCT only; one panel per trial size.
OBS panel      as the sweep panel, against the OBS sample size N; one panel per (design, trial size) of OBS_SWEEP.
Absolute       the same panels with the mean squared error and the mean risk themselves (*_abs).
Zoom           --zoom: the ratio panels leave out the naive pool, so that the intervals of the other curves show;
               the bias sweep keeps it, since its crossing is the point of that figure.
Table          table_ate.csv (every trial size) and table_ate.tex (n = 300, 900, 3000): bias, standard deviation,
               root mean squared error, mean standard error, coverage and mean width of the 95% Wald interval.
Failure        only with --failure: the designs with a useless OBS (SCM 3, lung cancer); not part of the paper.
Intervals are 95% percentile intervals of a bootstrap that resamples whole replications (2,000 draws), of the
mean squared error, of the mean risk, or of their ratio to RCT only.
"""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FormatStrFormatter

LAW, N_MAIN = "strong_shift", 15_000              # the OBS covariate law of the shifted panels; the default OBS size
MAIN = (("scm1", "same", "Synthetic, $r_0{=}1$", "scm1_same"), ("scm1", LAW, r"Synthetic, $r_0{\neq}1$", "scm1_shift"),
        ("ihdp", LAW, "IHDP", "ihdp"), ("acic", LAW, "ACIC", "acic"))         # (source, law, title, file tag)
OBS_SWEEP = (("scm1", 300, "Synthetic, $n{=}300$", "scm1_n300"), ("scm1", 900, "Synthetic, $n{=}900$", "scm1_n900"),
             ("ihdp", 300, "IHDP, $n{=}300$", "ihdp_n300"), ("ihdp", 900, "IHDP, $n{=}900$", "ihdp_n900"))
FAILURE = (("scm3", LAW, "SCM 3", "scm3"), ("lung", "native", "Lung cancer", "lung"))
SWEEP, SCM1_BIAS = "scm1_bias", 1.5               # sources of the bias sweep; "scm1" itself has the bias size 1.5
# Curves: (learner, method, label, colour, line, marker).  ATE curves are relative to RCT only; CATE curves to
# the RCT-only DR-learner, and a dashed CATE line is the version with the pseudo-outcome of the R-learner.
RCT = ("", "rct_only", "RCT only", "#52514e", "-", "o")
ATE_CURVES = (("", "omega_only", "PPI++", "#eb6834", "-", "o"),
              ("", "shrinkage", "Shrinkage", "#d62728", "-", "s"),
              ("", "pretest", "Pretest", "#8c564b", "-", "^"),
              ("", "fusion", "Fusion", "#2a78d6", "-", "o"),
              ("", "naive_pool", "Naive pool", "#1baf7a", "-", "o"))
RCT_DR = ("DRF", "rct_only", "RCT (DR)", "#52514e", "-", "o")
CATE_CURVES = (("RF", "rct_only", "RCT (R)", "#52514e", "--", "o"),
               ("DRF", "omega_only", "PPI++ DRF", "#eb6834", "-", "o"),
               ("RF", "omega_only", "PPI++ RF", "#eb6834", "--", "o"),
               ("DRF", "two_step", "2-step DR", "#d39a00", "-", "s"),
               ("RF", "two_step", "2-step R", "#d39a00", "--", "s"),
               ("RF", "integrative", "IR", "#e377c2", "-.", "D"),
               ("DRF", "fusion", "DRF", "#2a78d6", "-", "o"),
               ("RF", "fusion", "RF", "#9b4dca", "--", "o"),
               ("pool", "naive_pool", "Naive pool", "#1baf7a", "-", "o"))
PANEL, BOOTSTRAP = 1.9, 2000                      # panel width in inches; bootstrap draws
H2 = PANEL * 1.15 + 0.3                           # height of a two-row panel (ATE, sweep, obs): 2.49 in
H1 = PANEL * 0.70 + 0.3                           # height of a one-row panel (CATE, sweep ATE or CATE): 1.63 in
STYLE = {"font.size": 12, "axes.titlesize": 12.5, "axes.labelsize": 12.5, "xtick.labelsize": 10.5,
         "ytick.labelsize": 10.5, "legend.fontsize": 12, "lines.linewidth": 2.4, "lines.markersize": 4.5,
         "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 1.0, "pdf.fonttype": 42}


# ---------------------------------------------------------------------------- results
def load(folder):
    """values[(source, law, N)][n_RCT][(estimand, learner, method)] = array over replications, in replication
    order; ("SE", learner, method) holds the estimated standard errors of the ATE."""
    values, failed = {}, 0
    for path in sorted(Path(folder).glob("*.csv")):
        with open(path, newline="") as handle:
            for r in csv.DictReader(handle):
                if r["estimand"] == "FAILED":
                    failed += 1
                    continue
                N = int(float(r.get("N") or N_MAIN))
                cell = values.setdefault((r["source"], r["shift"], N), {}).setdefault(3 * int(r["m"]), {})
                rep = int(r["rep"])
                cell.setdefault((r["estimand"], r["learner"], r["method"]), {})[rep] = float(r["value"])
                if r["se"]:
                    cell.setdefault(("SE", r["learner"], r["method"]), {})[rep] = float(r["se"])
    if failed:
        print(f"{failed} replications failed and are left out")
    for cells in values.values():
        for cell in cells.values():
            reps = sorted(set.intersection(*(set(v) for v in cell.values())))
            for key in cell:
                cell[key] = np.array([cell[key][rep] for rep in reps])
    return values


def ratio_with_interval(numerator, denominator, rng):
    """mean(numerator) / mean(denominator) and its 95% bootstrap interval, resampling replications."""
    draws = rng.integers(0, len(numerator), size=(BOOTSTRAP, len(numerator)))
    ratios = numerator[draws].mean(1) / denominator[draws].mean(1)
    return numerator.mean() / denominator.mean(), np.quantile(ratios, 0.025), np.quantile(ratios, 0.975)


# ---------------------------------------------------------------------------- curves and axes
def curve(ax, x, numerators, denominators, rng, spec, tag):
    """One curve with its 95% interval as a band and whiskers: the ratio of means, or the mean when the
    denominators are ones."""
    _, _, label, color, line, marker = spec
    ratio, low, high = (np.array(v) for v in zip(*(ratio_with_interval(a, b, rng)
                                                  for a, b in zip(numerators, denominators))))
    ax.fill_between(x, low, high, color=color, alpha=0.30, lw=0)
    ax.errorbar(x, ratio, yerr=[ratio - low, high - ratio], fmt=marker + line, color=color, label=label,
                capsize=3.5, capthick=1.3, elinewidth=1.3)
    for at, value, lo, hi in zip(x, ratio, low, high):
        print(f"{tag} {at:g}: {label:12s} {value:9.4f}  [{lo:.4f}, {hi:.4f}]")


def ate_curves(ax, x, cells, rng, tag, absolute=False, zoom=False):
    """Mean squared error of the ATE relative to RCT only, or the mean squared error itself (absolute).
    zoom leaves the naive pool out of a ratio panel."""
    squared = lambda method: [c[("ATE", "", method)] ** 2 for c in cells]
    if absolute:
        base, specs = [np.ones(len(c[("ATE", "", RCT[1])])) for c in cells], (RCT,) + ATE_CURVES
    else:
        base, specs = squared(RCT[1]), ATE_CURVES
        ax.axhline(1.0, color=RCT[3], lw=1.2, ls=(0, (4, 2.5)))
    for spec in specs:
        if zoom and not absolute and spec[1] == "naive_pool":
            continue
        curve(ax, x, squared(spec[1]), base, rng, spec, f"ATE {tag}")


def cate_curves(ax, x, cells, rng, tag, absolute=False, zoom=False):
    """CATE risks relative to the RCT-only DR-learner, or the mean risks themselves (absolute).
    zoom leaves the naive pool out of a ratio panel."""
    key = ("CATE",) + RCT_DR[:2]
    if absolute:
        base, specs = [np.ones(len(c[key])) for c in cells], (RCT_DR,) + CATE_CURVES
    else:
        base, specs = [c[key] for c in cells], CATE_CURVES
        ax.axhline(1.0, color=RCT_DR[3], lw=2.6, label=RCT_DR[2])
    for spec in specs:
        if zoom and not absolute and spec[1] == "naive_pool":
            continue
        curve(ax, x, [c[("CATE", spec[0], spec[1])] for c in cells], base, rng, spec, f"CATE {tag}")


def style(ax, ticks, labels=None, log_x=True, ratio=False):
    """Log axes, the given horizontal ticks, and on a ratio axis plain ticks with the line at one in view."""
    if log_x:
        ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels or [f"{t:g}" for t in ticks])
    ax.minorticks_off()
    if ratio:
        low, high = ax.get_ylim()
        high = max(high, 1.08)
        ax.set_ylim(low, high)
        ax.yaxis.set_major_locator(FixedLocator(ratio_ticks(low, high)))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.grid(axis="y", color="#e3e2de", lw=0.5)
    ax.set_axisbelow(True)


def ratio_ticks(low, high):
    """At most four plain ticks inside [low, high] on a ratio axis, one of them at one."""
    wide = high / low > 3                            # 0.5, 1, 2 rather than crowded ticks below one
    candidates = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50) if wide else \
                 (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1, 1.5, 2, 3)
    ticks = [t for t in candidates if low <= t <= high]
    while len(ticks) > 4:                            # keep every second tick, counted from one
        ticks = [t for i, t in enumerate(ticks) if (i - ticks.index(1)) % 2 == 0]
    return ticks


def ends(sizes):
    """The smallest, the middle and the largest trial size: the ticks of the horizontal axis."""
    return [sizes[0], sizes[len(sizes) // 2], sizes[-1]]


# ---------------------------------------------------------------------------- panels
def ate_panel(top, bottom, by_size, title, rng, absolute, zoom=False):
    """One ATE column: the mean squared error, and below it its ratio to RCT only (no second row if absolute)."""
    sizes = sorted(by_size)
    cells = [by_size[n] for n in sizes]
    ate_curves(top, sizes, cells, rng, title, absolute=True)
    top.set_title(title)
    top.set_ylabel("MSE")
    style(top, ends(sizes))
    if bottom is None:
        top.set_xlabel("RCT size $n$")
        return
    top.tick_params(labelbottom=False)
    ate_curves(bottom, sizes, cells, rng, title, zoom=zoom)
    bottom.set_ylabel("MSE / RCT")
    style(bottom, ends(sizes), ratio=True)
    bottom.set_xlabel("RCT size $n$")


def cate_panel(ax, by_size, title, rng, absolute, zoom=False):
    sizes = sorted(by_size)
    cate_curves(ax, sizes, [by_size[n] for n in sizes], rng, title, absolute=absolute, zoom=zoom)
    ax.set_title(title)
    ax.set_ylabel("risk" if absolute else "risk / RCT")
    style(ax, ends(sizes), ratio=not absolute)
    ax.set_xlabel("RCT size $n$")


def x_panel(top, bottom, by_x, title, axis, rng, absolute, zoom=False):
    """One column of a sweep: ATE (top) and CATE (bottom) against the swept quantity.
    axis = (ticks, tick labels, axis label, log scale)."""
    ticks, labels, name, log_x = axis
    xs = sorted(by_x)
    cells = [by_x[x] for x in xs]
    ate_curves(top, xs, cells, rng, f"{title}:", absolute=absolute, zoom=zoom)
    top.set_title(title)
    top.set_ylabel("MSE" if absolute else "MSE / RCT")
    style(top, ticks, labels, log_x=log_x, ratio=not absolute)
    top.tick_params(labelbottom=False)
    cate_curves(bottom, xs, cells, rng, f"{title}:", absolute=absolute, zoom=zoom)
    bottom.set_ylabel("risk" if absolute else "risk / RCT")
    style(bottom, ticks, labels, log_x=log_x, ratio=not absolute)
    bottom.set_xlabel(name)


# ---------------------------------------------------------------------------- figures
def save(fig, path):
    fig.tight_layout(pad=0.3)
    fig.savefig(path)
    plt.close(fig)


def present(values, columns):
    """The columns (source, law, title, tag) whose results are there, at the default OBS size."""
    found = [(values[(s, law, N_MAIN)], title, tag) for s, law, title, tag in columns if (s, law, N_MAIN) in values]
    missing = [f"{s} under {law}" for s, law, _, _ in columns if (s, law, N_MAIN) not in values]
    if missing:
        print("not in the results:", ", ".join(missing))
    return found


def draw_columns(found, rows, height, draw, out, name):
    """Each column on its own (<name>_<tag>) and all columns side by side (<name>)."""
    for by, title, tag in found:
        fig, axes = plt.subplots(rows, 1, figsize=(PANEL, height), sharex=True, squeeze=False)
        draw(axes[:, 0], by, title)
        save(fig, out(f"{name}_{tag}"))
    fig, axes = plt.subplots(rows, len(found), figsize=(PANEL * len(found), height), sharex="col", squeeze=False)
    for column, (by, title, tag) in enumerate(found):
        draw(axes[:, column], by, title)
    for ax in axes[:, 1:].ravel():
        ax.set_ylabel("")
    save(fig, out(name))


def draw_ate(values, columns, out, absolute=False, zoom=False):
    rows, height = (1, H1) if absolute else (2, H2)
    draw = lambda axes, by, title: ate_panel(axes[0], None if absolute else axes[1], by, title,
                                             np.random.default_rng(0), absolute, zoom)
    draw_columns(present(values, columns), rows, height, draw, out, "ate")


def draw_cate(values, columns, out, absolute=False, zoom=False):
    draw = lambda axes, by, title: cate_panel(axes[0], by, title, np.random.default_rng(0), absolute, zoom)
    draw_columns(present(values, columns), 1, H1, draw, out, "cate")


def draw_sweep(values, law, out, absolute=False):
    """The bias sweep, if its sources are among the results.  It keeps the naive pool in every mode."""
    designs = {float(s[len(SWEEP):]): s for s, at, N in values if s.startswith(SWEEP) and at == law and N == N_MAIN}
    if not designs or ("scm1", law, N_MAIN) not in values:
        return False
    designs[SCM1_BIAS] = "scm1"
    sizes = ends(sorted(set.intersection(*(set(values[(s, law, N_MAIN)]) for s in designs.values()))))
    found = [({bias: values[(s, law, N_MAIN)][n] for bias, s in designs.items()}, f"$n = {n}$", f"n{n}")
             for n in sizes]
    axis = ([0, 0.5, 1, 1.5], None, r"bias / sd $\tau_0$", False)
    draw = lambda axes, by, title: x_panel(axes[0], axes[1], by, title, axis, np.random.default_rng(0), absolute)
    draw_columns(found, 2, H2, draw, out, "sweep")
    by_bias, n = found[-1][0], sizes[-1]              # sweep_ate_n<n> and sweep_cate_n<n> at the largest n
    biases = sorted(by_bias)
    for name, curves, label in (("ate", ate_curves, "MSE" if absolute else "MSE / RCT"),
                                ("cate", cate_curves, "risk" if absolute else "risk / RCT")):
        fig, ax = plt.subplots(figsize=(PANEL, H1))
        curves(ax, biases, [by_bias[b] for b in biases], np.random.default_rng(0), f"sweep {name} n={n}:", absolute)
        ax.set_title(f"$n = {n}$")
        ax.set_ylabel(label)
        style(ax, axis[0], log_x=False, ratio=not absolute)
        ax.set_xlabel(axis[2])
        save(fig, out(f"sweep_{name}_n{n}"))
    return True


def draw_obs(values, law, out, absolute=False, zoom=False):
    """The OBS-size sweep, if its results are there: one column per (design, trial size) of OBS_SWEEP."""
    found = []
    for source, n, title, tag in OBS_SWEEP:
        by_N = {N: cells[n] for (s, at, N), cells in values.items() if s == source and at == law and n in cells}
        if len(by_N) > 1:
            found.append((by_N, title, tag))
    if not found:
        return False
    Ns = sorted(set.union(*(set(by) for by, _, _ in found)))
    axis = (Ns, [f"{N / 1000:g}k" for N in Ns], "OBS size $N$", True)
    draw = lambda axes, by, title: x_panel(axes[0], axes[1], by, title, axis, np.random.default_rng(0), absolute,
                                           zoom)
    draw_columns(found, 2, H2, draw, out, "obs")
    return True


def draw_failure(values, columns, out, absolute=False):
    """The designs with a useless OBS, if they are among the results."""
    found = present(values, columns)
    if not found:
        return False
    draw = lambda axes, by, title: x_panel(axes[0], axes[1], by, title,
                                           (ends(sorted(by)), None, "RCT size $n$", True),
                                           np.random.default_rng(0), absolute)
    draw_columns(found, 2, H2, draw, out, "failure")
    return True


def draw_legends(out, zoom=False):
    """The legends as their own figures, in two rows when they are long: legend_ate and legend_cate.  The CATE
    legend keeps each DR curve above its R version and puts IR and the naive pool in the last column."""
    ate = (RCT,) + ATE_CURVES
    cate = (RCT_DR,) + tuple(sorted((spec for spec in CATE_CURVES if not (zoom and spec[1] == "naive_pool")),
                                    key=lambda spec: spec[1] in ("integrative", "naive_pool")))
    for name, specs in (("ate", ate), ("cate", cate)):
        columns = len(specs) if len(specs) <= 5 else -(-len(specs) // 2)
        fig = plt.figure(figsize=(1.6 * columns, 0.55))
        handles = [plt.Line2D([], [], color=color, ls=line, marker=marker) for _, _, _, color, line, marker in specs]
        fig.legend(handles, [spec[2] for spec in specs], loc="center", ncol=columns, frameon=False,
                   handlelength=2.2, columnspacing=0.8, handletextpad=0.4, borderpad=0.1)
        fig.savefig(out(f"legend_{name}"), bbox_inches="tight", pad_inches=0.02)
        plt.close(fig)


# ---------------------------------------------------------------------------- table
TABLE_METHODS = (("rct_only", "RCT only"), ("omega_only", "PPI++"), ("fusion", "Fusion"),
                 ("shrinkage", "Shrinkage"), ("pretest", "Pretest"), ("naive_pool", "Naive pool"))


def table_rows(values, columns):
    """Per design, trial size and ATE method: bias, SD, RMSE, and for the methods with a Wald interval the mean
    standard error, the coverage and the mean width of the 95% interval."""
    rows = []
    for source, law, title, tag in columns:
        for n, cell in sorted(values.get((source, law, N_MAIN), {}).items()):
            for method, name in TABLE_METHODS:
                if ("ATE", "", method) not in cell:
                    continue
                error = cell[("ATE", "", method)]
                row = {"design": tag, "n": n, "method": name, "bias": error.mean(), "sd": error.std(ddof=1),
                       "rmse": np.sqrt(np.mean(error ** 2))}
                if ("SE", "", method) in cell:
                    se = cell[("SE", "", method)]
                    row.update(se=se.mean(), coverage=np.mean(np.abs(error) <= 1.96 * se), width=np.mean(2 * 1.96 * se))
                rows.append(row)
    return rows


def write_table(values, columns, folder):
    rows = table_rows(values, columns)
    fields = ["design", "n", "method", "bias", "sd", "rmse", "se", "coverage", "width"]
    with open(folder / "table_ate.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, restval="")
        writer.writeheader()
        writer.writerows(rows)
    titles = {tag: title for _, _, title, tag in columns}
    number = lambda row, key, digits=3: f"{row[key]:.{digits}f}" if key in row else "--"
    lines = [r"\begin{tabular}{llrrrrrrr}", r"\toprule",
             r"Design & $n$ & Method & Bias & SD & RMSE & SE & Cov. & Width \\", r"\midrule"]
    previous = None
    for row in (r for r in rows if r["n"] in (300, 900, 3000)):
        block = (row["design"], row["n"])
        if previous is not None and block != previous:
            lines.append(r"\addlinespace" if row["design"] == previous[0] else r"\midrule")
        design = titles[row["design"]] if previous is None or row["design"] != previous[0] else ""
        size = str(row["n"]) if block != previous else ""
        lines.append(" & ".join([design, size, row["method"], number(row, "bias"), number(row, "sd"),
                                 number(row, "rmse"), number(row, "se"), number(row, "coverage", 2),
                                 number(row, "width")]) + r" \\")
        previous = block
    lines += [r"\bottomrule", r"\end{tabular}"]
    (folder / "table_ate.tex").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", default="results")
    parser.add_argument("--law", default=LAW, choices=("weak_shift", "strong_shift"),
                        help="the OBS covariate law of the shifted panels and of the sweeps")
    parser.add_argument("--out", default="figures")
    parser.add_argument("--format", default="pdf", choices=("pdf", "png"))
    parser.add_argument("--zoom", action="store_true", help="ratio panels without the naive pool")
    parser.add_argument("--failure", action="store_true", help="also draw the designs with a useless OBS")
    parser.add_argument("--baseline", default="cf", choices=("cf", "full"),
                        help="trial-only reference: cf = the cross-fitted trial-only estimators (nuisance on one "
                             "block), full = the practitioner's ones (nuisance on two blocks, no tuning sample)")
    a = parser.parse_args()
    if a.baseline == "full":                         # same curves, labels and legend; other reference rows
        RCT = RCT[:1] + ("rct_only_full",) + RCT[2:]
        RCT_DR = RCT_DR[:1] + ("rct_full",) + RCT_DR[2:]
        CATE_CURVES = (CATE_CURVES[0][:1] + ("rct_full",) + CATE_CURVES[0][2:],) + CATE_CURVES[1:]
    main = tuple((source, a.law if law == LAW else law, title, tag) for source, law, title, tag in MAIN)
    values = load(a.results)
    if not values:
        raise SystemExit(f"no results in {a.results}")
    folder = Path(a.out)
    folder.mkdir(parents=True, exist_ok=True)
    out = lambda name: folder / f"{name}.{a.format}"
    plt.rcParams.update(STYLE)
    drawn = []
    for absolute in (False, True):
        suffix = "_abs" if absolute else ""
        name = lambda stem, suffix=suffix: out(stem + suffix)
        draw_ate(values, main, name, absolute, a.zoom)
        draw_cate(values, main, name, absolute, a.zoom)
        drawn += [f"ate{suffix}", f"cate{suffix}"]
        if draw_sweep(values, a.law, name, absolute):
            drawn.append(f"sweep{suffix}")
        if draw_obs(values, a.law, name, absolute, a.zoom):
            drawn.append(f"obs{suffix}")
        if a.failure and draw_failure(values, FAILURE, name, absolute):
            drawn.append(f"failure{suffix}")
    draw_legends(out, a.zoom)
    write_table(values, main, folder)
    print(f"written to {folder}/ with the law {a.law}: " + ", ".join(drawn) + ", legends, table_ate")
