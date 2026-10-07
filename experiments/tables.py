"""Writes the LaTeX tables of Appendix C that visualization.py does not write, from the CSV files of simulation.py.

    python3 tables.py --results results --misspec results/misspec --hethi results/hethi --out tables

    table_ate_ratio.tex       ATE mean squared error relative to the cross-fitted trial-only AIPW
    table_cate_ratio.tex      CATE risk relative to the cross-fitted trial-only DR-learner
    table_ate_extra.tex       lambda only, PPI++, fusion, and fusion with the rule omega = D / B
    table_crossfit.tex        fusion relative to the practitioner's trial-only estimators (nuisance on 2/3)
    table_cate_cf.tex         DRF and RF relative to both trial-only DR-learners
    table_cate_samesieve.tex  the trial loss on the sieve with g, and DRF and RF relative to it
    table_misspec.tex         the misspecified shift (results of --shifts misspec_shift in the folder --misspec)
    table_hethi.tex           the design with a large effect heterogeneity (results of scm1_hethi in --hethi)

A ratio is a ratio of means over the replications of a cell: of squared ATE errors, or of CATE risks.  The
replications of a cell are read as in visualization.py (a failed replication is left out of every method).
The table of bias, SD, RMSE, SE, coverage and width of the ATE is table_ate.tex of visualization.py.
A table whose results are missing is skipped.
"""
import argparse
from pathlib import Path

import numpy as np

from visualization import N_MAIN, load

SIZES = (300, 900, 3000)                          # the trial sizes n = 3m shown in the tables
COLUMNS = tuple((source, law) for source in ("scm1", "ihdp", "acic") for law in ("same", "weak_shift", "strong_shift"))
PPI = r"PPI\textsuperscript{++}"
HEADER = [r"\begin{tabular}{@{}llccccccccc@{}}", r"\toprule",
          r"& & \multicolumn{3}{c}{Synthetic} & \multicolumn{3}{c}{IHDP} & \multicolumn{3}{c}{ACIC~2016}\\",
          r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}\cmidrule(lr){9-11}",
          r"$n$ & Method & $r_0{=}1$ & weak & strong & $r_0{=}1$ & weak & strong & $r_0{=}1$ & weak & strong\\",
          r"\midrule"]
FOOTER = [r"\bottomrule", r"\end{tabular}", r"\end{center}", r"\end{table}"]


def ate(method):
    return ("ATE", "", method)


def cate(learner, method):
    return ("CATE", learner, method)


# Ratio tables: (file, caption, label, rows), each row (name, numerator, denominator)
RATIO_TABLES = (
    ("table_ate_ratio", "ATE mean squared error relative to trial-only AIPW (cross-fitted), by trial size $n$ and "
     "covariate law of the OBS.", "tab:app-ate-ratio",
     ((PPI, ate("omega_only"), ate("rct_only")), ("Fusion", ate("fusion"), ate("rct_only")),
      ("Shrinkage", ate("shrinkage"), ate("rct_only")), ("Pretest", ate("pretest"), ate("rct_only")),
      ("Naive pool", ate("naive_pool"), ate("rct_only")))),
    ("table_cate_ratio", "CATE risk relative to the trial-only DR-learner (cross-fitted); otherwise as "
     r"Table~\ref{tab:app-ate-ratio}. RCT (R) is the trial-only R-learner.", "tab:app-cate-ratio",
     tuple((name, cate(*key), cate("DRF", "rct_only")) for name, key in (
         ("DRF", ("DRF", "fusion")), ("RF", ("RF", "fusion")), (PPI + " DRF", ("DRF", "omega_only")),
         (PPI + " RF", ("RF", "omega_only")), ("2-step DR", ("DRF", "two_step")), ("2-step R", ("RF", "two_step")),
         ("IR", ("RF", "integrative")), ("RCT (R)", ("RF", "rct_only")), ("Naive pool", ("pool", "naive_pool"))))),
    ("table_ate_extra", r"ATE mean squared error relative to trial-only AIPW: $\lambda$ only, PPI\textsuperscript{++}, "
     r"fusion ($\hat\omega=\hat D/(\hat B+\hat B_{\mathrm{cal}})$), and fusion with the rule $\hat\omega=\hat D/\hat B$.",
     "tab:ate-extra",
     ((r"$\lambda$ only", ate("lambda_only"), ate("rct_only")), (PPI, ate("omega_only"), ate("rct_only")),
      ("Fusion", ate("fusion"), ate("rct_only")),
      (r"Fusion, $\omega=\hat D/\hat B$", ate("fusion_nocal"), ate("rct_only")))),
    ("table_crossfit", "Cross-fitted fusion relative to the practitioner's trial-only estimators (3-fold, nuisance on "
     "two thirds).", "tab:crossfit",
     (("Fusion / trial-only (2/3)", ate("fusion"), ate("rct_only_full")),
      ("Trial-only (1/3) / (2/3)", ate("rct_only"), ate("rct_only_full")),
      ("DRF / DR (2/3)", cate("DRF", "fusion"), cate("DRF", "rct_full")),
      ("RF / DR (2/3)", cate("RF", "fusion"), cate("DRF", "rct_full")),
      ("RF / R (2/3)", cate("RF", "fusion"), cate("RF", "rct_full")),
      ("DR (1/3) / DR (2/3)", cate("DRF", "rct_only"), cate("DRF", "rct_full")))),
    ("table_cate_cf", "Cross-fitted DRF and RF relative to the two trial-only DR-learners.", "tab:cate-cf",
     (("DRF / DR (1/3)", cate("DRF", "fusion"), cate("DRF", "rct_only")),
      ("RF / DR (1/3)", cate("RF", "fusion"), cate("DRF", "rct_only")),
      ("DRF / DR (2/3)", cate("DRF", "fusion"), cate("DRF", "rct_full")),
      ("RF / DR (2/3)", cate("RF", "fusion"), cate("DRF", "rct_full")))),
    ("table_cate_samesieve", r"The trial loss on the sieve with $\hat g$ relative to the trial-only DR-learner, and DRF "
     "and RF relative to it.", "tab:samesieve",
     ((r"Trial loss, sieve with $\hat g$ (DR)", cate("DRF", "trial_sieve_g"), cate("DRF", "rct_only")),
      (r"Trial loss, sieve with $\hat g$ (R)", cate("RF", "trial_sieve_g"), cate("DRF", "rct_only")),
      ("DRF / same sieve (DR)", cate("DRF", "fusion"), cate("DRF", "trial_sieve_g")),
      ("RF / same sieve (R)", cate("RF", "fusion"), cate("RF", "trial_sieve_g")))),
)

# Error tables: (file, caption, blocks (label, source, law), methods (name, method))
MISSPEC = ("table_misspec", r"Misspecified shift (tilt $\exp\{a(h^2-1)/\sqrt2\}$, effective sample fraction 0.40), "
           r"cross-fitted: AIPWF with the calibrated ratio $\hat r_{\mathrm{BAL}}$ (Fusion) and with the uncalibrated "
           r"logistic ratio $\hat r_{\mathrm{CLS}}$. MSE ratio: relative to RCT only.",
           (("Synthetic", "scm1", "misspec_shift"), ("IHDP", "ihdp", "misspec_shift")),
           (("RCT only", "rct_only"), ("Fusion", "fusion"),
            (r"Fusion, uncalibrated $\hat r_{\mathrm{CLS}}$", "fusion_cls")))
HETHI = ("table_hethi", r"Synthetic design with large effect heterogeneity (outcome noise 0.2, $\mathrm{Var}_R\tau_0=5.88$, "
         r"confounding bias $0.75\,\mathrm{sd}_R\tau_0$), cross-fitted. MSE ratio: relative to RCT only.",
         ((r"$r_0{=}1$", "scm1_hethi", "same"), ("strong", "scm1_hethi", "strong_shift")),
         (("RCT only", "rct_only"), (PPI, "omega_only"), (r"$\lambda$ only", "lambda_only"), ("Fusion", "fusion"),
          (r"Fusion, $\hat\omega=\hat D/\hat B$", "fusion_nocal")))


def size_label(n):
    """3000 -> 3{,}000."""
    return f"{n // 1000}{{,}}{n % 1000:03d}" if n >= 1000 else str(n)


def mean_loss(cell, key):
    """The mean squared ATE error, or the mean CATE risk, of one method over the replications of a cell."""
    return float(np.mean(cell[key] ** 2)) if key[0] == "ATE" else float(np.mean(cell[key]))


def ratio_table(values, caption, label, rows):
    lines = [r"\begin{table}[h]", rf"\caption{{{caption}}}", rf"\label{{{label}}}", r"\begin{center}\footnotesize"]
    lines += HEADER
    for i, n in enumerate(SIZES):
        if i:
            lines.append(r"\addlinespace")
        cells = [values[(source, law, N_MAIN)][n] for source, law in COLUMNS]
        for j, (name, numerator, denominator) in enumerate(rows):
            ratios = [mean_loss(cell, numerator) / mean_loss(cell, denominator) for cell in cells]
            lines.append(f"{size_label(n) if j == 0 else ''} & {name} & "
                         + " & ".join(f"{ratio:.2f}" for ratio in ratios) + r" \\")
    return lines + FOOTER


def error_table(values, caption, blocks, methods):
    """Bias, SD, RMSE, coverage of the 95% Wald interval, and the MSE relative to RCT only, of the ATE errors."""
    lines = [r"\begin{table}[h]", rf"\caption{{{caption}}}", r"\begin{center}\footnotesize",
             r"\begin{tabular}{@{}lllrrrrr@{}}", r"\toprule",
             r"Design & $n$ & Method & Bias & SD & RMSE & Cov. & MSE ratio \\", r"\midrule"]
    for b, (label, source, law) in enumerate(blocks):
        if b:
            lines.append(r"\midrule")
        for i, n in enumerate(SIZES):
            if i:
                lines.append(r"\addlinespace")
            cell = values[(source, law, N_MAIN)][n]
            base = mean_loss(cell, ate("rct_only"))
            for j, (name, method) in enumerate(methods):
                error, se = cell[ate(method)], cell[("SE", "", method)]
                fields = [label if i == 0 and j == 0 else "", size_label(n) if j == 0 else "", name,
                          f"{error.mean():.3f}", f"{error.std(ddof=1):.3f}", f"{np.sqrt(np.mean(error ** 2)):.3f}",
                          f"{np.mean(np.abs(error) <= 1.96 * se):.2f}", f"{np.mean(error ** 2) / base:.2f}"]
                lines.append(" & ".join(fields) + r" \\")
    return lines + FOOTER


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", default="results", help="results of the main runs (scm1, ihdp, acic)")
    parser.add_argument("--misspec", default="results/misspec", help="results of the misspecified shift")
    parser.add_argument("--hethi", default="results/hethi", help="results of scm1_hethi")
    parser.add_argument("--out", default="tables")
    a = parser.parse_args()
    folder = Path(a.out)
    folder.mkdir(parents=True, exist_ok=True)
    main, misspec, hethi = load(a.results), load(a.misspec), load(a.hethi)
    jobs = [(name, lambda caption=caption, label=label, rows=rows: ratio_table(main, caption, label, rows))
            for name, caption, label, rows in RATIO_TABLES]
    jobs += [(name, lambda values=values, caption=caption, blocks=blocks, methods=methods:
              error_table(values, caption, blocks, methods))
             for values, (name, caption, blocks, methods) in ((misspec, MISSPEC), (hethi, HETHI))]
    for name, make in jobs:
        try:
            lines = make()
        except KeyError as missing:
            print(f"{name}: skipped, results missing ({missing})")
            continue
        (folder / f"{name}.tex").write_text("\n".join(lines) + "\n")
        print(f"{name}: written to {folder / (name + '.tex')}")
