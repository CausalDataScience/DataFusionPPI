# DataFusionPPI

This repository is a working research package for **prediction-powered data
fusion for treatment-effect estimation**. The randomized controlled trial (RCT)
remains the source of causal identification. Observational-study (OBS) outcome
models and covariates enter through two separately interpretable coefficients:

- \(\lambda\): outcome-model fusion inside an RCT-valid AIPW pseudo-outcome;
- \(\omega\): a prediction-powered covariate correction centered across the
  RCT and OBS samples.

The manuscript derives exact fixed-coefficient ATE results and a
prediction-powered squared-loss construction for CATE. It does **not** claim that
estimated coefficients always improve finite-sample performance. The completed
experiments concern ATE estimation; CATE and real-outcome studies are planned.

## Repository structure

- `manuscript/`: ITBound v2-based LaTeX manuscript and bibliography.
- `memo/PPI-Fusion.md`: canonical detailed theory memo.
- `code/ssem_ate_pilot.py`: initial SSEM ATE pilot.
- `code/total_budget_nested_benchmark.py`: nested four-method benchmark driver.
- `materials/`: protocols, replication-level results, summaries, and plots.

## Compile the manuscript

From the repository root:

```bash
cd manuscript
latexmk -pdf main.tex
```

The resulting PDF is `manuscript/main.pdf`. Generated LaTeX auxiliaries are
ignored, while the compiled manuscript PDF is intentionally trackable.

## Reproduce the existing ATE benchmarks

The benchmark script refuses to overwrite existing outputs unless explicitly
asked. Use a new output stem for a fresh run.

```bash
python3 code/total_budget_nested_benchmark.py \
  --replications 20 \
  --n-rct-values 20,30 \
  --output-stem reproduced_fixed_budget

python3 code/total_budget_nested_benchmark.py \
  --replications 20 \
  --n-rct-values 20,30,50,100 \
  --obs-multiplier 100 \
  --output-stem reproduced_proportional_budget
```

The script writes outputs under `materials/`. The fixed-budget default uses
\(N_O=5000\); the proportional design uses the exact mapping
\(N_O=100n_R\).

## Estimator labels

User-facing labels and retained implementation codes are:

| User-facing method | Internal code |
|---|---|
| RCT-only AIPW | `rct_aipw` |
| \(\lambda\)-only outcome-model fusion | `haipw_only` |
| \(\omega\)-only PPI correction | `ppi_power_tuned_only` |
| Joint \((\lambda,\omega)\) fusion | `full_datafusionppi` |

`haipw_only` is a stale code label and does not associate this estimator with
the separate HAIPW literature.

## Current empirical evidence

Each reported cell has 20 Monte Carlo repetitions.

| Benchmark | Cells | Best RMSE: joint / \(\lambda\)-only / \(\omega\)-only / AIPW | Joint vs AIPW |
|---|---:|---:|---|
| Fixed OBS budget, \(N_O=5000\) | 16 | 11 / 5 / 0 / 0 | Lower RMSE in 16/16 cells; mean RMSE ratio 0.648 |
| Proportional, \(N_O=100n_R\) | 32 | 18 / 10 / 1 / 3 | Lower RMSE in 27/32 cells; mean RMSE ratio 0.760 |

These results support joint tuning in many designs, but not a universal
finite-sample dominance claim. The \(\lambda\)-only channel is the more stable
single-channel estimator in the current runs.

The STAR experiment in this repository is **real-\(X\), semi-synthetic**:
empirical STAR covariates are reused, while treatment, confounding, potential
outcomes, and effect heterogeneity are generated. It is not an analysis of the
original STAR outcomes. No WHI participant-level data are included.

## Status

This is a research draft. Exact density-ratio transport results distinguish the
common covariate marginal \(P_R^X=P_O^X\) from covariate shift
\(P_R^X\neq P_O^X\). Estimated-ratio CATE inference, adaptive balancing,
end-to-end neural theory, the proposed R-loss extension, and real-outcome
validation remain open work.
