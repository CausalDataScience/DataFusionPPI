# DataFusionPPI ATE implementation and results

This note is the ATE-only publication record for the `rev4-2026-09-12` results.
The statistical specification remains `2026-09-09-implementation-handoff-JA.md`.
The independent replication counts below count generated data sets, not estimator,
ratio-route, metric, or long-format rows.

## Reproducible surface

- `exp_ate1.py` generates ATE-1. It uses 44 design cells and 100 independent
  replications per cell, for 4,400 generated data sets.
- `exp_ate2.py` generates ATE-2. The synthetic panel has 1,200 generated data
  sets, paired across five ratio routes. The NSW panel has 600 generated data
  sets, paired across four feasible routes. The CSV contains 42,000 estimator
  rows because route and estimator outputs repeat each generated data set.
- `exp_conditional_variance.py` freezes 15 fitted nuisance/coefficient states and
  generates 400 independent evaluation draws per state, for 6,000 evaluation
  draws and 42,000 estimator rows.
- `exp_ate_star_real.py` reproduces the ATE projection of the STAR real-outcome
  stress test: 8 cells and 100 independent replications per cell, for 800 data
  sets and 8,000 estimator rows.
- `analysis_variance.py` computes cellwise ATE bias, empirical variance, mean
  squared error (MSE), root MSE (RMSE), interval coverage, and paired-bootstrap
  variance reductions. `make_ate_figures.py` reads only frozen ATE CSVs.

## Aggregation estimands

ATE-1 tables give each of the 44 prespecified design cells equal weight. A pooled
variance-reduction rate is the mean of the 44 cell-specific rates, recomputed in
each paired bootstrap draw. It is not the ratio of two variances after pooling all
replications. RMSE ratios are computed within a cell and then averaged. Coverage
is the equal-cell mean of the empirical 95 percent Wald-interval coverage.

ATE-2 synthetic bias in Figure 4 is the mean across six cells of the absolute
cell bias, where a cell is one family and one trial size at a fixed shift and
ratio route. Its interval resamples replication indices inside every cell, keeps
all five routes paired, and recomputes the full statistic. It is not the absolute
value of a pooled signed bias.

The conditional-variance experiment conditions on each frozen fitted state.
Its target is the exact evaluation-sample variance in Theorem 1. The STAR ATE
stress test uses the law-weighted within-pattern contrast as its reference.

## Main results

For ATE-1, the mean cellwise RMSE ratio to trial-only is 0.748 for the joint
estimator, 0.759 for lambda-only, 0.994 for omega-only, and 0.707 for the oracle
joint coefficients. With lambda held fixed, estimated omega reduces whole-
procedure variance by 2.34 percent, with paired-bootstrap 95 percent interval
[1.04, 3.19]. Conditional on frozen fitted objects, the corresponding reduction
is 4.11 percent [3.06, 5.02]. The joint estimator has pooled coverage 0.9423,
against 0.9416 for trial-only.

In synthetic ATE-2, the joint estimator's mean absolute cell bias under mild
shift is 0.0290 with the oracle ratio, 0.0295 with the classifier, 0.0288 with
BAL-X, 0.0302 with BAL-X+g, and 0.0626 when the ratio is set to one. Under
strong shift, the corresponding values are 0.0261, 0.0305, 0.0283, 0.0324,
and 0.0945. Figure 4 reports the same statistic with paired within-cell
bootstrap intervals. The result isolates transport: ignoring the ratio adds
substantial bias, while the three estimated weighted routes remain close to
the oracle route.

Across the 105 frozen-fit and coefficient-rule comparisons, every empirical
variance lies within three Monte Carlo standard errors of the Theorem 1 value;
the largest absolute standardized gap is 1.679. The 400-draw expansion is an
exploratory variance diagnostic.

In the STAR real-outcome ATE stress test, the joint estimator has equal-cell mean
bias -0.00006, RMSE 0.07862, and coverage 0.9613. Trial-only has bias 0.00008,
RMSE 0.07918, and coverage 0.9588. This is a stress-test comparison against the
law-weighted reference, not identification of an individual-level conditional
effect.

## Limits and provenance

The ATE-1 Wald intervals are asymptotic. The cross-fitted joint row is a
reference-only panel because overlapping folds require a separate variance
argument. NSW results compare routes against an experimental benchmark and do
not supply a known population truth. Oracle coefficients and the conditional-
variance expansion are simulation diagnostics.

The run used base seed 190602 with deterministic derived seeds. The environment
is recorded in `requirements-ate.txt`. Source data are hash-checked by
`fusion_data.py`. `ate-artifact-manifest.sha256` binds the published code,
results, figures, report, and retained historical ATE pilot/benchmark surface.
The manifest intentionally does not hash itself.

## Commands

```bash
python3 -m pytest -q -p no:cacheprovider code/test_ate_acceptance.py
python3 code/exp_ate1.py --replications 100 --stem ate1_when_fusion_helps
python3 code/exp_ate2.py --replications 100 --stem ate2_shift
python3 code/exp_conditional_variance.py --training-draws 5 --evaluations 400 --stem conditional_variance_400
python3 code/exp_ate_star_real.py --replications 100 --stem cate2_star_real_ate
python3 code/analysis_variance.py --stem ate1_when_fusion_helps
python3 code/analysis_variance.py --stem conditional_variance_400 --group family,training_draw
python3 code/make_ate_figures.py
```
