# Manuscript-minimal CATE experiment: result

Plan `manuscript-minimal-v1`, run 2026-09-13.  200 replications per cell, 400
cell records, 33,000 raw rows, 0 failures.  Verification 16 of 16, run checks
6 of 6.

## What was asked

Implement the manuscript's doubly robust fusion and R-Fusion learners exactly,
then compare risk, squared bias and variance for the same trial target under two
covariate conditions.  The trial population and the true conditional effect are
identical in both cells; only the observational covariate law moves.

## Risk relative to the trial-only learner

Each number is a mean risk over one fixed grid of 50,000 trial-population points,
divided by the same learner's trial-only risk on the same replications.  Below
one means the fusion rule did better.  Intervals come from a paired bootstrap of
2,000 resamples that moves whole seed bundles at once.

| condition | learner | lambda only | omega only | joint |
| --- | --- | ---: | ---: | ---: |
| shared, r0 = 1 exactly | DRF | 0.792 | 0.936 | **0.754** [0.724, 0.785] |
| shared, r0 = 1 exactly | RF | 0.821 | 0.897 | **0.749** [0.717, 0.779] |
| shifted, true ratio | DRF | 0.795 | 0.981 | 0.800 [0.762, 0.839] |
| shifted, true ratio | RF | 0.824 | 0.948 | 0.803 [0.767, 0.841] |
| shifted, estimated ratio | DRF | 0.795 | 0.997 | 0.806 [0.770, 0.845] |
| shifted, estimated ratio | RF | 0.824 | 0.967 | 0.801 [0.768, 0.835] |
| shifted, ratio ignored | DRF | 0.795 | 1.820 | **1.737** [1.493, 2.020] |
| shifted, ratio ignored | RF | 0.824 | 1.874 | **1.693** [1.461, 1.968] |

## Four readings

**Under a common covariate law both coefficients earn their place.**  Each
channel alone lowers the risk and the two together lower it more than either.
The paired increment from adding omega on top of the selected lambda is -0.027
[-0.039, -0.015] for the doubly robust learner and -0.051 [-0.067, -0.037] for
the R-Fusion learner.  Both intervals exclude zero.

**Under a covariate shift, the incremental gain from adding omega on top of
lambda is not detected.**  With the true ratio, the paired increment is +0.004
[-0.016, +0.026] for the doubly robust learner and -0.015 [-0.034, +0.009] for
the R-Fusion learner.  Neither interval excludes zero.  This result does not
show that covariate shift disables the omega channel.  For example, the
omega-only R-Fusion risk ratio is 0.948.  It shows only that this design does
not establish an additional average improvement from omega after lambda has
already been selected.  The lambda-only ratios are 0.795 and 0.824; lambda does
not multiply the transport weight.

**No ratio-estimation cost is detected in this well-specified Gaussian shift.**
The classifier route gives 0.806 and 0.801 against the true ratio's 0.800 and
0.803.  The classifier converged in every replication, and the effective sample
size fraction it produces is 0.379 against the true ratio's 0.377, with the
analytic value for a unit Mahalanobis shift being 1/e = 0.368.  Equal covariance
Gaussian mean shift has a linear log density ratio, which matches the logistic
ratio model used here.  This experiment therefore does not show that estimating
the ratio is generally cost-free.

**Ignoring the shift is worse than not fusing at all.**  Setting the ratio to one
when it is not one raises the risk to 1.74 and 1.69 times the trial-only risk.
This is the only condition in which the method is actively harmful, and the
decomposition says why: the squared bias share rises from about 0.06 of the
trial-only risk in every other panel to 0.28 here, while variance rises too.
The failure distribution is heavy-tailed.  For the joint rule, the doubly robust
learner has median risk 0.770, 90th percentile 2.389, 99th percentile 6.765 and
maximum 11.191.  The corresponding R-Fusion values are 0.736, 2.312, 5.863 and
9.698.  The joint rule is worse than its paired trial-only fit in 51% and 52% of
replications, respectively.  Thus, the elevated mean is produced by both common
negative transfer and a small number of very large failures.

## The diagnostic the plan asked for

At the fixed candidate (0.5, 0.5) the selection score difference against the
trial-only candidate minus the population risk difference cancels the common
noise floor.  With an exact ratio the residual should carry no extra term; with a
wrong ratio the manuscript's `omega E_O[(r - r0)(ghat - zeta)^2]` remains.  The
measured residuals are -0.0259 for the shared cell, -0.0262 for the shifted cell
with the true ratio, -0.0288 with the estimated ratio and -0.1425 with the ratio
ignored.  The first three agree to within 0.003; the fourth is about five times
larger.  The common negative residual can include selection-sample variation as
well as fixed-grid integration error, so this diagnostic does not identify its
source.  It is a diagnostic on a fixed 50,000-point grid, not a proof that the
expected residual is zero.

## What was checked before the numbers were read

| check | result |
| --- | --- |
| pseudo-outcome against the written-out score | 0 to machine precision |
| both losses against their direct derivative | largest relative gradient 4.3e-15 |
| the two learners coincide at a trial propensity of one half | 4.5e-15 |
| dropping the R-Fusion weight is detected at 0.35 | detected |
| removing the observational correction from the score is detected | detected |
| using a separate validation omega is detected | detected |
| inserting the R-Fusion weight into the score is detected | detected |
| the exact ratio against the Gaussian log-density difference | 1e-10 |
| a zero shift gives a ratio of exactly one | exact |
| the trial arrays are identical across the two cells | exact |
| ties break lexicographically | verified both ways |
| the selection score takes no truth argument | verified by signature |
| the decomposition identity | 7.9e-16 |
| unpaired resampling is detected | detected |
| scheduled equals recorded | 200 per group, 0 failures |
| omega zero agrees across routes inside a cell | exactly 0 |

## Settings

Trial propensity 0.35, chosen because at one half the two learners coincide
observation by observation.  Trial sizes 100 for each of nuisance, tuning and
selection; observational sizes 3,000, 1,000 and 1,000.  One cubic spline basis
with three knots, fixed ridge 0.01 and no jitter.  Nine
candidates on the product grid of lambda and omega in {0, 0.5, 1}.  Shift set so
its Mahalanobis length is exactly one.  Truth grid of 50,000 points generated
once from its own seed.

There is no clipping of fitted outcome predictions or density ratios.  The
spline basis evaluation does clip each covariate to the boundary knots learned
from the trial nuisance sample before evaluating the B-splines.  The recorded
main-run time is 26.14 seconds.

## What this does not show

One structural model, one budget, one basis.  The exact ratio restores the target
of the expected loss; it does not promise the finite-sample risk of no shift, a
free ratio estimate, or an improvement in every replication.  Coverage is not
reported because no interval estimator was run.  The bootstrap interval is a
Monte Carlo interval conditional on the fixed truth grid and carries no
simultaneous guarantee across the panels.

## Preservation

The original raw rows, sufficient statistics, run configuration, run checks,
failure log, environment, run-time record, run code manifest, numerical
analysis, table and figures are preserved byte-for-byte.  The post-run analysis
code and its inputs and outputs are recorded in `analysis_provenance.json`.
`tail_diagnostics.csv` reports the joint-rule tail summaries for all eight
cell-by-learner panels.  `SHA256SUMS` covers every regular file in this stage
except itself.  These additions do not rerun the simulation or alter the
manuscript, mirror or Git history.
