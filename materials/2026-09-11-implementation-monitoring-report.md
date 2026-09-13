# Implementation Monitoring Report: DataFusionPPI Experiments

For the monitoring agent.  Revision 2, 2026-09-11, after the feedback that
$\omega$ must be judged on variance directly.  Every number came out of the code
in `code/` and is reproducible from the seed stored in each replication file.
Nothing has been committed or pushed.

## 1. Status

| gate | content | status |
| --- | --- | --- |
| G0 code | shared library, data loaders, baselines, five drivers, figures | done |
| G1 validation | eight acceptance tests, two algebraic identities | done, all pass |
| G2 studies | ATE-1, ATE-2, CATE-1, CATE-2, conditional-variance | done |
| G3 write-up | four variance figures, two supporting figures, three table files | done |
| G4 publication | mirror sync, commit, push | not started, awaiting instruction |

## 2. What the feedback changed

The earlier revision judged $\omega$ on root-mean-square error, which is the
wrong scale.  The coefficient is defined to minimize variance, and with an
unbiased estimator a 20 percent variance reduction is only a 10.6 percent
reduction in root-mean-square error, so a ratio near 0.99 can hide a real
variance effect.  This revision adds four things.

1. A **conditional-variance experiment** that freezes the fitted objects and the
   coefficients and regenerates only the evaluation data, which is the state the
   exact variance of Theorem 1 describes.
2. **Paired same-lambda comparisons** for $\omega$, with a paired bootstrap that
   resamples replication indices once and reuses them for both arms.
3. A **full reporting table** per cell and estimator with the mean estimate,
   bias, empirical variance, empirical standard deviation, mean estimated
   variance, mean squared error, root-mean-square error, coverage, interval
   length, effective replications, failure rate, and the Monte Carlo uncertainty
   of each.
4. A **variance and squared-bias split of the conditional risk** in the CATE
   studies, at a fixed anchor of covariates shared by every replication.

Two algebraic identities are checked automatically and hold to machine
precision.  For the average effect,
$\widehat{\mathrm{MSE}}=\frac{K-1}{K}S^2+\widehat{\mathrm{Bias}}^2$, largest gap
$3.6\times10^{-15}$.  For the conditional effect, risk on the anchor equals
$\frac{K-1}{K}$ times the prediction variance plus the squared bias, largest gap
$7.1\times10^{-14}$.  Both confirm that variance, bias, and error were
aggregated over the same replications with the same filter.

## 3. The omega question, answered on variance

### 3.1 Does the exact variance formula hold?

The conditional experiment fixes the nuisances, the transport weight and the
coefficients once per training draw, then draws only fresh evaluation data.
With five training draws for each of three families, seven coefficient rules,
and 400 evaluation replications, the replicated spread matches Theorem 1 in
**105 of 105 frozen fits** within three Monte Carlo standard errors, with mean
$z$ between $-0.16$ and $-0.04$.  Examples: the joint rule has empirical
variance 0.1883 against a formula value of 0.1880, and the trial-only rule
0.2957 against 0.2956.

At 100 replications the same check passes in 100 of 105 fits.  The five misses
sit in one training draw whose pseudo-outcome has excess kurtosis 5.9 and a
maximum absolute value of 46.6 against a standard deviation of 2.8.  A direct
follow-up on that draw shows the gap closing as replications grow: $z=-0.58$ at
100 replications, $+0.03$ at 400, $+0.22$ at 1600, $-0.61$ at 6400.  The
expansion to 400 was decided after seeing the 100-replication result and is
labelled as such.  Sample variances of heavy-tailed means are noisy, so the
Monte Carlo standard error of every variance in this report is a bootstrap
standard error rather than the normal-theory $S^2\sqrt{2/(K-1)}$.

### 3.2 How much variance does omega remove?

Every comparison holds $\lambda$ fixed across the two arms inside the same
replication.  Rates are relative to the variance with $\omega$ switched off.

| scope | comparison | mean rate | cells whose interval is above zero | oracle mean rate |
| --- | --- | ---: | ---: | ---: |
| conditional | omega alone | 2.7% | 60% of 15 | 3.8% |
| conditional | omega given lambda | 4.1% | 53% of 15 | 5.9% |
| conditional | whole method | 35.3% | 100% of 15 | |
| whole procedure | omega alone | 1.2% | 4.5% of 44 | 2.3% |
| whole procedure | omega given lambda | 2.3% | 6.8% of 44 | 3.2% |
| whole procedure | whole method | 38.0% | 84.1% of 44 | |

Conditionally the reduction is real: the paired interval sits above zero in more
than half the frozen fits, far above the 5 percent a 95 percent interval
produces by chance.  Once the whole procedure is refitted in every replication
the signal disappears into the extra variability of fitting, with 4.5 and 6.8
percent of cells above zero, which is what a 95 percent interval gives when the
true effect is small relative to the noise.

The size is the story.  **The oracle coefficient reaches only 2.3 to 5.9
percent, so the reduction available in these designs is itself small.**  Where
the oracle gain is separated from zero, which happens in 53 to 67 percent of
frozen fits and in 7 to 14 percent of whole-procedure cells, the estimated
$\widehat\omega$ recovers 60 to 71 percent of it.  Elsewhere the recovery rate
is reported as not applicable, because dividing by a gain that is itself
indistinguishable from zero is unstable.  The answer to the diagnostic question
is therefore the first branch: there is little variance to remove here, not a
channel that fails to remove it.

For scale, 4 percent of variance is 2 percent of root-mean-square error, which
is exactly the 0.98 to 0.99 ratio the earlier table showed for the
$\omega$-only estimator.  The two scales agree once both are computed.

### 3.3 What does omega cost?

| scope | comparison | bias change | RMSE change | coverage change | interval length change |
| --- | --- | ---: | ---: | ---: | ---: |
| conditional | omega alone | -0.0009 | -0.0072 | 0.0000 | -0.028 |
| conditional | omega given lambda | -0.0014 | -0.0090 | +0.0007 | -0.035 |
| whole procedure | omega alone | +0.0013 | -0.0030 | +0.0007 | -0.014 |
| whole procedure | omega given lambda | -0.0004 | -0.0059 | -0.0016 | -0.018 |

Nothing is paid.  Absolute bias moves by at most 0.0014, coverage by at most
0.0016, the interval shortens, and the root-mean-square error improves slightly.

### 3.4 Judgment on the three axes

| axis | verdict |
| --- | --- |
| evidence of variance reduction | present conditionally, absent for the whole procedure at 100 replications |
| practically large | no: 1 to 4 percent, against an oracle ceiling of 2 to 6 percent |
| cost in accuracy or inference | none detected |

The honest summary is that $\omega$ works as the theory says and the available
reduction in these designs is small.  Whether 4 percent is worth reporting as a
contribution is a separate decision about practical importance, and the minimum
that counts should be fixed before any new design is run.

## 4. The lambda channel and the whole method

The variance reduction of the whole method is 35 percent conditionally and 38
percent for the whole procedure, with the interval above zero in 100 and 84
percent of cells.  Almost all of it comes from $\lambda$.  In root-mean-square
error, averaged over four families:

| estimator | 50 | 100 | 200 | 400 |
| --- | ---: | ---: | ---: | ---: |
| trial only | 1.000 | 1.000 | 1.000 | 1.000 |
| lambda only | 0.690 | 0.719 | 0.755 | 0.790 |
| omega only | 0.994 | 0.984 | 0.998 | 0.992 |
| joint, Algorithm 1 | 0.685 | 0.703 | 0.746 | 0.783 |
| joint at oracle coefficients | 0.664 | 0.673 | 0.733 | 0.725 |

Coverage stays at nominal for every fusion estimator, 0.905 to 0.962, and
collapses for every estimator that borrows the observational effect: shrinkage
0.285 to 0.585, adaptive 0.092 to 0.348, naive pooling 0.033 to 0.165,
observational AIPW 0.013 to 0.030.

Figure 3 separates their two failure modes.  The ratio of the mean estimated
variance to the empirical variance is 0.99 to 1.03 for every fusion estimator,
0.38 for shrinkage and 0.24 for adaptive, and 0.97 for the observational AIPW.
So shrinkage and adaptive fail because their variance estimate ignores the
selection they perform, while the observational AIPW estimates its variance
correctly and fails on bias.

## 5. The conditional effect, split into variance and bias

The conditional risk is measured on a fixed anchor of 500 covariates shared by
every replication, so it splits exactly into the variance of the fitted function
across replications and the squared bias of the mean fitted function.  Values
below are for the doubly robust learner, averaged over four families.

| sieve | rule | quantity | 100 | 200 | 400 |
| --- | --- | --- | ---: | ---: | ---: |
| spline, 3 knots | trial only | prediction variance | 19.47 | 30.79 | 29.87 |
| spline, 3 knots | trial only | squared bias | 0.23 | 0.59 | 0.18 |
| spline, 3 knots | selected | prediction variance | 2.86 | 1.81 | 1.29 |
| spline, 3 knots | selected | squared bias | 0.06 | 0.04 | 0.04 |
| neural | trial only | prediction variance | 16.11 | 10.78 | 13.94 |
| neural | trial only | squared bias | 0.20 | 0.09 | 0.11 |
| neural | selected | prediction variance | 1.98 | 0.97 | 0.68 |
| neural | selected | squared bias | 0.13 | 0.14 | 0.12 |

Prediction variance is 83 to 100 percent of the risk in every rule, so the
conditional gain is a variance gain.  Selection cuts the prediction variance by
66 to 96 percent.  The squared bias falls too on the spline sieve, by 48 to 93
percent, and moves little on the neural sieve, where it is already small.  So
the observational channel mostly stabilizes repeated fitting rather than moving
the average prediction toward the truth, and the paper should say that.

Selection improves on the trial-only learner in 95.2 percent of replications.
The median per-replication risk ratio is 0.115 to 0.202 for the neural sieve and
0.198 to 0.620 for the splines, shrinking as the trial grows.  The selected
learner sits 5 to 17 percent above the best candidate on the grid.

Theorem 6 holds in 100 percent of replications and is loose by four orders of
magnitude: median regret 0.07 to 0.12 against a median radius near 1,300.

## 6. CATE-2, real treatment and real outcomes

The construction induces a naive observational contrast of 0.58 at the weak
level and 1.04 at the strong one against an experimental 0.174, with stratum
shares matched to about $2\times10^{-3}$.

The conditional effect is unknown, so the report separates prediction
variability from the score, and makes no squared-bias claim.

| sieve | learner | prediction variance, trial only | prediction variance, selected | score gap | z |
| --- | --- | ---: | ---: | ---: | ---: |
| spline, 3 knots | DRF | 0.99 | 7.52 | +0.215 | 2.5 |
| spline, 3 knots | RF | 1.09 | 7.30 | +0.187 | 2.3 |
| neural | DRF | 25.35 | 19.54 | -1.067 | -8.0 |
| neural | RF | 25.83 | 19.23 | -1.354 | -10.7 |

On the neural basis selection both stabilizes the fitted function, a 23 to 26
percent cut in prediction variance, and improves the score.  On the coarse
spline basis over eight one-hot covariates the trial-only learner is nearly
constant, with a mean prediction range of 0.63, so it is stable and
uninformative; adding the channels multiplies its prediction variance sevenfold
without improving the score.  Stability and correctness are reported separately
here, and neither is claimed from the other.

For the average effect the fusion estimators match the trial-only estimator,
bias 0.023 against 0.024 and coverage 0.958 against 0.960, while naive pooling
carries bias 0.190 with coverage 0.379 and the observational AIPW bias 0.652
with coverage 0.110.

## 7. Specification changes forced by running

The first three change what the experiments mean and belong in the handoff.

1. **Propensity trimming applies to the trial only.**  Applied to the
   observational sample it distorts that sample's covariate law and breaks the
   transport identity $\mathbb E_O(r_0\widehat g)=\mathbb E_R(\widehat g)$ that
   the $\omega$ channel rests on.
2. **Every candidate, every fitted prediction, and the pseudo-outcome are
   clipped to a prespecified bound.**  Theorem 6 assumes bounded scores; without
   the bound an exploding candidate drove one STAR score to $-1.4\times10^{13}$
   against a typical 5 and captured the selector.  The bound is three times the
   largest absolute pseudo-outcome on the tuning sample.  After the fix the
   worst of 6,400 STAR rows is $-47.3$, which is the bound itself.
3. **The bound is read off the pseudo-outcome alone**, because including the
   fitted prediction lets a diverged nuisance fit inflate the bound.
4. The sieve basis is built from the nuisance sample, the shift levels are 0.5
   and 1.0 of the existing shift vector because larger levels break the overlap
   cap, the sieve solver carries a trace-scaled jitter, the STAR real-outcome
   study uses a linear outcome map on its one-hot covariates, and the Monte
   Carlo sizes are $2\times10^4$ for the oracle and $5\times10^4$ for the test
   sample.

One bug was found: in the exponential tilt of Definition 12 the overflow shift
was subtracted from the exponent but never added back to the objective value, so
the optimizer minimized a function whose gradient belonged to a different
problem.  The balance residual fell from $7\times10^{-4}$ to $3\times10^{-11}$.

## 8. Artifacts

Code: `fusion_core.py`, `fusion_data.py`, `fusion_baselines.py`, `fusion_ate.py`,
`fusion_cate.py`, `fusion_io.py`, `analysis_variance.py`, `exp_ate1.py`,
`exp_ate2.py`, `exp_cate1.py`, `exp_cate2.py`, `exp_conditional_variance.py`,
`make_figures.py`, `make_figures_variance.py`, `test_acceptance.py`,
`requirements.txt`.

Figures: `figure1_omega_variance_reduction.png`,
`figure2_variance_bias_rmse.png`, `figure3_variance_and_coverage.png`,
`figure4_cate_variance_bias.png`, with `figure1_ate_rmse.png` and
`figure2_cate_risk.png` kept as supporting views.

Tables: `judgment_axes.csv`, `paper_tables.csv`, `acceptance_tests.csv`, and per
study a replication file, a summary file, a reporting file and a variance file.

## 9. Compute

| study | cells | replications | wall clock |
| --- | ---: | ---: | ---: |
| ATE-1 | 44 | 4,400 | 260 s |
| ATE-2 | 18 | 1,800 | 194 s |
| conditional variance, 100 and 400 | 15 frozen fits | 3,500 and 14,000 | 26 s and 66 s |
| CATE-1, splines | 24 | 4,800 | 1,891 s |
| CATE-1, neural | 24 | 2,400 | 1,840 s |
| CATE-2 | 8 | 800 | 631 s |

## 10. What the monitoring agent should check next

1. Whether the paper reports the $\omega$ result as it stands: the channel works,
   the reduction is 1 to 4 percent, and the ceiling in these designs is 2 to 6
   percent.  A design with a stronger observational predictor would raise the
   ceiling, and that is a design choice to make before running, not after.
2. Whether the minimum reduction that counts as practically important is fixed
   before any new design, with the denominator stated as the $\omega$-off
   variance.
3. Whether items 1 to 3 of Section 7 are folded into the handoff and into the
   manuscript as stated conditions.
4. Whether the paper reports the two loose findings: the Theorem 6 radius is
   about $10^4$ times the realized regret, and the conditional risk is heavy
   tailed so medians carry the message.
5. Nothing has been committed or pushed.
