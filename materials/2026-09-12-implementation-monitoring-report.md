# Implementation Monitoring Report: DataFusionPPI Experiments

Revision 4, 2026-09-12, after the second review that returned FAIL.  Run
identifier `rev4-2026-09-12`.  Every number below comes from this run.  All five
studies were re-executed from one code state after the last fix landed, so no
table mixes a corrected module with an uncorrected result.  Results from `rev1`,
`rev2` and `rev3` are never mixed into the same table.  Nothing has been
committed or pushed.

The short version: every finding in the second review was confirmed against the
code, all of them are fixed, and fixing them changed two more headline
conclusions.  One of those two is unfavourable to the paper.  A comparator that
the first review asked for, and that revision 3 computed but did not report,
beats the conditional-effect method in every synthetic design.  Section 7
explains why, and the reason is a property of the designs rather than a coding
error.

## 1. What the second review found, what the code showed, what changed

| finding | confirmed? | what the fix changed |
| --- | --- | --- |
| P0-1: the real-outcome driver kept three errors that revision 3 reported as fixed | yes, all three | revision 3 fixed `fusion_cate.py` and left the parallel code path in `exp_cate2.py` untouched.  The evaluation pseudo-outcome was still clipped on lines 87, 92 and 118, the trial-only rule was still pinned to the smallest ridge on line 97, and the source indicator was still added to every coefficient on line 117.  All three now call the same shared functions the other driver calls |
| P0-2: the bounded design returned the pre-clipping treatment effect | yes | the outcome is built inside its budget instead of being clipped afterwards, so the returned effect is the effect of the sampler.  Test 16 checks this exactly |
| P0-3: the Theorem 6 bound of 4.333 was not guaranteed | yes.  The review's counterexample gives 4.667, and the correct worst case is 8.667 | the design now declares 9 |
| P0-4: the STAR reference was the raw cohort difference | yes.  The raw difference is 0.1736 for mathematics, the law-weighted value is 0.1678 | the reference is the law-weighted contrast, which is the estimand the two samplers share |
| P0-5: the average-effect comparators differed from the approved specification | yes.  The covariance term between the two channels was missing, and pooling used the wrong denominator | all four comparators are rewritten to the handoff formulas in sections 4.2 to 4.4 |
| P0-6: four acceptance tests did not test what they claimed | yes, all four | rewritten with error injection; see Section 9 |
| P0-2: the bounded confounding knob changed the seed, and the heterogeneity knob was unused | yes | both knobs now enter the sampler and nothing else |

Four of my own tests were written so that they could not fail.  Test 10 read
`not same or True`.  Test 9 compared a quantity against itself.  Test 11 rebuilt
the correct formula inside the test instead of calling the production function.
Test 12 only checked that a number existed.  Each now injects the specific error
it is meant to catch and checks that the check fires.

## 2. Two conclusions that changed

**The neural conditional-effect gain on real outcomes disappears.**  Revision 3
reported that on the STAR data the selected learner improved the held-out score
on the neural basis by 0.683 with a z statistic of minus 21.0.  That comparison
used a trial-only baseline pinned to the smallest ridge, because the real-outcome
driver never received the fix.  With the same ridge choice available to both
rules, the pooled gap is minus 0.008 with z of minus 0.32 for the doubly robust
loss and minus 0.018 with z of minus 0.78 for the residual loss.  Neither
excludes zero.  Splitting by outcome shows the two outcomes disagree: reading
improves with z of minus 3.07, mathematics does not move.  Almost all of the
revision 3 neural result was the unequal comparison.

**A naive pooled learner beats the method on every synthetic conditional
design.**  The source-indicator pooled learner, which stacks the trial
pseudo-outcomes and the observational pseudo-outcomes with one indicator column,
reaches a median true risk between 5 and 9 times lower than the selected fusion
rule, at every confounding level, in all five families.  Section 7 gives the
measurement and the reason.

The three conclusions that changed in revision 3 survive this revision.  The
omega channel reduces variance in the whole procedure, the synthetic
conditional-effect gain over a trial-only baseline is real but much smaller than
revision 2 claimed, and the STAR spline result holds.

## 3. The omega question, stated carefully

| scope | comparison | rate | 95 percent interval | oracle rate | oracle interval |
| --- | --- | ---: | --- | ---: | --- |
| conditional | omega alone | 2.73% | 1.91 to 3.56% | 3.80% | 2.88 to 4.70% |
| conditional | omega given lambda | 4.11% | 3.06 to 5.02% | 5.85% | 4.70 to 6.91% |
| conditional | whole method | 35.28% | 33.27 to 36.98% | | |
| whole procedure | omega alone | 1.18% | 0.24 to 1.98% | 2.32% | 1.42 to 3.18% |
| whole procedure | omega given lambda | 2.34% | 1.04 to 3.19% | 3.17% | 2.11 to 4.10% |
| whole procedure | whole method | 37.99% | 23.07 to 43.93% | | |

Every arm of a comparison is paired on lambda inside the same replication, and
each interval comes from one resample of replications within cells, with cells
held as the fixed design they are and the whole weighted average recomputed each
time.  These six rows are identical to revision 3 to the reported precision,
which is the expected result: the second review's fixes touched the comparators,
the bounded design and the real-outcome driver, none of which enter this table.

Per family, the whole-procedure reduction with lambda held fixed is 4.82 percent
on the first Gaussian family with an interval of 2.27 to 7.12, 2.75 percent on
the second with 0.34 to 4.99, 1.10 percent on the third with minus 0.44 to 2.50,
and 0.69 percent on the STAR-covariate family with minus 2.53 to 1.89.  The
pooled signal therefore comes from two of the four families and is not visible in
the other two.

On the cost side, absolute bias moves by at most 0.0014, pooled coverage by at
most 0.002, and the interval shortens.  These are the measured differences in
these cells at 100 replications, which is weaker than a claim that nothing is
paid anywhere.

## 4. The exact variance of Theorem 1

Fifteen frozen fits, seven coefficient rules each, 400 evaluation draws per fit,
42,000 rows.  The replicated spread matches the formula in all 105 fit and rule
combinations within three bootstrap standard errors.  The z statistics run from
minus 1.54 to 1.68, and the mean z by rule lies between minus 0.156 and minus
0.039.  The formula value is itself a numerical integral over 200,000 draws, so
it carries its own approximation error of about half a percent.  The expansion
from 100 to 400 evaluations was decided after seeing the 100-replication result
and is recorded as an exploratory diagnostic.

## 5. Average effect, with the specified comparators

Root-mean-square error as a ratio to the trial-only estimator, averaged over the
four families, trial-size axis:

| estimator | 50 | 100 | 200 | 400 |
| --- | ---: | ---: | ---: | ---: |
| trial only | 1.000 | 1.000 | 1.000 | 1.000 |
| lambda only | 0.690 | 0.719 | 0.755 | 0.790 |
| omega only | 0.994 | 0.984 | 0.998 | 0.992 |
| joint | 0.685 | 0.703 | 0.746 | 0.783 |
| joint at oracle coefficients | 0.664 | 0.673 | 0.733 | 0.725 |
| fixed half weight | 0.766 | 0.772 | 0.816 | 0.826 |
| semi-supervised | 1.009 | 1.014 | 1.011 | 1.005 |
| shrinkage | 0.585 | 0.719 | 0.977 | 1.057 |
| adaptive | 0.580 | 0.740 | 1.109 | 1.258 |
| naive pooling | 1.460 | 1.274 | 1.323 | 1.300 |
| observational AIPW, transported | 0.564 | 0.722 | 1.074 | 1.427 |

The five fusion rows are unchanged from revision 3, and the five comparator rows
moved, which is the expected signature of a change confined to the comparators.
Shrinkage and the adaptive rule beat the joint estimator at the smallest trial
and lose to it from 200 upward, because their squared bias is fixed while the
trial-only variance they are competing against keeps falling.

Pooled coverage of the nominal 95 percent interval over all 44 cells, with the
point and the interval recomputed inside one resample:

| estimator | coverage | interval |
| --- | ---: | --- |
| trial only | 0.9416 | 0.9345 to 0.9484 |
| omega only | 0.9423 | 0.9355 to 0.9493 |
| lambda only | 0.9439 | 0.9368 to 0.9502 |
| joint | 0.9423 | 0.9350 to 0.9491 |
| shrinkage | 0.5798 | 0.5661 to 0.5934 |
| adaptive | 0.4289 | 0.4164 to 0.4418 |
| naive pooling | 0.2218 | 0.2141 to 0.2296 |
| observational AIPW, transported | 0.2698 | 0.2611 to 0.2791 |

The four fusion estimators sit about one point below the nominal level.  Of the
four, only the lambda-only interval reaches 0.95; the upper ends for the joint,
omega-only and trial-only rules are 0.9491, 0.9493 and 0.9484, so those three
intervals stop short of the nominal level.  The shortfall is shared with the
trial-only estimator, which uses no observational data at all, so it is a
property of the normal approximation at these trial sizes rather than a cost of
fusion.

Every reporting row has 100 effective replications and a failure rate of zero.
The identity relating mean squared error, variance and squared bias holds to 3.6
times ten to the minus fifteen.

## 6. Covariate shift, with paired draws

This is the only study in which the two covariate laws differ.  In the other
four studies the observational and trial covariate laws are equal by
construction, so the transport ratio is exactly one and the pipeline is given
the constant one rather than an estimate.  The synthetic families are drawn in
the `shared` regime, where the oracle ratio returns one by definition, and the
real-outcome study draws both sources from a single pattern law.  Everything
reported in Sections 3, 5, 7 and 8 therefore measures the two fusion channels
with the transport weight switched off, and the cost of estimating that weight
appears only here.

All five ratio routes see the same draw in a replication.  Joint estimator,
three synthetic families and two trial sizes:

| shift | ratio | bias | z | coverage |
| --- | --- | ---: | ---: | ---: |
| mild | oracle | -0.0016 | -0.11 | 0.947 |
| mild | unweighted | +0.0331 | +2.34 | 0.945 |
| mild | classifier | -0.0004 | -0.03 | 0.948 |
| mild | balanced on covariates | -0.0026 | -0.18 | 0.948 |
| mild | balanced on covariates and prediction | -0.0002 | -0.01 | 0.947 |
| strong | oracle | +0.0255 | +1.75 | 0.930 |
| strong | unweighted | +0.0943 | +6.52 | 0.912 |
| strong | classifier | +0.0298 | +2.04 | 0.930 |
| strong | balanced on covariates | +0.0275 | +1.88 | 0.930 |
| strong | balanced on covariates and prediction | +0.0303 | +2.06 | 0.928 |

At the strong shift the trial-only estimator, which uses no observational
channel at all, shows the same offset of 0.0299 with z of 1.71.  The offset is
therefore common to the design rather than caused by the fusion channels or the
weight.  What the weight removes is the extra 0.065 that the unweighted route
carries on top of it.

On the NSW panels the comparison is against the experimental benchmark rather
than a known truth, so the level of the bias is not interpretable on its own.
What is comparable is the gap between routes at a fixed panel and trial size.
At the CPS panel with a trial of 100 the three weighted routes carry bias of
minus 0.41, minus 0.44 and minus 0.46 with coverage 0.94 to 0.95, while the
unweighted route carries minus 1.16 with coverage 0.91.  The trial-only estimator at the same cell carries
minus 0.08 with coverage 0.94, so part of the weighted routes' offset is the
observational channel and part is the panel.  These panels are reported as a
stress test of the weighting routes, not as evidence that the method recovers
the NSW effect.

## 7. Conditional effect on synthetic designs

Theorem 6 is checked only on the bounded design, whose outcome, propensity and
ratio are bounded before any data is drawn.  The declared bound is now 9 rather
than the 4.333 of revision 3, which the review showed was not guaranteed.  The
radius scales with the square of the bound, so it moves from 42.4 to 183.1.  The
median regret is 0.0014 and the bound holds in every replication, which makes the
statement true and very loose.  The Gaussian and STAR families report the radius
as not applicable, because their outcome is unbounded and the range condition of
the theorem does not hold.

Median per-replication risk ratio of the selected learner to a trial-only learner
that chooses its ridge on the same score, pooled over the doubly robust and
residual losses:

| sieve | 100 | 200 | 400 |
| --- | ---: | ---: | ---: |
| spline, 3 knots | 0.231 | 0.408 | 0.625 |
| spline, 5 knots | 0.182 | 0.368 | 0.586 |
| neural | 0.645 | 0.754 | 0.834 |

Selection beats the trial-only learner in 88.2 percent of replications.  Split
by loss, the ratio pooled over all cells is 0.438 for the doubly robust loss and
0.411 for the residual loss on the three-knot spline, and 0.765 and 0.756 on the
neural basis.  Wherever a single number for the neural sieve appears in the
manuscript it must say which loss it refers to.

The median risks of the feasible single-channel rules sit between the trial-only
and the selected rule, and below each sits its own true-risk oracle, as it
should.  For the three-knot spline the medians are 2.923 for trial only, 2.098
for the lambda channel with 1.951 for its oracle, 1.749 for the omega channel
with 1.634 for its oracle, 1.299 for the selected rule and 1.118 for the grid
oracle.

The risk on the fixed anchor splits into prediction variance and squared bias,
and the split is exact to 7.1 times ten to the minus fourteen.  For the doubly
robust loss on the three-knot spline the prediction variance falls from 15.7,
24.7 and 23.9 to 2.3, 1.5 and 1.0 across the three trial sizes, while the
squared bias falls from 0.18, 0.47 and 0.14 to 0.05, 0.04 and 0.03.  Prediction
variance is the larger part of the risk in every rule, so the gain is mostly a
stabilisation of repeated fitting.

### 7.1 The comparator that wins

The first review asked for a source-indicator pooled learner as the conditional
analogue of naive pooling.  It was implemented and computed in revision 3 and
was not reported.  It should have been.  Median true risk on the anchor, doubly
robust loss, by family and confounding strength:

| family | confounding | fusion | trial only | pooled | pooled / fusion |
| --- | ---: | ---: | ---: | ---: | ---: |
| bounded | 0 | 0.043 | 0.102 | 0.006 | 0.14 |
| bounded | 2 | 0.045 | 0.098 | 0.005 | 0.11 |
| Gaussian 1 | 0 | 1.282 | 3.035 | 0.186 | 0.15 |
| Gaussian 1 | 2 | 1.401 | 3.363 | 0.187 | 0.13 |
| Gaussian 2 | 2 | 1.744 | 4.223 | 0.318 | 0.18 |
| Gaussian 3 | 2 | 2.064 | 5.172 | 0.366 | 0.18 |
| STAR covariates | 2 | 0.995 | 1.650 | 0.504 | 0.51 |

The advantage does not shrink as confounding grows, which is the clue.  Measured
directly on 200,000 observational draws, the confounding bias of the
observational contrast is large in level and nearly flat in the covariates.  On
the first Gaussian family at confounding 2 the mean bias is 0.557 while its
spread across quintiles of the leading covariate is 0.030, against a true
conditional effect with standard deviation 0.588.  The other two families behave
the same way, with quintile spreads of 0.070 and 0.096.

The confounder enters these outcome equations additively and independently of
the covariates, so it shifts the observational contrast by an almost constant
amount.  A source indicator is exactly a constant, so it absorbs almost all of
that shift, and the pooled learner then keeps the precision of the large
observational sample with very little of its bias.

Two things follow.  The synthetic conditional designs cannot discriminate
between the fusion rule and naive pooling, because they confound the level of
the contrast and not its shape.  A design that supports such a claim has to make
the confounding depend on the covariates.  Until that design exists, the
conditional-effect claim the experiments support is the narrower one: against a
trial-only learner with a fair choice of ridge, selection lowers the risk in
every family and at every trial size measured here.

The other specified comparator, the observational fit with a trial-fitted
correction, is worse than the trial-only learner on the neural basis, with a
median risk ratio of 5.5, and indistinguishable from it on both splines, with
1.002.

Both halves of this argument are in `materials/figure5_pooled_comparator.png`.
The flatness measurement is produced by `code/probe_confounding_shape.py`, which
writes `materials/confounding_shape.csv`.

### 7.2 The conditional effect when the two covariate laws differ

Everything above draws both sources from one covariate law, so the transport
ratio is one and the pipeline is handed the constant one.  This study moves the
trial covariate mean away from the observational mean and lets the same five
ratio routes compete, which is the conditional counterpart of Section 6.  Under
the strong shift the true ratio ranges over 0.105 to 13.5 on the first Gaussian
family, 0.199 to 7.1 on the second and 0.510 to 1.95 on the third.

Median risk ratio of the selected learner to the trial-only learner, three
Gaussian families, strong shift, trial size 200:

| route | spline, 3 knots | neural |
| --- | ---: | ---: |
| true ratio | 0.373 | 0.758 |
| fitted by a classifier | 0.384 | 0.789 |
| fitted by calibrated balancing | 0.384 | 0.789 |
| ignored, ratio set to one | 0.405 | 0.806 |
| the same cells with no shift at all | 0.375 | 0.764 |

Three readings follow.  A known transport weight fully restores the no-shift
result, 0.373 against 0.375 on the spline and 0.758 against 0.764 on the neural
basis, so the shift itself costs nothing once the weight is right.  A fitted
weight recovers part of the gap between ignoring the weight and knowing it, 66
percent on the spline and 35 percent on the neural basis.  Ignoring the weight
is the worst route in both, and the whole procedure still beats the trial-only
learner even then.

The mechanism is visible in the channels.  The lambda channel never touches the
ratio, and its risk ratio is identical to three decimals across all four routes,
0.794 on the spline and 0.821 on the neural basis.  The entire effect of the
ratio sits in the omega channel, where the spline moves from 0.504 at the true
ratio to 0.529 when the ratio is ignored, and the neural basis moves from 0.958
to 1.000, which is the omega channel losing its whole value.

At the mild shift the five routes are indistinguishable, with the pooled risk
ratio between 0.381 and 0.396 and the ignoring route at 0.390, sitting inside
that range rather than below it.

The design is `code/exp_cate3.py`, run at 100 replications with the spline over
both shift levels and both trial sizes and with the neural basis at the strong
shift and trial size 200.  The risk decomposition is exact to 3.6 times ten to
the minus fifteen.  Figure 6 shows this study beside Section 6.

## 8. Real outcomes

The STAR construction draws both sources from one law over 27 covariate patterns
covering 98.3 percent of the retained cohort.  The reference effect is the
law-weighted contrast, 0.1678 for mathematics and 0.1585 for reading, rather
than the raw cohort difference of 0.1736 and 0.1723.  The naive observational
contrast runs from 0.59 to 0.99 depending on the tilt, against an experimental
0.17, so the confounding is large.

Held-out score gap of the selected rule against the trial-only rule, where a
negative number favours selection:

| outcome | sieve | loss | gap | z |
| --- | --- | --- | ---: | ---: |
| mathematics | spline, 3 knots | doubly robust | -0.064 | -10.8 |
| mathematics | spline, 3 knots | residual | -0.085 | -11.5 |
| mathematics | neural | doubly robust | +0.043 | +1.00 |
| mathematics | neural | residual | +0.037 | +0.88 |
| reading | spline, 3 knots | doubly robust | -0.077 | -11.0 |
| reading | spline, 3 knots | residual | -0.097 | -9.6 |
| reading | neural | doubly robust | -0.058 | -3.07 |
| reading | neural | residual | -0.073 | -3.72 |

On the coarse spline selection improves the score for both outcomes and both
losses.  On the neural basis the two outcomes disagree, and pooling them gives
minus 0.008 with z of minus 0.32 for the doubly robust loss, which is the
cancellation of a real reading effect against no mathematics effect rather than
a common null.

The prediction variance, weighted by the trial population law, moves in the same
direction as the score on the spline and against it on the neural basis.  On the
three-knot spline it falls from 0.308 to 0.242 for the doubly robust loss and
from 0.330 to 0.239 for the residual loss.  On the neural basis, averaged over
the two outcomes, it rises from 0.301 to 0.349 and from 0.299 to 0.318; the one
cell that moves the other way is reading under the residual loss, where it falls
from 0.304 to 0.284.  Stability and correctness are reported
separately; neither is inferred from the other.

The pooled comparator behaves here as it does on the synthetic designs, but not
uniformly.  On the spline it beats the fusion rule, with a gap of minus 0.200
against minus 0.064.  On the neural basis it is worse than the trial-only rule,
with a gap of plus 0.042.  The real cohort therefore does not reproduce the
clean dominance the synthetic designs show, which is consistent with the
explanation in Section 7.1: real confounding is not a constant shift.

For the average effect the fusion estimators match the trial-only estimator.
Against the law-weighted reference, pooled over the eight cells, the joint
estimator carries bias of minus 0.0001 with coverage 0.961, the trial-only
estimator 0.0001 with 0.959, shrinkage 0.018 with 0.934, the adaptive rule 0.016
with 0.923, naive pooling 0.163 with 0.410, and the transported observational
AIPW 0.624 with 0.071.

## 9. Acceptance tests

Fifteen tests, all passing.  The four that the second review found vacuous are
rewritten so that each injects the error it is meant to catch:

- Test 9 now perturbs one pattern of the observational law and checks that the
  agreement test fails on the perturbed law while passing on the real one.  The
  largest difference between the two samplers over 27 patterns is 2.1 times ten
  to the minus three.
- Test 10 perturbs only the reporting half and checks that the selection is
  unchanged in 12 of 12 replications while the reported score moves in 12 of 12.
- Test 11 calls the production function rather than rebuilding the formula, and
  also evaluates the rule that revision 3 shipped.  The production rule matches
  the direct design-matrix product to 8.9 times ten to the minus sixteen, and
  the earlier rule differs by 2.419, so the check would have caught it.
- Test 12 perturbs the evaluation outcomes and checks that the weight and the
  detection decision do not move while the estimate does.

Test 16 is new and checks the bounded design against the assumptions it exists
to satisfy.  The outcome reaches 0.618 against its bound of 1, the trial
pseudo-outcome reaches 1.403 against the declared 9, the prediction reaches
0.328 against its bound of 2, and the returned treatment effect equals the
effect of the sampler exactly.

## 10. What is still open

1. **The synthetic conditional designs do not challenge naive pooling.**  This is
   the largest open item and it is a design question, not a coding one.  Section
   7.1 gives the measurement and what a discriminating design would need.
2. **The manuscript still describes the old study.**  Section 5, the abstract and
   the simulation appendix say the average-effect experiments used 20
   replications and that the conditional-effect experiments are planned.  None of
   this revision is in the manuscript.
3. **Submission form.**  The source is a working draft with the author's name and
   no AI use statement, and the official 2027 style file was not published at the
   time of checking.
4. **Novelty against current literature.**  One bibliographic record, predictions
   as surrogates, is still unverified, and no assumption-by-assumption comparison
   table against the closest prior work exists yet.
5. **The transport ratio is exercised on synthetic designs only.**  Four of the
   six studies hold the two covariate laws equal.  Section 6 covers a shifted
   law for the average effect and Section 7.2 now covers it for the conditional
   effect, both on Gaussian families.  Neither real-data study varies the two
   laws, and the conditional shift study covers one sieve at the neural basis
   and one trial size.
6. **The bounded design is loose.**  Its declared bound of 9 is correct and
   conservative, and the observed pseudo-outcome only reaches 1.4, so the
   Theorem 6 radius it produces is 130 times the observed regret.  The statement
   is honest and carries little information at these sample sizes.

The traceability map is in `materials/2026-09-12-traceability-map.md`; it lists,
for every claim, the manuscript location, the implementing function, the test,
the result file and the figure, and classifies each run rule as manuscript,
extension or exploratory.
