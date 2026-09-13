# Traceability Map: Manuscript, Implementation, Test, Result

Answers P0-A of the revised plan.  Every claim the experiments support is traced
from the manuscript to the function that implements it, the test that checks it,
the file that holds the numbers, and the figure that shows them.  Each run rule
is classified as `manuscript` when it implements the stated algorithm exactly,
`extension` when it adds a practical device the manuscript does not state, and
`exploratory` when it is a diagnostic rather than evidence.

Run identifier of this build: `rev4-2026-09-12`.  The earlier builds are
`rev1-2026-09-11` (before the variance analysis), `rev2-2026-09-11` (before the
corrections of this map) and `rev3-2026-09-12` (before the second review).
Results from different identifiers are never placed in one table.  In `rev4`
every one of the five studies was re-executed from one code state after the last
fix landed, so no table mixes a corrected module with an uncorrected result.

## 1. Average-effect claims

| claim | manuscript | implementation | test | result file | figure |
| --- | --- | --- | --- | --- | --- |
| AIPW score is trial-valid for any fixed regression | Equation (5), Proposition 1 | `fusion_core.aipw_score` | test 2 | `ate1_*_replications.csv` | |
| fused estimator at fixed coefficients | Definition 1, Equation (15) | `fusion_core.ate_estimate` | test 2 | same | |
| exact conditional variance | Theorem 1, Equation (13) | `fusion_core.theorem1_variance` | test 3, conditional experiment | `conditional_variance_400_summary.csv` | |
| oracle coefficients | Theorem 2, Equation (22) | `fusion_ate.oracle_coefficients` | test 4 | `ate1_*_variance.csv` | figure 1 |
| ordering of the three variances | Theorem 3 | comparison rows `whole_method` | test 4 | `ate1_*_variance.csv` | figure 1 |
| plug-in selection rule | Algorithm 1, Equation (39) | `fusion_core.algorithm1` | test 2 | `ate1_*_reporting.csv` | figure 2 |
| evaluation-sample variance and interval | Equation (43), Theorem 4 | `fusion_core.evaluation_variance` | conditional experiment | `ate1_*_pooled.csv` | figure 3 |
| drift of an estimated weight | Theorem 1, Equation (12) | `exp_ate2.drift_reference` | test 5 | `ate2_shift_*.csv` | |
| calibrated balancing | Definition 12, Proposition 6 | `fusion_core.ratio_balanced` | tests 5 and 13 | `ate2_shift_summary.csv` | |
| shrinkage, adaptive, naive pooling and transported observational comparators | handoff sections 4.2 to 4.4, not manuscript claims | `fusion_baselines.shrinkage_weight`, `adaptive_weight`, `combine`, `naive_pool_spec`, `transported_obs_estimate` | tests 12 and 14 | `ate1_*_reporting.csv`, `ate1_*_pooled.csv` | figure 3 |

## 2. Conditional-effect claims

| claim | manuscript | implementation | test | result file | figure |
| --- | --- | --- | --- | --- | --- |
| the two fusion losses | Definition 15 | `fusion_cate.cate_replication` | test 6 | `cate1_*_replications.csv` | |
| sieve closed form and lambda path | Proposition 10 | `fusion_core.sieve_solver` | test 6 | same | |
| validation score | Definition 18, Equation (63) | `fusion_core.validation_score` | test 7 | same | |
| honest validation bound | Theorem 6, Equation (64) | `fusion_core.eval_radius` | bounded design only | `cate1_*_replications.csv`, rows with `radius_applicable = 1` | |
| risk of a sieve candidate | Theorem 7 | `fusion_cate._oracle_calculus` | decomposition identity | `cate1_*_decomposition.csv` | figure 4 |
| oracle sieve coefficients | Corollary 7.1 | `fusion_cate._oracle_calculus` | | `cate1_*_summary.csv` | |
| source-indicator pooled comparator | handoff section 4.6, not a manuscript claim | `fusion_baselines.domain_indicator_pool`, `source_indicator_coefficients` | test 11 | `cate1_*_replications.csv`, `cate2_star_real_replications.csv`, rows with `learner = domain_indicator_pool` | figure 5 |
| shape of the observational confounding bias | diagnostic for the comparator result, not a manuscript claim | `probe_confounding_shape.py` | | `confounding_shape.csv` | figure 5 |
| the conditional effect when the two covariate laws differ | Definition 12 and Theorem 1 applied to the conditional rule | `exp_cate3.py` with the five routes of `exp_ate2.ratio_specs` | tests 5 and 13 | `cate3_shift_spline_*.csv`, `cate3_shift_neural_*.csv` | figure 6 |
| observational fit with a trial-fitted correction | handoff section 4.5, not a manuscript claim | `fusion_baselines.experimental_grounding` | | same, rows with `learner = experimental_grounding` | |
| bounded design that satisfies the Theorem 6 range condition | supports Theorem 6 | `fusion_data.BoundedFamily` | test 16 | `cate1_*` rows with `family = bounded` | |
| one covariate law shared by both real-data samplers | supports the ratio of one | `fusion_data.star_cohort`, `star_cate2_replication` | tests 8 and 9 | `cate2_star_real_*.csv` | |

## 3. Classification of every run rule

| rule | class | note |
| --- | --- | --- |
| AIPW score, fused estimator, plug-in selection | manuscript | implemented as stated |
| exact variance and evaluation variance | manuscript | implemented as stated |
| rectangle `[0,1]^2` and thresholds `1e-3 * Var(Z0)` | extension | the manuscript leaves the rectangle and the thresholds to the user |
| trial-only propensity trimming at `[0.15, 0.85]` | extension | needed for the AIPW and R-Fusion denominators; never applied to the observational sample, because that would break the transport identity |
| clipping of the candidate and of the fitted prediction | extension | a prespecified deterministic transformation of the candidate class; Theorem 6 then applies to the clipped class |
| clipping of the evaluation pseudo-outcome | removed in `rev3` | it changed the conditional mean and so changed the target; the score is bounded below without it, because the pseudo-outcome enters only through a nonnegative square |
| ridge `rho` inside the candidate grid | extension | enlarges the candidate set, which Theorem 6 allows |
| trace-scaled jitter in the sieve solver | extension | numerical; inert once the grid ridge is present |
| radius of Theorem 6 | manuscript on the bounded design, `N/A` elsewhere | the range condition holds a priori only where the design bounds the outcome, the propensity and the ratio |
| oracle coefficients from fresh draws | exploratory | a ceiling for comparison, not an estimator |
| grid oracle chosen with the true effect | exploratory | a ceiling; the feasible single-channel rules select on the score |
| conditional-variance experiment at 400 evaluations | exploratory expansion | the expansion from 100 was decided after seeing the 100-replication result and is labelled as a diagnostic |
| declared bound `B = 9` on the bounded design | manuscript | the worst case of the sampler is 8.667, so 9 is valid and conservative; the 4.333 used in `rev3` was not attainable as a bound |
| law-weighted STAR reference | manuscript | the estimand both samplers share is the contrast weighted by the pattern law, not the raw cohort difference |
| comparator weights fixed on the tuning samples | extension | the handoff leaves the split unspecified; deciding the weight on the evaluation sample would invalidate the interval |
| transport ratio fixed to one in ATE-1, CATE-1, CATE-2 and the conditional-variance study | manuscript | those designs draw both sources from one covariate law, so the ratio is one before any data is drawn; the shift studies of Section 6 and Section 7.2 are where the ratio has to be estimated |
| conditional shift study at one sieve and one trial size for the neural basis | exploratory extension | added after the second review to close the gap the review process surfaced; the spline covers both shift levels and both trial sizes |

## 4. Provenance of each run

| item | value |
| --- | --- |
| base seed | 190602, with per-cell seeds from `stable_seed` |
| environment | Python 3.14, numpy 2.4.6, pandas 3.0.3, scikit-learn 1.9.0, torch 2.13.0, pinned in `code/requirements.txt` |
| data sources | STAR csv sha256 `0e8b179e...3ae6`; NSW, CPS and PSID sha256 recorded in `fusion_data.NSW_FILES`; every download is hash-checked at run time |
| role sizes | recorded per replication as `n_rct_eval` and `n_obs_eval`, and the attempted and effective replication counts appear in every reporting file |
| determinism | acceptance test 1 reruns a cell and compares the replication file byte for byte |

## 5. What this map does not claim

The map records which code implements which statement.  It is not a proof check.
The manuscript's proofs were not re-derived here, and the tests verify numerical
agreement and structural properties rather than the theorems themselves.
