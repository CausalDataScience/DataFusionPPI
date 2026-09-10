# AISTATS Experiment Protocol: Four Studies for the ATE and CATE Results

## Status and scope

This document specifies the four experiments that will support the manuscript
sections 3 and 4 as of 2026-09-09 (mirror commit `d8bfb51` plus the removal
of Section 4.5).  Every experiment is tied to the theorems it checks, uses
100 Monte Carlo repetitions per design cell, and reuses the nested honest
machinery of `code/total_budget_nested_benchmark.py` wherever the estimator
is an ATE estimator.  The CATE studies need a new script that shares the
sampling, nuisance, and density-ratio code with the ATE script.

Decisions fixed with the author on 2026-09-09:

1. ATE-2 includes a real covariate-shift panel built from the NSW experiment
   and the CPS comparison group.
2. CATE-1 includes a frozen neural basis in addition to additive cubic
   B-splines.
3. Repetitions are 100 per cell in every study.

Notation follows the manuscript: $n_S^{\mathrm{nuis}},n_S^{\mathrm{tune}},n_S^{\mathrm{eval}}$
are the role-specific sample sizes of source $S\in\{R,O\}$, $Z_\lambda=Z_0+\lambda\Delta$
is the fused pseudo-outcome, $\widehat g$ is the OBS treatment-effect
prediction, $r_0$ is the exact density ratio, $(A,B,C,D)$ are the ATE variance
coefficients of Definition 5, $\mathcal R$ is the CATE risk of Definition 11,
$\widehat C_j$ is the validation score of Definition 18, and
$\varepsilon_{\mathrm{eval}}(\alpha)$ is the radius of Theorem 6.

## Zotero and source preflight

Baseline methods are named by type here.  Before any baseline is implemented
or cited, its bibliographic record must be located in the Zotero local library
(first source) and then, if absent, in Semantic Scholar or arXiv, following
`research/AGENT.md` rule 11.  Items already in `manuscript/reference.bib`:
shrinkage combination (`rosenman2020shrinkage`), adaptive combination
(`cheng2021adaptive`), experimental grounding (`kallus2018removing`), direct
multi-source CATE learning (`li2022directfusion`), PPI and PPI++
(`angelopoulos2023ppi`, `angelopoulos2024ppipp`).  Items to be verified before
use: semi-supervised ATE estimation with unlabeled covariates, prognostic
covariate adjustment from historical data, prediction-powered generalization
of causal inferences, and predictions-as-surrogates.

Data sources and their access status:

- Arbitrary nonlinear SCM families: `code/` (already implemented).
- Tennessee STAR: `https://vincentarelbundock.github.io/Rdatasets/csv/AER/STAR.csv`
  (already used for covariates; sha256 pinned in the script).  The same file
  contains the kindergarten class-type assignment and the kindergarten reading
  and mathematics scores, which CATE-2 uses as real treatment and real outcomes.
- NSW experimental sample (Dehejia and Wahba subset, 445 units) and the CPS
  comparison sample (15,992 units): public files distributed by Dehejia.  The
  download is not yet done; it requires the author's approval and a pinned
  sha256 before use.

## Common conventions

- Honest roles.  Each replication draws or splits data into nuisance, tuning,
  and evaluation roles for both sources, as in Assumption 1.  Default split
  fractions are $0.4/0.3/0.3$ of the RCT sample and $0.6/0.2/0.2$ of the OBS
  sample.  A cross-fitted variant (five folds, same total budget) is run only
  where stated.
- Nuisances.  $\widehat\mu_R$, $\widehat\mu_O$, $\widehat g$, and the density
  ratio are fitted on the nuisance samples with the existing `fit_outcome` and
  `fit_ratio` code.  The trial propensity $e$ is the design propensity.
- Coefficients.  ATE coefficients use `tune_all` (Algorithm 1 with the box
  $[0,1]^2$).  CATE coefficients use the grid
  $\mathcal G=\{0,0.25,0.5,0.75,1\}^2$ and Algorithm 2.
- Repetitions and seeds.  100 repetitions per cell, seed derived from the
  existing `stable_seed` scheme with base seed 190602.
- Reporting.  Every table reports the mean and a Monte Carlo standard error
  over the 100 repetitions.  Every RMSE is reported as a ratio to the RCT-only
  estimator in the same cell.

## ATE-1: when fusion helps, and whether Algorithm 1 finds it

Theorems checked: Theorem 1 (exact variance), Theorem 2 and Corollary 2.1
(oracle coefficients), Theorem 3 (three-way comparison), Theorem 4
(normal approximation after honest evaluation).

Design factors (full factorial unless stated):

- RCT size $n_R\in\{50,100,200,400\}$ (total, before the honest split).
- OBS size $N_O\in\{1000,5000,20000\}$.
- OBS prediction quality, two knobs: unmeasured-confounding strength in the
  OBS treatment assignment $\{0,\text{moderate},\text{strong}\}$, and OBS
  outcome-model specification $\{\text{correct family},\text{linear
  misspecification}\}$.  These knobs move $D=\operatorname{Cov}_R(Z_0,\widehat g)$
  and $C=\operatorname{Cov}_R(Z_0,\Delta)$.
- Data families: the three SCM families and the STAR real-$X$ family, common
  covariate marginal only ($r_0=1$).

Estimators: RCT-only AIPW; $\lambda$-only; $\omega$-only (power-tuned PPI++
form); joint (Algorithm 1); oracle-coefficient joint, where $(A,B,C,D)$ are
computed by Monte Carlo from $10^5$ fresh draws given the fitted nuisances and
$(\lambda^\star,\omega^\star)=(-C/A,D/B)$ is used directly; and the four
literature baselines listed in the preflight (semi-supervised ATE, shrinkage,
adaptive combination, naive pooling).

Metrics:

1. RMSE ratio to AIPW, bias, and empirical variance ratio.
2. Exact-variance check: for the fixed-coefficient estimator at the oracle
   coefficients, the Monte Carlo variance against Equation (13) of Theorem 1,
   evaluated with the same fitted nuisances.
3. Coefficient recovery: scatter of $(\widehat\lambda,\widehat\omega)$ from
   Algorithm 1 against $(-C/A,D/B)$, and the boundary-hit rate.
4. Inference: coverage and mean length of the 95% interval of Theorem 4 for
   the joint estimator, and of the AIPW interval.
5. Honest split against cross-fitting: RMSE ratio and coverage of the joint
   estimator under both schemes, one panel per data family.

Expected outputs: one figure with RMSE ratio against $n_R$, panels by
prediction quality and data family; one table with coverage; one scatter
figure for coefficient recovery.

## ATE-2: covariate shift and the density ratio

Theorems checked: Theorem 1 (drift identity), Definition 12 and Proposition 6
(calibrated balancing ratio and automatic balance), Theorem 5 and
Corollary 5.1 (feasible normal approximation with an estimated ratio).

Synthetic panel.  The shifted regime is regenerated so that $r_0$ has a
closed form: the OBS covariate law is the RCT law tilted by
$\exp(\gamma^\top x)$ with normalizing constant computed numerically, and the
shift strength is $\|\gamma\|\in\{\text{mild},\text{strong}\}$ with overlap
kept above $0.05$.  Sizes $n_R\in\{100,400\}$, $N_O=5000$, three SCM
families.  Ratios compared:

1. oracle $r_0$;
2. unweighted $r\equiv1$, which ignores the shift;
3. classifier odds (current `fit_ratio`);
4. calibrated balancing ratio with $\phi$ equal to the raw covariates;
5. calibrated balancing ratio with $\phi$ equal to the raw covariates plus
   $\widehat g$ (loss-targeted balancing).

Metrics: bias of the joint estimator against the drift formula
$\omega\,\mathbb E_O[\{\widehat r(X)-r_0(X)\}\widehat g(X)]$ computed by Monte
Carlo; RMSE; coverage of the feasible interval of Theorem 5; the balancing
residual $\mathbb P_O^{\mathrm{nuis}}(\widehat r\phi)-\mathbb P_R^{\mathrm{nuis}}(\phi)$
as a diagnostic.

NSW and CPS panel.  The RCT is the NSW experimental sample (185 treated,
260 control); the OBS sample is built from a random half of the NSW treated
units together with the CPS controls, so that no unit appears in both
sources.  The remaining NSW units form the RCT, subsampled to
$n_R\in\{100,200,\text{all}\}$.  The reference value is the full-experiment
difference in means of 1978 earnings.  Because the reference value carries its
own sampling error, the panel reports point estimates, interval lengths, and
the drift diagnostic across the 100 subsamples rather than an RMSE against
a known truth.  Ratios 3 to 5 are compared; ratio 1 is unavailable.

## CATE-1: honest validation on the sieve and the oracle calculus

Theorems checked: Theorem 6 (honest validation), Theorem 7 (risk of a sieve
candidate), Corollary 7.1 (oracle sieve coefficients).

Data: the three SCM families with a heterogeneous effect $\tau(x)$ of two
strengths, and the STAR real-$X$ family with a synthetic $\tau$.  Common
covariate marginal.  Factors: $n_R\in\{100,200,400\}$, $N_O=5000$, the two
prediction-quality knobs of ATE-1, heterogeneity strength.

Learners: DRF and RF, each on two sieves.

- Additive cubic B-splines with $K\in\{3,5\}$ interior knots per continuous
  covariate and dummies for categorical covariates.
- Frozen neural basis: a two-hidden-layer network ($64$ units, ReLU) with a
  shared trunk $h(x)\in\mathbb R^{64}$ and two treatment heads, trained on the
  OBS nuisance sample; $b(x)=(1,h(x))$ and $\widehat g=\widehat\mu_{O,1}-\widehat\mu_{O,0}$
  from the heads.  The ridge $\rho$ enters the candidate set as an extra grid
  axis $\{10^{-2},10^{-1},1\}$, which Theorem 6 allows.

Candidates and comparators: the grid $\mathcal G$ selected by Algorithm 2
(joint); RCT-only $(0,0)$; $\lambda$-only; $\omega$-only; grid oracle chosen
with the true $\tau$; plug-in of Corollary 7.1 with $A_p,\ldots,D_p$
estimated on the tuning samples; experimental-grounding baseline and
domain-indicator pooling (after Zotero preflight).

Metrics:

1. RCT-population CATE risk $\mathcal R$ on an independent test sample of
   $10^5$ draws with the true $\tau$.
2. Regret $\mathcal R(\text{selected})-\min_{\mathcal G}\mathcal R$ against
   $2\varepsilon_{\mathrm{eval}}(0.05)$, with $B$ and $\bar r$ replaced by the
   observed maxima on the evaluation sample; report the fraction of
   repetitions with regret below the radius and the ratio of the two.
3. Leading-risk check: Monte Carlo mean of $\mathcal R(\widehat\zeta_j)-a_p^2$
   against $\mathcal J_j(\lambda,\omega)$ of Theorem 7 at four fixed grid
   points.
4. Selected coefficients against the Corollary 7.1 oracle, for DRF.

Expected outputs: risk against $n_R$ by learner and sieve; a regret histogram
with the radius marked; a coefficient scatter.

## CATE-2: real outcomes with a constructed confounded observational sample

Purpose: credibility of the full pipeline on real treatment, real outcomes,
and a representation learned from the OBS sample.

Data: STAR kindergarten file.  Treatment is small class against regular class
(regular-with-aide is dropped).  Outcomes are the kindergarten mathematics
score and, as a second outcome, the reading score.  Covariates are the five
columns already used plus school.  The design propensity is the within-school
share of small classes.

Construction per replication: split students into an RCT sample of
$n_R\in\{200,400\}$, an evaluation sample of 2000, and an OBS pool from the
rest.  The OBS sample is drawn by outcome-dependent selection: treated
students are included with probability increasing in their score and
controls with probability decreasing in it, so that the naive OBS contrast is
biased upward and the confounding is real rather than simulated.

Learners and comparators as in CATE-1, with the neural basis trained on the
OBS nuisance part.  Metrics: for the ATE, the full-file experimental estimate
is the reference value; for the CATE, the difference of validation scores
$\widehat C_j$ on the evaluation sample, which by Theorem 6 item 1 estimates
the risk difference between candidates unbiasedly when $r_0=1$.

## Outputs, expected shapes, and QA

- Replication-level CSV, summary CSV, and composite PNG per study, following
  the `output_paths` convention in `materials/`.
- Each summary row carries the design cell, estimator, repetitions, mean,
  Monte Carlo standard error, and the metric name.
- QA: seeds reproduce the CSV byte for byte; the RCT-only estimator matches
  the existing benchmark on the shared cells; the oracle-coefficient joint
  estimator has variance no larger than the RCT-only estimator in every cell
  (Corollary 2.1), which is a check of the implementation rather than of the
  method.

## Reuse of the existing script

Reused as is: `synthetic_sample`, `load_star_support`, `star_sample`,
`fit_outcome`, `fit_ratio`, `pseudo`, `tuning_moments`, `tune_all`,
`folds`, `inner_scores`, `evaluate`, `summarize`, `validate`, `make_figure`,
`write_outputs`.  New code: the closed-form tilted shift and its exact ratio;
the calibrated balancing ratio; the four ATE baselines; the NSW and CPS
loader; the sieve fitter (splines and neural basis), the grid validation of
Algorithm 2, and the risk, regret, and $\mathcal J_j$ computations; the STAR
real-outcome constructor.
