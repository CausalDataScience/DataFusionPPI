# Total-Budget Nested Data-Fusion Benchmark and WHI Protocol

## Status and scope

This document specifies an exploratory debugging benchmark, not a theorem
verification or a resource-neutral superiority claim.  The total sample budget
in every Monte Carlo replication is exactly

$$
n_R\in\{20,30\},\qquad N_O=5000.
$$

There is no separate training or calibration sample.  The same unique
observations rotate through training, coefficient-tuning, and evaluation roles
by honest nested sample splitting.  The four primary estimators are

1. cross-fitted RCT-only AIPW with $(\lambda,\omega)=(0,0)$;
2. HAIPW-only, which fixes $\omega=0$ and re-estimates $\lambda$ in that
   restricted problem;
3. power-tuned PPI++-only, which fixes $\lambda=0$ and re-estimates $\omega$
   in that restricted problem; and
4. Full DataFusionPPI, which jointly selects $(\lambda,\omega)$.

Here `PPI++-only` means the power-tuned coefficient described below.  The
classical fixed choice $(\lambda,\omega)=(0,1)$ is not an additional method in
this experiment.

The run uses seed `190602`, $B=20$ outer Monte Carlo replications, two
covariate regimes (`shared`, `shifted`), three fixed illustrative random SCM
draws, and one Tennessee STAR real-covariate semi-synthetic benchmark.  Twenty
replications give noisy bias and variance estimates; all conclusions must be
reported with Monte Carlo uncertainty.

## Zotero and source preflight

The Zotero-first search found the Project STAR paper by Athey, Chetty, and
Imbens under item key `3LUSRKB9` (PDF attachment `JEBP5BRP`).  It did not find a
dataset-specific STAR item.  The benchmark therefore uses the `STAR` data in
the CRAN `AER` package, whose manual calls it the Project STAR public-access
data set.  The current data have 11,598 rows and 47 variables.  The current
`AER` package declares `GPL-2 | GPL-3`; the upstream data-specific license is
not separately stated in the package documentation, so no raw data are
redistributed by this project.

- [CRAN AER package description](https://cran.r-project.org/web/packages/AER/DESCRIPTION)
- [CRAN AER reference manual](https://cran.r-project.org/web/packages/AER/AER.pdf)
- [Rdatasets transport copy of AER/STAR](https://vincentarelbundock.github.io/Rdatasets/csv/AER/STAR.csv)

The runtime transport copy is downloaded into memory, not the workspace.  The
preflight artifact is 2,124,167 bytes with SHA-256
`0e8b179ea3d883730b25008ca1293c4c8f886aa8d5d299c03ab6ff05d0123ae6`.
Any mismatch fails closed and requires source review.

The Zotero WHI preflight found item `IYTFFVSP` (Rosenman, Dominici, and
Miratrix; PDF `R3T7AHMF`) and item `W7AYUBNK` (Cheng and Cai).  WHI-CTOS itself
is gated through BioLINCC, so this run neither requests nor analyzes WHI data.

## Honest total-budget role rotation

For each source, split the one generated sample into $K=5$ outer folds.  For
outer fold $k$, the fold is evaluation data and its complement is development
data.  Thus the RCT evaluation/development sizes are $4/16$ for $n_R=20$ and
$6/24$ for $n_R=30$; the OBS sizes are $1000/4000$.

Within each outer development complement, make two inner folds.  For each
inner fold, fit $\widehat\mu_R,\widehat\mu_O,$ and $\widehat r$ on the other
inner fold, then predict the held-out inner fold.  Concatenating the inner
out-of-fold predictions gives development-only values

$$
Z_0,\qquad \Delta=Z_1-Z_0,\qquad g_R,\qquad \widehat r g_O.
$$

No evaluation-fold observation enters a nuisance fit or coefficient choice.
After tuning, refit the three nuisance objects on the full outer development
complement and score the untouched outer fold.  Each original observation is
evaluated exactly once; repeated appearances in training folds do not increase
the unique participant budget.

The training sets overlap across outer folds.  Consequently the independent
three-sample conditional variance theorem in the memo is not an exact
finite-sample variance theorem for the final algorithm.  We report empirical
Monte Carlo variance of the entire refitted algorithm.  Any score-based
plug-in variance is labeled diagnostic, never exact.

## Nuisances and four coefficient restrictions

The outcome regression is a deterministic joint ridge regression of $Y$ on a
generic feature map, $A$, and feature-by-$A$ interactions.  The ridge remains
defined even if a tiny inner RCT training fold contains only one arm; this is a
declared algebraic fallback, not a redraw.  The RCT propensity is known from
the generating design.

The density ratio is

$$
r(x)=\frac{dP_R^X}{dP_O^X}(x).
$$

It is estimated by standardized linear ridge logistic domain classification
with sample-prior correction.  Raw ratios are clipped to $[0.05,20]$ for the
estimator.  Evaluation ratios are not self-normalized; doing so would change
the estimator and conceal normalization error.

On inner out-of-fold development scores, define

$$
\begin{aligned}
A&=\widehat{\operatorname{Var}}(\Delta),
&C&=\widehat{\operatorname{Cov}}(Z_0,\Delta),\\
H&=\widehat{\operatorname{Cov}}(\Delta,g_R),
&D&=\widehat{\operatorname{Cov}}(Z_0,g_R),\\
B&=\widehat{\operatorname{Var}}(g_R)
 +\frac{n_{R,\mathrm{dev}}}{N_{O,\mathrm{dev}}}
 \widehat{\operatorname{Var}}(\widehat r g_O).
\end{aligned}
$$

Up to a coefficient-independent constant, the empirical calibration objective
is

$$
Q(\lambda,\omega)
=2\lambda C+\lambda^2A-2\omega D-2\lambda\omega H+\omega^2B.
$$

All methods use the common box $[-2,2]^2$.  The two one-channel estimators solve
their own restricted problems:

$$
\widehat\lambda_H
=\operatorname{clip}_{[-2,2]}(-C/A),
\qquad
\widehat\omega_H=0,
$$

and

$$
\widehat\lambda_P=0,
\qquad
\widehat\omega_P
=\operatorname{clip}_{[-2,2]}(D/B),
$$

with the corresponding coefficient set to zero when its denominator is zero.
Neither restricted estimator is made by taking the Full solution and zeroing
one coordinate.

For Full DataFusionPPI, the unconstrained stationary system is

$$
\begin{pmatrix}A&-H\\-H&B\end{pmatrix}
\begin{pmatrix}\lambda\\\omega\end{pmatrix}
=
\begin{pmatrix}-C\\D\end{pmatrix}.
$$

Coordinatewise clipping of this solution is not generally the box-constrained
minimum when $H\ne0$.  The implementation therefore evaluates the feasible
interior or pseudoinverse candidate, each edge's exact one-dimensional
constrained optimum, all four corners, the two restricted-axis solutions, and
the origin.  It deterministically selects the candidate with smallest value of
the original $Q$; no ridge-altered objective is used.  Every calibration fold
asserts

$$
Q(\widehat\lambda_F,\widehat\omega_F)
\leq
\min\{Q(0,0),Q(\widehat\lambda_H,0),Q(0,\widehat\omega_P)\}
$$

up to numerical tolerance.  Recording $H$ is important: population separation
gives $H=0$ under honest fixed nuisances, but a tiny empirical tuning sample
need not.

For outer fold $k$,

$$
\widehat\theta_{R,k}=\mathbb P_{R,k}Z_{0,k},
$$

and

$$
\widehat\theta_{F,k}
=\mathbb P_{R,k}(Z_{0,k}+\widehat\lambda_k\Delta_k)
+\widehat\omega_k\left\{
\mathbb P_{O,k}(\widehat r_k g_{O,k})-
\mathbb P_{R,k}g_{R,k}\right\}.
$$

The other two fold estimators use the same expression with
$(\lambda,\omega)=(\widehat\lambda_{H,k},0)$ and
$(0,\widehat\omega_{P,k})$, respectively.  All four estimates share the same
generated samples, outer and inner folds, nuisance fits, ratio fit, and
untouched evaluation observations.

The reported estimator averages the five equal-size fold estimates.  Columns
named `mean_*_used`, constrained-from-raw flags, actual boundary-hit flags,
fallback, restriction, and
`uses_ratio_channel` are method-specific.  Columns prefixed `common_*` describe
the shared nuisance, ratio, arm-count, and calibration moments and are
intentionally repeated across the four paired rows.  A `paired_design_id`
verifies that the rows came from one common split and fit.  Estimated
$r$ creates a finite-sample drift
$P_O[(\widehat r-r)g_O]$; cross-fitting prevents evaluation leakage but does
not make that drift vanish.  The output therefore records normalization error,
oracle product-bias diagnostics, ESS, maximum weight, clipping, arm counts, and
raw/clipped coefficient behavior.

## Arbitrary nonlinear SCM benchmark

The three SCMs are the same fixed random draws from the previously approved
generator, not seed-searched realizations.  For $d=4$ and scalar
$U,\varepsilon_Y\stackrel{\mathrm{iid}}\sim N(0,1)$,

$$
X=b_s+\ell U+\varepsilon_X,
\qquad
\varepsilon_X\sim N\{0,\operatorname{diag}(\sigma_X^2)\},
$$

with $b_R=b_O$ in the shared regime and $b_R=b_O+\delta$ in the shifted
regime.  The nonlinear baseline, treatment effect, RCT propensity, OBS
propensity, and coefficient-generation law are exactly those disclosed in the
preserved pilot protocol.  In particular,

$$
Y=m(X)+\gamma_UU+A\tau(X)+\sigma_Y\varepsilon_Y,
$$

the RCT treatment depends only on $X$, and OBS treatment depends on both $X$
and $U$.  The target $E_R\{\tau(X)\}$ and Gaussian density ratio are analytic.
Positive definiteness, propensity bounds, analytic target checks, and ratio
change-of-measure checks are smoke gates.

## Tennessee STAR real-$X$ semi-synthetic benchmark

The original STAR treatment variables (`stark`, `star1`, ...) and original
test outcomes (`read*`, `math*`) are never used.  The real covariate support
uses only kindergarten-era fields

`gender`, `ethnicity`, `birth`, `lunchk`, and `schoolk`.

Categorical values receive deterministic one-hot encoding; numeric values get
median imputation; missingness indicators are retained.  Constant columns are
removed and the finite support is standardized once.  Let its rows be
$x_1,\ldots,x_M$.  Independent fixed-seed draws attach $U_i$ and
$\varepsilon_i$ to every row and generate fixed potential outcomes

$$
Y_i(0)=m(x_i)+0.7U_i+\varepsilon_i,
\qquad
Y_i(1)=Y_i(0)+\tau(x_i),
$$

where $m$ and $\tau$ are fixed random nonlinear functions drawn once from a
prespecified Gaussian coefficient law.  They are not tuned using results.
Specifically, after centering and scaling the dictionary
$\Phi(x)=(x,x^2,\sin x,$ and the six pairwise products among the first four
coordinates$)$, with $p=\dim\Phi$,

$$
\beta_m\sim N(0,0.55^2I_p/p),\qquad
\beta_\tau\sim N(0,0.22^2I_p/p),
$$

and $m(x)=\Phi(x)^\top\beta_m$ and
$\tau(x)=1+\Phi(x)^\top\beta_\tau$.  Independently,
$\beta_R\sim N(0,0.25^2I_d/d)$ and
$\beta_O\sim N(0,0.35^2I_d/d)$, with

$$
e_R(x)=\operatorname{clip}\{\operatorname{expit}(0.1+x^\top\beta_R),0.2,0.8\},
$$

$$
e_O(x,u)=\operatorname{clip}\{\operatorname{expit}(-0.2+x^\top\beta_O+0.9u),0.05,0.95\}.
$$

All support-level random quantities are drawn exactly once from seed `190602`.

In the shared regime both studies sample the finite support uniformly with
replacement.  In the shifted regime, OBS remains uniform while

$$
p_R(i)=\frac{\exp\{0.35s_i\}}
{\sum_{j=1}^M\exp\{0.35s_j\}},
\qquad
s_i=\frac{v^\top x_i}{\sqrt d},\quad v\sim N(0,I_d)
$$

is drawn once from the canonical seed.  Thus

$$
r_i=M p_R(i),\qquad
\theta_R=\sum_{i=1}^M p_R(i)\{Y_i(1)-Y_i(0)\}
$$

are exact finite sums.  RCT treatment is newly randomized with a known bounded
propensity depending on $X$; OBS treatment depends on $X$ and fixed latent
$U_i$.  Samples use replacement, and exactly $n_R$ RCT and $N_O$ OBS draws are
made per replication.

This benchmark validates behavior on one empirical covariate support.  It is
not an analysis of Project STAR's original intervention, and it is not evidence
about education policy effects.

## WHI CT+OS gated-access execution specification

The official source is the [BioLINCC WHI-CTOS study
page](https://biolincc.nhlbi.nih.gov/studies/whi_ctos/) and its [data preparation
guide](https://biolincc.nhlbi.nih.gov/media/studies/whi_ctos/WHI_CTOS_Data_Guide.pdf).
BioLINCC labels WHI-CTOS an Open BioLINCC Study, but downloading it requires a
registered account and an approved request.  The current run performs no login,
request, download, or analysis.

The script exposes a fail-closed `whi-check` command.  It accepts a
user-supplied, locally authorized extracted data directory outside the
workspace, inventories file names and column headers, and checks a user-supplied
semantic crosswalk.  It writes no participant-level data.  Without gated data
and a completed crosswalk it exits as blocked, which is the expected current
state.

The primary candidate analysis is:

- population: women with an intact uterus who satisfy harmonized baseline
  eligibility;
- RCT: the estrogen-plus-progestin hormone-therapy trial;
- treatment: randomized active therapy versus placebo, analyzed by intention
  to treat;
- OBS index: baseline eligible current combined estrogen-plus-progestin users
  versus nonusers, with an explicit new-user/index-date sensitivity analysis;
- outcome: centrally adjudicated coronary heart disease by five years;
- target: the five-year ITT risk difference in the RCT-eligible CT population;
- baseline covariates: age, race/ethnicity, region, education, income,
  hysterectomy/uterus status, smoking, BMI, blood pressure, diabetes, prior
  cardiovascular history, medication use, reproductive history, and baseline
  hormone history, all measured before the index date.

Execution gates are: verified CT/OS membership; intact-uterus and E+P arm
mapping; common time zero; harmonized inclusion/exclusion criteria; endpoint
adjudication and event-date mapping; loss-to-follow-up/death accounting;
positivity and risk-set overlap; absence of direct identifiers; and a signed
semantic crosswalk.  A five-year binary endpoint with censoring is not covered
by the current complete-outcome ATE theorem.  Analysis remains blocked until an
IPCW/survival extension or a defensible complete-outcome estimand is approved.

## Outputs, expected shapes, and QA

At $B=20$, the replication CSV has

$$
3\times2\times2\times20\times4
+1\times2\times2\times20\times4=1280
$$

rows.  The summary CSV has 64 rows.  It reports bias, empirical variance, RMSE,
AIPW- and Full-relative ratios, paired mean and squared-error differences, and
the paired centered-squared-contribution diagnostic whose mean equals the Full
minus method empirical variance difference.  The figure has benchmark rows and
columns for bias with $\pm2$ Monte Carlo standard errors, empirical variance
relative to AIPW, and RMSE relative to AIPW.

Required gates are:

1. Python syntax compilation and a $B=1$ smoke run;
2. exact 1280/64 row, four-method, cell, and per-cell replication counts;
3. exact unique sample accounting and no redraw or dropped replication;
4. deterministic rerun hashes;
5. all numeric CSV values finite and no missing values;
6. analytic/finite-sum target validation and known-r normalization checks;
7. restricted and Full calibration-objective assertions, coefficient,
   density-ratio, arm-count, and fallback diagnostics present;
8. visual inspection of the PNG;
9. scoped `git diff --check` and confirmation that only the five approved new
   paths changed.

No dependency is added, no raw third-party data are stored, and no file is
staged or committed.  The canonical memo and all earlier pilot artifacts remain
untouched.
