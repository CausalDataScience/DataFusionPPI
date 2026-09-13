# Implementation Handoff: Four Experiments for DataFusionPPI

Version 3, 2026-09-10.  This replaces version 2 of 2026-09-10; the changes are
listed in Section 10.  This document is the single authority for the
experiments.  It is self-contained: it restates the parts of the manuscript
that the code must implement, describes the existing code and where it departs
from the manuscript, gives every estimator as a formula, specifies the data
with pinned hashes, sizes the designs and the compute, and lists the
acceptance tests.  No full run may start before the Agile gate in Section 5 has
a recorded decision.

Reading order: Section 1 (what the paper claims), Section 2 (repository and
conventions), Section 3 (building blocks), Section 4 (baselines), Section 5
(Agile gate), Section 6 (the four experiments), Section 7 (data), Section 8
(compute budget and schedule), Section 9 (QA, deliverables, and how the results
enter the paper), Section 10 (revision history).  The manuscript is
`manuscript/main.pdf`;
theorem and equation numbers below refer to the build of 2026-09-09.

## 1. What the paper claims and what the experiments must show

### 1.1 Setting

Two data sources.  A randomized controlled trial (RCT, source $R$) with
covariates $X$, a binary treatment $A$ assigned with a known propensity
$e(X)=P_R(A=1\mid X)$ bounded away from 0 and 1, and an outcome $Y$.  An
observational study (OBS, source $O$) with the same variables, whose treatment
assignment may depend on unmeasured confounders.  The target is the RCT
population: the conditional average treatment effect
$\tau(x)=\mathbb E_R\{Y(1)-Y(0)\mid X=x\}$ and its mean
$\theta=\mathbb E_R\{\tau(X)\}$.

ATE uses three independent roles per source (Assumption 1): a nuisance sample
$D_S^{\mathrm{nuis}}$ for outcome models and the density ratio, a tuning sample
$D_S^{\mathrm{tune}}$ for coefficients, and an evaluation sample
$D_S^{\mathrm{eval}}$ for the final estimate.  CATE uses four independent roles:
nuisance, tuning, selection $D_S^{\mathrm{select}}$, and reporting
$D_S^{\mathrm{report}}$.  The candidate is chosen only on the selection role;
its final validation score and all reported performance comparisons are
computed only on the reporting role.  The manuscript's $n,N$ denote the role
used by the displayed estimator or score.  We write $\mathbb P_S^{\mathrm{role}}$
when the distinction matters.

Covariate shift.  Either the two covariate laws coincide (Assumption 3, ratio
$r_0\equiv1$) or they differ with $r_0=dP_R^X/dP_O^X$ (Assumption 4).  A
candidate ratio $r$ transports OBS covariate averages to the RCT population;
$r=r_0$ is exact.

### 1.2 The two channels

Fitted outcome regressions $\widehat\mu_R=(\widehat\mu_{R,0},\widehat\mu_{R,1})$
from the RCT nuisance sample and
$\widehat\mu_O=(\widehat\mu_{O,0},\widehat\mu_{O,1})$ from the OBS nuisance
sample.  The OBS treatment-effect prediction is
$\widehat g(x)=\widehat\mu_{O,1}(x)-\widehat\mu_{O,0}(x)$.  For a pair
$q=(q_0,q_1)$ the AIPW score is

$$\varphi(V;q)=q_1(X)-q_0(X)+\frac{A}{e(X)}\{Y-q_1(X)\}-\frac{1-A}{1-e(X)}\{Y-q_0(X)\},$$

which satisfies $\mathbb E_R\{\varphi(V;q)\mid X\}=\tau(X)$ for every fixed $q$
(Equation 6).  The fused regression is
$\widehat\mu_\lambda=(1-\lambda)\widehat\mu_R+\lambda\widehat\mu_O$, the fused
pseudo-outcome is $Z_\lambda=\varphi(V;\widehat\mu_\lambda)=Z_0+\lambda\Delta$
with $\Delta=Z_1-Z_0$.

- Channel $\lambda$ changes the regression inside the pseudo-outcome.
- Channel $\omega$ adds the control variate
  $\omega\{\mathbb P_{O,N}(r\widehat g)-\mathbb P_{R,n}\widehat g\}$, whose mean
  is zero when $r=r_0$.

### 1.3 Results the experiments check

ATE (Section 3 of the manuscript).

- Theorem 1: for fixed $(\lambda,\omega)$ and $r=r_0$ the estimator is
  conditionally unbiased; for a fixed candidate $\widetilde r$ the mean drift
  is $\omega\,\mathbb E_O[\{\widetilde r(X)-r_0(X)\}\widehat g(X)]$; and

  $$\operatorname{Var}\{\widehat\theta_{r_0}(\lambda,\omega)\}=\frac1n\operatorname{Var}_R\{Z_\lambda-\omega\widehat g(X)\}+\frac{\omega^2}{N}\operatorname{Var}_O\{r_0(X)\widehat g(X)\}.$$

- Definition 5 and Theorem 2: with $A=\operatorname{Var}_R(\Delta)$,
  $B=\operatorname{Var}_R(\widehat g)+\frac nN\operatorname{Var}_O(r_0\widehat g)$,
  $C=\operatorname{Cov}_R(Z_0,\Delta)$, $D=\operatorname{Cov}_R(Z_0,\widehat g)$,
  the variance is
  $\frac1n\{\operatorname{Var}_R(Z_0)+2C\lambda-2D\omega+A\lambda^2+B\omega^2\}$
  and the unconstrained oracle is $\lambda^\star=-C/A$, $\omega^\star=D/B$.
- Theorem 3: $n\operatorname{Var}$ equals $\sigma_0^2$ for RCT-only,
  $\sigma_0^2-D^2/B$ for $\omega$-only, and $\sigma_0^2-D^2/B-C^2/A$ for the
  joint oracle, with $\sigma_0^2=\operatorname{Var}_R(Z_0)$.
- Algorithm 1 and Theorem 4: the plug-in rule with the evaluation-sample
  variance estimator gives a valid normal approximation under independent
  roles.
- Theorem 5 and Corollary 5.1: the same with an estimated ratio and an
  explicit drift.

CATE (Section 4).

The manuscript calls the sample used by Algorithm 2 the evaluation sample.
Version 3 calls that role the **selection** sample and adds a new independent
**reporting** sample; this prevents performance reporting from reusing the
sample that selected the candidate.

- Definition 15: the DRF and RF losses.
- Algorithm 2, Definition 18, Theorem 6: fit candidates on the tuning samples
  for every grid pair, score them on the evaluation samples, select the
  minimizer.  The selected candidate is within
  $2\varepsilon_{\mathrm{eval}}(\alpha)$ of the best candidate on the grid plus
  twice the largest transport error of $r$.
- Proposition 10, Theorem 7, Corollary 7.1: for a linear sieve the candidates
  have closed form, the risk splits into approximation plus a variance term
  $\mathcal J_j(\lambda,\omega)$, and for DRF that term is a separated
  quadratic in $(\lambda,\omega)$ with the same shape as the ATE variance.

Every experiment reports the RCT-only estimator as the reference; report gains
as ratios when the reference metric has a stable, nonzero denominator and as
paired differences otherwise.

The empirical program has four labels, which must not be conflated in tables
or prose.

- **Theorem check:** a prespecified DGP satisfying the theorem's assumptions,
  with the theorem's claimed identity or bound as the target.
- **Diagnostic:** a numerical check of a mechanism, moment, overlap quantity,
  or implementation identity; it does not validate a theorem.
- **Empirical extension:** a useful estimator or learner outside the exact
  scope of the theorem being cited.
- **Real-data stress test:** an analysis of STAR or NSW that probes behavior
  under an empirical source law.  It does not establish causal identification
  unless the separate causal assumptions are justified.

## 2. Repository, environment, and conventions

Canonical project root: `research/papers/DataFusionPPI/` in the PAIOS
workspace.  The publication mirror is the private GitHub repository
`CausalDataScience/DataFusionPPI`, checked out at
`/Users/yonghanjung/Dropbox/Personal/Research/Code/DataFusionPPI`.  Work in the
canonical root; the author syncs the mirror.

Layout: `code/` for scripts, `materials/` for protocols, replication CSV files,
summary CSV files, and PNG figures, `manuscript/` for LaTeX, `memo/` for notes.

Environment used for the existing artifacts: Python 3 with numpy 2.4,
pandas 3.0, scikit-learn 1.9, matplotlib; torch 2.13 is available for the
neural basis.  Pin what you use in `code/requirements.txt`.

Existing scripts.

- `code/ssem_ate_pilot.py`: the SCM families (`SCMParameters`,
  `generate_scm_parameters(scm_id)` for `scm_id` in 1, 2, 3), the true
  $\tau(x)$ (`tau_function`), the true $\theta$ (`true_theta`), the RCT and OBS
  propensities, the sampling (`sample_dataset`), and the oracle ratio for the
  shifted regime (`oracle_ratio`, a Gaussian likelihood ratio because the RCT
  covariate law is the OBS law with a mean shift).
- `code/total_budget_nested_benchmark.py`: the current ATE benchmark.  It
  imports the pilot module as `legacy`, adds the STAR real-$X$ family
  (`load_star_support`, `star_sample`), fits ridge outcome models on a fixed
  feature map (`fit_outcome`), fits a logistic domain classifier (`fit_ratio`),
  computes pseudo-outcomes (`pseudo`), tuning moments (`tuning_moments`) and
  coefficients (`tune_all`), and writes the three output artifacts
  (`write_outputs`, `summarize`, `make_figure`).

Implementation status as of 2026-09-10: these are legacy/pilot scripts.  They
do not yet implement the v3 gate, four-CATE-role split, bounded theorem-check
DGP, frozen CATE-2 source law, covariance-aware comparators, adaptive oracle
MCSE rule, or coverage extension.  No existing CSV or PNG is v3 evidence.

Seeds: base seed 190602, per-cell seeds from `stable_seed(...)` and
`rng_for(...)`.  Keep this scheme so every replication is reproducible from the
design cell and the replication index.

Outputs: one replication-level CSV with one row per design cell, replication,
and estimator; one summary CSV with the mean and the Monte Carlo standard error
over replications; one composite PNG per study.  All go under `materials/` with
a descriptive stem.  Every row carries the design cell fields, the estimator
name, the metric name, and the seed.

### 2.1 Where the existing script departs from the manuscript

These differences matter because the manuscript's guarantees require
independent roles and Algorithm 1 exactly as stated.

1. Sample reuse.  The script runs five outer folds and tunes coefficients by a
   two-fold inner cross-fit inside each fold.  Theorem 4 assumes independent
   nuisance, tuning, and evaluation samples and states that cross-fitted reuse
   is not covered.  The new ATE experiments use the honest three-way split;
   CATE additionally separates selection from reporting.  Cross-fitting is
   retained only as a comparison panel in ATE-1.
2. Coefficient rule.  The script minimizes the full empirical quadratic,
   including the cross term, over the box $[-2,2]^2$.  Algorithm 1 uses the
   separated plug-in rule with thresholds and projection onto a prespecified
   rectangle; implement Algorithm 1 exactly, with the rectangle $[0,1]^2$.
3. Variance.  The script records a non-exact variance proxy.  The experiments
   compute $\widehat V$ on the evaluation samples and report interval coverage.
4. Ratio clipping.  Clipping to $[0.05,20]$ changes the ratio target.  Report
   the clip fraction, run the unclipped ratio where it is finite, and treat
   clipping as a labeled variant.
5. Tuning-size factor.  The script uses the training-size ratio inside
   $\widehat B$.  The manuscript uses $n/N$ with the evaluation sizes, because
   the quantity being minimized is the evaluation-sample variance.
6. Version-3 protocol.  The scripts do not enforce the gate lock, independent
   CATE selection/reporting, fixed empirical CATE-2 source law, or the expanded
   assertions of Section 9.1.  Implement and pass G0 before treating any pilot
   artifact as evidence for the v3 experiment.

## 3. Building blocks

All quantities are computed conditional on the fitted objects, which is how the
theorems are stated.

### 3.1 Roles, nuisances, and pseudo-outcomes

1. For ATE, split each source into nuisance, tuning, and evaluation roles with
   fractions $0.4/0.3/0.3$ for the RCT and $0.6/0.2/0.2$ for the OBS sample,
   unless a study says otherwise.  For CATE, construct mutually independent
   nuisance, tuning, selection, and reporting samples as specified by the
   study.  Use distinct, named RNG namespaces for every role; never obtain a
   role by subdividing, recycling, or conditionally redrawing another role.
   Required namespaces are `design`, `R_nuis`, `O_nuis`, `R_tune`, `O_tune`,
   `R_select`, `O_select`, `R_report`, `O_report`, `truth`, `oracle`, and
   `gate_bootstrap`; a study may omit an unused namespace but may not alias two.
2. Fit $\widehat\mu_R$ on $D_R^{\mathrm{nuis}}$ and $\widehat\mu_O$ on
   $D_O^{\mathrm{nuis}}$ with `fit_outcome`; the misspecification knob in ATE-1
   replaces the feature map by the raw covariates.  Fit the density ratio on
   the two nuisance samples.
3. On any sample with known $e$, compute $Z_0=\varphi(V;\widehat\mu_R)$,
   $Z_1=\varphi(V;\widehat\mu_O)$, $\Delta=Z_1-Z_0$, and $\widehat g$.

Propensity trimming.  Both the AIPW score and the RF learner divide by factors
that vanish as $e$ approaches 0 or 1.  Restrict every analysis to units with
$e(X)\in[0.15,0.85]$ and report the share dropped.  The SCM families already
clip the RCT propensity to $[0.2,0.8]$, so nothing is dropped there.  For the
STAR data see Section 7.2.

### 3.2 ATE estimator, Algorithm 1, and inference

Estimator for fixed $(\lambda,\omega)$ and ratio $r$ on the evaluation samples:

$$\widehat\theta_r(\lambda,\omega)=\mathbb P_{R,n}Z_\lambda+\omega\{\mathbb P_{O,N}(r\widehat g)-\mathbb P_{R,n}\widehat g\}.$$

Tuning moments on the tuning samples:
$\widehat A=\widehat{\operatorname{Var}}_R(\Delta)$,
$\widehat B=\widehat{\operatorname{Var}}_R(\widehat g)+\frac nN\widehat{\operatorname{Var}}_O(\widehat r\widehat g)$
with $n=n_R^{\mathrm{eval}}$ and $N=n_O^{\mathrm{eval}}$,
$\widehat C=\widehat{\operatorname{Cov}}_R(Z_0,\Delta)$,
$\widehat D=\widehat{\operatorname{Cov}}_R(Z_0,\widehat g)$.

Plug-in rule with rectangle $[0,1]^2$ and thresholds
$t_A=t_B=10^{-3}\widehat{\operatorname{Var}}_R(Z_0)$:

$$\widehat\lambda=\begin{cases}\operatorname{Proj}_{[0,1]}(-\widehat C/\widehat A)&\widehat A>t_A\\0&\text{otherwise}\end{cases},\qquad
\widehat\omega=\begin{cases}\operatorname{Proj}_{[0,1]}(\widehat D/\widehat B)&\widehat B>t_B\\0&\text{otherwise.}\end{cases}$$

Restricted estimators: RCT-only fixes $(0,0)$; $\lambda$-only fixes $\omega=0$;
$\omega$-only fixes $\lambda=0$; joint uses both.

Variance estimator and interval on the evaluation samples:
$\psi_{R,i}=Z_{\widehat\lambda}(V_{R,i})-\widehat\omega\widehat g(X_{R,i})$,
$\psi_{O,j}=\widehat\omega\,\widehat r(X_{O,j})\widehat g(X_{O,j})$,
$\widehat V=\widehat{\operatorname{Var}}_{R,\mathrm{eval}}(\psi_R)/n+\widehat{\operatorname{Var}}_{O,\mathrm{eval}}(\psi_O)/N$
with unbiased sample variances, and the interval
$\widehat\theta\pm z_{0.975}\sqrt{\widehat V}$.  Coverage is the share of
replications whose interval contains $\theta$.

Oracle coefficients, computed once per replication.  Hold that replication's
fitted nuisances fixed and start with $2\times10^4$ fresh RCT and OBS draws
from the design.  Compute $(A,B,C,D)$ of Definition 5 with the same $n/N$, set
$(\lambda^\star,\omega^\star)=(-C/A,D/B)$ projected onto $[0,1]^2$, and compute
the right-hand side of the Theorem 1 variance formula.  For every moment report
its Monte Carlo standard error (MCSE), not a generic relative-error claim.  Let
$s_0^2$ be the concurrently simulated variance of $Z_0$.  The absolute stopping
tolerances are $0.01s_0^2$ for $A$ and $B$, and $0.005s_0^2$ for covariances
$C$ and $D$; the latter remain meaningful when a covariance is near zero and a
relative MCSE is undefined.  Double the draw count to $4\times10^4$ and then
$8\times10^4$ while any tolerance fails; use a final increment to at most
$10^5$ draws.  A still-unstable oracle after $10^5$ is reported as Monte
Carlo-inconclusive rather than silently redrawn.

### 3.3 Density-ratio estimators

Classifier (Definition 11): pool the nuisance covariates, label RCT as 1, fit a
logistic regression on the standardized covariates, and set

$$\widehat r_{\mathrm{CLS}}(x)=\frac{1-\widehat\pi}{\widehat\pi}\cdot\frac{\widehat s(x)}{1-\widehat s(x)},\qquad \widehat\pi=\frac{n_R^{\mathrm{nuis}}}{n_R^{\mathrm{nuis}}+n_O^{\mathrm{nuis}}}.$$

Report the normalization error $\mathbb P_O^{\mathrm{nuis}}(\widehat r)-1$, the
effective sample size $(\sum\widehat r)^2/\sum\widehat r^2$ on the OBS
evaluation sample, and the clip fraction when clipping is used.

Calibrated balancing ratio (Definition 12): with base ratio
$\widehat r_{\mathrm{base}}=\widehat r_{\mathrm{CLS}}$ and a feature map $\phi$
without intercept, solve the convex problem

$$\widehat\xi\in\arg\min_{\xi}\left[\log\mathbb P_O^{\mathrm{nuis}}\{\widehat r_{\mathrm{base}}(X)e^{\xi^\top\phi(X)}\}-\xi^\top\mathbb P_R^{\mathrm{nuis}}\{\phi(X)\}\right]$$

by L-BFGS.  The gradient is the balance residual
$\mathbb P_O^{\mathrm{nuis}}(\widehat r_{\mathrm{BAL}}\phi)-\mathbb P_R^{\mathrm{nuis}}\phi$,
where

$$\widehat r_{\mathrm{BAL}}=\frac{\widehat r_{\mathrm{base}}e^{\widehat\xi^\top\phi}}{\mathbb P_O^{\mathrm{nuis}}\{\widehat r_{\mathrm{base}}e^{\widehat\xi^\top\phi}\}}.$$

Verify Proposition 6 numerically: every component of the balance residual is
below $10^{-6}$.  Two feature maps are used: $\phi$ equal to the standardized
raw covariates, and $\phi$ equal to those plus $\widehat g$.  When the RCT
feature mean lies outside the convex hull of the OBS feature vectors there is
no interior minimizer; detect this by a residual that does not vanish and
record the cell as infeasible.

Oracle ratio: `oracle_ratio` for the SCM families, the support probabilities
for the STAR real-$X$ family, and $r_0\equiv1$ by construction for CATE-2.

### 3.4 CATE losses, sieve fit, validation, and metrics

Losses (Definition 15) on the tuning samples, with
$G_\zeta=(\widehat g-\zeta)^2$ and
$\chi(A,X)=\{A-e(X)\}^2/[e(X)\{1-e(X)\}]$:

$$\widehat L_{\mathrm{DRF}}(\zeta;\lambda,\omega,r)=\mathbb P_R^{\mathrm{tune}}\{Z_\lambda-\zeta(X)\}^2+\omega\left[\mathbb P_O^{\mathrm{tune}}\{r(X)G_\zeta(X)\}-\mathbb P_R^{\mathrm{tune}}\{G_\zeta(X)\}\right],$$

$$\widehat L_{\mathrm{RF}}(\zeta;\lambda,\omega,r)=\mathbb P_R^{\mathrm{tune}}\left[\frac{[Y-m_\lambda(X)-\{A-e(X)\}\zeta(X)]^2}{e(X)\{1-e(X)\}}\right]+\omega\left[\mathbb P_O^{\mathrm{tune}}\{r(X)G_\zeta(X)\}-\mathbb P_R^{\mathrm{tune}}\{\chi(A,X)G_\zeta(X)\}\right],$$

where $m_S(x)=e(x)\widehat\mu_{S,1}(x)+\{1-e(x)\}\widehat\mu_{S,0}(x)$ and
$m_\lambda=(1-\lambda)m_R+\lambda m_O$.

Sieve fit (Proposition 10).  For a feature map $b(x)\in\mathbb R^p$ and ridge
$\rho\geq0$, with $w_{\mathrm{DRF}}=1$,
$\widetilde Z_{\mathrm{DRF},\lambda}=Z_\lambda$, $w_{\mathrm{RF}}=\chi(A,X)$,
$\widetilde Z_{\mathrm{RF},\lambda}=\{Y-m_\lambda(X)\}/\{A-e(X)\}$, solve

$$\{\widehat H_j(\omega)+\rho I\}\beta=\mathbb P_R^{\mathrm{tune}}\left[w_jb\{\widetilde Z_{j,\lambda}-\omega\widehat g\}\right]+\omega\,\mathbb P_O^{\mathrm{tune}}(r\,b\,\widehat g),$$

$$\widehat H_j(\omega)=(1-\omega)\mathbb P_R^{\mathrm{tune}}(w_jbb^\top)+\omega\,\mathbb P_O^{\mathrm{tune}}(r\,bb^\top),$$

and set $\widehat\zeta_j(x;\lambda,\omega)=b(x)^\top\widehat\beta_j$.  For fixed
$\omega$ factorize once and reuse it for every $\lambda$, because the
right-hand side is affine in $\lambda$.  The RF weight
$\widetilde Z_{\mathrm{RF},\lambda}$ divides by $A-e(X)$, which the trimming
rule of Section 3.1 keeps at least $0.15$ in absolute value.

Feature maps.

- Splines: additive cubic B-splines with $K\in\{3,5\}$ interior knots per
  continuous covariate at empirical quantiles of the RCT tuning sample, dummies
  for categorical covariates, and an intercept.
- Frozen neural basis: on the OBS nuisance sample train a network with a shared
  trunk $h(x)\in\mathbb R^{64}$ of two hidden layers with 64 ReLU units and two
  linear heads $\widehat\mu_{O,a}(x)=\theta_a^\top h(x)+c_a$, mean squared
  error, Adam, early stopping on a 10 percent holdout of that nuisance sample.
  Freeze the trunk, set $b(x)=(1,h(x))$, and take $\widehat g$ from the heads.
  Estimate the centering and scaling constants for $h$ from nuisance samples
  only and freeze them before any tuning, selection, or reporting observation
  is used.  Train once per replication and reuse the frozen trunk and scaling
  constants for every learner and every grid point; the basis must not depend
  on the tuning, selection, or reporting samples.

Candidate grid: $\mathcal G=\{0,0.25,0.5,0.75,1\}^2$, so $M=25$.  For the
neural basis the candidate is the triple
$s=(\lambda,\omega,\rho)$ with
$\rho\in\{10^{-2},10^{-1},1\}$, so the joint selector minimizes over all 75
candidates.  Define its RCT-only comparator by

$$\widehat\rho_R\in\arg\min_{\rho}
\widehat C^{\mathrm{select}}(0,0,\rho),$$

breaking ties by the smallest ridge-grid index, and freeze $\widehat\rho_R$
before reporting.  For splines take $s=(\lambda,\omega,0)$ and suppress
$\rho=0$ from the notation.
The fixed-ridge neural system is an empirical candidate for the Theorem 6
finite-grid selection rule and the displayed normal equation.  It is not a
Theorem 7 or Corollary 7.1 check.  Those checks use only $\rho=0$, a
well-conditioned spline basis, and $p<n_R^{\mathrm{tune}}$.

Selection score (Definition 18) on the selection samples.  Write $\omega_s$ for
the $\omega$ component of $s$ and
$G_s=\{\widehat g(X)-\widehat\zeta(X;s)\}^2$:

$$\widehat C^{\mathrm{select}}(s)=
\mathbb P_R^{\mathrm{select}}\left[\{Z_0-\widehat\zeta(X;s)\}^2-
\omega_sG_s\right]+
\mathbb P_O^{\mathrm{select}}\{\omega_s r(X)G_s\}.$$

Select the minimizer over the candidate set, breaking ties by the smallest
index.  Freeze that choice, then evaluate its score on the independent
reporting samples.  The minimum selection score is an optimization statistic,
not an unbiased performance estimate for the selected candidate.

Radius (Theorem 6).  Test the Hoeffding radius only in a prespecified bounded
theorem-check DGP whose deterministic bounds $B_0$ and $\bar r_0$ are part of
the DGP before any replication:

$$\varepsilon_{\mathrm{select}}(\alpha)=2B_0^2\sqrt{2\log(4M/\alpha)}\left\{\frac{2}{\sqrt{n_R^{\mathrm{select}}}}+\frac{\bar r_0}{\sqrt{n_O^{\mathrm{select}}}}\right\},\qquad \alpha=0.05.$$

Evaluation-, selection-, or reporting-sample maxima may not be substituted for
$B_0$ or $\bar r_0$.  The radius is `N/A` for Gaussian DGPs and STAR; a
descriptive empirical range may be reported separately without theorem-check
language.

Metrics.

- CATE risk $\mathcal R(\zeta)=\mathbb E_R\{\zeta(X)-\tau(X)\}^2$ on an
  independent test sample of $10^5$ RCT-population covariates with the true
  $\tau$.  When $\tau$ is unknown (CATE-2), freeze the selected candidate and
  use the reporting-score difference.  For splines this is
  $\widehat C^{\mathrm{report}}(\lambda,\omega)-\widehat C^{\mathrm{report}}(0,0)$;
  for neural candidates it is
  $\widehat C^{\mathrm{report}}(\widehat s)-
  \widehat C^{\mathrm{report}}(0,0,\widehat\rho_R)$.
  Conditional on a fixed candidate and the source-law assumptions, this is an
  unbiased risk-difference score; the selected minimum itself is not.
- Regret $\mathcal R(\widehat s)-\min_{s\in\mathcal S}\mathcal R(s)$.  Only in
  the bounded theorem-check DGP, also report its ratio to
  $2\varepsilon_{\mathrm{select}}(0.05)$ and the share of replications with
  regret at most that radius.  Also report

  $$\texttt{bound\_vacuity\_ratio}=
  \frac{2\varepsilon_{\mathrm{select}}(0.05)}
  {\max_s\mathcal R(s)-\min_s\mathcal R(s)},$$

  using `Inf` when the denominator is zero.  A value at least one means the
  valid bound covers the entire observed candidate-risk range; satisfying such
  a radius is not substantive evidence of good performance.
- Leading risk: $\Gamma_R=\mathbb E_R(bb^\top)$,
  $\beta_p=\Gamma_R^{-1}\mathbb E_R(b\tau)$, $a_p^2=\mathcal R(b^\top\beta_p)$ on
  the test sample;
  $\psi_{R,j}(\lambda,\omega)=w_jb\{\widetilde Z_{j,\lambda}-\omega\widehat g-(1-\omega)\zeta_p\}$
  and $\psi_O(\omega)=\omega r_0b(\widehat g-\zeta_p)$ on fresh draws;

  $$\mathcal J_j=\frac{\operatorname{tr}[\Gamma_R^{-1}\operatorname{Var}_R(\psi_{R,j})]}{n_R^{\mathrm{tune}}}+\frac{\operatorname{tr}[\Gamma_R^{-1}\operatorname{Var}_O(\psi_O)]}{n_O^{\mathrm{tune}}};$$

  compare the Monte Carlo mean of $\mathcal R(\widehat\zeta_j)-a_p^2$ with
  $\mathcal J_j$ at $(\lambda,\omega)\in\{(0,0),(0.5,0),(0,0.5),(0.5,0.5)\}$.
- Oracle sieve coefficients (Corollary 7.1, DRF): $A_p,B_p,C_p,D_p$ from the
  same fresh draws with $\zeta_p$ known, giving
  $(\lambda_p^\star,\omega_p^\star)=(-C_p/A_p,D_p/B_p)$.  A feasible version
  replaces $\zeta_p$ by the RCT-only fit and the expectations by tuning-sample
  averages.  Report both against the selected coefficients.

## 4. Baselines

Each literature-inspired baseline below is a mechanism-matched heuristic in
this paper's notation, not the original authors' exact estimator, unless the
subsection explicitly says otherwise.  Preserve that label in every result
table, figure, and memo.  All baselines use the same roles, nuisance fits, and
evaluation or reporting samples as the fusion estimators so that the comparison
isolates the combination rule.

### 4.1 Semi-supervised ATE (ATE-1)

The idea is to use the OBS covariates only, as an auxiliary sample for the
covariate distribution, with no OBS outcome information.  Let
$\widehat g_R=\widehat\mu_{R,1}-\widehat\mu_{R,0}$ be the RCT-fitted contrast.
Define

$$\widehat\theta_{\mathrm{SS}}=\mathbb P_{R,n}Z_0+\widehat\eta\left\{\mathbb P_{O,N}(r\,\widehat g_R)-\mathbb P_{R,n}\widehat g_R\right\},\qquad
\widehat\eta=\frac{\widehat{\operatorname{Cov}}_R(Z_0,\widehat g_R)}{\widehat{\operatorname{Var}}_R(\widehat g_R)+\frac nN\widehat{\operatorname{Var}}_O(r\widehat g_R)},$$

with $\widehat\eta$ computed on the tuning samples and projected onto $[0,1]$.
This is the $\omega$ channel with $\widehat g_R$ in place of $\widehat g$, which
is exactly what distinguishes semi-supervised variance reduction from
prediction-powered fusion: it cannot import OBS outcome information.

### 4.2 Shrinkage combination (ATE-1)

On the tuning roles let $\widehat\theta_R^{\mathrm{tune}}=
\mathbb P_R^{\mathrm{tune}}Z_0$ and let

$$\widehat\theta_O^{\mathrm{tune}}=\mathbb P_R^{\mathrm{tune}}\widehat g+
\mathbb P_O^{\mathrm{tune}}\left[r\left\{\frac{A}{\widehat e_O}(Y-\widehat\mu_{O,1})-
\frac{1-A}{1-\widehat e_O}(Y-\widehat\mu_{O,0})\right\}\right]$$

be the transported OBS AIPW estimator that treats the OBS study as
unconfounded.  Because both estimates use the same RCT covariate sample, let
$u_O=A(Y-\widehat\mu_{O,1})/\widehat e_O-
(1-A)(Y-\widehat\mu_{O,0})/(1-\widehat e_O)$.  Conditional on nuisances, use

$$\widehat V_R=\frac{\widehat{\operatorname{Var}}_R(Z_0)}{n_R},\quad
\widehat V_O=\frac{\widehat{\operatorname{Var}}_R(\widehat g)}{n_R}+
\frac{\widehat{\operatorname{Var}}_O(r u_O)}{n_O},\quad
\widehat C_{RO}=\frac{\widehat{\operatorname{Cov}}_R(Z_0,\widehat g)}{n_R}.$$

The last term is the shared-RCT covariate covariance.  Estimate all three on
the tuning roles using the final evaluation-role values of $n_R,n_O$.  Every
quantity below, including both $\widehat\theta$ values, comes from the tuning
roles:

$$\widehat b^2=\max\{(\widehat\theta_R^{\mathrm{tune}}-
\widehat\theta_O^{\mathrm{tune}})^2-
(\widehat V_R+\widehat V_O-2\widehat C_{RO}),\;0\},$$

$$\widehat\gamma=\operatorname{Proj}_{[0,1]}
\frac{\widehat V_R-\widehat C_{RO}}
{\widehat V_R+\widehat V_O-2\widehat C_{RO}+\widehat b^2},\qquad
\widehat\theta_{\mathrm{SHR}}^{\mathrm{eval}}=
(1-\widehat\gamma)\widehat\theta_R^{\mathrm{eval}}+
\widehat\gamma\widehat\theta_O^{\mathrm{eval}}.$$

Freeze $\widehat\gamma$, then recompute only $\widehat\theta_R$ and
$\widehat\theta_O$ on the ATE evaluation roles for the displayed final
combination.

When the two estimates disagree by much more than their noise, $\widehat b^2$
is large, $\widehat\gamma$ is near zero, and the rule returns the RCT-only
estimator.

### 4.3 Adaptive combination (ATE-1)

The same two ingredients with a hard test in place of soft thresholding:

$$\widehat\gamma=\begin{cases}
\operatorname{Proj}_{[0,1]}\dfrac{\widehat V_R-\widehat C_{RO}}
{\widehat V_R+\widehat V_O-2\widehat C_{RO}}
&\text{if }|\widehat\theta_R^{\mathrm{tune}}-
\widehat\theta_O^{\mathrm{tune}}|\leq
z_{0.975}\sqrt{\widehat V_R+\widehat V_O-2\widehat C_{RO}},\\[2mm]
0&\text{otherwise,}
\end{cases}$$

and $\widehat\theta_{\mathrm{ADA}}^{\mathrm{eval}}=
(1-\widehat\gamma)\widehat\theta_R^{\mathrm{eval}}+
\widehat\gamma\widehat\theta_O^{\mathrm{eval}}$.
The test decision and $\widehat\gamma$ are computed on tuning data and frozen
before the ATE evaluation sample is used; only the two final component
estimates are recomputed on evaluation data.  When no incompatibility is detected,
the rule uses the correlated-estimator minimum-variance weight; otherwise it
falls back to the trial.

### 4.4 Naive pooling (ATE-1 and CATE-1)

Fit the pooled outcome model and source-specific propensities on nuisance data.
On ATE evaluation data or CATE reporting data form the RCT AIPW pseudo-outcome
$Z_R$ using known $e$ and an
OBS pseudo-outcome $Z_O$ using fitted $\widehat e_O$, transporting the latter
by $r$.  The operational estimator is the literal pooled average

$$\widehat\theta_{\mathrm{NP}}=
\frac{n_R\mathbb P_R Z_R+n_O\mathbb P_O(rZ_O)}{n_R+n_O}.$$

Thus OBS outcome information enters the final average directly.  This is true
naive causal pooling: it treats the OBS study as a second trial and can be
biased under unmeasured confounding.  For CATE, fit the same sieve to the
concatenated RCT pseudo-outcomes and transported OBS pseudo-outcomes.

### 4.5 Experimental grounding (CATE-1 and CATE-2)

Two stages.  First fit an OBS CATE estimate $\widehat\tau_O$ on
$D_O^{\mathrm{tune}}$ by regressing the OBS AIPW pseudo-outcome on $X$ in the
same sieve class.  Then fit a correction $\widehat q$ in the same class by
regressing $Z_0-\widehat\tau_O(X)$ on $X$ using $D_R^{\mathrm{tune}}$.  Output
$\widehat\tau_O+\widehat q$.  The correction is what the trial contributes.

### 4.6 Domain-indicator pooling (CATE-1 and CATE-2)

Fit the DRF loss on the concatenation of the RCT and OBS pseudo-outcomes in the
same sieve class, with the basis augmented by a source indicator and its
interaction with the intercept.  This is the CATE analogue of naive pooling.

### 4.7 Bibliographic status

Zotero, the author's local library of 590 items, was checked on 2026-09-09.
Present: Hatt and others 2022 on representation learning with both sources,
Athey, Chetty, and Imbens 2020 on long-term outcomes, Oberst and others 2022 on
bias-robust integration.  Absent: semi-supervised average treatment effect
estimation, prognostic covariate adjustment, prediction-powered generalization,
and predictions as surrogates.  Semantic Scholar supplied the records below.

| baseline or related work | record | status |
| --- | --- | --- |
| Semi-supervised ATE | Cheng, Ananthakrishnan, and Cai, Biometrics, DOI 10.1111/biom.13298 | verified |
| Shrinkage combination | `rosenman2020shrinkage` | in `reference.bib` |
| Adaptive combination | `cheng2021adaptive` | in `reference.bib` |
| Experimental grounding | `kallus2018removing` | in `reference.bib` |
| Naive pooling, domain-indicator pooling | no citation needed | not applicable |
| Prediction-powered generalization of causal inferences | Demirel, Alaa, Philippakis, and others, ICML 2024, arXiv 2406.02873 | verified |
| Prediction-powered causal inferences | Cadei, Demirel, De Bartolomeis, and others, NeurIPS 2025, arXiv 2502.06343 | verified; read before writing Section 1 |
| Cross-prediction-powered inference | Zrnic and Candes, PNAS, DOI 10.1073/pnas.2322083121 | verified |
| Prognostic covariate adjustment | Schuler, Walsh, Hall, and others, The International Journal of Biostatistics, DOI 10.1515/ijb-2021-0072, arXiv 2012.09935 | verified |
| Predictions as surrogates | see Section 10 | pending |

Implement a baseline only after its record is verified, and let the author
confirm the final list before it enters the manuscript.

## 5. Agile staged gate before the full run

The immediate goal is to learn whether the implementation is valid and whether
the two-channel scientific story merits the full compute budget.  Run this gate
before any full study.  The panels are deliberately narrow:

| panel | cells | cap per cell | maximum replications |
| --- | ---: | ---: | ---: |
| ATE-1 | SCM 1 reference; SCM 1 with confounding 0/correct outcome model; SCM 1 with confounding 2/linear outcome model; STAR real-$X$ reference | 100 | 400 |
| ATE-2 synthetic | SCM 1, $n_R=100$, at the prespecified mild and strong shifts | 50 | 100 |
| CATE-1 | DRF, spline $K=3$, SCM-B reference and SCM-B confounding-2 cells | 50 | 100 |
| total | eight | | 600 |

Exclude the neural basis, NSW/PSID, CATE-2, and the full comparator set from
G1/G2 decision runs.  The gate estimators are RCT-only, $\lambda$-only,
$\omega$-only, joint, and oracle where defined.  G0 runs one replication per
listed cell as a smoke test.  G0 also runs two non-replication hard smoke tests:

1. execute the neural path end to end on the CATE-1 SCM-B reference cell; and
2. execute the CATE-2 source-law, DRF spline-$K=3$, selection, and reporting
   pipeline for primary mathematics, $\alpha=0.8$, and $n_R=200$.

These two checks emit assertion records but no performance-replication rows,
are not reusable prefixes, and therefore do not alter the 600-replication cap.
The CATE-2 hard assertions are the exact source probability, exact-cell
propensity, and role/RNG independence checks.  Directional movement of its
naive OBS contrast is a soft diagnostic, never an engineering failure.

After all G0 assertions pass, G1 runs ten replications per cell to check direction,
finite outputs, and runtime.  G2 extends only the cells needed for a decision,
incrementally, up to the table's caps.  G0, G1, and G2 counts are cumulative;
replication IDs are prefixes of $0,\ldots,\text{cap}-1$.  Do not automatically
fill every cap.

The prespecified scientific gate uses the ATE-1 confounding-0/correct-model cell
and the CATE-1 SCM-B reference cell as the expected-useful regimes.  For each,
let $L$ be MSE for ATE and true CATE risk for CATE and define the paired,
RCT-normalized incremental gain

$$G_\omega=\frac{\min\{\widehat L_{\mathrm{RCT}},
\widehat L_{\lambda\text{-only}}\}-\widehat L_{\mathrm{joint}}}
{\widehat L_{\mathrm{RCT}}}.$$

Use a replication-level paired bootstrap with namespace `gate_bootstrap` to
form a 90 percent interval.  Denote the minimum scientifically relevant gain by
$\delta=0.05$.  An ATE coverage loss is
unacceptable when the lower endpoint of the 90 percent interval for
$\operatorname{Coverage}_{\mathrm{RCT}}-
\operatorname{Coverage}_{\mathrm{joint}}$ exceeds 0.05.  These definitions,
bootstrap seed, and interval convention are frozen before G1.

The decision must be recorded in `materials/gate_v3_results.md` before a
full-run command is permitted, using exactly one of these outcomes:

- `GO_FULL_TWO_CHANNEL`: engineering and theorem checks pass, at least one
  expected-useful regime has a 90 percent lower endpoint strictly above
  $\delta=0.05$, and the coverage guard does not fail.
- `REFRAME`: validity checks pass, the relevant cells have reached their caps,
  and neither `GO_FULL_TWO_CHANNEL` nor `REDESIGN` applies.  This means there is
  insufficient evidence for the full two-channel claim; it does not assert
  that the incremental effect is zero.
- `REDESIGN`: an engineering assertion or theorem-check identity fails under
  its stated assumptions.  Diagnose and repair before further performance runs.
- `INCONCLUSIVE`: validity checks pass, `GO_FULL_TWO_CHANNEL` does not apply,
  and at least one relevant cell remains below its cap.  Extend only those
  inconclusive cells incrementally.  `INCONCLUSIVE` is not a terminal state at
  the cap; at the cap, valid runs that do not satisfy
  `GO_FULL_TWO_CHANNEL` terminate as `REFRAME`.

These cases are exhaustive: a hard validity/theorem/engineering failure gives
`REDESIGN`; absent such a failure, a satisfied scientific and coverage gate
gives `GO_FULL_TWO_CHANNEL`, any other valid state below the cap gives
`INCONCLUSIVE`, and every other valid state at the cap gives `REFRAME`.

The gate criterion and its numerical tolerance must be written to the results
memo before G1 outcomes are inspected.  A joint estimate that does not beat
RCT-only is never, by itself, evidence of a code failure.

Gate replications may become the prefix of the full study only when the cell is
already in the Section 6 design and code, protocol version, and seed mapping are
unchanged.  A stress corner outside the full grid is non-reusable and must be
labeled and counted separately.  Version 3 includes no such extra corner and
does not change the 94-cell full plan.  Any later interaction diagnostic is a
gate-only, non-reusable addition unless the protocol is explicitly revised.

## 6. The four experiments

The full plan is conditional on a Section 5 gate decision.  Performance cells
use 100 replications as screening, the RCT-only estimator as the reference, and
report the mean and Monte Carlo standard error of every metric.  The four
sentinel coverage cells specified below receive 400 additional replications
each after screening, for 500 total per sentinel cell.

The designs below vary one axis at a time around a reference configuration
rather than crossing every factor.  Version 1 crossed them, which produced 288
cells for ATE-1 alone; one axis at a time answers the same questions with 44.

### 6.1 ATE-1: when fusion helps, and whether Algorithm 1 finds it

Classification: SCM cells check the fixed-coefficient mean/variance identities
and oracle ordering of Theorems 1, 2, 2.1, and 3.  Interval behavior for
Algorithm 1 and Theorem 4 is a theorem check only where its regularity
conditions are verified; otherwise it is a calibration diagnostic.  STAR
real-$X$ is an empirical extension.

Data families: SCM 1, 2, 3 and STAR real-$X$; common covariate marginal only,
so $r_0\equiv1$ and no ratio is estimated.

Reference configuration: $n_R=100$, $N_O=5000$, moderate confounding, correct
OBS outcome model.  Three axes, 11 configurations per family, 44 cells.

| axis | configurations | new cells per family |
| --- | --- | ---: |
| RCT size | $n_R\in\{50,100,200,400\}$ at $N_O=5000$ | 4 |
| OBS size | $N_O\in\{1000,20000\}$ at $n_R=100$ | 2 |
| Prediction quality | confounding in $\{0,1,2\}$ crossed with outcome model in {correct, linear}, minus the reference | 5 |

The confounding level scales `obs_latent_coefficient` by 0, 1, or 2; level 0
removes the unmeasured confounder from the OBS assignment.  The outcome-model
knob replaces the existing feature map by a linear model on the raw covariates.

Estimators: RCT-only, $\lambda$-only, $\omega$-only, joint by Algorithm 1,
oracle-coefficient joint, and the four baselines of Sections 4.1 to 4.4.  Add
one panel that runs the joint estimator with the existing cross-fitted scheme
at the reference configuration only.

Metrics: RMSE ratio to RCT-only, bias, variance ratio; the exact-variance
check; coefficient recovery, meaning a scatter of
$(\widehat\lambda,\widehat\omega)$ against $(\lambda^\star,\omega^\star)$ with
the boundary-hit rate and the threshold-fallback rate; coverage and mean length
of the 95 percent interval for joint and RCT-only; the honest versus
cross-fitted panel.

Coverage sentinels are the reference configuration in each of SCM 1, SCM 2,
SCM 3, and STAR real-$X$.  Run 100 screening replications and then 400
additional replications in each sentinel, for 500 total.  At nominal coverage
0.95, the binomial MCSE is
$\sqrt{0.95(0.05)/500}=0.00975$.  For an exact-ratio, correctly specified
sentinel, the three-MCSE acceptance interval is
$0.95\pm3\sqrt{0.95(0.05)/500}=[0.921,0.979]$ after rounding.  Coverage outside
that interval is a no-go signal requiring diagnosis; it is not repaired by
dropping or redrawing replications.
Use replication IDs 0--99 for screening and 100--499 for the extension.

Outputs: `materials/ate1_when_fusion_helps_replications.csv`,
`materials/ate1_when_fusion_helps_summary.csv`,
`materials/ate1_when_fusion_helps.png`,
`materials/ate1_coefficient_recovery.png`.

### 6.2 ATE-2: covariate shift and the density ratio

Classification: the synthetic panel checks the drift and balance identities of
Theorem 1, Definition 12, and Proposition 6 under their applicable assumptions.
Gaussian mean-shift overlap, estimated-ratio coverage, and clipping behavior are
diagnostics rather than bounded-ratio theorem checks.  NSW/PSID is a real-data
stress test.

Synthetic panel.  In the Gaussian mean-shift regime the two laws have analytic
common support, but the density ratio is unbounded.  Let $\delta$ be the mean
shift and $\Sigma$ the common covariance.  Prespecify mild and strong levels by
their squared Mahalanobis distances
$d^2=\delta^\top\Sigma^{-1}\delta$, equivalently by the population oracle-weight
ESS fraction $\exp(-d^2)=\{\mathbb E_O(r_0^2)\}^{-1}$.  Use the existing shift
direction $v$ but rescale it as
$\delta=v\sqrt{d^2/(v^\top\Sigma^{-1}v)}$.  Fix mild at
$d^2=\log 2$ (ESS fraction $0.50$) and strong at $d^2=\log 5$ (ESS fraction
$0.20$) before outcomes are simulated; do not tune them using realized extrema.
Design:
3 SCM families $\times$ $n_R\in\{100,400\}$ $\times$ 2 shift levels, so 12 cells,
with $N_O=5000$ throughout.

Five ratios per replication: oracle $r_0$; unweighted $r\equiv1$; classifier;
calibrated balancing on the raw covariates; calibrated balancing on the raw
covariates plus $\widehat g$.

Metrics: bias of the joint estimator against the drift
$\omega\,\mathbb E_O[\{\widehat r(X)-r_0(X)\}\widehat g(X)]$ computed on $10^5$
draws; RMSE; coverage of the feasible interval of Theorem 5; the balance
residual; the effective sample size and the normalization error of each ratio.
Also report the analytic and empirical $\mathbb E_O(r_0^2)$, ESS fraction,
log-ratio distribution, and fixed ratio quantiles.  Sample minima or maxima are
descriptive only and may not support a global-bound claim.  The unclipped
oracle $r_0$ is primary.  Any clipped ratio is a separately labeled candidate
whose estimand drift is reported; it cannot replace the oracle silently.

NSW and CPS panel.  Per replication, split the 185 NSW treated units at random
into two halves.  The OBS sample is one half together with the comparison
group; the RCT is the other half together with the 260 NSW controls, subsampled
to $n_R\in\{100,200,\text{all}\}$.  The RCT design propensity is the treated
share of that RCT sample.  Design: 3 RCT sizes $\times$ 2 comparison groups
(CPS with 15,992 rows, PSID with 2,490 rows), so 6 cells.

The reference value is the full-experiment difference in means of `re78`,
namely $6349.1-4554.8=1794.3$.  It carries its own sampling error, so report
point estimates, interval lengths, the drift diagnostic
$\omega\{\mathbb P_O(\widehat r\widehat g)-\mathbb P_R\widehat g\}$, and the
effective sample size of each ratio, rather than an RMSE against a known truth.
Ratios 3 to 5 are available here; the oracle ratio is not.

Outputs: `materials/ate2_shift_replications.csv`,
`materials/ate2_shift_summary.csv`, `materials/ate2_shift.png`,
`materials/ate2_nsw_replications.csv`, `materials/ate2_nsw_summary.csv`,
`materials/ate2_nsw.png`.

### 6.3 CATE-1: honest selection, reporting, and the oracle calculus

Classification: the prespecified bounded DGP below is the Theorem 6 check.
Theorem 7 and Corollary 7.1 are checked only with an unregularized spline sieve
satisfying $p<n_R^{\mathrm{tune}}$ and condition number at most $10^8$.
Fixed-ridge neural results
are empirical extensions.  Theorem 7 rows use the unclipped linear-sieve output
and are distinct from the clipped SCM-B Theorem 6 rows.

Data: the family axis is SCM-B, Gaussian SCM 2, Gaussian SCM 3, and STAR
real-$X$, all with a common covariate marginal.  This replaces SCM 1 only in
CATE-1 and leaves the four-family, 24-cell count unchanged.  SCM-B is the
prespecified bounded theorem-check DGP:

$$X_j,U,\epsilon\stackrel{\mathrm{iid}}{\sim}\operatorname{Unif}[-1,1],\quad
e_R(X)=0.5+0.2X_1,$$

$$\tau(X)=0.15+0.10X_1-0.05X_2^2,\qquad
Y(a)=0.15\sin(\pi X_1)+0.05X_2X_3+0.10U+0.05\epsilon+a\tau(X).$$

The OBS treatment log-odds add the designated confounding multiplier times
$U$.  The strong-heterogeneity variant replaces $\tau$ by $2\tau$; still
$|Y|<1$.  Clip nuisance outcome predictions to $[-1,1]$, their contrast to
$[-2,2]$, and CATE candidates to $[-2,2]$.  Therefore the deterministic
Theorem 6 constants are the conservative bounds $B_0=9$ and $\bar r_0=1$,
specified before simulation.
Gaussian and STAR-real-$X$ cells keep their existing synthetic $\tau$ and have
radius status `N/A`.

Reference configuration: $n_R=200$, $N_O=5000$, moderate confounding, base
heterogeneity.  Three axes, 6 configurations per family, 24 cells.

For CATE-1, each displayed $n_R$ and $N_O$ is the size of each independently
generated nuisance, tuning, selection, and reporting role.  The $10^5$ truth
sample is a fifth independent RNG namespace and is never used for selection.

| axis | configurations | new cells per family |
| --- | --- | ---: |
| RCT size | $n_R\in\{100,200,400\}$ | 3 |
| Prediction quality | confounding in $\{0,2\}$ | 2 |
| Heterogeneity | strong variant | 1 |

Learners: DRF and RF.  Sieves: splines with $K=3$, splines with $K=5$, and the
frozen neural basis.  Candidates and comparators: joint by Algorithm 2 over
$\mathcal G$; RCT-only $(0,0)$ for splines and
$(0,0,\widehat\rho_R)$ for neural; $\lambda$-only; $\omega$-only; the grid
oracle chosen with the true $\tau$; the Corollary 7.1 oracle calculus and its
feasible plug-in for DRF; experimental grounding; domain-indicator pooling.
Every neural restriction and oracle ranges over its allowed $\rho$ values as
well as the allowed $\lambda,\omega$ values.

Metrics: CATE risk on the $10^5$ test sample; regret against
$2\varepsilon_{\mathrm{select}}(0.05)$ only for the bounded theorem-check DGP;
the leading-risk check at the four grid points only on Theorem-7-eligible rows;
selected coefficients against the oracle calculus.

Outputs: `materials/cate1_sieve_validation_replications.csv`,
`materials/cate1_sieve_validation_summary.csv`,
`materials/cate1_sieve_validation.png`, `materials/cate1_regret.png`,
`materials/cate1_coefficients.png`.

### 6.4 CATE-2: STAR real-data source-law stress test

Purpose and classification: this is a real-data stress test and an exact
source-marginal check.  It asks whether the pipeline behaves coherently when
treatment, outcomes, covariates, and a learned representation come from STAR.
It is not a theorem check.  Historical within-school assignment rates describe
the STAR randomization context only; they are not the propensity used by the
CATE-2 estimators, and this construction does not by itself establish causal
identification.

Before any replication, construct and freeze one empirical source-law family
for each outcome and $\alpha$ level.  Treatment is the small class; the control
pools regular and regular-with-aide classes.  The exact learner covariate tuple is

$$X=(\text{gender},\text{ethnicity3},\text{free lunch},\text{school type},
\text{school ID},\text{all corresponding missingness indicators}).$$

Mode-impute a categorical value before forming its tuple and retain its
missingness indicator.  No learner covariate may be omitted from or added to
this tuple without a protocol revision.  Let $H_0$ be the rows with treatment
and both outcomes observed.  On $H_0$, compute each exact cell's treated share
$\widetilde e_x$.  Retain entire cells for which
$\widetilde e_x\in[0.15,0.85]$, call the resulting frozen cohort $H$, and let
$P_R$ be uniform on $H$.  Then the exact-cell empirical propensity is

$$e_x=P_R(A=1\mid X=x)=
\frac{\#\{i\in H:X_i=x,A_i=1\}}{\#\{i\in H:X_i=x\}}.$$

Thus $e_x=\widetilde e_x\in[0.15,0.85]$.  Freeze the row IDs, preprocessing,
exact-cell probability table, and exact-cell propensity table.  Standardize
kindergarten mathematics (primary) and reading (secondary) using constants
computed once on $H$.  Every CATE-2 AIPW score and RF learner uses
$e(X)=e_x$; the within-school rate is retained only as descriptive STAR
context.  Record the eligible-row count and the number and size distribution
of exact $X$ cells before G0; these are frozen design facts, not replication
outcomes.

Using fixed folds and a design-only RNG namespace, cross-fit a pooled outcome
regression on the eligible cohort.  Compute out-of-fold residuals, their ranks
within treatment arm, and positive source weights

$$w_i(\alpha)=\begin{cases}
1+\alpha(2R_i-1),&A_i=1,\\
1-\alpha(2R_i-1),&A_i=0,
\end{cases}\qquad \alpha\in\{0.4,0.8\}.$$

Freeze these outcome-specific ranks and weights before replication.  A
whole-cohort in-sample residual is prohibited.  The ranks and weights are
design metadata and are never supplied as learner features.

Generate an OBS row from $P_O$ in three steps: draw an exact tuple
$X=x$ from $P_R^X$; draw $A=a$ from the frozen empirical $P_R(A=a\mid X=x)$;
then, among eligible rows with that same exact $(X=x,A=a)$, draw a row with
probability proportional to its frozen $w_i(\alpha)$.  Consequently
$P_O^X=P_R^X$ exactly at the source-law level and $r_0(x)=1$ by construction,
while arm-specific outcome ranks induce confounding.  Realized independent
samples need not have identical empirical $X$ frequencies.

Within each replication, draw every role independently with replacement from
its source law using disjoint RNG namespaces.  Draw RCT and OBS nuisance
samples, RCT and OBS tuning samples, RCT and OBS selection samples, and RCT and
OBS reporting samples separately; never conditionally redraw a failed fit.
Use $n_R\in\{200,400\}$ for each RCT nuisance and tuning role, $N_O=5000$ for
each OBS nuisance and tuning role, and 1000 observations per source for each of
selection and reporting.  The neural basis and its standardization use only
the nuisance roles.  The candidate is fit on tuning, chosen on selection, then
frozen and assessed on reporting.

Design: 2 RCT sizes $\times$ 2 confounding levels $\times$ 2 outcomes, so 8
cells.  Learners and comparators as in CATE-1, with the neural basis trained on
the OBS nuisance part.

Metrics: after freezing eligibility, recompute the primary descriptive ATE
benchmark as the treated-minus-control mathematics contrast in that changed
eligible cohort.  The previously reported 8.3-point full-cohort contrast is
context only.  For the CATE, the reporting-score difference is, for splines,
$\widehat C^{\mathrm{report}}(\text{selected})-
\widehat C^{\mathrm{report}}(0,0)$ and, for neural candidates,
$\widehat C^{\mathrm{report}}(\widehat s)-
\widehat C^{\mathrm{report}}(0,0,\widehat\rho_R)$.  Negative means improvement
over the basis-matched RCT-only candidate under the frozen empirical source
law.  Do not call the selected
selection-score minimum unbiased performance.  Diagnostics include the naive
OBS contrast, the source-law $X$ probability-vector difference (which must be
zero to machine precision), and realized-sample $X$ discrepancies (descriptive
sampling noise only).

Outputs: `materials/cate2_star_real_replications.csv`,
`materials/cate2_star_real_summary.csv`, `materials/cate2_star_real.png`.

## 7. Data

### 7.1 SCM families and STAR real-$X$

Implemented in the two existing scripts.  Four covariates from a one-factor
Gaussian model, an unmeasured confounder $U$ entering the OBS assignment and
the outcome, a known RCT propensity clipped to $[0.2,0.8]$, nonlinear outcome
bases, and a mean-shift regime with a closed-form ratio.  The STAR real-$X$
family uses five STAR covariates with synthetic $U$, outcomes, and $\tau$.

### 7.2 STAR kindergarten cohort, real treatment and outcomes

Source: `https://vincentarelbundock.github.io/Rdatasets/csv/AER/STAR.csv`,
sha256 `0e8b179ea3d883730b25008ca1293c4c8f886aa8d5d299c03ab6ff05d0123ae6`,
11,598 rows.  Counts verified on 2026-09-09 and 2026-09-10:

- 5,786 students have a kindergarten class type and both kindergarten scores.
- With the pooled control, 1,738 are in small classes and 4,048 are not.
- The 79 schools have within-school small-class shares from 0.16 to 0.43.  This
  describes the historical randomization context; CATE-2 instead constructs
  the exact-cell $e_x$ of Section 6.4 and retains only cells in
  $[0.15,0.85]$.
- The sensitivity subset that uses regular classes only as the control has
  3,743 students; historical within-school shares there reach 1.00.  The former
  school-level trimming diagnostic kept 78 schools and 3,730 students, but it
  does not define the v3 CATE-2 propensity or eligible cohort.
- The naive small-versus-pooled-control difference in kindergarten mathematics
  is 8.3 points, which is a descriptive benchmark for the ATE part of CATE-2,
  not an identified causal truth.  It is context only after the eligible cohort
  changes; the primary descriptive contrast is recomputed on that cohort.
- Missing values in the covariates number at most nine per column; impute the
  mode and add a missingness indicator.

The earlier 23-stratum summary (median 120 students, residual stratum 202) is a
descriptive coarse-stratification diagnostic only.  CATE-2 does not use it to
claim exact balance; Section 6.4 defines exact equality at the frozen
source-law level using the complete learner-$X$ tuple, including school ID and
missingness indicators.

### 7.3 NSW experimental sample and the CPS and PSID comparison samples

Files from Dehejia's distribution, fetched and hashed on 2026-09-09:

| file | URL | sha256 | rows |
| --- | --- | --- | --- |
| `nsw_dw.dta` | `https://users.nber.org/~rdehejia/data/nsw_dw.dta` | `d1bd2680a1c6f799f1c6d2455bf29633fdf19be01cb19490621c20a560b4e072` | 445, of which 185 treated |
| `cps_controls.dta` | `https://users.nber.org/~rdehejia/data/cps_controls.dta` | `80f0123eaf723870bd060ba9f8569617a7cb519a462dc20badbfb3443cd8374b` | 15,992 |
| `psid_controls.dta` | `https://users.nber.org/~rdehejia/data/psid_controls.dta` | `7beebae8928035d6abf662b994ce77a4fbea7ae6575037ff94e92fe099b43e14` | 2,490 |

Columns: `treat, age, education, black, hispanic, married, nodegree, re74,
re75, re78`.  Covariates for all models: the eight pre-treatment columns plus
the indicators `re74 == 0` and `re75 == 0`.  Load with `pandas.read_stata`,
verify the sha256 and the row count before use, and do not commit the files;
download them at run time as the STAR file is downloaded.

## 8. Compute budget and schedule

Cell and replication counts under the conditional full designs of Section 6:

| study | cells | replications | dominant cost per replication |
| --- | ---: | ---: | --- |
| ATE-1 | 44 | 4,400 | two ridge fits, ten estimators, adaptive oracle draws starting at $2\times10^4$ |
| ATE-2 synthetic | 12 | 1,200 | five ratios, each a logistic fit and a calibration solve |
| ATE-2 NSW and PSID | 6 | 600 | three ratios on small data |
| CATE-1 | 24 | 2,400 | one neural training, six sieve fits over 25 or 75 candidates |
| CATE-2 | 8 | 800 | one neural training, six sieve fits |
| performance total | 94 | 9,400 | |
| four ATE-1 coverage extensions | 4 sentinel cells | 1,600 additional | 400 beyond screening per cell |
| unique full-plan total | 94 | 11,000 | performance plus coverage extension |

Thus the accounting is exactly $9{,}400$ screening/performance replications plus
$1{,}600$ additional coverage replications, or $11{,}000$ unique full-plan
replications.  The Section 5 gate has a maximum of 600 replications and is not
added when its runs are reusable prefixes; any explicitly approved non-reusable
stress corner is counted separately.

The version-2 estimate was about 8 single-core hours for 9,400 replications.
Budget 9--10 single-core hours for 11,000 plus adaptive oracle draws, subject to
measurement at G1; do not promise an eight-core time before parallel scaling is
measured.  Neural training is expected to dominate CATE at about 5 seconds per
replication.  Report measured role-level and replication-level costs in every
results memo.

Required order: Section 5 G0, G1, and any decision-directed G2; record the gate
decision; only `GO_FULL_TWO_CHANNEL` authorizes the full sequence ATE-1, ATE-2
synthetic, CATE-1 with splines, CATE-1 with the neural basis, ATE-2 with NSW,
then CATE-2.  `REFRAME` and `REDESIGN` require a revised protocol; `INCONCLUSIVE`
permits only incremental extensions of inconclusive gate cells that remain
below their caps.

## 9. QA, deliverables, and how the results enter the paper

### 9.1 Acceptance tests

Automate these checks in the scripts.  G0 must pass the applicable engineering
checks before G1; every applicable gate check and a recorded decision must pass
before a full-run command is enabled.

1. Determinism: rerunning a cell with the same seed reproduces the replication
   CSV byte for byte.
2. Role independence: nuisance, tuning, selection, and reporting row draws are
   disjoint where sampling without replacement is used; for CATE-2's
   with-replacement design, the random streams are independent and have
   distinct recorded RNG namespaces.  Selection and reporting samples are
   never the same draw.
3. RCT-only agreement: on shared cells the RCT-only estimator agrees with the
   existing benchmark up to the difference in the split scheme.
4. Exact variance: the Monte Carlo variance of the fixed-coefficient estimator
   is within three Monte Carlo standard errors of the Theorem 1 formula in
   every prespecified theorem-check ATE-1 cell.  The code also tests the
   algebraic variance identity on fixed arrays, including covariance terms.
5. Oracle ordering: in every exact-ratio ATE theorem-check cell for which the
   oracle is defined, the oracle-coefficient joint estimator has variance no
   larger than the RCT-only estimator, which is Corollary 2.1.
6. Balance: the calibrated ratio satisfies the balance equations to $10^{-6}$
   wherever the problem is feasible.
7. Sieve fit: the normal-equation solution matches direct minimization of the
   loss to $10^{-8}$ on one cell per learner and sieve.
8. Score identity: on synthetic cells with $r=r_0$, for candidates fixed before
   reporting, the mean reporting-score difference is within three Monte Carlo
   standard errors of the corresponding risk difference.  Unit tests also
   verify the score identity exactly on fixed arrays.
9. Selection/reporting honesty: changing reporting outcomes cannot change the
   selected candidate; changing selection outcomes while holding reporting
   fixed can change the candidate but not the stored per-candidate reporting
   scores.  Neural joint selection enumerates all 75 triples, and its RCT-only
   comparator freezes $\widehat\rho_R$ selected from exactly the three
   $(0,0,\rho)$ candidates.
10. CATE-2 source law: the frozen probability vectors satisfy
    $P_O^X=P_R^X$ and $r_0=1$ to machine precision over the complete exact-$X$
    support, and

    $$\max_x\left|e(x)-\sum_{a\in\{0,1\}}aP_R(A=a\mid X=x)\right|\leq10^{-12}.$$

    This is not replaced by equality of realized sample shares.  The source
    probabilities, exact-cell propensities, and role dependencies are hard G0
    assertions; the direction of induced bias is not.
11. Prespecified bounds: no code path computes a Theorem 6 radius from a sample
    maximum.  Gaussian and STAR result rows carry radius status `N/A`; Gaussian
    tables and prose cannot claim a finite global ratio bound.  Bounded-DGP rows
    include `bound_vacuity_ratio`; passing a vacuous bound is not recorded as a
    substantive performance success.
12. Theorem-7 scope: every Theorem 7 or Corollary 7.1 table row has spline
    basis, $\rho=0$, $p<n_R^{\mathrm{tune}}$, and condition number at most
    $10^8$ (with the realized value logged).
    It also uses the unclipped linear-sieve output.  Neural and clipped SCM-B
    rows are labeled separately and excluded from those tables.
13. Comparator freeze: shrinkage and adaptive weights and decisions use tuning
    data only, include $\widehat C_{RO}$, and are unchanged by perturbing
    ATE evaluation outcomes.  Only $\widehat\theta_R^{\mathrm{eval}}$ and
    $\widehat\theta_O^{\mathrm{eval}}$ change under that perturbation.  Naive
    pooling has nonzero direct OBS pseudo-outcome weight in its final average.
14. Failure accounting: every scheduled seed produces either one result row or
    one explicit failure row.  No replication is dropped, replaced, or redrawn;
    adaptive oracle sample-size increases retain the same replication ID.
15. Gate lock: the full-run entry point fails closed unless a gate memo records
    one valid decision token and only `GO_FULL_TWO_CHANNEL` unlocks the v3
    94-cell plan.

For CATE-2 additionally verify that frozen source weights are derived from OOF
residuals, that no learner receives ranks or source weights, and that the naive
OBS contrast moves in the intended direction.  Failure of the last directional
diagnostic is reported; it is not repaired by resampling.

### 9.2 Deliverables

- Scripts under `code/`, each with a one-line description, a
  `code/requirements.txt`, and the acceptance tests runnable with one command.
- A gate results memo recording G0/G1/G2 counts, prespecified tolerances, failed
  assertions, measured runtime, and exactly one Section 5 decision token.
- A replication CSV, a summary CSV, and a PNG per study under `materials/`.
- A results memo per study under `materials/` in Markdown: the design, what was
  run, the figure, the table, the measured compute cost, and three to five
  sentences of findings tied to the labels and theorem numbers of Section 1.3.
  Mark every result as theorem check, diagnostic, empirical extension, or
  real-data stress test, and retain the mechanism-matched-heuristic label for
  relevant comparators.  Report
  negative and mixed findings plainly; the manuscript already says that
  estimated coefficients carry no universal no-harm guarantee, so a mixed
  result is a finding rather than a failure.

### 9.3 What goes into Section 5 of the paper

Deciding this now avoids redrawing figures later.  The body of the paper has
room for two figures and two tables.

| item | content | source |
| --- | --- | --- |
| Figure 1 | RMSE ratio to RCT-only against $n_R$, two panels for the correct and the misspecified OBS outcome model, one line per estimator | ATE-1 |
| Figure 2 | CATE risk against $n_R$ for DRF and RF, splines and the neural basis | CATE-1 |
| Table 1 | Interval coverage and mean length for the ATE, and bounded-DGP regret against the radius for the CATE | ATE-1 and CATE-1 |
| Table 2 | Real-data results: the NSW panel and the STAR real-outcome panel | ATE-2 and CATE-2 |

Everything else goes to the appendix: the ratio comparison of ATE-2, the
coefficient-recovery scatter, the exact-variance check, the leading-risk check,
and the honest versus cross-fitted panel.  Produce every figure at a width of
one column and a height of about 55 millimetres, with font sizes readable at
that size, so that no figure needs to be redrawn for the submission.

### 9.4 Questions to raise before deviating

The split sizes, the rectangle $[0,1]^2$, the thresholds, the candidate grid,
the network architecture, the trimming rule, the STAR control definition, the
exact learner-$X$ tuple, source-law construction, gate criteria, and seed/RNG
mapping are fixed by this document.  Raise any change with the author first;
changing one invalidates reuse of earlier gate replications and may change what
a theorem check means.

## 10. Revision history

### 10.1 Version 1 to Version 2

1. ATE-1 was a full factorial with 288 cells and an oracle step of $10^5$ draws
   per replication.  It is now three axes around a reference configuration with
   44 cells and $2\times10^4$ oracle draws.  CATE-1 was 72 cells and is now 24.
   The total workload falls from about 38,000 replications to 9,400.
2. CATE-2 selected the OBS sample on the outcome, which shifts the covariate
   distribution and breaks the $r_0\equiv1$ claim that its CATE metric rests
   on.  The construction now selects on the within-arm residual rank and then
   rebalances the covariate strata exactly, so $r_0\equiv1$ holds by
   construction.  The role split of the constructed OBS sample is now
   specified.
3. The baselines were named but not specified.  Sections 4.1 to 4.6 now give
   each one as a formula.
4. A propensity trimming rule is stated in Section 3.1, and Section 7.2 reports
   the STAR propensity range that makes the RF denominator safe.
5. Section 8 sizes the compute, and Section 9.3 fixes which figures and tables
   go into the body of the paper.

One bibliographic record is still pending, predictions as surrogates, because
Semantic Scholar rate-limited that query on both 2026-09-09 and 2026-09-10.
Prognostic covariate adjustment was verified on 2026-09-10 and is now in the
table of Section 4.7.  Neither item is a baseline to implement; both belong to
related work, so no experiment is blocked.  Retry before the related work
section is written.

### 10.2 Version 2 to Version 3

1. Section 5 adds the maximum-600-replication Agile gate and makes the 94-cell
   full plan conditional on a recorded scientific/engineering decision.
2. Results are classified as theorem checks, diagnostics, empirical extensions,
   or real-data stress tests.  Theorem 6's Hoeffding radius is now restricted to
   a prespecified bounded DGP; Gaussian and STAR radii are `N/A`.
3. CATE selection and reporting use independent samples.  Neural
   standardization is nuisance-only, and fixed-ridge neural results are excluded
   from Theorem 7 and Corollary 7.1 checks.
4. CATE-2 now freezes an OOF-residual-weighted empirical source law before
   replication and enforces $P_O^X=P_R^X$ over the complete learner covariate
   tuple.  Its interpretation is limited to a real-data stress test and
   source-marginal check.
5. Gaussian shift levels use prespecified Mahalanobis distance and ESS rather
   than realized ratio extrema.  Oracle moment simulation now reports
   moment-specific MCSE and adapts from $2\times10^4$ up to $10^5$ draws.
6. Shrinkage and adaptive comparators include their shared-RCT covariance and
   freeze all decisions on tuning data.  Naive pooling now has a literal OBS
   pseudo-outcome contribution, and literature-inspired comparators retain the
   mechanism-matched-heuristic label.
7. The 9,400 performance replications are screening runs.  Four ATE-1 coverage
   sentinels add 1,600 replications, yielding 11,000 unique full-plan
   replications, with expanded honesty, identity, scope, failure-accounting,
   and gate-lock assertions in Section 9.1.
8. Final v3 review made the CATE-2 propensity exact-cell-specific, matched the
   neural RCT-only comparator on selected ridge, made the gate terminal states
   exhaustive, added neural and CATE-2 G0 smoke paths, replaced the coverage
   tolerance by three MCSEs, and added `bound_vacuity_ratio`.
