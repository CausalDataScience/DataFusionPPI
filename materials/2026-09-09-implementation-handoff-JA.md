# Implementation Handoff for JA: Four Experiments for DataFusionPPI

Prepared 2026-09-09 for the person implementing the experiments.  This
document is self-contained: it restates the parts of the manuscript that the
code must implement, describes the existing code and where it departs from
the manuscript, specifies the data, the estimators, the four experiments, the
outputs, and the acceptance tests.  The companion protocol
`materials/2026-09-09-aistats-experiment-protocol.md` is the shorter
statement of the same plan; where the two differ, this document wins.

Reading order: Section 1 (what the paper claims), Section 2 (repository and
conventions), Section 3 (building blocks with formulas), Section 4 (the four
experiments), Section 5 (data), Section 6 (baselines), Section 7 (QA and
deliverables).  The manuscript is `manuscript/main.pdf`; equation and theorem
numbers below refer to the build of 2026-09-09 (mirror commit `d8bfb51` plus
the removal of Section 4.5).

## 1. What the paper claims and what the experiments must show

### 1.1 Setting

Two data sources.  A randomized controlled trial (RCT, source $R$) with
covariates $X$, a binary treatment $A$ assigned with a known propensity
$e(X)=P_R(A=1\mid X)$ bounded away from 0 and 1, and an outcome $Y$.  An
observational study (OBS, source $O$) with the same variables but with
treatment assignment that may depend on unmeasured confounders.  The target is
the RCT population: the conditional average treatment effect
$\tau(x)=\mathbb E_R\{Y(1)-Y(0)\mid X=x\}$ and its mean $\theta=\mathbb E_R\{\tau(X)\}$.
Each source is split into three roles that must be independent (manuscript
Assumption 1): a nuisance sample $D_S^{\mathrm{nuis}}$ for fitting outcome
models and the density ratio, a tuning sample $D_S^{\mathrm{tune}}$ for
choosing coefficients and fitting CATE candidates, and an evaluation sample
$D_S^{\mathrm{eval}}$ for the final estimate.  Sizes are $n_S^{\mathrm{nuis}}$,
$n_S^{\mathrm{tune}}$, $n_S^{\mathrm{eval}}$; the manuscript writes
$n=n_R^{\mathrm{eval}}$ and $N=n_O^{\mathrm{eval}}$, and
$\mathbb P_{R,n}$, $\mathbb P_{O,N}$ for averages over the evaluation samples,
$\mathbb P_S^{\mathrm{tune}}$ and $\mathbb P_S^{\mathrm{nuis}}$ for averages over the
other roles.

Covariate shift.  Either the two covariate laws coincide (Assumption 3,
density ratio $r_0\equiv1$) or they differ with $r_0=dP_R^X/dP_O^X$
(Assumption 4).  A candidate ratio $r$ transports OBS covariate averages to
the RCT population; $r=r_0$ is exact.

### 1.2 The two channels

Fitted outcome regressions $\widehat\mu_R=(\widehat\mu_{R,0},\widehat\mu_{R,1})$
(from the RCT nuisance sample) and $\widehat\mu_O=(\widehat\mu_{O,0},\widehat\mu_{O,1})$
(from the OBS nuisance sample).  The OBS treatment-effect prediction is
$\widehat g(x)=\widehat\mu_{O,1}(x)-\widehat\mu_{O,0}(x)$.  For a pair
$q=(q_0,q_1)$ the AIPW score is

$$\varphi(V;q)=q_1(X)-q_0(X)+\frac{A}{e(X)}\{Y-q_1(X)\}-\frac{1-A}{1-e(X)}\{Y-q_0(X)\},$$

which satisfies $\mathbb E_R\{\varphi(V;q)\mid X\}=\tau(X)$ for every fixed $q$
(Equation 6).  The fused regression is
$\widehat\mu_\lambda=(1-\lambda)\widehat\mu_R+\lambda\widehat\mu_O$, the fused
pseudo-outcome is $Z_\lambda=\varphi(V;\widehat\mu_\lambda)$, and
$Z_\lambda=Z_0+\lambda\Delta$ with $\Delta=Z_1-Z_0$.

- Channel $\lambda$ changes the regression inside the pseudo-outcome.
- Channel $\omega$ adds the control variate
  $\omega\{\mathbb P_{O,N}(r\widehat g)-\mathbb P_{R,n}\widehat g\}$, whose mean is
  zero when $r=r_0$.

### 1.3 Results the experiments check

ATE (Section 3).

- Theorem 1: for fixed $(\lambda,\omega)$ and $r=r_0$,
  $\mathbb E\{\widehat\theta_{r_0}(\lambda,\omega)\}=\theta$; for a fixed
  candidate $\widetilde r$ the mean drift is
  $\omega\,\mathbb E_O[\{\widetilde r(X)-r_0(X)\}\widehat g(X)]$; and

  $$\operatorname{Var}\{\widehat\theta_{r_0}(\lambda,\omega)\}=\frac1n\operatorname{Var}_R\{Z_\lambda-\omega\widehat g(X)\}+\frac{\omega^2}{N}\operatorname{Var}_O\{r_0(X)\widehat g(X)\}.$$

- Definition 5 and Theorem 2: with $A=\operatorname{Var}_R(\Delta)$,
  $B=\operatorname{Var}_R(\widehat g)+\frac nN\operatorname{Var}_O(r_0\widehat g)$,
  $C=\operatorname{Cov}_R(Z_0,\Delta)$, $D=\operatorname{Cov}_R(Z_0,\widehat g)$,
  the variance is $\frac1n\{\operatorname{Var}_R(Z_0)+2C\lambda-2D\omega+A\lambda^2+B\omega^2\}$
  and the unconstrained oracle is $\lambda^\star=-C/A$, $\omega^\star=D/B$.
- Theorem 3: $n\operatorname{Var}$ equals $\sigma_0^2$ for RCT-only,
  $\sigma_0^2-D^2/B$ for $\omega$-only at $\omega^\star$, and
  $\sigma_0^2-D^2/B-C^2/A$ for the joint oracle, where $\sigma_0^2=\operatorname{Var}_R(Z_0)$.
- Algorithm 1 and Theorem 4: the plug-in rule of Section 3.2.1 with the
  evaluation-sample variance estimator $\widehat V$ gives a valid normal
  approximation when the three roles are independent.
- Theorem 5 and Corollary 5.1: the same with an estimated ratio, with an
  explicit drift.

CATE (Section 4).

- Definition 15: the DRF and RF losses (Section 3.4 below).
- Algorithm 2, Definition 18, Theorem 6: fit candidates on the tuning samples
  for every grid pair $(\lambda,\omega)\in\mathcal G$, compute the validation
  score on the evaluation samples, select the minimizer.  The selected
  candidate's CATE risk is within $2\varepsilon_{\mathrm{eval}}(\alpha)$ of the
  best candidate on the grid plus twice the largest transport error of $r$.
- Proposition 10, Theorem 7, Corollary 7.1: for a linear sieve, closed-form
  candidates; the risk decomposes as approximation plus a variance term
  $\mathcal J_j(\lambda,\omega)$ that is a separated quadratic in
  $(\lambda,\omega)$ for DRF with the same shape as the ATE variance.

Every experiment reports the RCT-only estimator as the reference and all
gains as ratios to it.

## 2. Repository, environment, and conventions

Canonical project root: `research/papers/DataFusionPPI/` in the PAIOS
workspace; the publication mirror is the private GitHub repository
`CausalDataScience/DataFusionPPI` (checked out at
`/Users/yonghanjung/Dropbox/Personal/Research/Code/DataFusionPPI`).  Work in
the canonical root; the author syncs the mirror.

Layout: `code/` (scripts), `materials/` (protocols, CSV replications,
summaries, PNG figures), `manuscript/` (LaTeX), `memo/`.

Environment used for the existing artifacts: Python 3 with numpy 2.4,
pandas 3.0, scikit-learn 1.9, matplotlib; torch 2.13 is available for the
neural basis.  Pin the versions you use in a `code/requirements.txt`.

Existing scripts.

- `code/ssem_ate_pilot.py`: the SCM families (`SCMParameters`,
  `generate_scm_parameters(scm_id)` for `scm_id` in {1,2,3}), the true
  $\tau(x)$ (`tau_function`), the true $\theta$ (`true_theta`), the RCT and
  OBS propensities, the sampling (`sample_dataset`), the oracle density ratio
  for the shifted regime (`oracle_ratio`; the RCT covariate law is the OBS law
  with a mean shift, so $r_0$ is a Gaussian likelihood ratio), and smoke
  tests.
- `code/total_budget_nested_benchmark.py`: the current ATE benchmark.  It
  imports the pilot module as `legacy`, adds the STAR real-$X$ family
  (`load_star_support`, `star_sample`), fits ridge outcome models on a fixed
  feature map (`fit_outcome`), fits a logistic domain classifier
  (`fit_ratio`, clipped to $[0.05,20]$), computes pseudo-outcomes (`pseudo`),
  tuning moments (`tuning_moments`), coefficients (`tune_all`), and writes
  `materials/<stem>_replications.csv`, `materials/<stem>_summary.csv`, and
  `materials/<stem>.png` (`write_outputs`, `summarize`, `make_figure`).

Seeds: base seed 190602; per-cell seeds from `stable_seed(...)` and
`rng_for(...)`.  Keep this scheme so that every replication is reproducible
from the design cell and the replication index.

Output conventions: one replication-level CSV (one row per design cell,
replication, and estimator), one summary CSV (mean and Monte Carlo standard
error over replications), and one composite PNG per study, all under
`materials/` with a descriptive stem.  Every row carries the design cell
fields, the estimator name, the metric name, and the seed.

### 2.1 Where the existing script departs from the manuscript

These differences matter because the manuscript's guarantees hold for the
honest three-way split and for Algorithm 1 exactly as stated.

1. Sample reuse.  The script runs five outer folds; within each fold it
   tunes coefficients by a two-fold inner cross-fit on the training part and
   evaluates on the held-out fold, then averages the five fold estimates.
   The manuscript's Theorem 4 assumes independent nuisance, tuning, and
   evaluation samples and states that cross-fitted reuse is not covered.  The
   new experiments use the honest three-way split as the primary scheme and
   keep the cross-fitted scheme only as a comparison panel in ATE-1.
2. Coefficient rule.  The script minimizes the full empirical quadratic
   (including the cross term `h = cov(delta, gr)`) over the box $[-2,2]^2$
   (`exact_box_minimizer`).  Algorithm 1 uses the separated plug-in rule
   with thresholds and projection onto a prespecified rectangle; the
   experiments implement Algorithm 1 exactly (Section 3.2 below) with the
   rectangle $[0,1]^2$.
3. Variance.  The script records a `score_variance_proxy_nonexact`.  The
   experiments compute $\widehat V$ of Equation (ate-evaluation-variance) on
   the evaluation samples and report interval coverage.
4. Ratio clipping.  Clipping to $[0.05,20]$ changes the ratio target
   (manuscript remark after Definition 11).  Report the clip fraction and run
   the unclipped ratio where it is finite; use clipping only as a labeled
   variant.
5. Tuning-size factor.  The script uses the training-size ratio in
   $\widehat B$.  The manuscript's $\widehat B$ uses $n/N$ with the evaluation
   sizes, because the quantity being minimized is the evaluation-sample
   variance.

## 3. Building blocks

All quantities below are computed conditional on the fitted objects, which
is exactly how the theorems are stated.

### 3.1 Roles, nuisances, and pseudo-outcomes

1. Split each source into nuisance, tuning, and evaluation roles with the
   fractions $0.4/0.3/0.3$ (RCT) and $0.6/0.2/0.2$ (OBS) unless a study
   says otherwise.  Draw the split with the cell seed.
2. Fit $\widehat\mu_R$ on $D_R^{\mathrm{nuis}}$ and $\widehat\mu_O$ on
   $D_O^{\mathrm{nuis}}$ with the existing `fit_outcome` (ridge on the fixed
   feature map); in ATE-1 the misspecification knob replaces the feature map
   by the raw covariates.  Fit the density ratio on the two nuisance samples.
3. On any sample with known $e$, compute $Z_0=\varphi(V;\widehat\mu_R)$,
   $Z_1=\varphi(V;\widehat\mu_O)$, $\Delta=Z_1-Z_0$, and $\widehat g$.

### 3.2 ATE estimator, Algorithm 1, and inference

Estimator for fixed $(\lambda,\omega)$ and ratio $r$ on the evaluation samples:

$$\widehat\theta_r(\lambda,\omega)=\mathbb P_{R,n}Z_\lambda+\omega\{\mathbb P_{O,N}(r\widehat g)-\mathbb P_{R,n}\widehat g\}.$$

Tuning moments on the tuning samples (Equation ate-empirical-moments):
$\widehat A=\widehat{\operatorname{Var}}_R(\Delta)$,
$\widehat B=\widehat{\operatorname{Var}}_R(\widehat g)+\frac nN\widehat{\operatorname{Var}}_O(\widehat r\widehat g)$
with $n=n_R^{\mathrm{eval}}$, $N=n_O^{\mathrm{eval}}$,
$\widehat C=\widehat{\operatorname{Cov}}_R(Z_0,\Delta)$,
$\widehat D=\widehat{\operatorname{Cov}}_R(Z_0,\widehat g)$.

Plug-in rule (Equation ate-separated-plugin) with rectangle
$[L_\lambda,U_\lambda]\times[L_\omega,U_\omega]=[0,1]^2$ and thresholds
$t_A=10^{-3}\widehat{\operatorname{Var}}_R(Z_0)$, $t_B=10^{-3}\widehat{\operatorname{Var}}_R(Z_0)$:

$$\widehat\lambda=\begin{cases}\operatorname{Proj}_{[0,1]}(-\widehat C/\widehat A)&\widehat A>t_A\\0&\text{otherwise}\end{cases},\qquad
\widehat\omega=\begin{cases}\operatorname{Proj}_{[0,1]}(\widehat D/\widehat B)&\widehat B>t_B\\0&\text{otherwise}\end{cases}.$$

Restricted estimators: RCT-only fixes $(0,0)$; $\lambda$-only fixes
$\omega=0$ and uses $\widehat\lambda$; $\omega$-only fixes $\lambda=0$ and
uses $\widehat\omega$; joint uses both.

Variance estimator and interval on the evaluation samples:
$\psi_{R,i}=Z_{\widehat\lambda}(V_{R,i})-\widehat\omega\widehat g(X_{R,i})$,
$\psi_{O,j}=\widehat\omega\,\widehat r(X_{O,j})\widehat g(X_{O,j})$,
$\widehat V=\widehat{\operatorname{Var}}_{R,\mathrm{eval}}(\psi_R)/n+\widehat{\operatorname{Var}}_{O,\mathrm{eval}}(\psi_O)/N$
(unbiased sample variances), interval
$\widehat\theta\pm z_{0.975}\sqrt{\widehat V}$.  Coverage is the fraction of
replications whose interval contains $\theta$.

Oracle coefficients for the ceiling estimator: given the fitted nuisances of
a replication, draw $10^5$ fresh RCT and OBS observations from the design,
compute $(A,B,C,D)$ from Definition 5 with the same $n/N$, and use
$(\lambda^\star,\omega^\star)=(-C/A,D/B)$ projected onto $[0,1]^2$.

Exact-variance check (ATE-1 metric 2): with the oracle coefficients held
fixed across replications of the same cell, compare the Monte Carlo variance
of $\widehat\theta_{r_0}(\lambda^\star,\omega^\star)$ with the right-hand side
of the Theorem 1 variance formula evaluated on the $10^5$ draws.

### 3.3 Density-ratio estimators

Classifier (Definition 11): pool the nuisance covariates, label RCT as 1,
fit a logistic regression on the standardized covariates, and set
$\widehat r_{\mathrm{CLS}}(x)=\frac{1-\widehat\pi}{\widehat\pi}\frac{\widehat s(x)}{1-\widehat s(x)}$
with $\widehat\pi=n_R^{\mathrm{nuis}}/(n_R^{\mathrm{nuis}}+n_O^{\mathrm{nuis}})$.
Report the normalization error $\mathbb P_O^{\mathrm{nuis}}(\widehat r)-1$, the
effective sample size $(\sum\widehat r)^2/\sum\widehat r^2$ on the OBS
evaluation sample, and the clip fraction if clipping is used.

Calibrated balancing ratio (Definition 12): with base ratio
$\widehat r_{\mathrm{base}}=\widehat r_{\mathrm{CLS}}$ and a feature map
$\phi$ without intercept, solve the convex problem

$$\widehat\xi\in\arg\min_{\xi}\left[\log\mathbb P_O^{\mathrm{nuis}}\{\widehat r_{\mathrm{base}}(X)e^{\xi^\top\phi(X)}\}-\xi^\top\mathbb P_R^{\mathrm{nuis}}\{\phi(X)\}\right]$$

by L-BFGS with an analytic gradient (the gradient is the balance residual
$\mathbb P_O^{\mathrm{nuis}}(\widehat r_{\mathrm{BAL}}\phi)-\mathbb P_R^{\mathrm{nuis}}\phi$),
and set $\widehat r_{\mathrm{BAL}}=\widehat r_{\mathrm{base}}e^{\widehat\xi^\top\phi}/\mathbb P_O^{\mathrm{nuis}}\{\widehat r_{\mathrm{base}}e^{\widehat\xi^\top\phi}\}$.
Verify Proposition 6 numerically: the balance residual is below $10^{-6}$ in
every component.  Two feature maps are used: $\phi=$ standardized raw
covariates, and $\phi=$ raw covariates plus $\widehat g$.  If the RCT feature
mean lies outside the convex hull of the OBS feature vectors the problem has
no interior minimizer; detect this by a residual that does not vanish and
record the cell as infeasible.

Oracle ratio: `oracle_ratio` for the SCM families and the STAR support
probabilities for the STAR family.

### 3.4 CATE losses, sieve fit, validation, and metrics

Losses (Definition 15) on the tuning samples, with $G_\zeta=(\widehat g-\zeta)^2$
and $\chi(A,X)=\{A-e(X)\}^2/[e(X)\{1-e(X)\}]$:

$$\widehat L_{\mathrm{DRF}}(\zeta;\lambda,\omega,r)=\mathbb P_R^{\mathrm{tune}}\{Z_\lambda-\zeta(X)\}^2+\omega\left[\mathbb P_O^{\mathrm{tune}}\{r(X)G_\zeta(X)\}-\mathbb P_R^{\mathrm{tune}}\{G_\zeta(X)\}\right],$$

$$\widehat L_{\mathrm{RF}}(\zeta;\lambda,\omega,r)=\mathbb P_R^{\mathrm{tune}}\left[\frac{[Y-m_\lambda(X)-\{A-e(X)\}\zeta(X)]^2}{e(X)\{1-e(X)\}}\right]+\omega\left[\mathbb P_O^{\mathrm{tune}}\{r(X)G_\zeta(X)\}-\mathbb P_R^{\mathrm{tune}}\{\chi(A,X)G_\zeta(X)\}\right],$$

where $m_S(x)=e(x)\widehat\mu_{S,1}(x)+\{1-e(x)\}\widehat\mu_{S,0}(x)$ and
$m_\lambda=(1-\lambda)m_R+\lambda m_O$.

Sieve fit (Proposition 10).  For a feature map $b(x)\in\mathbb R^p$ and ridge
$\rho\geq0$, with $w_{\mathrm{DRF}}=1$, $\widetilde Z_{\mathrm{DRF},\lambda}=Z_\lambda$,
$w_{\mathrm{RF}}=\chi(A,X)$, $\widetilde Z_{\mathrm{RF},\lambda}=\{Y-m_\lambda(X)\}/\{A-e(X)\}$,
solve

$$\{\widehat H_j(\omega)+\rho I\}\beta=\mathbb P_R^{\mathrm{tune}}\left[w_jb\{\widetilde Z_{j,\lambda}-\omega\widehat g\}\right]+\omega\,\mathbb P_O^{\mathrm{tune}}(r\,b\,\widehat g),\qquad \widehat H_j(\omega)=(1-\omega)\mathbb P_R^{\mathrm{tune}}(w_jbb^\top)+\omega\,\mathbb P_O^{\mathrm{tune}}(r\,bb^\top),$$

and set $\widehat\zeta_j(x;\lambda,\omega)=b(x)^\top\widehat\beta_j$.  For fixed
$\omega$ solve once and reuse the factorization for all $\lambda$ (the
right-hand side is affine in $\lambda$).  Check the fit against a direct
minimization of the loss on one cell.

Feature maps.

- Splines: additive cubic B-splines with $K\in\{3,5\}$ interior knots per
  continuous covariate at empirical quantiles of the RCT tuning sample,
  dummies for categorical covariates, and an intercept.
- Frozen neural basis: on the OBS nuisance sample train a network with a
  shared trunk $h(x)\in\mathbb R^{64}$ (two hidden layers of 64 units, ReLU)
  and two linear heads $\widehat\mu_{O,a}(x)=\theta_a^\top h(x)+c_a$, mean
  squared error, Adam, early stopping on a 10% holdout of the nuisance
  sample; freeze the trunk and set $b(x)=(1,h(x))$; take $\widehat g$ from
  the heads.  Standardize $h$ on the RCT tuning sample before the sieve fit.

Candidate grid: $\mathcal G=\{0,0.25,0.5,0.75,1\}^2$; for the neural basis
add the ridge axis $\rho\in\{10^{-2},10^{-1},1\}$ so that the candidate set
has $M=75$ elements.  Theorem 6 applies to any finite candidate set.

Validation score (Definition 18) on the evaluation samples, with
$G=\{\widehat g(X)-\widehat\zeta(X)\}^2$:

$$\widehat C(\lambda,\omega)=\mathbb P_{R,n}\left[\{Z_0-\widehat\zeta(X;\lambda,\omega)\}^2-\omega G\right]+\mathbb P_{O,N}\{\omega\,r(X)G\}.$$

Select the minimizer over the candidate set with the lexicographically
smallest index on ties.

Radius (Theorem 6): $\varepsilon_{\mathrm{eval}}(\alpha)=2B^2\sqrt{2\log(4M/\alpha)}\{2/\sqrt n+\bar r/\sqrt N\}$
with $\alpha=0.05$, $B=\max$ of $|Z_0|$, $|\widehat g|$, and $|\widehat\zeta|$
over the evaluation samples and candidates, $\bar r=\max r$ on the OBS
evaluation sample.

Metrics.

- CATE risk $\mathcal R(\zeta)=\mathbb E_R\{\zeta(X)-\tau(X)\}^2$ estimated on an
  independent test sample of $10^5$ RCT-population covariates with the true
  $\tau$ (synthetic families) or replaced by the validation-score difference
  $\widehat C(\lambda,\omega)-\widehat C(0,0)$ on the evaluation sample when
  $\tau$ is unknown (CATE-2).
- Regret $\mathcal R(\text{selected})-\min_{\mathcal G}\mathcal R$, its ratio to
  $2\varepsilon_{\mathrm{eval}}(0.05)$, and the fraction of replications with
  regret at most $2\varepsilon_{\mathrm{eval}}(0.05)$.
- Leading risk: $\Gamma_R=\mathbb E_R(bb^\top)$, $\beta_p=\Gamma_R^{-1}\mathbb E_R(b\tau)$,
  $a_p^2=\mathcal R(b^\top\beta_p)$ on the test sample;
  $\psi_{R,j}(\lambda,\omega)=w_jb\{\widetilde Z_{j,\lambda}-\omega\widehat g-(1-\omega)\zeta_p\}$
  and $\psi_O(\omega)=\omega r_0 b(\widehat g-\zeta_p)$ on fresh draws;
  $\mathcal J_j=\operatorname{tr}[\Gamma_R^{-1}\operatorname{Var}_R(\psi_{R,j})]/n_R^{\mathrm{tune}}+\operatorname{tr}[\Gamma_R^{-1}\operatorname{Var}_O(\psi_O)]/n_O^{\mathrm{tune}}$;
  compare the Monte Carlo mean of $\mathcal R(\widehat\zeta_j)-a_p^2$ with
  $\mathcal J_j$ at $(\lambda,\omega)\in\{(0,0),(0.5,0),(0,0.5),(0.5,0.5)\}$.
- Oracle sieve coefficients (Corollary 7.1, DRF): $A_p,B_p,C_p,D_p$ from the
  same fresh draws with $\zeta_p$ known ("oracle calculus"), giving
  $(\lambda_p^\star,\omega_p^\star)=(-C_p/A_p,D_p/B_p)$; a feasible version
  replaces $\zeta_p$ by the RCT-only fit $\widehat\zeta_j(\cdot;0,0)$ and the
  expectations by tuning-sample averages.  Report both against the selected
  coefficients.

## 4. The four experiments

Every study: 100 replications per design cell, honest three-way split,
RCT-only estimator as the reference, mean and Monte Carlo standard error for
every metric.

### 4.1 ATE-1: when fusion helps, and whether Algorithm 1 finds it

Data families: SCM 1, 2, 3 (`generate_scm_parameters`) and STAR real-$X$;
common covariate marginal only (`regime="shared"`, $r_0\equiv1$, no ratio
estimation).

Factors.

- $n_R\in\{50,100,200,400\}$ total RCT size, $N_O\in\{1000,5000,20000\}$.
- Confounding strength: scale `obs_latent_coefficient` by
  $\{0,1,2\}$ (0 removes the unmeasured confounder from the OBS assignment).
- OBS outcome-model specification: the existing feature map (`fit_outcome`)
  versus a linear model on the raw covariates.

Estimators: RCT-only, $\lambda$-only, $\omega$-only, joint (Algorithm 1),
oracle-coefficient joint, and the baselines of Section 6 (semi-supervised
ATE, shrinkage, adaptive combination, naive pooling).  Add one panel that
runs the joint estimator with the existing cross-fitted scheme.

Metrics: RMSE ratio to RCT-only, bias, variance ratio; the exact-variance
check; coefficient recovery (scatter of $(\widehat\lambda,\widehat\omega)$
against $(\lambda^\star,\omega^\star)$, boundary-hit rate, threshold-fallback
rate); coverage and mean length of the 95% interval for joint and RCT-only;
the honest-versus-cross-fitted panel.

Outputs: `materials/ate1_when_fusion_helps_{replications,summary}.csv`,
`materials/ate1_when_fusion_helps.png` (RMSE ratio against $n_R$, panels by
confounding strength and family), a coverage table, a coefficient scatter.

### 4.2 ATE-2: covariate shift and the density ratio

Synthetic panel.  `regime="shifted"` with the existing mean shift as the
mild level and twice the shift vector as the strong level (verify overlap:
$\min r_0\geq0.05$ and $\max r_0\leq20$ on a $10^5$ draw; reduce the strong
level until this holds).  $n_R\in\{100,400\}$, $N_O=5000$, SCM families 1 to 3.

Ratios: oracle $r_0$; unweighted $r\equiv1$; classifier; calibrated
balancing with raw covariates; calibrated balancing with raw covariates plus
$\widehat g$.

Metrics: bias of the joint estimator against the drift
$\omega\,\mathbb E_O[\{\widehat r(X)-r_0(X)\}\widehat g(X)]$ computed on $10^5$
draws; RMSE; coverage of the feasible interval (Theorem 5); balance
residual; effective sample size and normalization error of each ratio.

NSW and CPS panel (Section 5.3 for the data).  Per replication: split the
185 NSW treated units at random into two halves; the OBS sample is one half
together with the 15,992 CPS controls; the RCT is the other half together
with the 260 NSW controls, subsampled to $n_R\in\{100,200,\text{all}\}$.  The
design propensity of the RCT is the treated share of the RCT sample.
Reference value: the full-experiment difference in means of `re78`,
$6349.1-4554.8=1794.3$.  Ratios 3 to 5 above.  Report point estimates,
interval lengths, the drift diagnostic
$\omega\{\mathbb P_O(\widehat r\widehat g)-\mathbb P_R\widehat g\}$, and the
effective sample size of each ratio; do not report an RMSE against a known
truth.  Repeat with the 2,490 PSID controls as a second OBS variant.

Outputs: `materials/ate2_shift_{replications,summary}.csv`,
`materials/ate2_shift.png`, `materials/ate2_nsw_cps_{replications,summary}.csv`,
`materials/ate2_nsw_cps.png`.

### 4.3 CATE-1: honest validation on the sieve and the oracle calculus

Data: SCM families 1 to 3 with the existing $\tau(x)$ and a strong variant
(treatment coefficients multiplied by 2), and STAR real-$X$ with the
existing synthetic $\tau$.  Common covariate marginal.

Factors: $n_R\in\{100,200,400\}$, $N_O=5000$, confounding strength
$\{0,1,2\}$, heterogeneity $\{$existing, strong$\}$.

Learners: DRF and RF; sieves: splines ($K=3$ and $K=5$) and the frozen neural
basis.  Candidates and comparators: joint (Algorithm 2 over $\mathcal G$),
RCT-only $(0,0)$, $\lambda$-only, $\omega$-only, grid oracle (true $\tau$),
Corollary 7.1 oracle calculus and feasible plug-in (DRF), experimental
grounding, and domain-indicator pooling (Section 6).

Metrics: CATE risk on the $10^5$ test sample; regret against
$2\varepsilon_{\mathrm{eval}}(0.05)$; the leading-risk check at four grid
points; selected coefficients against the oracle calculus.

Outputs: `materials/cate1_sieve_validation_{replications,summary}.csv`,
`materials/cate1_sieve_validation.png` (risk against $n_R$ by learner and
sieve), `materials/cate1_regret.png` (histogram with the radius),
`materials/cate1_coefficients.png`.

### 4.4 CATE-2: STAR with real treatment, real outcomes, and a constructed confounded OBS sample

Data: the STAR kindergarten cohort (Section 5.2).  Treatment is the small
class; the primary control pools regular and regular-with-aide classes
(5,786 students with both kindergarten scores), and the sensitivity analysis
uses regular only (3,743 students).  Outcome: `mathk` (primary) and `readk`
(secondary), standardized by the cohort mean and standard deviation.
Covariates: gender, ethnicity, birth quarter, free-lunch status, school
type.  Design propensity: the within-school share of small classes among the
analysis subset, by `schoolidk`; drop schools whose share is 0 or 1.

Construction per replication: from the cohort draw an RCT sample of
$n_R\in\{200,400\}$ and an evaluation pool of 2,000; the remaining students
form the OBS pool.  Include an OBS-pool student in the OBS sample with
probability $\sigma(\alpha_0+\alpha_1\tilde y)$ if treated and
$\sigma(\alpha_0-\alpha_1\tilde y)$ if control, where $\tilde y$ is the
standardized outcome, $\alpha_0=0$, and $\alpha_1\in\{0.5,1.0\}$ sets the
confounding strength; record the naive OBS difference in means as a
diagnostic of the induced bias.  The evaluation pool plays the role of
$D_R^{\mathrm{eval}}$ with the design propensity; the RCT sample supplies
$D_R^{\mathrm{nuis}}$ and $D_R^{\mathrm{tune}}$.

Learners and comparators as in CATE-1; the neural basis is trained on the
OBS nuisance part.  Metrics: ATE against the full-cohort experimental
estimate; CATE by the validation-score difference
$\widehat C(\text{selected})-\widehat C(0,0)$ on the evaluation pool, negative
values meaning improvement over RCT-only (Theorem 6 item 1).

Outputs: `materials/cate2_star_real_{replications,summary}.csv`,
`materials/cate2_star_real.png`.

## 5. Data

### 5.1 SCM families and STAR real-$X$

Implemented in `code/ssem_ate_pilot.py` and `code/total_budget_nested_benchmark.py`.
Four covariates from a one-factor Gaussian model, an unmeasured confounder
$U$ entering the OBS assignment (`obs_latent_coefficient`) and the outcome
(`outcome_latent_coefficient`), a known RCT propensity clipped to
$[0.2,0.8]$, nonlinear outcome bases, and a mean-shift regime with a
closed-form ratio.  The STAR real-$X$ family uses the five STAR covariates
with synthetic $U$, outcomes, and $\tau$ (`load_star_support`).

### 5.2 STAR kindergarten cohort (real treatment and outcomes)

Source: `https://vincentarelbundock.github.io/Rdatasets/csv/AER/STAR.csv`,
sha256 `0e8b179ea3d883730b25008ca1293c4c8f886aa8d5d299c03ab6ff05d0123ae6`,
11,598 rows.  Kindergarten fields: `stark` (small, regular, regular+aide),
`readk`, `mathk`, `schoolidk` (79 schools), `schoolk` (rural, suburban,
urban, inner-city), `gender`, `ethnicity`, `birth`, `lunchk`.  Counts
checked on 2026-09-09: 6,325 students with `stark`; 5,786 with both
kindergarten scores; among them 1,738 small and 2,005 regular; within-school
small-class shares range from 0.23 to 1.00 in the small-versus-regular
subset, which is why schools with share 0 or 1 are dropped.  Missing values
in the covariates are at most nine per column; impute the median or the
mode and add a missingness indicator for `birth`.

### 5.3 NSW experimental sample and CPS and PSID comparison samples

Source files (Dehejia's distribution), fetched and hashed on 2026-09-09:

| file | URL | sha256 | rows |
| --- | --- | --- | --- |
| `nsw_dw.dta` | `https://users.nber.org/~rdehejia/data/nsw_dw.dta` | `d1bd2680a1c6f799f1c6d2455bf29633fdf19be01cb19490621c20a560b4e072` | 445 (185 treated, 260 control) |
| `cps_controls.dta` | `https://users.nber.org/~rdehejia/data/cps_controls.dta` | `80f0123eaf723870bd060ba9f8569617a7cb519a462dc20badbfb3443cd8374b` | 15,992 |
| `psid_controls.dta` | `https://users.nber.org/~rdehejia/data/psid_controls.dta` | `7beebae8928035d6abf662b994ce77a4fbea7ae6575037ff94e92fe099b43e14` | 2,490 |

Columns: `treat, age, education, black, hispanic, married, nodegree, re74,
re75, re78`.  Covariates for all models: the eight pre-treatment columns
plus the indicators `re74 == 0` and `re75 == 0`.  Load with
`pandas.read_stata`, verify the sha256 and the row counts before use, and do
not commit the data files; download at run time like the STAR file.

## 6. Baselines and their bibliographic status

Zotero (the author's local library, 590 items) was checked on 2026-09-09.
Present: Hatt et al. 2022 (representation learning with RCT and OBS),
Athey, Chetty, and Imbens 2020 (long-term outcomes), Oberst et al. 2022
(bias-robust integration).  Absent: semi-supervised ATE estimation,
prognostic covariate adjustment, prediction-powered generalization, and
predictions-as-surrogates.  Semantic Scholar was queried for the absent
items; verified records are listed with a DOI, the rest are marked pending.

| baseline | role in the experiments | status |
| --- | --- | --- |
| Semi-supervised ATE (Cheng, Ananthakrishnan, Cai; Biometrics; DOI 10.1111/biom.13298) | ATE-1: RCT outcome model imputed on the OBS covariates, the closest relative of the $\omega$ channel | verified |
| Shrinkage combination (`rosenman2020shrinkage`) | ATE-1: shrink the OBS estimate toward the RCT estimate | in reference.bib |
| Adaptive combination (`cheng2021adaptive`) | ATE-1: weight by an estimated OBS bias | in reference.bib |
| Naive pooling | ATE-1 and CATE-1: fit on RCT and OBS together with a source indicator | no citation needed |
| Experimental grounding (`kallus2018removing`) | CATE-1: OBS CATE plus an RCT-fitted linear correction | in reference.bib |
| Cross-prediction-powered inference (Zrnic and Candès; PNAS; DOI 10.1073/pnas.2322083121) | discussion only | verified |
| Prediction-powered generalization of causal inferences (Demirel, Alaa, Philippakis, and others; ICML 2024; arXiv 2406.02873) | related work: PPI for transporting trial effects | verified |
| Prediction-powered causal inferences (Cadei, Demirel, Bartolomeis, and others; NeurIPS 2025; arXiv 2502.06343) | related work: closest prior use of PPI for causal targets; read before writing Section 1 | verified |
| Prognostic covariate adjustment, predictions-as-surrogates | related work | pending (Semantic Scholar rate-limited on 2026-09-09; retry) |

Implement a baseline only after its record is verified; the author confirms
the final list.

## 7. QA, acceptance, and deliverables

Acceptance tests, all automated in the scripts.

1. Determinism: rerunning a cell with the same seed reproduces the
   replication CSV byte for byte.
2. RCT-only agreement: on the shared cells the RCT-only estimator agrees with
   the existing benchmark up to the split scheme.
3. Exact variance: the Monte Carlo variance of the fixed-coefficient
   estimator is within three Monte Carlo standard errors of the Theorem 1
   formula in every ATE-1 cell.
4. Oracle ordering: the oracle-coefficient joint estimator has variance no
   larger than the RCT-only estimator in every cell (Corollary 2.1).
5. Balance: the calibrated ratio satisfies the balance equations to
   $10^{-6}$ where the problem is feasible.
6. Sieve fit: the normal-equation solution matches direct minimization of
   the loss to $10^{-8}$ on one cell per learner and sieve.
7. Score identity: on synthetic cells with $r=r_0$, the mean of
   $\widehat C(\lambda,\omega)-\widehat C(0,0)$ over replications is within
   three Monte Carlo standard errors of
   $\mathcal R\{\widehat\zeta(\cdot;\lambda,\omega)\}-\mathcal R\{\widehat\zeta(\cdot;0,0)\}$.

Deliverables.

- Scripts under `code/` with a `README` line per script, a
  `requirements.txt`, and the acceptance tests runnable with one command.
- Replication CSV, summary CSV, and PNG per study under `materials/`.
- A results memo per study under `materials/` (Markdown): design, what was
  run, the figure, the table, and three to five sentences of findings tied
  to the theorem numbers above.  Report negative or mixed findings plainly.

Suggested order: ATE-1 (extends the existing code), ATE-2 synthetic, CATE-1
with splines, CATE-1 with the neural basis, ATE-2 NSW and CPS, CATE-2.

Questions to raise with the author before deviating: any change to the
split fractions, the rectangle, the thresholds, the grid, the network
architecture, or the STAR control definition.
