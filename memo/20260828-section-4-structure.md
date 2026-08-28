# Section 4 Structure: Prediction-Powered Fusion Learning for CATE

Date: 2026-08-28

## Purpose

Section 4 should explain, in dependency order, how an RCT-valid CATE loss is
augmented by a large observational covariate sample without treating the OBS
prediction as causal truth. The primary method is the Doubly Robust Fusion
(DRF) Learner based on an AIPW pseudo-outcome. The secondary method is the
normalized R-Fusion (RF) Learner based on a curvature-normalized R-loss.

The reader should encounter each object only when it becomes necessary. In
particular, the candidate CATE function is denoted by $\zeta$, the ideal risk
is $\mathcal R(\zeta)$, and the fitted observational treatment-effect
prediction is $\widehat g_O$. Its loss feature is

$$
\widehat G_\zeta(X)
\triangleq
\{\widehat g_O(X)-\zeta(X)\}^2.
$$

The coefficient $\omega$ is not part of $\widehat G_\zeta$. It appears only
when the transport imbalance is scaled to form the prediction-powered
correction.

## Opening

The opening paragraph should state that the RCT identifies the CATE and that
the OBS supplies a potentially accurate, but not necessarily causal,
treatment-effect prediction. The construction preserves an RCT-valid target
and uses the OBS prediction only through a transported loss correction.

The opening should distinguish three roles in prose:

1. **Target:** $\tau$ and $\mathcal R(\zeta)$.
2. **Nuisances and design objects:** The nuisance functions are the RCT and
   OBS outcome regressions and $\widehat g_O$. The propensity $e_0$ is a known
   RCT design object, not a nuisance. Under covariate shift, the density ratio
   is an additional nuisance.
3. **Tuning:** $\lambda$ controls nuisance fusion and $\omega$ controls the
   prediction-powered correction.

## 4.1 CATE targets and RCT-valid learning losses

### Reader question

How can the unobserved CATE be learned from an RCT, and what are the two
RCT-valid baseline losses?

### Objects and results

Define

$$
\tau(x)
\triangleq
\mathbb E_R\{Y(1)-Y(0)\mid X=x\},
\qquad
\mathcal R(\zeta)
\triangleq
\mathbb E_R\{\zeta(X)-\tau(X)\}^2.
$$

State the honest-fitting convention before using any learned nuisance. Define

$$
q_R\triangleq\widehat\mu_R,
\qquad
q_O\triangleq\widehat\mu_O,
\qquad
q_\lambda\triangleq(1-\lambda)q_R+\lambda q_O,
$$

and

$$
Z_\lambda\triangleq\varphi(V;q_\lambda).
$$

Known RCT propensity and randomization imply

$$
\mathbb E_R(Z_\lambda\mid X)=\tau(X)
$$

for every fixed $\lambda$. Hence

$$
\mathbb E_R\{Z_\lambda-\zeta(X)\}^2
=
Q_{\mathrm{DRF}}(\lambda)+\mathcal R(\zeta),
$$

where

$$
Q_{\mathrm{DRF}}(\lambda)
\triangleq
\mathbb E_R\{Z_\lambda-\tau(X)\}^2.
$$

$Q_{\mathrm{DRF}}(\lambda)$ is the integrated pseudo-outcome noise. It is not
the final risk of the fitted CATE learner.

For the residual route, define

$$
v_0(x)
\triangleq
e_0(x)\{1-e_0(x)\}
$$

and the fused marginal-outcome nuisance $m_\lambda$. Derive

$$
\mathbb E_R[\{A-e_0(X)\}\{Y-m(X)\}\mid X]
=v_0(X)\tau(X)
$$

and

$$
\mathbb E_R[\{A-e_0(X)\}^2\mid X]=v_0(X).
$$

The raw R-loss therefore induces a $v_0(X)$-weighted projection in a
restricted class. Define the normalized loss

$$
\ell_R^{\mathrm{norm}}(\zeta;m)
\triangleq
\frac{[Y-m(X)-\{A-e_0(X)\}\zeta(X)]^2}{v_0(X)}.
$$

Then

$$
\mathbb E_R\{\ell_R^{\mathrm{norm}}(\zeta;m_\lambda)\}
=
Q_{\mathrm{RF}}(\lambda)+\mathcal R(\zeta),
$$

where

$$
Q_{\mathrm{RF}}(\lambda)
\triangleq
\mathbb E_R\{\ell_R^{\mathrm{norm}}(\tau;m_\lambda)\}.
$$

$Q_{\mathrm{DRF}}$ and $Q_{\mathrm{RF}}$ are method-specific noise floors.
They need not have the same minimizer, and neither is generally the final
learner risk.

## 4.2 Prediction-powered common-marginal fusion and curvature matching

### Reader question

How is the large OBS prediction sample added when the RCT and OBS have the
same covariate marginal, and why does RF use $\kappa$?

Define the fitted OBS prediction by reference to Section 2 and define

$$
\widehat G_\zeta(X)
\triangleq
\{\widehat g_O(X)-\zeta(X)\}^2.
$$

Every definition of $\widehat G_\zeta$ should remind the reader that
$\widehat g_O$ is the fitted OBS treatment-effect prediction from Section 2.
It need not equal $\tau$.

Under the common covariate marginal, the OBS and RCT expectations of
$\widehat G_\zeta$ agree. The DRF empirical correction therefore subtracts
the RCT loss feature from its OBS average.

For RF, define

$$
\kappa(A,X)
\triangleq
\frac{\{A-e_0(X)\}^2}{v_0(X)}.
$$

It has two roles. First,

$$
\mathbb E_R\{\kappa(A,X)\widehat G_\zeta(X)\mid X\}
=\widehat G_\zeta(X).
$$

Second, it matches the empirical curvature of the normalized R-loss. For a
linear sieve and $0\leq\omega\leq1$, the RF quadratic coefficient becomes

$$
(1-\omega)\widehat G_{R,\kappa}
+\omega\widehat G_{O,r},
$$

a positive-semidefinite convex blend when the weights are nonnegative. The
factor $\kappa$ is not necessary for population mean cancellation alone. It
is required for the proposed curvature-matched RF construction.

## 4.3 Exact transport imbalance and method-specific empirical corrections

### Reader question

What is transported under covariate shift, and why do DRF and RF have the
same population correction but different empirical subtractions?

Define the unscaled transport imbalance

$$
B(r,\zeta)
\triangleq
\mathbb E_O\{r(X)\widehat G_\zeta(X)\}
-\mathbb E_R\{\widehat G_\zeta(X)\}
$$

and its scaled correction

$$
M_\omega(r,\zeta)
\triangleq
\omega B(r,\zeta).
$$

Under the common covariate marginal, $r_0=1$. Under covariate shift, assume
$P_R^X\ll P_O^X$ and define

$$
r_0(x)\triangleq\frac{dP_R^X}{dP_O^X}(x).
$$

In both regimes, $B(r_0,\zeta)=M_\omega(r_0,\zeta)=0$.

Use $M_\omega$ directly in the two population losses rather than repeatedly
expanding the transport difference. At the empirical level, define
method-specific corrections:

$$
\widehat M_{\mathrm{DRF}}
\triangleq
\omega[\mathbb P_O(r\widehat G_\zeta)
-\mathbb P_R(\widehat G_\zeta)]
$$

and

$$
\widehat M_{\mathrm{RF}}
\triangleq
\omega[\mathbb P_O(r\widehat G_\zeta)
-\mathbb P_R(\kappa\widehat G_\zeta)].
$$

Their population expectations agree because
$\mathbb E_R(\kappa\mid X)=1$, but their empirical curvatures differ.

## 4.4 Estimated-ratio drift and loss-targeted transport balancing

### Reader question

What remains when $r_0$ is replaced by $\widehat r$, and what should the
ratio learner balance?

The exact population drift is

$$
B(\widehat r,\zeta)
=
\mathbb E_O[\{\widehat r(X)-r_0(X)\}\widehat G_\zeta(X)],
$$

and the scaled drift is

$$
M_\omega(\widehat r,\zeta)
=
\omega\mathbb E_O[\{\widehat r-r_0\}\widehat G_\zeta].
$$

For a candidate class $\mathcal Z$, define the uniform loss-targeted
transport imbalance

$$
\Delta_{\mathcal Z}(\widehat r)
\triangleq
\sup_{\zeta\in\mathcal Z}
\left|
\mathbb E_O[\{\widehat r-r_0\}\widehat G_\zeta]
\right|.
$$

Classifier odds, Riesz-representer learning, and positive calibration from
Section 3 provide candidate ratio learners. For CATE, their required balance
must cover the generated loss-feature class
$\{\widehat G_\zeta:\zeta\in\mathcal Z\}$. ATE-specific balance conditions
do not automatically establish this CATE condition.

## 4.5 Orthogonality audit

### Reader question

Which nuisance directions are first-order orthogonal, and which require
direct control?

Define the Gateaux derivative before defining projection scores. The DRF
projection score is orthogonal in the outcome-regression nuisance when the
RCT propensity is known. The normalized RF projection score is orthogonal in
the marginal-outcome nuisance. These are score statements, not claims that
the numerical loss value is orthogonal.

The ratio direction satisfies

$$
D_rM_\omega(r_0,\zeta)[h]
=
\omega\mathbb E_O\{h(X)\widehat G_\zeta(X)\},
$$

which is generally nonzero. Loss-targeted transport balancing controls the
realized ratio-error pairing; it does not create Neyman orthogonality. The
section must also retain the caveat that shared OBS regressions can generate
both a score nuisance and $\widehat g_O$, so joint nuisance orthogonality is
not established.

## 4.6 Noise-floor coefficients and final learner-risk coefficients

### Reader question

Does minimizing $Q_j(\lambda)$ choose the coefficients that minimize the
error of the fitted CATE function?

Distinguish the two noise-floor choices

$$
\lambda_{Q,j}^\star
\in
\operatorname*{arg\,min}_\lambda Q_j(\lambda)
$$

from the algorithm-induced final learner-risk oracle

$$
\mathfrak R_j(\lambda,\omega)
\triangleq
\mathbb E_{\mathrm{fit}}
\mathcal R(\widehat\zeta_{\lambda,\omega}^{,j}).
$$

The final coefficients minimize $\mathfrak R_j$ over a prespecified domain.
They generally have no universal closed form and need not coincide with the
noise-floor choices. Honest validation should compare one finite candidate
set containing the identical RCT-only learner. Candidate inclusion gives an
oracle comparator, not universal finite-sample dominance for the selected
learner.

## 4.7 Linear-sieve finite-sample theory

### Reader question

What can be established when the CATE learner is restricted to a fixed
linear or spline sieve?

Define

$$
\zeta_{\boldsymbol\beta}(x)
\triangleq
b_p(x)^\top\boldsymbol\beta.
$$

For DRF and RF, define the unpenalized empirical Hessian and score before the
penalized normal equation. The finite-sample sequence should contain:

1. exact coefficient and projection-risk identities;
2. deterministic Gram-stability and ridge-bias bounds;
3. source-centered RCT and OBS score concentration;
4. estimated-ratio score and Hessian drifts;
5. a finite-candidate honest-validation inequality;
6. an operational additive-spline corollary with deterministic complexity or
   a finite grid, uniform eigenvalue and sub-Gaussian conditions, explicit
   ridge order, and the validation remainder.

The short rate may omit the ridge term only when the ridge bias is zero or of
the displayed stochastic order. The result should explain when the small RCT
sample remains the bottleneck.

## 4.8 Scope and open theory

The closing subsection should state exactly what is established:

- exact-ratio DRF and normalized RF target identities;
- method-specific score orthogonality;
- explicit ratio nonorthogonality and loss-targeted drift control;
- honest finite-candidate selection;
- fixed-basis linear-sieve and conditional additive-spline results.

It should also state what remains open:

- joint orthogonality with shared OBS nuisances;
- unrestricted end-to-end neural learner theory;
- a universal closed form for final learner-risk-optimal coefficients;
- universal finite-sample no-harm with estimated coefficients or ratios;
- a general claim that OBS data change the RCT-driven convergence exponent.

## Narrative dependency

The final reader path is

$$
\text{RCT-valid target}
\longrightarrow
\text{normalized RF curvature}
\longrightarrow
\text{exact transport correction}
\longrightarrow
\text{estimated-ratio balance}
\longrightarrow
\text{orthogonality audit}
\longrightarrow
\text{final-risk tuning}
\longrightarrow
\text{finite-sample sieve theory}.
$$
