# Random Nonlinear SCM ATE-Fusion Debugging Pilot: PERR Plan

**Date:** 2026-08-22

**Master seed:** 190602

**Status:** manager-approved execution specification

## Necessity Gate

This experiment asks whether fully feasible prediction-powered fusion has similar bias and lower variance than RCT-only AIPW when the data-generating mechanism is not hand-tuned to favor a particular fusion direction. Three nonlinear SCMs are generated once from a prespecified fixed-seed law. Every nuisance function, density ratio, and tuning coefficient is refitted in each outer replication.

This is explicitly a **debugging pilot**. The three SCM draws are illustrative fixed realizations, and each cell has only 20 Monte Carlo replications. Results do not establish a theorem, a random-SCM-family claim, calibrated coverage, or general efficiency.

## Random SCM generator

There are $K=3$ SCM draws. SCM $k\in\{0,1,2\}$ uses

```text
SeedSequence([190602, 7000, k])
```

exactly once. The seed is not searched or changed in response to results. All quantities below are mutually independent unless an equation links them.

Let $U\sim N(0,1)$, $\varepsilon_X\sim N(0,I_4)$, and $\varepsilon_Y\sim N(0,1)$. Draw

$$
b\sim N(0,0.15^2I_4),
\qquad
\ell\sim N(0,0.30^2I_4),
$$

$$
d_j\stackrel{\mathrm{iid}}{\sim}\operatorname{Unif}(0.75,1.05),
\qquad
\delta\sim N(0,0.18^2I_4),
$$

and define $D=\operatorname{diag}(d_1,\ldots,d_4)$. The OBS covariates and shared-regime RCT covariates follow

$$
X=b+\ell U+D\varepsilon_X.
$$

The shifted-regime RCT covariates follow

$$
X=b+\delta+\ell U+D\varepsilon_X.
$$

Thus

$$
X_O\sim N(b,\Sigma),
\qquad
\Sigma=\ell\ell^\top+D^2,
$$

while $X_R\sim N(b,\Sigma)$ in shared and $X_R\sim N(b+\delta,\Sigma)$ in shifted. The shared regime has $P_R^X=P_O^X$. The shifted regime has $P_R^X\ne P_O^X$ almost surely. The common covariance makes the shifted Gaussian log density ratio exactly linear in $x$ and gives common support and finite fixed-order ratio moments.

## Nonlinear treatment and outcome mechanisms

Define

$$
\phi_m(x)=
(x_1,x_2,x_3,x_4,\sin x_1,\cos x_2,x_1x_2,x_3x_4,x_1^2)^\top,
$$

$$
\phi_\tau(x)=
(x_1,x_2,\sin x_3,x_1x_4,x_2^2)^\top,
$$

and

$$
\psi_A(x)=
(x_1,x_2,\tanh x_3,\sin x_4,x_1x_2)^\top.
$$

Draw

$$
\alpha_m\sim N(0,0.2^2),
\qquad
\beta_m\sim N(0,0.30^2I_9),
$$

$$
\eta_Y=S_YV_Y,
\qquad
P(S_Y=-1)=P(S_Y=1)=\frac12,
\qquad
V_Y\sim\operatorname{Unif}(0.4,0.8),
$$

$$
\alpha_\tau\sim\operatorname{Unif}(0.8,1.2),
\qquad
\beta_\tau\sim N(0,0.22^2I_5),
\qquad
\sigma_Y\sim\operatorname{Unif}(0.8,1.2).
$$

The common outcome mechanism is

$$
Y
=
\alpha_m+\beta_m^\top\phi_m(X)+\eta_YU
+A\tau(X)+\sigma_Y\varepsilon_Y,
$$

$$
\tau(X)=\alpha_\tau+\beta_\tau^\top\phi_\tau(X).
$$

For the RCT, draw

$$
\alpha_R\sim N(0,0.1^2),
\qquad
\gamma_R\sim N(0,0.18^2I_5),
$$

and use the known propensity

$$
e_R(x)
=
\operatorname{clip}
\left[
\operatorname{expit}\{\alpha_R+\gamma_R^\top\psi_A(x)\},
0.20,
0.80
\right].
$$

For OBS, draw

$$
\alpha_O\sim N(-0.2,0.1^2),
\qquad
\gamma_O\sim N(0,0.30^2I_5),
$$

$$
\eta_A=S_AV_A,
\qquad
P(S_A=-1)=P(S_A=1)=\frac12,
\qquad
V_A\sim\operatorname{Unif}(0.7,1.1),
$$

and

$$
P(A=1\mid X,U,S=O)
=
\operatorname{clip}
\left[
\operatorname{expit}
\{\alpha_O+\gamma_O^\top\psi_A(X)+\eta_AU\},
0.05,
0.95
\right].
$$

The OBS treatment is confounded because both treatment and outcome depend on $U$. Random signs prevent fixing the confounding direction in advance.

## Exact RCT target

Let $\mu_R=b$ in shared and $\mu_R=b+\delta$ in shifted. Write $\beta_\tau=(q_1,\ldots,q_5)^\top$. Since $X_R\sim N(\mu_R,\Sigma)$,

$$
\begin{aligned}
\theta_R
={}&\mathbb E_R\{\tau(X)\}\\
={}&\alpha_\tau+q_1\mu_{R,1}+q_2\mu_{R,2}
+q_3e^{-\Sigma_{33}/2}\sin(\mu_{R,3})\\
&+q_4(\Sigma_{14}+\mu_{R,1}\mu_{R,4})
+q_5(\Sigma_{22}+\mu_{R,2}^2).
\end{aligned}
$$

For every SCM and regime, 250,000 independent Gaussian covariates check this analytic target. The smoke gate requires absolute error at most

$$
\max\{0.01,5\operatorname{MCSE}\}.
$$

## Honest learning with small RCT samples

The fitted outcome feature map is generic and deliberately not identical to the SCM bases. It contains all four linear terms, all four squares, all six pairwise interactions, four sine terms, and four cosine terms: 22 raw features. Training-sample means and standard deviations standardize them. A single joint ridge model uses

$$
Y\approx\beta_0^\top\Phi(X)+A\beta_A^\top\Phi(X),
$$

where $\Phi$ adds an intercept. The joint design has $2(22+1)=46$ parameters and fixed ridge penalty $\alpha=5$; only the global baseline intercept is unpenalized. RCT and OBS use the same fitting procedure. This is stable when $n<p$ and does not invert separate treatment-arm matrices.

RCT nuisance training requires at least three observations in each arm. OBS nuisance training requires at least 30. Fixed-seed preflight found no RCT failures and minimum arm count four. If any rerun fails, the program aborts and reports it; it never drops a replication or redraws a seed.

Because the shifted Gaussian log ratio is exactly linear, the estimated ratio uses only the four raw covariates in a standardized linear ridge logistic domain classifier. It has five parameters including its intercept, fixed $C=1$, prior correction,

$$
\widehat r(x)
=
\frac{\widehat P(S=R\mid x)}{1-\widehat P(S=R\mid x)}
\frac{1-\widehat\pi}{\widehat\pi},
$$

and clipping to $[0.05,20]$. Evaluation normalization is forbidden. Shared also estimates $r$ rather than being told that $r=1$, so it measures generic pipeline overhead.

## Honest tuning and primary estimators

Training fits $\widehat\mu_R$, $\widehat\mu_O$, and $\widehat r$. Independent calibration data select

$$
\widehat\lambda
=
\Pi_{[-2,2]}
\left{
-\frac{\widehat C}{\widehat A}
\right},
$$

$$
\widehat\omega
=
\Pi_{[-2,2]}
\left{
\frac{\widehat D}{\widehat B}
\right},
$$

where

$$
\widehat A=\widehat{\operatorname{Var}}_R(\Delta),
\qquad
\widehat C=\widehat{\operatorname{Cov}}_R(Z_0,\Delta),
$$

$$
\widehat D=\widehat{\operatorname{Cov}}_R(Z_0,g),
$$

$$
\widehat B
=
\widehat{\operatorname{Var}}_R(g)
+\frac{n}{5000}\widehat{\operatorname{Var}}_O(\widehat r g).
$$

Denominators at or below $10^{-10}$ fall back to zero. Record coefficient means, standard deviations, clipping counts and rates, nominal unconstrained gains, and gains attained by clipped coefficients.

Only two primary estimators are reported:

$$
\widehat\theta_{\mathrm{AIPW}}=\mathbb P_{R,\mathrm{eval}}Z_0,
$$

and

$$
\widehat\theta_{\mathrm{fusion}}
=
\mathbb P_{R,\mathrm{eval}}Z_{\widehat\lambda}
+\widehat\omega
\left\{
\mathbb P_{O,\mathrm{eval}}(\widehat r g)
-\mathbb P_{R,\mathrm{eval}}g
\right\}.
$$

Oracle-r, ignore-r, and lambda-only estimators are not primary methods.

## Ratio debugging diagnostics

The known Gaussian oracle ratio is used only for diagnostics. Every evaluation pair records

$$
\mathbb P_O(\widehat r)-1,
\qquad
\operatorname{ESS}(\widehat r),
\qquad
\max_j\widehat r_j,
$$

the fraction clipped to either boundary,

$$
\left\{\mathbb P_O(\widehat r-r)^2\right\}^{1/2},
$$

$$
\mathbb P_O\{(\widehat r-r)\widehat g\},
$$

and the mechanically available feasible-minus-oracle ratio-channel difference

$$
\widehat\omega\mathbb P_O\{(\widehat r-r)\widehat g\}.
$$

These diagnostics do not add an oracle estimator to the primary comparison. Low ESS or heavy clipping is reported rather than filtered.

## Repetitions and sample accounting

For every $(\text{SCM},\text{regime},n,\text{outer replication})$, independently generate a training pair, calibration pair, and evaluation pair. Set

$$
n\in\{20,30\},
\qquad
N=5000,
\qquad
B=20.
$$

The values $n$ and $N$ are **per-stage evaluation sizes and per-stage auxiliary sizes**, not total sample sizes.

- For per-stage $n=20$, one outer replication generates 60 RCT and 15,000 OBS observations.
- For per-stage $n=30$, one outer replication generates 90 RCT and 15,000 OBS observations.
- AIPW uses RCT training and evaluation data. Fusion additionally uses the RCT calibration sample and all three OBS samples. This is a pipeline diagnostic, not a resource-matched contest.

## Outputs

The replication CSV has

$$
3\times2\times2\times20\times2=480
$$

rows. The summary CSV, grouped by SCM, regime, per-stage $n$, and estimator, has

$$
3\times2\times2\times2=24
$$

rows. Each group reports bias and MCSE, empirical unconditional pipeline variance and normal-theory MCSE, variance ratio, RMSE, coverage, paired fusion-minus-AIPW difference and MCSE, conditional variance estimates, coefficient instability, channel gains, arm counts, and ratio diagnostics.

Across-SCM statements use medians, ranges, and counts across the 12 SCM-regime-$n$ cells. Three fixed draws are descriptive and do not estimate a random-SCM population.

The PNG has two regime rows and four columns: per-stage $n=20$ bias, $n=20$ variance ratio, $n=30$ bias, and $n=30$ variance ratio. SCM ID is the horizontal axis. The title states per-stage and total-generated sample accounting.

## Realized fixed SCM draws

The fixed generator produced the following coefficients. Values below are rounded to six decimals for display; the seed law above and the script deterministically reproduce the exact floating-point values.

### SCM 0

```text
b       = [ 0.126877,  0.237062, -0.230769, -0.153203]
ell     = [ 0.271265, -0.455727,  0.264416,  0.258733]
d       = [ 0.800852,  0.863037,  0.930385,  0.853494]
delta   = [-0.305772,  0.202848,  0.366962,  0.142662]
alpha_m =  0.132129
beta_m  = [-0.189500, 0.206827, 0.477616, 0.282286, 0.075600,
           0.116060, 0.089714, 0.708373, 0.325870]
eta_Y   = -0.654198; sigma_Y = 0.936249
alpha_tau = 1.159079
beta_tau  = [0.012746, 0.022896, -0.107276, -0.292286, -0.299325]
alpha_R = -0.044606
gamma_R = [0.028347, -0.208019, -0.162193, 0.088927, 0.092225]
alpha_O = -0.196877
gamma_O = [0.054635, 0.374653, 0.249265, -0.005720, 0.173085]
eta_A   = 0.968690
theta_shared = 0.864726; theta_shifted = 0.793643
```

### SCM 1

```text
b       = [ 0.060589, -0.087839, -0.133907, -0.158678]
ell     = [-0.073173,  0.239811,  0.087905,  0.448555]
d       = [ 0.777056,  0.860516,  0.961163,  0.805616]
delta   = [-0.392535,  0.009757, -0.118062, -0.096632]
alpha_m =  0.171123
beta_m  = [0.114991, 0.247140, -0.481117, -0.182443, 0.193259,
          -0.129625, -0.115778, -0.005690, 0.168853]
eta_Y   = 0.650849; sigma_Y = 0.842962
alpha_tau = 0.876012
beta_tau  = [0.392121, -0.429120, 0.381514, -0.205015, -0.119533]
alpha_R = 0.099314
gamma_R = [0.135293, -0.266865, -0.172982, -0.243322, 0.000399]
alpha_O = -0.250966
gamma_O = [0.233774, -0.231026, 0.387778, -0.317368, 0.190309]
eta_A   = 0.833731
theta_shared = 0.817886; theta_shifted = 0.612896
```

### SCM 2

```text
b       = [ 0.066236, -0.200684, -0.032132, -0.109828]
ell     = [ 0.371747, -0.469922, -0.406955,  0.224583]
d       = [ 1.002944,  0.940977,  1.012009,  0.890738]
delta   = [ 0.133660, -0.227001,  0.256139,  0.175928]
alpha_m = -0.054976
beta_m  = [-0.075851, -0.175808, 0.306504, -0.091342, 0.073839,
           -0.113565, -0.317101, -0.140155, -0.030472]
eta_Y   = -0.624802; sigma_Y = 0.999818
alpha_tau = 0.926866
beta_tau  = [-0.353099, -0.221503, -0.206916, 0.107358, -0.385504]
alpha_R = 0.175014
gamma_R = [0.289507, -0.027857, 0.258837, 0.148741, -0.449125]
alpha_O = -0.086524
gamma_O = [0.101010, -0.089341, -0.103314, 0.250734, -0.070630]
eta_A   = 0.912618
theta_shared = 0.517784; theta_shifted = 0.439059
```

## Review gates

The run passes only if:

- all SCM covariances are positive definite and all mean shifts are nonzero;
- shared oracle ratios equal one exactly;
- shifted oracle ratio normalization and weighted first/second moments pass Monte Carlo smoke checks;
- analytic targets pass the 250,000-draw checks;
- the generic outcome map has 22 raw features, joint outcome fit has 46 parameters, and ratio logit has five parameters;
- all 240 RCT nuisance-training cells pass the prespecified arm gate without redraw or drop;
- compact and expanded variance calculations agree;
- repeated miniature and full runs are exactly deterministic;
- outputs contain exactly `rct_aipw` and `estimated_r_datafusion`, with 20 replications per group;
- row counts are 480 and 24, and all values are finite and nonmissing;
- coefficient clipping and ratio instability diagnostics are present and reported prominently;
- the figure passes visual QA;
- scoped `git diff --check` passes and only the five approved paths change.

## Retrospect

Interpret empirical variance as unconditional over the three-sample refitting pipeline but conditional on each illustrative fixed SCM. Interpret row-level variance estimates as conditional on fitted nuisances, ratio, and tuning coefficients. Report Monte Carlo uncertainty, coefficient clipping, ratio $L^2$ error, bias functional, ESS, and feasible-minus-oracle channel discrepancy. Do not convert three SCMs and 20 repetitions into a general validity or efficiency claim.

In the completed debugging run, fusion had lower empirical variance in 10 of 12 SCM-regime-$n$ cells. Its median variance ratio was $0.624$, with range $[0.291,1.148]$. All six shifted cells but only four of six shared cells had ratios below one. The two losses occurred for shared SCM 1, showing that the feasible pipeline does not mechanically improve every fixed SCM.

Observed absolute bias was smaller for fusion in 8 of 12 cells. Every paired fusion-minus-AIPW mean lay within two of its own Monte Carlo standard errors, but only six lay within one. With $B=20$, this does not establish equal bias.

Fusion's $\widehat\omega$ hit a coefficient boundary 45 times among 240 cell replications (18.75%); $\widehat\lambda$ hit a boundary once. By regime and per-stage size, omega clipping counts were 14/60 for shared $n=20$, 7/60 for shared $n=30$, 17/60 for shifted $n=20$, and 7/60 for shifted $n=30$. This is material tuning instability and must accompany any variance result.

Across the 12 cells, mean fitted-ratio ESS ranged from 3120 to 4448 of 5000, mean maximum weight from 3.45 to 10.88, mean clipping fraction from 0 to 0.0022, and mean empirical $L^2(\widehat r-r)$ error from 0.335 to 0.573. The mean ratio bias functional ranged from $-0.0495$ to $0.0230$, while the mean feasible-minus-oracle ratio-channel difference ranged from $-0.0488$ to $0.0150$. Ratio error is therefore non-negligible even though clipping is rare.
