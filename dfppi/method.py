"""What AIPWF, DRF and RF share: the outcome regressions, the density ratio and the AIPW score.

Data.  A sample is a dict {"x": covariates, "a": treatment in {0, 1}, "y": outcome}.  R is the trial (RCT)
and O the observational source (OBS).  Each source is split into three independent samples (Sec. 2):

    nuis   fits the regressions, the density ratio and the sieve    (line 1 of both algorithms)
    tune   selects (lambda, omega) for the ATE and fits the CATE candidates
    eval   computes the ATE estimate and scores the CATE candidates

The estimators

    aipwf.py   Algorithm 1: the average treatment effect
    drf.py     Algorithm 2 with the doubly robust pseudo-outcome: the CATE
    rf.py      Algorithm 2 with the pseudo-outcome of the R-learner: the CATE

Where each shared object of the paper is computed

    Lemma 1   AIPW score phi(V; q), Eq. (1)                  aipw_score
    Table 1   blended regression mu_lambda, score Z^lambda   blended_regression, fused_score
    Sec. 3.2  classifier ratio r_CLS                         classifier_ratio
    Def. 2    calibrated balancing ratio, Eq. (6)            balancing_ratio
    line 1    regressions, g, components, classifier ratio   fit_nuisances
    Prop. 6   the density ratio of one estimator             density_ratio, called by aipwf.py and drf.py
"""
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from scipy.optimize import minimize, root
from sklearn.linear_model import LogisticRegression


def var(v):
    """Unbiased sample variance."""
    return np.var(v, ddof=1)


def remember(function):
    """The same function, with its values kept per input array, so that every sample is evaluated once."""
    seen = {}

    def remembered(x):
        if id(x) not in seen:
            seen[id(x)] = (x, function(x))          # holding x keeps its id from being reused
        return seen[id(x)][1]
    return remembered


# ---------------------------------------------------------------------------- outcome regressions
def flexible_features(x):
    """x, x^2, sin(x) and the pairwise products of the first six covariates."""
    d = min(x.shape[1], 6)
    pairs = [x[:, i] * x[:, j] for i in range(d) for j in range(i + 1, d)]
    return np.column_stack([x, x ** 2, np.sin(x)] + pairs)


def fit_outcome_regression(sample, ridge=5.0, clip_quantile=None):
    """Returns mu with mu(x) = (mu(x, 0), mu(x, 1)), an estimate of E[Y | X = x, A = a].

    Ridge regression of Y on (f, A, A f) with an unpenalised intercept, f = standardised flexible features.
    With clip_quantile = c, each arm's prediction is clipped to the [c, 1 - c] quantile range of that arm's
    predictions on the fitting sample.
    """
    f = flexible_features(sample["x"])
    mean, sd = f.mean(0), f.std(0)
    sd[sd < 1e-8] = 1.0

    def design(x, a):
        z = (flexible_features(x) - mean) / sd
        return np.column_stack([z, a, a[:, None] * z])

    D, y = design(sample["x"], sample["a"]), sample["y"]
    D_mean, y_mean = D.mean(0), y.mean()
    Dc = D - D_mean
    coef = np.linalg.solve(Dc.T @ Dc + ridge * np.eye(D.shape[1]), Dc.T @ (y - y_mean))

    def mu(x):
        return tuple(y_mean + (design(x, np.full(len(x), a)) - D_mean) @ coef for a in (0.0, 1.0))

    if clip_quantile is None:
        return mu
    bounds = [np.quantile(m, [clip_quantile, 1 - clip_quantile]) for m in mu(sample["x"])]
    return lambda x: tuple(np.clip(m, lo, hi) for m, (lo, hi) in zip(mu(x), bounds))


def fit_propensity(sample, C=1.0):
    """Returns e with e(x) an estimate of P(A = 1 | X = x) in the sample's source: logistic regression on the
    standardised flexible features, clipped to [0.01, 0.99].  Only the comparison learner RF.integrative uses it."""
    f = flexible_features(sample["x"])
    mean, sd = f.mean(0), f.std(0)
    sd[sd < 1e-8] = 1.0
    model = LogisticRegression(C=C, max_iter=5000).fit((f - mean) / sd, sample["a"])
    return lambda x: np.clip(model.predict_proba((flexible_features(x) - mean) / sd)[:, 1], 0.01, 0.99)


# ---------------------------------------------------------------------------- density ratio (Sec. 3.2)
def classifier_ratio(x_R, x_O):
    """Sec. 3.2: r_CLS = {(1 - pi) / pi} s / (1 - s).

    s(x) estimates P(S = 1 | X = x) by a logistic regression on the pooled inputs (S = 1 for the trial),
    and pi is the trial share.  s / (1 - s) is computed as exp(log-odds).
    """
    x = np.vstack([x_R, x_O])
    mean, sd = x.mean(0), x.std(0)
    sd[sd < 1e-8] = 1.0
    source = np.r_[np.ones(len(x_R)), np.zeros(len(x_O))]
    model = LogisticRegression(C=1.0, max_iter=5000).fit((x - mean) / sd, source)
    if model.n_iter_[0] >= model.max_iter:
        raise RuntimeError("the classifier of the density ratio did not converge")
    pi = len(x_R) / len(x)
    return lambda v: np.exp(model.decision_function((v - mean) / sd) + np.log((1 - pi) / pi))


def balancing_ratio(r_base, f, x_R, x_O, within_standard_error):
    """Def. 2: r_BAL = r_base e^{xi'f} / P_O{r_base e^{xi'f}}, with xi chosen so that the balance equations
    P_O(r_BAL f_k) = P_R(f_k), Eq. (6), hold on the nuisance samples x_R, x_O.

    within_standard_error = False   exact balance: xi is the root of Eq. (6).            (ATE, f = g)
    within_standard_error = True    each equation holds up to delta_k, the standard error of its two sides,
                                    delta_k^2 = Var_O(r_base f_k) / N + Var_R(f_k) / n:  (sieve CATE)
                                    xi minimises the convex function
                                        log P_O{r_base e^{xi'f}} - xi'P_R(f) + sum_k delta_k |xi_k|,
                                    whose optimality condition is |P_O(r_BAL f_k) - P_R(f_k)| <= delta_k.
    """
    F_O, F_R, w = f(x_O), f(x_R), r_base(x_O)
    scale = F_O.std(0)                      # dividing a feature by a constant does not change r_BAL
    scale[scale < 1e-8] = 1.0
    F_O, F_R = F_O / scale, F_R / scale
    target, K, N, n = F_R.mean(0), F_O.shape[1], len(F_O), len(F_R)

    def tilt(xi):
        """log P_O{r_base e^{xi'f}}, and r_BAL on the OBS nuisance rows."""
        u = np.log(w) + F_O @ xi
        log_mean = u.max() + np.log(np.mean(np.exp(u - u.max())))
        return log_mean, np.exp(u - log_mean)

    def imbalance(xi):
        """P_O(r_BAL f_k) - P_R(f_k) for every k."""
        return tilt(xi)[1] @ F_O / N - target

    if not within_standard_error:
        delta = np.zeros(K)
        xi = root(imbalance, np.zeros(K), tol=1e-12).x
    else:
        delta = np.sqrt(np.var(w[:, None] * F_O, axis=0, ddof=1) / N + np.var(F_R, axis=0, ddof=1) / n)

        def objective(v):                   # v = (xi+, xi-) >= 0 and xi = xi+ - xi-, so |xi_k| = xi+_k + xi-_k
            xi = v[:K] - v[K:]
            gap = imbalance(xi)
            return tilt(xi)[0] - xi @ target + delta @ (v[:K] + v[K:]), np.r_[gap + delta, delta - gap]

        v = minimize(objective, np.zeros(2 * K), jac=True, method="L-BFGS-B", bounds=[(0, None)] * (2 * K),
                     options={"maxiter": 5000, "maxfun": 50000, "ftol": 1e-15, "gtol": 1e-9}).x
        xi = v[:K] - v[K:]

    if np.max(np.abs(imbalance(xi)) - delta) > 1e-6:
        raise RuntimeError("the calibration of the density ratio did not reach its balance equations")
    log_mean = tilt(xi)[0]
    return lambda x: r_base(x) * np.exp(f(x) / scale @ xi - log_mean)


# ---------------------------------------------------------------------------- line 1 of Algorithms 1 and 2
def principal_components(x_nuis, n_components):
    """x -> (1, first principal components of the standardised covariates, each with unit standard deviation),
    computed from the trial nuisance covariates only."""
    mean, sd = x_nuis.mean(0), x_nuis.std(0)
    sd[sd < 1e-8] = 1.0
    z = (x_nuis - mean) / sd
    directions = np.linalg.svd(z, full_matrices=False)[2][:n_components]
    pc_sd = (z @ directions.T).std(0)
    return lambda x: np.column_stack([np.ones(len(x)), ((x - mean) / sd) @ directions.T / pc_sd])


@dataclass
class Nuisances:
    """What line 1 of Algorithms 1 and 2 fits before any density ratio is calibrated."""
    e: float                     # known trial propensity P_R(A = 1 | X)
    mu_R: Callable               # trial outcome regression, x -> (mu_R(x, 0), mu_R(x, 1))
    mu_O: Callable               # OBS outcome regression
    g: Callable                  # OBS effect prediction g = mu_O(., 1) - mu_O(., 0)
    components: Callable         # x -> (1, principal components of the trial covariates)
    r_base: Optional[Callable]   # classifier ratio (Sec. 3.2); None when r_0 = 1 is known
    x_R: np.ndarray              # covariates of the trial nuisance sample, on which a ratio is calibrated
    x_O: np.ndarray              # covariates of the OBS nuisance sample


def fit_nuisances(R_nuis, O_nuis, e, same_domain, n_components=5, clip_quantile=0.005, mu_O=None):
    """Line 1 of Algorithms 1 and 2, up to the density ratio.  Only the two nuisance samples enter.

    same_domain = True    Assumption 3(i): r_0 = 1 is known, and no ratio is estimated.
    same_domain = False   Assumption 3(ii): r_0 is estimated.  The classifier ratio of Sec. 3.2 is fitted
                          here, and each estimator calibrates it to the functions it needs (density_ratio).

    Both outcome regressions are clipped, so that neither predicts on new rows beyond what it predicted on
    its own sample.  mu_O, if given, is the OBS regression
    already fitted on O_nuis: the cross-fitted replication rotates the trial roles only, so it fits mu_O once.
    """
    mu_R = remember(fit_outcome_regression(R_nuis, clip_quantile=clip_quantile))
    mu_O = mu_O or remember(fit_outcome_regression(O_nuis, clip_quantile=clip_quantile))
    components = remember(principal_components(R_nuis["x"], n_components))
    x_R, x_O = R_nuis["x"], O_nuis["x"]

    def g(x):
        mu0, mu1 = mu_O(x)
        return mu1 - mu0

    if same_domain:
        return Nuisances(e, mu_R, mu_O, g, components, None, x_R, x_O)

    # Sec. 3.2 leaves the classifier free.  Its inputs are the principal components, not the raw covariates:
    # with one coefficient per covariate and a small trial sample, a rare covariate value can give a single
    # row an enormous ratio.
    def scores(x):
        return components(x)[:, 1:]

    r_scores = classifier_ratio(scores(x_R), scores(x_O))

    def r_base(x):
        return r_scores(scores(x))

    return Nuisances(e, mu_R, mu_O, g, components, r_base, x_R, x_O)


def density_ratio(nuis, f, within_standard_error):
    """The density ratio r of an estimator whose target needs the functions f balanced (Prop. 6).

    r_0 = 1 known (Assumption 3(i))       r = 1.
    r_0 estimated (Assumption 3(ii))      the classifier ratio, calibrated on the nuisance samples so that it
                                          balances f (Def. 2).
    """
    if nuis.r_base is None:
        return lambda x: np.ones(len(x))
    return remember(balancing_ratio(nuis.r_base, f, nuis.x_R, nuis.x_O, within_standard_error))


# ---------------------------------------------------------------------------- scores
def aipw_score(sample, mu, e):
    """Lemma 1, Eq. (1): phi(V; q) = q(X, 1) - q(X, 0) + A {Y - q(X, 1)} / e - (1 - A) {Y - q(X, 0)} / (1 - e)."""
    q0, q1 = mu(sample["x"])
    a, y = sample["a"], sample["y"]
    return q1 - q0 + a * (y - q1) / e - (1 - a) * (y - q0) / (1 - e)


def blended_regression(nuis, lam):
    """Table 1: mu_lambda = (1 - lambda) mu_R + lambda mu_O."""
    def mu(x):
        (r0, r1), (o0, o1) = nuis.mu_R(x), nuis.mu_O(x)
        return (1 - lam) * r0 + lam * o0, (1 - lam) * r1 + lam * o1
    return mu


def fused_score(sample, nuis, lam):
    """Table 1: Z^lambda = phi(V; mu_lambda).  Z^0 is the trial score and Z^1 - Z^0 is Delta."""
    return aipw_score(sample, blended_regression(nuis, lam), nuis.e)
