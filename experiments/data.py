"""The simulation designs: the synthetic designs SCM 1 and SCM 3, the benchmarks IHDP and ACIC 2016, and the lung
cancer RCT-EHR pair of the intFRT package.  The paper uses SCM 1 (with its bias sweep and its large-heterogeneity
variant), IHDP and ACIC 2016; SCM 3 and the lung cancer pair are not used in the paper.

In the first three designs a latent U ~ N(0, 1) enters the outcome in both sources and the treatment in the OBS only:

    trial   A ~ Bernoulli(0.35)
    OBS     A ~ Bernoulli(clip(expit(eta(X) + nu(X) U), 0.05, 0.95)),      nu(x) = nu0 {1 + eta_std(x)}

eta_std is eta standardised over the trial covariates: the confounding is stronger for units that are more
likely to be treated, and it changes sign where eta_std < -1.  The OBS contrast is biased by

    delta(x) = gamma {E[U | A = 1, x] - E[U | A = 0, x]},

which varies with x, and nu0 sets its size (see each design).  The OBS covariate law equals the trial law
("same") or is an exponential tilt of it ("weak_shift", "strong_shift") with effective sample fraction
1 / E_O[r_0^2] = 0.70 or 0.40.  A tilt changes the law of X only: (U, A, Y) given X stays the same.

Every design offers   laws (its OBS covariate laws),   draw(rows, rng) -> random numbers,   trial(numbers),
obs(numbers, law),   and the truth   theta = E_R tau_0(X),   truth_x, truth_tau = trial covariates and tau_0 at them.
"""
import hashlib
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import brentq
from scipy.special import expit

E_TRIAL = 0.35                                    # known trial propensity
EFFECTIVE_SAMPLE_FRACTION = {"same": 1.0, "weak_shift": 0.70, "strong_shift": 0.40}
# The misspecified shift: the tilt exp{a q(X)} with q = (h^2 - 1) / sqrt(2), a > 0, of the same
# standardised direction h as the shifts above.  It widens the law of h, so the OBS covers the trial tails, and a
# logistic classifier on the principal components (linear in h) cannot represent it.  Effective sample fraction
# 0.40.  For a standard normal h the tilted h is N(0, s^2) with 1 / E_O[r_0^2] = sqrt(2 s^2 - 1) / s^2.
MISSPEC_FRACTION = 0.40
MISSPEC_SD = float(np.sqrt((1.0 + np.sqrt(1.0 - MISSPEC_FRACTION ** 2)) / MISSPEC_FRACTION ** 2))


def obs_propensity(eta, nu, u):
    """P_O(A = 1 | x, U = u) for propensity index eta and confounding strength nu."""
    return np.clip(expit(eta + nu * u), 0.05, 0.95)


def confounding_bias(eta, nu, u_mean, u_var, gamma):
    """delta = gamma {E[U | A = 1] - E[U | A = 0]} under the OBS law, for units with propensity index eta,
    confounding strength nu and U ~ N(u_mean, u_var).  Gauss-Hermite quadrature over U, one row per unit."""
    nodes, weights = hermgauss(96)
    u = u_mean[:, None] + np.sqrt(2.0 * u_var) * nodes
    p = obs_propensity(eta[:, None], nu[:, None], u)
    treated, control = p * weights, (1.0 - p) * weights
    return gamma * ((treated * u).sum(1) / treated.sum(1) - (control * u).sum(1) / control.sum(1))


# ============================================================================ synthetic SCM 1
#   X = l U + 0.9 N in R^10, N ~ N(0, I); l loads U on the first two covariates (+-0.5), so Var(U | X) = 0.62
#   Y = m(X) + 2 U + A tau(X) + 0.7 eps;  tau has mean 0.8 and variance 3 * 0.7^2 over the trial covariates
#   nu0 makes the root mean square of delta(X) equal to 1.5 times the standard deviation of tau(X);
#        the bias sweep keeps everything else and lowers this size down to 0, where the OBS is unconfounded
#   shift: the OBS covariates are N(a d, Sigma) instead of N(0, Sigma), Sigma = Cov(X), d' Sigma^{-1} d = 1;
#          this is the tilt exp(a d' Sigma^{-1} x), and its effective sample fraction is exp(-a^2)
DIM, LOADING, NOISE_SD = 10, 0.5, 0.9
GAMMA_U, SIGMA = 2.0, 0.7
TAU_MEAN, TAU_VARIANCE = 0.8, 3.0 * SIGMA ** 2
ETA_COEF = np.array([0.5, -0.4, 0.3, 0.3, 0.2, -0.2, 0.0, 0.0, 0.0, 0.0])
RMS_BIAS_IN_SD_TAU = 1.5                          # SCM 1
NU0 = {0.0: 0.0, 0.25: 0.18024124161053823, 0.5: 0.3942782215652453, 0.75: 0.6762819794571086,
       1.0: 1.0876807492904734, 1.25: 1.7948417294819745, 1.5: 3.5061666324764635}
#                                                   nu0 for each size of the bias: SyntheticSCM().calibrate_nu0(size)


def tau_shape(x):
    return (1.2 * x[:, 0] - 0.8 * x[:, 1] + 0.9 * np.sin(x[:, 2]) + 0.7 * x[:, 0] * x[:, 3]
            + 0.5 * x[:, 1] ** 2 + 0.6 * x[:, 4] - 0.5 * x[:, 3] * x[:, 5])


class SyntheticSCM:
    laws = tuple(EFFECTIVE_SAMPLE_FRACTION) + ("misspec_shift",)

    def __init__(self, bias=RMS_BIAS_IN_SD_TAU, noise_sd=SIGMA, tau_variance=TAU_VARIANCE):
        """bias: the root mean square of delta(X) in units of sd tau(X).  1.5 is SCM 1; the other sizes of NU0
        are the bias sweep.  Only the OBS treatment depends on it.  noise_sd and tau_variance set the outcome
        noise and Var_R tau(X) apart (scm1_hethi: 0.2 and 4 * 1.47); nu0 stays NU0[bias] in both cases, so a
        larger tau_variance keeps the absolute size of delta and lowers its size in units of sd tau(X)."""
        self.noise_sd = noise_sd
        rng = np.random.default_rng(3026)                       # the coefficients of the design
        self.load = np.zeros(DIM)
        self.load[:2] = LOADING * np.where(rng.random(DIM) < 0.5, -1.0, 1.0)[:2]
        self.precision = np.linalg.inv(np.outer(self.load, self.load) + NOISE_SD ** 2 * np.eye(DIM))
        self.pairs = [(i, j) for i in range(6) for j in range(i + 1, 6)]
        self.b_linear = rng.normal(0, 0.35, DIM)
        self.b_pair = rng.normal(0, 0.35, len(self.pairs))
        self.b_sin = rng.normal(0, 0.4, DIM)
        self.b_square = rng.normal(0, 0.2, DIM)
        direction = rng.normal(size=DIM)                        # direction d of the covariate shift
        self.shift_direction = direction / np.sqrt(direction @ self.precision @ direction)
        # scales fixed on a large reference sample of trial covariates: tau = t0 + t1 tau_shape, and eta_std
        reference = self.covariates(self.draw(400_000, np.random.default_rng(1)))
        shape, eta = tau_shape(reference), self.eta(reference)
        self.t1 = float(np.sqrt(tau_variance / shape.var()))
        self.t0 = float(TAU_MEAN - self.t1 * shape.mean())
        self.eta_mean, self.eta_sd = float(eta.mean()), float(eta.std())
        self.nu0 = NU0[bias]
        # truth: theta_0 = E_R tau(X) in closed form (E tau_shape = 0.5 Var of the second covariate),
        # and a grid of trial covariates for the CATE risk
        self.theta = self.t0 + self.t1 * 0.5 * (LOADING ** 2 + NOISE_SD ** 2)
        self.truth_x = self.covariates(self.draw(50_000, np.random.default_rng(7001)))
        self.truth_tau = self.tau(self.truth_x)

    # ---- structural functions
    def baseline(self, x):
        pair_terms = np.column_stack([x[:, i] * x[:, j] for i, j in self.pairs])
        return x @ self.b_linear + np.sin(x) @ self.b_sin + (x ** 2 - 1.0) @ self.b_square + pair_terms @ self.b_pair

    def tau(self, x):
        return self.t0 + self.t1 * tau_shape(x)

    def eta(self, x):
        z = x.copy()
        z[:, 2] = np.tanh(z[:, 2])
        return -0.3 + z @ ETA_COEF

    def nu(self, x, nu0=None):
        """Confounding strength nu(x) = nu0 {1 + eta_std(x)}."""
        return (self.nu0 if nu0 is None else nu0) * (1.0 + (self.eta(x) - self.eta_mean) / self.eta_sd)

    # ---- sampling
    @staticmethod
    def draw(rows, rng):
        """The random numbers of `rows` units.  One draw can be realised under several covariate laws."""
        return {"u": rng.normal(size=rows), "noise": rng.normal(size=(rows, DIM)),
                "a_uniform": rng.random(rows), "eps": rng.normal(size=rows)}

    def covariates(self, numbers, shift="same"):
        return self._u_and_x(numbers, shift)[1]

    def _u_and_x(self, numbers, shift):
        if shift == "misspec_shift":
            # h = d' Sigma^{-1} x = alpha u + beta' noise with alpha^2 + |beta|^2 = 1; the tilt exp{a q(h)} acts on
            # (u, noise) through h alone, so it stretches (u, noise) along (alpha, beta) by the factor MISSPEC_SD
            alpha = self.shift_direction @ self.precision @ self.load
            beta = NOISE_SD * self.precision @ self.shift_direction
            h = alpha * numbers["u"] + numbers["noise"] @ beta
            u = numbers["u"] + (MISSPEC_SD - 1.0) * h * alpha
            noise = numbers["noise"] + (MISSPEC_SD - 1.0) * np.outer(h, beta)
            return u, np.outer(u, self.load) + NOISE_SD * noise
        a = np.sqrt(-np.log(EFFECTIVE_SAMPLE_FRACTION[shift]))
        c = a * self.precision @ self.shift_direction           # the tilt exp(c'x) moves E[X] to a d
        u = numbers["u"] + self.load @ c
        noise = numbers["noise"] + NOISE_SD * c
        return u, np.outer(u, self.load) + NOISE_SD * noise

    def _sample(self, numbers, u, x, propensity):
        a = (numbers["a_uniform"] < propensity).astype(float)
        y = self.baseline(x) + GAMMA_U * u + a * self.tau(x) + self.noise_sd * numbers["eps"]
        return {"x": x, "a": a, "y": y}

    def trial(self, numbers):
        u, x = self._u_and_x(numbers, "same")
        return self._sample(numbers, u, x, E_TRIAL)

    def obs(self, numbers, shift="same"):
        u, x = self._u_and_x(numbers, shift)
        return self._sample(numbers, u, x, obs_propensity(self.eta(x), self.nu(x), u))

    # ---- confounding
    def delta(self, x, nu0=None):
        """delta(x), with U | X = x ~ N(x' Sigma^{-1} l, 1 - l' Sigma^{-1} l)."""
        u_mean, u_var = x @ self.precision @ self.load, 1.0 - self.load @ self.precision @ self.load
        return confounding_bias(self.eta(x), self.nu(x, nu0), u_mean, u_var, GAMMA_U)

    def calibrate_nu0(self, bias=RMS_BIAS_IN_SD_TAU):
        """nu0 with sqrt(E_R delta(X)^2) = bias * sd_R tau(X), on a fixed sample of 400,000 trial covariates."""
        blocks = [self.covariates(self.draw(20_000, np.random.default_rng([9000, 1, j]))) for j in range(20)]
        x = np.vstack(blocks)
        sd_tau = self.tau(x).std()
        return brentq(lambda nu0: np.sqrt(np.mean(self.delta(x, nu0) ** 2)) / sd_tau - bias, 0.0, 40.0, xtol=1e-10)


# ============================================================================ synthetic SCM 3: a useless OBS
#   the trial of SCM 1, and an OBS whose outcome and treatment mechanisms are those of SCM 1 applied to the
#   covariates in reverse order (x_10, ..., x_1).  The OBS still has unmeasured confounding through U, but its
#   baseline, its treatment effect and its treatment selection depend on other covariates than the trial's, so
#   its outcome regression and its effect prediction carry no information about the trial: the two channels
#   of the fusion (lambda and omega, or g in the sieve) have nothing to use, and a naive pool is misled.
class UnrelatedSCM(SyntheticSCM):
    def obs(self, numbers, shift="same"):
        u, x = self._u_and_x(numbers, shift)
        mirrored = x[:, ::-1]
        a = (numbers["a_uniform"] < obs_propensity(self.eta(mirrored), self.nu(mirrored), u)).astype(float)
        y = self.baseline(mirrored) + GAMMA_U * u + a * self.tau(mirrored) + self.noise_sd * numbers["eps"]
        return {"x": x, "a": a, "y": y}


# ============================================================================ benchmarks: IHDP, ACIC 2016
#   real covariates and the published noiseless response surfaces mu_0, mu_1; U and the noise are ours
#   trial law of X: uniform on the rows of the data set;  U ~ N(0, 1) independent of X
#   Y = mu_A(X) + U + eps, eps ~ N(0, 1);  mu is scaled so that Var_R mu_A(X) = 1, hence Var_R(Y) = 3
#   nu0 makes |E_R delta(X)| = 0.50 sd_R(Y)
#   shift: the OBS draws row i with probability proportional to exp(a h_i), where h is the standardised
#          first principal component of the covariates
DATASETS = Path(__file__).parent / ".cache"       # downloaded on first use
URLS = {"ihdp": "https://www.fredjo.com/files/",
        "acic2016": "https://raw.githubusercontent.com/IBM/causallib/master/causallib/datasets/data/"
                    "acic_challenge_2016/"}
SHA256 = {
    "ihdp_npci_1-100.train.npz": "750697c71b4f8d7a3aafff771b56a4ac4cd83ec649bf69afb04f8a5aee41a240",
    "ihdp_npci_1-100.test.npz": "a70a8acbcc4e8deb677cc9bf9e9dabeb17caaa37cdbb1d7ba06be7ffb929c41c",
    "acic2016/x.csv": "0d6387ad45d23e54b11cbf967248fb68ad5e578db80d193208b8b4c98b1ed0c1",
    "acic2016/zymu_1.csv": "002fad96bb55d54edec08a0af46ef4f8ad69e0cfe6294ab48f351f22976a6e76",
    "acic2016/zymu_2.csv": "bc3a03bb461bb1629add09c01dfc77053cf36e8faf8e48a7645b35ddb4005d1c",
    "acic2016/zymu_3.csv": "8090562d3fcf5a3b02ff5493592a2880b95b106e1c12e2c3a6e6c28ae8468ff6",
    "acic2016/zymu_4.csv": "61b405c545a451b7e2ce31614150834f94aa4b27069e2755876ed6e11cfee5fd",
    "acic2016/zymu_5.csv": "e350ac76d0b8e83f45f7f04a4b76b9a32a0d53098452075617e61ce62ad0463a",
    "acic2016/zymu_6.csv": "a059950f0b51056660364eca14524a3d3254f93c608d8a1df9240c2c4179a7b5",
    "acic2016/zymu_7.csv": "ac102c3f34680ff901c30fb95985297eea0f20efaf02108f287ceec60b6a2a1f",
    "acic2016/zymu_8.csv": "bbcde6ad666640adcd39fd50cab00865d22b5c463ce9db6d0940d8eccd5afcf3",
    "acic2016/zymu_9.csv": "02663b776fa87cd3d91ee7871e98c1b4d63579957b12ef9bb1ecf367b7afc366",
    "acic2016/zymu_10.csv": "c8ee8ca537ef4315a0ee5753bd01d792a1b3cdacb1fd69e1288b2bb89944d1b7",
    "intfrt/lungcancer.csv": "33469b9a7a6601a203eae7e2773d984a5603075e5fdace6924ffa2cfdd2ea0ee",
    "intfrt/lungcancer_truth.csv": "dc2c91aba1a78385af37f95bd004354c728b9610cdac2703cd7fe8695e4a8626",
}
MEAN_BIAS_IN_SD_Y = 0.50


def dataset(name):
    """Path of a data file.  It is downloaded with curl on first use, and it must match its registered SHA-256."""
    path = DATASETS / name
    if not path.exists() and name.startswith("intfrt/"):
        raise RuntimeError(f"{path} is missing; it is the package data of github.com/ke-zhu/intFRT written with R "
                           "(see the lung cancer section of data.py)")
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        url = URLS["acic2016" if name.startswith("acic2016/") else "ihdp"] + path.name
        subprocess.run(["curl", "-fsSL", "-o", str(path), url], check=True)
    if hashlib.sha256(path.read_bytes()).hexdigest() != SHA256[name]:
        raise RuntimeError(f"{path} does not match its registered SHA-256")
    return path


def load_ihdp(index):
    """IHDP (Hill 2011), realisation `index` in 0..99: 747 units, 25 covariates."""
    parts = [np.load(dataset(f"ihdp_npci_1-100.{part}.npz")) for part in ("train", "test")]
    x = np.concatenate([part["x"][:, :, index] for part in parts]).astype(float)
    return x, np.concatenate([p["mu0"][:, index] for p in parts]), np.concatenate([p["mu1"][:, index] for p in parts])


def load_acic(index):
    """ACIC 2016 (Dorie et al. 2019), setting `index` in 1..10: 4,802 units, 58 covariates, of which the
    three categorical ones are one-hot encoded."""
    x = pd.read_csv(dataset("acic2016/x.csv"))
    categorical = list(x.select_dtypes(exclude="number").columns)
    x = pd.get_dummies(x, columns=categorical, drop_first=True, dtype=float).to_numpy(float)
    surface = pd.read_csv(dataset(f"acic2016/zymu_{index}.csv"))
    return x, surface["mu0"].to_numpy(float), surface["mu1"].to_numpy(float)


class Benchmark:
    laws = tuple(EFFECTIVE_SAMPLE_FRACTION) + ("misspec_shift",)

    def __init__(self, name, index):
        x, mu0, mu1 = load_ihdp(index) if name == "ihdp" else load_acic(index)
        e = E_TRIAL
        # scale the surface so that the noiseless outcome mu_A(X) has variance one under the trial law
        scale = np.sqrt(np.mean(e * mu1 ** 2 + (1 - e) * mu0 ** 2) - np.mean(e * mu1 + (1 - e) * mu0) ** 2)
        self.x, self.mu0, self.mu1 = x, mu0 / scale, mu1 / scale
        z = (x[:, :6] - x[:, :6].mean(0)) / x[:, :6].std(0)
        self.eta = -0.2 + 0.5 * z[:, 0] - 0.4 * z[:, 1] + 0.3 * z[:, 2]

        # row laws: uniform for the trial, tilted by exp(a h) for a shifted OBS
        sd = x.std(0)
        sd[sd == 0] = 1.0
        standardised = (x - x.mean(0)) / sd
        h = standardised @ np.linalg.svd(standardised, full_matrices=False)[2][0]
        h = (h - h.mean()) / h.std()
        self.row_law = {"same": np.full(len(x), 1.0 / len(x))}
        for shift in ("weak_shift", "strong_shift"):
            # for the tilt exp(a h) of a uniform law, 1 / E_O[r_0^2] = 1 / {mean(e^{a h}) mean(e^{-a h})}
            fraction = lambda a: 1.0 / (np.mean(np.exp(a * h)) * np.mean(np.exp(-a * h)))
            a = brentq(lambda a: fraction(a) - EFFECTIVE_SAMPLE_FRACTION[shift], 0.0, 20.0, xtol=1e-14)
            self.row_law[shift] = np.exp(a * h) / np.exp(a * h).sum()
        q = (h ** 2 - 1.0) / np.sqrt(2.0)                       # the misspecified shift, a > 0
        fraction = lambda a: 1.0 / (np.mean(np.exp(a * q)) * np.mean(np.exp(-a * q)))
        a = brentq(lambda a: fraction(a) - MISSPEC_FRACTION, 0.0, 20.0, xtol=1e-14)
        self.row_law["misspec_shift"] = np.exp(a * q) / np.exp(a * q).sum()

        # confounding per row: nu = nu0 {1 + eta_std}, with nu0 such that |E_R delta| = 0.50 sd_R(Y),
        # where Var_R(Y) = Var mu_A(X) + Var(U) + Var(eps) = 3
        eta_std, zero = (self.eta - self.eta.mean()) / self.eta.std(), np.zeros(len(x))
        bias = lambda nu0: confounding_bias(self.eta, nu0 * (1.0 + eta_std), zero, 1.0, 1.0)
        nu0 = brentq(lambda nu0: abs(bias(nu0).mean()) / np.sqrt(3.0) - MEAN_BIAS_IN_SD_Y, 0.0, 80.0, xtol=1e-10)
        self.nu, self.delta = nu0 * (1.0 + eta_std), bias(nu0)

        # truth: the trial law is uniform on the rows, so averages over the rows are exact
        self.truth_x, self.truth_tau = self.x, self.mu1 - self.mu0
        self.theta = float(self.truth_tau.mean())

    @staticmethod
    def draw(rows, rng):
        return {"row_uniform": rng.random(rows), "u": rng.normal(size=rows),
                "a_uniform": rng.random(rows), "eps": rng.normal(size=rows)}

    def _sample(self, numbers, shift, confounded):
        cdf = np.cumsum(self.row_law[shift])
        cdf[-1] = 1.0
        rows = np.searchsorted(cdf, numbers["row_uniform"], side="right")
        u = numbers["u"]
        propensity = obs_propensity(self.eta[rows], self.nu[rows], u) if confounded else E_TRIAL
        a = (numbers["a_uniform"] < propensity).astype(float)
        y = np.where(a == 1, self.mu1[rows], self.mu0[rows]) + u + numbers["eps"]
        return {"x": self.x[rows], "a": a, "y": y}

    def trial(self, numbers):
        return self._sample(numbers, "same", confounded=False)

    def obs(self, numbers, shift="same"):
        return self._sample(numbers, shift, confounded=True)


# ============================================================================ lung cancer: the intFRT RCT-EHR pair
#   a semi-synthetic pair built on the CALGB 9633 trial and an EHR cohort (R package intFRT,
#   github.com/ke-zhu/intFRT, MIT licence): covariates by synthpop, survival times by survival random forests fitted within
#   each cohort and arm, and in lungcancer_truth the potential survival times T1, T0 (years) of every subject.
#   The two cohorts differ in covariate law, treatment mechanism and outcome surface as the generator made
#   them: nothing is added here, and the five covariates both cohorts share are used.
#   trial law:  uniform on the 335 C9633 rows; treatment randomised here; outcome the potential time T_A
#   OBS law:    uniform on the 16,217 EHR rows; treatment and outcome T_treat as in the data (no censoring)
#   truth:      theta_0 = mean of T1 - T0 over the trial rows; truth_tau = the individual effects T1 - T0,
#               so a CATE risk contains the variance of the individual effect around the CATE, the same for
#               every learner; compare learners by differences of risks
#   files:      the package's lungcancer.rda and lungcancer_truth.rda, written to .cache/intfrt/ with R 4.4:
#               load("lungcancer.rda"); load("lungcancer_truth.rda"); write.csv(lungcancer, "lungcancer.csv",
#               row.names = FALSE); write.csv(lungcancer_truth, "lungcancer_truth.csv", row.names = FALSE)
LUNG_COVARIATES = ["sex", "age", "race", "hist", "tsize"]


def load_lungcancer():
    """(trial rows, EHR rows) of the lung cancer pair, each with its potential survival times T0 and T1."""
    data = pd.read_csv(dataset("intfrt/lungcancer.csv")).merge(pd.read_csv(dataset("intfrt/lungcancer_truth.csv")),
                                                              on="patid", validate="one_to_one")
    return data[data.cohort == "C9633"], data[data.cohort == "EHR"]


class LungCancer:
    laws = ("native",)                            # the EHR covariate law as it is: r_0 != 1

    def __init__(self):
        trial, ehr = load_lungcancer()
        self.x_R = trial[LUNG_COVARIATES].to_numpy(float)
        self.t0_R, self.t1_R = trial["T0"].to_numpy(float), trial["T1"].to_numpy(float)
        self.x_O, self.a_O = ehr[LUNG_COVARIATES].to_numpy(float), ehr["treat"].to_numpy(float)
        self.y_O = np.where(self.a_O == 1, ehr["T1"].to_numpy(float), ehr["T0"].to_numpy(float))
        self.truth_x, self.truth_tau = self.x_R, self.t1_R - self.t0_R
        self.theta = float(self.truth_tau.mean())

    @staticmethod
    def draw(rows, rng):
        return {"row_uniform": rng.random(rows), "a_uniform": rng.random(rows)}

    def trial(self, numbers):
        rows = (numbers["row_uniform"] * len(self.x_R)).astype(int)
        a = (numbers["a_uniform"] < E_TRIAL).astype(float)
        return {"x": self.x_R[rows], "a": a, "y": np.where(a == 1, self.t1_R[rows], self.t0_R[rows])}

    def obs(self, numbers, law="native"):
        rows = (numbers["row_uniform"] * len(self.x_O)).astype(int)
        return {"x": self.x_O[rows], "a": self.a_O[rows], "y": self.y_O[rows]}
