"""Data sources for the DataFusionPPI experiments (handoff Section 6).

Three families:
  * synthetic structural causal models with knobs for confounding, covariate
    shift, and effect heterogeneity (reuses ssem_ate_pilot.py);
  * the Tennessee STAR file, used twice: as a real covariate support with a
    synthetic effect, and as a real trial with real outcomes for CATE-2;
  * the NSW experimental sample with the CPS and PSID comparison groups.
Every download is hash-checked; no data file is committed.
"""
from __future__ import annotations

import dataclasses
import hashlib
import io
import math
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

import ssem_ate_pilot as legacy
from fusion_core import BASE_SEED, rng_for, take

CACHE = Path("/private/tmp/claude-502/datafusionppi-data")
STAR_URL = "https://vincentarelbundock.github.io/Rdatasets/csv/AER/STAR.csv"
STAR_SHA = "0e8b179ea3d883730b25008ca1293c4c8f886aa8d5d299c03ab6ff05d0123ae6"
NSW_FILES = {
    "nsw_dw.dta": ("https://users.nber.org/~rdehejia/data/nsw_dw.dta",
                   "d1bd2680a1c6f799f1c6d2455bf29633fdf19be01cb19490621c20a560b4e072", 445),
    "cps_controls.dta": ("https://users.nber.org/~rdehejia/data/cps_controls.dta",
                         "80f0123eaf723870bd060ba9f8569617a7cb519a462dc20badbfb3443cd8374b", 15992),
    "psid_controls.dta": ("https://users.nber.org/~rdehejia/data/psid_controls.dta",
                          "7beebae8928035d6abf662b994ce77a4fbea7ae6575037ff94e92fe099b43e14", 2490),
}
STRATUM_MIN = 40          # handoff Section 5.4: merge rarer strata
NSW_REFERENCE = 1794.3    # full-experiment difference in means of re78


def _download(url: str, sha: str, name: str) -> bytes:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / name
    if path.exists():
        payload = path.read_bytes()
    else:
        payload = subprocess.run(["curl", "-fsSL", url], check=True, capture_output=True,
                                 timeout=180).stdout
        path.write_bytes(payload)
    digest = hashlib.sha256(payload).hexdigest()
    if digest != sha:
        raise RuntimeError(f"{name} hash mismatch: {digest}")
    return payload


# ------------------------------------------------------------------- SCM
def scm_family(scm_id: int, confounding: float = 1.0, shift: float = 1.0,
               heterogeneity: float = 1.0) -> legacy.SCMParameters:
    """Knobs of handoff Sections 5.1 to 5.3, applied to a base SCM."""
    scm = legacy.generate_scm_parameters(scm_id)
    return dataclasses.replace(
        scm,
        obs_latent_coefficient=scm.obs_latent_coefficient * confounding,
        rct_mean_shift=scm.rct_mean_shift * shift,
        treatment_coefficients=scm.treatment_coefficients * heterogeneity,
    )


def sample_scm(scm: legacy.SCMParameters, study: str, n: int, regime: str,
               rng: np.random.Generator) -> dict[str, np.ndarray]:
    data = legacy.sample_dataset(rng, n, study, regime, scm)
    data["true_ratio"] = legacy.oracle_ratio(data["x"], regime, scm)
    return data


def scm_oracle_ratio(scm: legacy.SCMParameters, regime: str):
    return lambda x: legacy.oracle_ratio(x, regime, scm)


def scm_true_theta(scm: legacy.SCMParameters, regime: str) -> float:
    return legacy.true_theta(scm, regime)


def scm_tau(x: np.ndarray, scm: legacy.SCMParameters) -> np.ndarray:
    return legacy.tau_function(x, scm)


def scm_test_sample(scm: legacy.SCMParameters, regime: str, n: int,
                    rng: np.random.Generator) -> dict[str, np.ndarray]:
    """RCT-population covariates with the true conditional effect."""
    u = rng.normal(size=n)
    mean = scm.rct_mean(regime)
    x = mean + np.outer(u, scm.latent_loading) + rng.normal(size=(n, 4)) * scm.covariate_noise_sd
    return {"x": x, "tau": legacy.tau_function(x, scm)}


def overlap_check(scm: legacy.SCMParameters, regime: str, n: int = 100000) -> dict[str, float]:
    rng = rng_for(99, scm.scm_id)
    x = scm.rct_mean(regime) + np.outer(rng.normal(size=n), scm.latent_loading) \
        + rng.normal(size=(n, 4)) * scm.covariate_noise_sd
    r = legacy.oracle_ratio(x, regime, scm)
    return {"min": float(r.min()), "max": float(r.max()), "mean": float(r.mean())}


# ------------------------------------------------------- STAR real covariates
def star_realx_support():
    """The synthetic-effect STAR family already used by the existing benchmark."""
    import total_budget_nested_benchmark as bench
    return bench.load_star_support()


def star_realx_sample(support, regime: str, n_rct: int, n_obs: int, rep: int):
    import total_budget_nested_benchmark as bench
    saved = bench.N_OBS_FIXED, bench.OBS_MULTIPLIER
    bench.N_OBS_FIXED, bench.OBS_MULTIPLIER = n_obs, None
    try:
        rct, obs, theta = bench.star_sample(support, regime, n_rct, rep)
    finally:
        bench.N_OBS_FIXED, bench.OBS_MULTIPLIER = saved
    return rct, obs, theta


def star_realx_test(support, regime: str, n: int, rep: int) -> dict[str, np.ndarray]:
    p_r, _ = support.probabilities(regime)
    rng = rng_for(6200, rep, n)
    idx = rng.choice(len(support.x), size=n, replace=True, p=p_r)
    return {"x": support.x[idx], "tau": support.tau[idx]}


# --------------------------------------------------- STAR real outcomes (CATE-2)
@dataclasses.dataclass(frozen=True)
class StarCohort:
    """The kindergarten cohort, organised by distinct covariate pattern.

    Both sources of CATE-2 are drawn from one explicitly specified law over the
    retained patterns, so the covariate distributions are equal by construction
    and the transport weight is identically one as a statement about
    populations, not about realised sample shares.
    """
    pattern_x: np.ndarray        # one row per retained covariate pattern
    pattern_law: np.ndarray      # the specified probability of each pattern
    pattern_propensity: np.ndarray   # within-pattern share of the small class
    members: list                # cohort row indices of each (pattern, arm) cell
    tilt: list                   # residual-rank tilt weights of each cell
    a: np.ndarray
    y: np.ndarray
    ate_reference: float
    outcome: str
    control: str
    retained_share: float
    raw_cohort_difference: float
    pattern_contrast: np.ndarray


STAR_MIN_PER_ARM = 5


def star_cohort(outcome: str = "mathk", control: str = "pooled") -> StarCohort:
    payload = _download(STAR_URL, STAR_SHA, "STAR.csv")
    frame = pd.read_csv(io.BytesIO(payload))
    if len(frame) != 11598:
        raise RuntimeError("STAR row count mismatch")
    keep = frame["stark"].notna() & frame["readk"].notna() & frame["mathk"].notna()
    if control == "regular":
        keep &= frame["stark"].isin(["small", "regular"])
    df = frame[keep].copy()
    df["treat"] = (df["stark"] == "small").astype(float)

    eth = df["ethnicity"].where(df["ethnicity"].isin(["cauc", "afam"]), "other")
    cats = pd.DataFrame({"gender": df["gender"].astype(str),
                         "ethnicity": eth.astype(str),
                         "lunchk": df["lunchk"].fillna("Missing").astype(str),
                         "schoolk": df["schoolk"].astype(str)})
    x = pd.get_dummies(cats, drop_first=True, dtype=float).to_numpy()
    y_raw = df[outcome].to_numpy(dtype=float)
    y = (y_raw - y_raw.mean()) / y_raw.std()
    a = df["treat"].to_numpy(dtype=float)

    patterns, inverse = np.unique(x, axis=0, return_inverse=True)
    retained, counts = [], []
    for k in range(len(patterns)):
        mask = inverse == k
        if int(a[mask].sum()) >= STAR_MIN_PER_ARM and int((1 - a[mask]).sum()) >= STAR_MIN_PER_ARM:
            retained.append(k)
            counts.append(int(mask.sum()))
    counts = np.array(counts, dtype=float)
    law = counts / counts.sum()

    # the design-side residual, used only to build the observational selection
    from sklearn.linear_model import Ridge
    design = np.column_stack([x, a, a[:, None] * x])
    residual = y - Ridge(alpha=1.0).fit(design, y).predict(design)

    members, tilt, propensity = [], [], []
    for k in retained:
        mask = inverse == k
        rows = np.where(mask)[0]
        propensity.append(float(a[rows].mean()))
        cell, weights = [], []
        for arm in (0.0, 1.0):
            idx = rows[a[rows] == arm]
            res = residual[idx]
            rank = (np.argsort(np.argsort(res)) + 0.5) / len(idx)
            cell.append(idx)
            weights.append(2 * rank - 1)          # in (-1, 1), mean zero by construction
        members.append(cell)
        tilt.append(weights)

    # The experiments draw both sources from `law`, so the target is the
    # law-weighted within-pattern contrast, not the raw cohort difference.  The
    # raw difference is kept as a separate field for reference only.
    per_pattern = np.array([y[members[i][1]].mean() - y[members[i][0]].mean()
                            for i in range(len(retained))])
    reference = float(np.sum(law * per_pattern))
    raw_difference = float(y[a == 1].mean() - y[a == 0].mean())
    return StarCohort(pattern_x=patterns[retained], pattern_law=law,
                      pattern_propensity=np.array(propensity), members=members, tilt=tilt,
                      a=a, y=y, ate_reference=reference, outcome=outcome, control=control,
                      retained_share=float(counts.sum() / len(x)),
                      raw_cohort_difference=raw_difference,
                      pattern_contrast=per_pattern)


def star_cate2_replication(cohort: StarCohort, n_rct: int, alpha: float, rep: int,
                           n_eval: int = 2000, n_obs: int = 1200, n_rct_ate: int | None = None):
    """Draw both sources from the same pattern law (handoff CATE-2, revised).

    Trial and evaluation students are drawn uniformly inside a pattern, so their
    treatment follows the real within-pattern randomisation.  Observational
    students are drawn inside a pattern and arm with a weight that rises with the
    outcome residual for treated units and falls for controls, which is the
    unmeasured confounding.  The selection weight has mean one inside every
    pattern and arm, so the pattern law is unchanged and the transport weight is
    identically one.
    """
    rng = rng_for(7100, n_rct, int(alpha * 100), rep)
    law, prop = cohort.pattern_law, cohort.pattern_propensity

    def draw_trial(size: int) -> dict[str, np.ndarray]:
        pats = rng.choice(len(law), size=size, p=law)
        rows = np.empty(size, dtype=int)
        for j, k in enumerate(pats):
            pool = np.concatenate(cohort.members[k])
            rows[j] = pool[rng.integers(len(pool))]
        return {"x": cohort.pattern_x[pats], "a": cohort.a[rows], "y": cohort.y[rows],
                "propensity": prop[pats], "pattern": pats,
                "true_ratio": np.ones(size)}

    def draw_observational(size: int) -> dict[str, np.ndarray]:
        pats = rng.choice(len(law), size=size, p=law)
        arms = (rng.random(size) < prop[pats]).astype(float)
        rows = np.empty(size, dtype=int)
        for j, (k, arm) in enumerate(zip(pats, arms)):
            slot = int(arm)
            idx = cohort.members[k][slot]
            signed = cohort.tilt[k][slot] * (1.0 if arm == 1.0 else -1.0)
            weight = 1.0 + alpha * signed
            rows[j] = idx[rng.choice(len(idx), p=weight / weight.sum())]
        return {"x": cohort.pattern_x[pats], "a": arms, "y": cohort.y[rows],
                "propensity": prop[pats], "pattern": pats,
                "true_ratio": np.ones(size)}

    rct, evaluation, obs = draw_trial(n_rct), draw_trial(n_eval), draw_observational(n_obs)
    ate_rct = draw_trial(n_rct_ate) if n_rct_ate else None
    naive = float(obs["y"][obs["a"] == 1].mean() - obs["y"][obs["a"] == 0].mean())
    diagnostics = {"obs_size": len(obs["x"]), "naive_obs_contrast": naive,
                   "law_weighted_reference": cohort.ate_reference,
                   "raw_cohort_difference": cohort.raw_cohort_difference,
                   "retained_cohort_share": cohort.retained_share,
                   "patterns": len(law),
                   "n_trial_total": len(rct["x"]), "n_eval_pool": len(evaluation["x"]),
                   "n_trial_for_ate": len(ate_rct["x"]) if ate_rct else 0}
    return rct, obs, {"eval": evaluation, "ate_rct": ate_rct, "diagnostics": diagnostics}



# ------------------------------------------------- bounded design for Theorem 6
BOUNDED_PROPENSITY = (0.3, 0.7)
BOUNDED_OUTCOME = 1.0
BOUNDED_MU = 0.20              # sup of the baseline mean
BOUNDED_TAU = 0.20             # sup of the effect at heterogeneity one
BOUNDED_U = 0.30               # sup of the latent term
BOUNDED_NOISE = 0.10           # sup of the noise
BOUNDED_PREDICTION = 1.0       # the outcome model is clipped to the outcome range
BOUNDED_CANDIDATE = 2.0        # the prediction and every candidate are clipped here


@dataclasses.dataclass(frozen=True)
class BoundedFamily:
    """A design whose scores are bounded before any data is drawn.

    Theorem 6 assumes the validation scores lie in intervals of fixed length,
    which a Gaussian outcome never satisfies.  Here every potential outcome
    lies in [-1, 1] by construction rather than by a clip applied afterwards,
    so the returned effect is the exact difference of potential outcomes:

        sup|mu| + sup|tau| + sup|u term| + sup|noise|
            = 0.20 + 0.20 h + 0.30 + 0.10 <= 1  for h <= 2.

    With the outcome model clipped to [-1, 1] and the prediction and every
    candidate clipped to [-2, 2], the pseudo-outcome obeys

        |Z0| <= 2 * 1 + (1 + 1) / 0.3 = 8.667 < 9,

    so the common bound is 9 and the ratio cap is one.
    """
    seed: int
    confounding: float = 1.0
    heterogeneity: float = 1.0
    dimension: int = 4

    def __post_init__(self):
        budget = BOUNDED_MU + BOUNDED_TAU * self.heterogeneity + BOUNDED_U + BOUNDED_NOISE
        if budget > BOUNDED_OUTCOME + 1e-12:
            raise ValueError(f"the outcome budget {budget:.3f} exceeds the bound")

    @property
    def score_bound(self) -> float:
        return 2 * BOUNDED_PREDICTION + (BOUNDED_OUTCOME + BOUNDED_PREDICTION) / BOUNDED_PROPENSITY[0]

    @property
    def score_range(self) -> tuple[float, float]:
        """The (B, ratio cap) pair of Theorem 6, rounded up to a round number."""
        return 9.0, 1.0

    def _coefficients(self):
        rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, 9100, self.seed]))
        d = self.dimension
        return (rng.normal(0, 1, d), rng.normal(0, 1, d), rng.normal(0, 1, d))

    def tau(self, x: np.ndarray) -> np.ndarray:
        _, effect, _ = self._coefficients()
        return BOUNDED_TAU * self.heterogeneity * np.tanh(x @ effect / math.sqrt(self.dimension))

    def sample(self, study: str, n: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
        base, effect, assign = self._coefficients()
        d = self.dimension
        x = rng.uniform(-1.0, 1.0, size=(n, d))
        u = rng.uniform(-1.0, 1.0, size=n)
        tau = self.tau(x)
        mu = BOUNDED_MU * np.tanh(x @ base / math.sqrt(d))
        lo, hi = BOUNDED_PROPENSITY
        score = x @ assign / math.sqrt(d)
        if study == "OBS":
            score = score + 2.0 * self.confounding * u
        e = lo + (hi - lo) * (0.5 + 0.5 * np.tanh(score))
        a = rng.binomial(1, e).astype(float)
        noise = BOUNDED_NOISE * rng.uniform(-1.0, 1.0, size=n)
        y = mu + a * tau + BOUNDED_U * u + noise
        return {"x": x, "a": a, "y": y, "propensity": e, "tau": tau,
                "true_ratio": np.ones(n)}

    def test_sample(self, n: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
        x = rng.uniform(-1.0, 1.0, size=(n, self.dimension))
        return {"x": x, "tau": self.tau(x)}


def bounded_family(seed: int, confounding: float = 1.0, heterogeneity: float = 1.0) -> BoundedFamily:
    return BoundedFamily(seed=seed, confounding=confounding, heterogeneity=heterogeneity)


# ------------------------------------------------------------------- NSW
NSW_COVARIATES = ("age", "education", "black", "hispanic", "married", "nodegree", "re74", "re75")


def _nsw_frame(name: str) -> pd.DataFrame:
    url, sha, rows = NSW_FILES[name]
    _download(url, sha, name)
    frame = pd.read_stata(CACHE / name)
    if len(frame) != rows:
        raise RuntimeError(f"{name} row count mismatch: {len(frame)}")
    return frame


def nsw_design(comparison: str = "cps") -> dict[str, pd.DataFrame]:
    trial = _nsw_frame("nsw_dw.dta")
    control = _nsw_frame("cps_controls.dta" if comparison == "cps" else "psid_controls.dta")
    return {"trial": trial, "comparison": control}


def _nsw_matrix(frame: pd.DataFrame) -> np.ndarray:
    base = frame[list(NSW_COVARIATES)].to_numpy(dtype=float)
    zeros = np.column_stack([(frame["re74"].to_numpy() == 0).astype(float),
                             (frame["re75"].to_numpy() == 0).astype(float)])
    return np.column_stack([base, zeros])


def nsw_replication(design: dict[str, pd.DataFrame], n_rct: int | None, rep: int):
    """Split the trial treated units so that no unit appears in both sources."""
    rng = rng_for(7200, n_rct or 0, rep)
    trial = design["trial"]
    treated = trial[trial["treat"] == 1].reset_index(drop=True)
    controls = trial[trial["treat"] == 0].reset_index(drop=True)
    perm = rng.permutation(len(treated))
    half = len(treated) // 2
    rct_treated, obs_treated = treated.iloc[perm[:half]], treated.iloc[perm[half:]]

    rct_frame = pd.concat([rct_treated, controls], ignore_index=True)
    if n_rct is not None and n_rct < len(rct_frame):
        idx = rng.choice(len(rct_frame), size=n_rct, replace=False)
        rct_frame = rct_frame.iloc[idx].reset_index(drop=True)
    obs_frame = pd.concat([obs_treated, design["comparison"]], ignore_index=True)

    def pack(frame: pd.DataFrame, propensity: float | np.ndarray) -> dict[str, np.ndarray]:
        a = frame["treat"].to_numpy(dtype=float)
        e = np.full(len(frame), propensity) if np.isscalar(propensity) else propensity
        return {"x": _nsw_matrix(frame), "a": a,
                "y": frame["re78"].to_numpy(dtype=float) / 1000.0,
                "propensity": np.clip(e, 0.05, 0.95), "true_ratio": np.ones(len(frame))}

    share = float(rct_frame["treat"].mean())
    obs_share = float(obs_frame["treat"].mean())
    return pack(rct_frame, share), pack(obs_frame, max(obs_share, 0.01)), NSW_REFERENCE / 1000.0
