"""One replication of the simulation of Sec. 5.  simulation.py repeats it; visualization.py draws the results.

Sources: scm1 (synthetic) and ihdp, acic (benchmarks), the bias sweep scm1_bias*, and scm1_hethi (a synthetic
design with a large effect heterogeneity).  One replication draws a trial of 3m rows, the blocks
B1, B2, B3 of m rows, and an OBS sample of N rows, 60% nuisance, 20% tuning and 20% evaluation rows (N = 15,000 by
default: 9000, 3000, 3000).  The same OBS random numbers are realised under every covariate law asked for, so the
laws are compared on the same draw.

Cross-fitting.  Every estimator that splits the trial runs three times, with the blocks rotated through the
roles (nuisance, tuning, evaluation) = (B1, B2, B3), (B2, B3, B1), (B3, B1, B2); the OBS roles stay fixed.  The ATE
is the mean of the three estimates, the CATE the mean of the three fitted functions (their sieves differ, so the
coefficients are not averaged).  The first rotation alone is the single-split version of each estimator.  The
naive pool uses no split; shrinkage and pretest pool it with the cross-fitted trial-only estimate.

Besides Algorithms 1 and 2 and their restrictions, every replication computes the comparison estimators RCT
only, naive pool, shrinkage and pretest (ATE), the two-step learner of Kallus et al. and the integrative
R-learner of Wu and Yang (CATE), and these additional estimators and diagnostics:
    lambda_only         ATE with omega = 0
    fusion_nocal        ATE with omega = D / B, which ignores the calibration variance; Algorithm 1 itself takes
                        omega = D / (B + B_cal) when the ratio is estimated.  PPI++ uses D / B.
    fusion_cls          ATE with the bare classifier ratio, under the misspecified shift only; it takes
                        omega = D / B and its standard error has no calibration term
    trial_sieve_g       CATE at the grid point (0, 0): the trial loss on the sieve with g
    rct_only_full,      the practitioner's trial-only ATE and CATE learners: 3-fold, mu_R fitted on two blocks,
    rct_full            no tuning sample
    DIAG rows           Var of Z^0 and of Z^lambda on the tuning samples, and Var_R tau_0 (diagnostics)

Covariate laws of the OBS (data.py):  "same" has r_0 = 1, which the estimators then use as known
(Assumption 3(i));  "weak_shift", "strong_shift" and "misspec_shift" have r_0 != 1, and every estimator estimates
its own balancing ratio (Assumption 3(ii)).
"""
import numpy as np

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # the package dfppi when it is not installed

from dfppi import aipwf as AIPWF, drf as DRF, rf as RF  # noqa: E402
from dfppi.crossfit import EVALUATED_IN, NUISANCE_IN, ROTATIONS, crossfit_ate  # noqa: E402,F401
from data import (E_TRIAL, NU0, RMS_BIAS_IN_SD_TAU, TAU_VARIANCE, Benchmark, LungCancer, SyntheticSCM,
                  UnrelatedSCM)
from dfppi.method import (aipw_score, fit_nuisances, fit_outcome_regression, fit_propensity, fused_score,
                    principal_components, remember, var)

SOURCES = ("scm1", "ihdp", "acic", "lung", "scm3", "scm1_hethi")
# The bias sweep: SCM 1 with a smaller confounding bias, "scm1_bias0" (no unmeasured confounding) to
# "scm1_bias1.25".  The number is the root mean square of delta in units of sd(tau_0); SCM 1 itself has 1.5.
SWEEP = {f"scm1_bias{bias:g}": bias for bias in NU0 if bias != RMS_BIAS_IN_SD_TAU}
HETHI = {"noise_sd": 0.2, "tau_variance": 4.0 * TAU_VARIANCE}   # scm1_hethi: nu0 of SCM 1 kept
N_OBS = 15_000                                   # OBS rows, split 60 / 20 / 20 into nuis, tune, eval
SHIFTS = ("same", "weak_shift", "strong_shift", "native")   # a design runs the laws it has
EXTRA_SHIFTS = ("misspec_shift",)                           # run only when asked for (misspecified shift)
ATE_METHODS = {"rct_only": (False, False), "omega_only": (False, True), "lambda_only": (True, False),
               "fusion": (True, True)}                              # (use_lambda, use_omega)
CATE_LEARNERS = {"DRF": DRF, "RF": RF}
CATE_GRIDS = {"omega_only": tuple(c for c in DRF.GRID if c[0] == 0.0), "fusion": DRF.GRID}
COLUMNS = ["source", "m", "N", "shift", "rep", "estimand", "learner", "method", "lambda", "omega", "value", "se",
           "covered"]
SEED = 20261005
DESIGNS = {}


def design_for(source, rep):
    """SCM 1, SCM 3, scm1_hethi and each design of the bias sweep are one fixed law.  The benchmarks rotate their
    published surfaces (Appendix C): replication rep uses IHDP realisation rep mod 100 and ACIC setting
    (rep mod 10) + 1."""
    key = {"ihdp": ("ihdp", rep % 100), "acic": ("acic", rep % 10 + 1)}.get(source, (source,))
    if key not in DESIGNS:
        if len(key) == 2:
            DESIGNS[key] = Benchmark(*key)
        elif source == "lung":
            DESIGNS[key] = LungCancer()
        elif source == "scm3":
            DESIGNS[key] = UnrelatedSCM()
        elif source == "scm1_hethi":
            DESIGNS[key] = SyntheticSCM(**HETHI)
        else:
            DESIGNS[key] = SyntheticSCM(SWEEP.get(source, RMS_BIAS_IN_SD_TAU))
    return DESIGNS[key]


def split(sample, sizes):
    cuts = np.cumsum([0] + list(sizes.values()))
    return {role: {k: v[cuts[i]:cuts[i + 1]] for k, v in sample.items()} for i, role in enumerate(sizes)}


def obs_sizes(N):
    nuis, tune = round(0.6 * N), round(0.2 * N)
    return {"nuis": nuis, "tune": tune, "eval": N - nuis - tune}


def rotation(blocks, roles, O, same, mu_O, e_O, misspec, truth_x):
    """Every split-based estimator with the trial blocks in the roles (nuisance, tuning, evaluation).  Returns
    the outputs of aipwf, the ratio of each ATE method at the OBS nuisance rows, the fitted CATE functions at
    truth_x with their (lambda, omega), and the diagnostic score variances on the tuning sample."""
    R_nuis, R_tune, R_eval = (blocks[b] for b in roles)
    nuis = fit_nuisances(R_nuis, O["nuis"], E_TRIAL, same_domain=same, mu_O=mu_O)   # line 1 of both algorithms
    r_ate, r_cate = AIPWF.density_ratio(nuis), DRF.density_ratio(nuis)
    samples = (R_tune, O["tune"], R_eval, O["eval"])
    x_On = O["nuis"]["x"]

    # ATE: Algorithm 1, its restrictions, and two variants
    ate = {name: AIPWF.aipwf(nuis, r_ate, *samples, use_lambda, use_omega, omega_cal=name == "fusion")
           for name, (use_lambda, use_omega) in ATE_METHODS.items()}
    ate["fusion_nocal"] = AIPWF.aipwf(nuis, r_ate, *samples, omega_cal=False)
    ratio_O = {name: r_ate(x_On) for name in ate}
    if misspec:
        ate["fusion_cls"] = AIPWF.aipwf(nuis, nuis.r_base, *samples, omega_cal=False, calibrated=False)
        ratio_O["fusion_cls"] = nuis.r_base(x_On)

    # CATE: Algorithm 2 with each pseudo-outcome, and the comparison learners, as functions at truth_x
    C, S, g_truth = nuis.components(truth_x), DRF.sieve(nuis, truth_x), nuis.g(truth_x)
    cate = {}
    for learner, estimator in CATE_LEARNERS.items():
        cate[(learner, "rct_only")] = (C @ estimator.rct_only(nuis, R_tune), 0.0, 0.0)   # no OBS input, no g
        for name, grid in CATE_GRIDS.items():
            lam, omega, beta = estimator.fit(nuis, r_cate, *samples, grid)
            cate[(learner, name)] = (S @ beta, lam, omega)
        cate[(learner, "two_step")] = (g_truth + C @ estimator.two_step(nuis, R_tune), np.nan, np.nan)
        beta = DRF.fit_candidate(0.0, 0.0, nuis, r_cate, R_tune, O["tune"], pseudo_outcome=estimator.pseudo_outcome)
        cate[(learner, "trial_sieve_g")] = (S @ beta, 0.0, 0.0)
    beta = RF.integrative(nuis, e_O, *samples, c_basis=RF.c_basis(nuis, x_On))   # Wu and Yang: tau in (1, Z)
    cate[("RF", "integrative")] = (C @ beta, np.nan, np.nan)

    diag = {"var_Z0": var(fused_score(R_tune, nuis, 0.0)),
            "var_Zlam": var(fused_score(R_tune, nuis, ate["fusion"]["lambda"]))}
    return {"ate": ate, "ratio_O": ratio_O, "cate": cate, "diag": diag}




def trial_only_full(trial, blocks, truth_x):
    """The practitioner's trial-only estimators, which need no tuning sample.  Each block is scored
    with mu_R fitted on the two other blocks; the ATE is the mean of all 3m scores, and each CATE learner regresses
    all 3m pseudo-outcomes on (1, five principal components of all trial covariates) with the ridge of Alg. 2."""
    components = principal_components(trial["x"], 5)
    scores, pseudo = [], {learner: [] for learner in CATE_LEARNERS}
    for b in range(3):
        train = {key: np.concatenate([blocks[j][key] for j in range(3) if j != b]) for key in ("x", "a", "y")}
        mu = fit_outcome_regression(train, clip_quantile=0.005)
        scores.append(aipw_score(blocks[b], mu, E_TRIAL))
        for learner, estimator in CATE_LEARNERS.items():
            pseudo[learner].append(estimator.pseudo_outcome(blocks[b], mu, E_TRIAL))
    scores = np.concatenate(scores)
    basis = np.vstack([components(block["x"]) for block in blocks])
    C = components(truth_x)
    cate = {}
    for learner, parts in pseudo.items():
        kappa, Z = np.concatenate([k for k, _ in parts]), np.concatenate([z for _, z in parts])
        H = (basis * kappa[:, None]).T @ basis / len(Z)
        cate[learner] = C @ np.linalg.solve(H + 0.01 * np.eye(basis.shape[1]), basis.T @ (kappa * Z) / len(Z))
    return float(scores.mean()), float(np.std(scores, ddof=1) / np.sqrt(len(scores))), cate


def replicate(source, m, rep, shifts=SHIFTS, N=N_OBS):
    """Rows of results, one per (covariate law, estimand, method).  `value` is the ATE error (estimate minus
    theta_0) or the CATE risk E_R (t - tau_0)^2 over the trial covariates."""
    design = design_for(source, rep)
    seeded_as = "scm1" if source in SWEEP or source == "scm3" else source   # the same random numbers as SCM 1
    rng = np.random.default_rng([SEED, SOURCES.index(seeded_as), m, rep])
    trial = design.trial(design.draw(3 * m, rng))
    blocks = list(split(trial, {"B1": m, "B2": m, "B3": m}).values())
    obs_numbers = design.draw(N, rng)
    truth_x, truth_tau = design.truth_x, design.truth_tau
    risk = lambda t: float(np.mean((t - truth_tau) ** 2))
    theta_full, se_full, cate_full = trial_only_full(trial, blocks, truth_x)   # the same under every law
    rows = []
    for shift in shifts:
        if shift not in design.laws:
            continue
        obs = design.obs(obs_numbers, shift)
        O = split(obs, obs_sizes(N))
        same = shift == "same"
        cell = {"source": source, "m": m, "N": N, "shift": shift, "rep": rep}
        mu_O = remember(fit_outcome_regression(O["nuis"], clip_quantile=0.005))   # the same in every rotation
        e_O = fit_propensity(O["nuis"])              # the OBS propensity of IR, the same in every rotation
        try:                                         # line 1 of both algorithms in each rotation
            parts = [rotation(blocks, roles, O, same, mu_O, e_O, shift == "misspec_shift", truth_x)
                     for roles in ROTATIONS]
        except RuntimeError as failure:              # recorded, never replaced silently
            rows.append({**cell, "estimand": "FAILED", "method": str(failure)})
            continue

        # ATE: the cross-fitted Algorithm 1, its restrictions and variants
        g = lambda x: mu_O(x)[1] - mu_O(x)[0]
        g_blocks, g_On = [g(block["x"]) for block in blocks], g(O["nuis"]["x"])
        ate = {}
        for name in parts[0]["ate"]:
            ate[name] = crossfit_ate([(p["ate"][name], p["ratio_O"][name]) for p in parts], g_blocks, g_On,
                                     calibrated_ratio=not same and name != "fusion_cls")
            error, se = ate[name][0] - design.theta, ate[name][1]
            rows.append({**cell, "estimand": "ATE", "method": name,
                         "lambda": np.mean([p["ate"][name]["lambda"] for p in parts]),
                         "omega": np.mean([p["ate"][name]["omega"] for p in parts]),
                         "value": error, "se": se, "covered": float(abs(error) <= 1.96 * se)})
        error = theta_full - design.theta
        rows.append({**cell, "estimand": "ATE", "method": "rct_only_full", "value": error, "se": se_full,
                     "covered": float(abs(error) <= 1.96 * se_full)})

        # CATE: the mean of the three fitted functions
        for (learner, name), _ in parts[0]["cate"].items():
            fits = [p["cate"][(learner, name)] for p in parts]
            row = {**cell, "estimand": "CATE", "learner": learner, "method": name,
                   "value": risk(np.mean([f[0] for f in fits], axis=0))}
            if not np.isnan(fits[0][1]):
                row["lambda"], row["omega"] = np.mean([f[1] for f in fits]), np.mean([f[2] for f in fits])
            rows.append(row)
        for learner, tau_hat in cate_full.items():
            rows.append({**cell, "estimand": "CATE", "learner": learner, "method": "rct_full", "value": risk(tau_hat)})

        # diagnostics: the score variances on the tuning samples, and Var_R tau_0
        for name in ("var_Z0", "var_Zlam"):
            rows.append({**cell, "estimand": "DIAG", "method": name, "value": np.mean([p["diag"][name] for p in parts])})
        rows.append({**cell, "estimand": "DIAG", "method": "var_tau", "value": float(np.var(truth_tau))})

        # Naive pool: one outcome regression on all trial and OBS rows, ignoring the source (the ridge regression
        # of mu_R and mu_O, without their clipping)
        pooled = {key: np.concatenate([trial[key], obs[key]]) for key in ("x", "a", "y")}
        mu = fit_outcome_regression(pooled)
        contrast = lambda x: mu(x)[1] - mu(x)[0]
        theta_pool = contrast(trial["x"]).mean()
        rows.append({**cell, "estimand": "ATE", "method": "naive_pool", "value": theta_pool - design.theta})
        # Shrinkage and pretest: the cross-fitted trial-only estimate pooled with the naive pool
        theta_R, se_R = ate["rct_only"]
        for name, estimate in (("shrinkage", AIPWF.shrinkage(theta_R, se_R, theta_pool)),
                               ("pretest", AIPWF.pretest(theta_R, se_R, theta_pool))):
            rows.append({**cell, "estimand": "ATE", "method": name, "value": estimate - design.theta})
        rows.append({**cell, "estimand": "CATE", "learner": "pool", "method": "naive_pool",
                     "value": np.mean((contrast(truth_x) - truth_tau) ** 2)})
    return rows
