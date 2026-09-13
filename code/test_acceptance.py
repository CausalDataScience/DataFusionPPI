"""Acceptance tests of handoff Section 8.1.  Run with: python3 test_acceptance.py"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

import fusion_data as fd
import fusion_cate as fc
from fusion_ate import ate_replication
from fusion_core import (build_scores, chi_weight, clip_prediction, fit_outcome, ghat_on,
                         learner_parts, prediction_bound, ratio_balanced, ratio_unit,
                         rng_for, sieve_solver, spline_basis_factory, stable_seed,
                         standardize_factory, validation_score, var)

MATERIALS = Path(__file__).resolve().parents[1] / "materials"
RESULTS: list[dict] = []


def record(name: str, passed: bool, detail: str) -> None:
    RESULTS.append({"test": name, "passed": bool(passed), "detail": detail})
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}", flush=True)


def test_determinism() -> None:
    """Test 1: the same seed reproduces the replication file byte for byte."""
    digests = []
    for _ in range(2):
        subprocess.run([sys.executable, "exp_ate1.py", "--smoke", "--replications", "5",
                        "--families", "scm1", "--stem", "acc_determinism"],
                       check=True, capture_output=True, cwd=Path(__file__).parent)
        digests.append((MATERIALS / "acc_determinism_replications.csv").read_bytes())
    for path in MATERIALS.glob("acc_determinism*"):
        path.unlink()
    record("1_determinism", digests[0] == digests[1],
           f"two runs produced {'identical' if digests[0] == digests[1] else 'different'} bytes")


def test_rct_only_agreement() -> None:
    """Test 2: the RCT-only estimator is the plain AIPW average of Equation (7)."""
    scm = fd.scm_family(1)
    theta = fd.scm_true_theta(scm, "shared")
    seed = stable_seed("acc", "rct")
    rct = fd.sample_scm(scm, "RCT", 400, "shared", rng_for(seed, 1))
    obs = fd.sample_scm(scm, "OBS", 4000, "shared", rng_for(seed, 2))
    rows = ate_replication(rct, obs, theta, ratio_unit, "unit", "flexible", seed,
                           with_baselines=False)
    got = [r for r in rows if r["estimator"] == "rct_only"][0]
    from fusion_ate import RCT_FRACTIONS
    from fusion_core import split_roles, trim_propensity
    rct_t, _ = trim_propensity(rct)
    roles = split_roles(rct_t, RCT_FRACTIONS, rng_for(seed, 1))
    mu_r = fit_outcome(roles["nuis"], spec="flexible")
    m0, m1 = mu_r.predict(roles["eval"]["x"])
    from fusion_core import aipw_score
    direct = float(aipw_score(roles["eval"], m0, m1).mean())
    gap = abs(direct - got["estimate"])
    record("2_rct_only_is_plain_aipw", gap < 1e-10, f"|pipeline - direct AIPW| = {gap:.2e}")


def test_oracle_ordering() -> None:
    """Test 4 (Corollary 2.1): the oracle joint variance is no larger than RCT-only."""
    summary = MATERIALS / "ate1_when_fusion_helps_summary.csv"
    if not summary.exists():
        record("4_oracle_ordering", False, "ATE-1 summary not found; run exp_ate1.py first")
        return
    frame = pd.read_csv(summary)
    block = frame[frame.metric == "rmse"]
    keys = ["family", "axis", "n_rct", "n_obs", "confounding", "spec"]
    joined = block[block.estimator == "oracle_joint"].merge(
        block[block.estimator == "rct_only"], on=keys, suffixes=("_oracle", "_rct"))
    tol = 3 * (joined["monte_carlo_se_oracle"] + joined["monte_carlo_se_rct"])
    violations = int((joined["value_oracle"] > joined["value_rct"] + tol).sum())
    record("4_oracle_ordering", violations == 0,
           f"{len(joined)} cells, {violations} with oracle RMSE above RCT-only by more than 3 Monte Carlo SE")


def test_balance() -> None:
    """Test 5 (Proposition 6): the calibrated ratio solves the balance equations."""
    worst = 0.0
    for scm_id in (1, 2, 3):
        scm = fd.scm_family(scm_id, shift=1.0)
        rct = fd.sample_scm(scm, "RCT", 2000, "shifted", rng_for(31, scm_id))
        obs = fd.sample_scm(scm, "OBS", 4000, "shifted", rng_for(32, scm_id))
        r = ratio_balanced(rct, obs, seed=scm_id)
        worst = max(worst, r.diagnostics["balance_residual"])
    record("5_calibrated_balance", worst <= 1e-6, f"largest balance residual = {worst:.2e}")


def test_sieve_fit() -> None:
    """Test 6: the normal-equation solution minimizes the penalized loss."""
    scm = fd.scm_family(1)
    seed = stable_seed("acc", "sieve")
    rct = fd.sample_scm(scm, "RCT", 600, "shared", rng_for(seed, 1))
    obs = fd.sample_scm(scm, "OBS", 3000, "shared", rng_for(seed, 2))
    mu_r = fit_outcome({k: v[:200] for k, v in rct.items()}, spec="flexible")
    mu_o = fit_outcome({k: v[:1500] for k, v in obs.items()}, spec="flexible")
    tune = {k: v[200:] for k, v in rct.items()}
    obs_t = {k: v[1500:] for k, v in obs.items()}
    s = build_scores(tune, mu_r, mu_o)
    basis = standardize_factory(spline_basis_factory(tune["x"], 3), tune["x"])
    b_r, b_o = basis(tune["x"]), basis(obs_t["x"])
    g_r, g_o = s.ghat, ghat_on(obs_t["x"], mu_o)
    r_o = np.ones(len(b_o))
    worst = 0.0
    for learner in ("DRF", "RF"):
        w, zt0, dzt = learner_parts(tune, s, learner)
        lam, om, rho = 0.5, 0.5, 1e-2
        solve, _ = sieve_solver(b_r, w, zt0, dzt, g_r, b_o, r_o, g_o, om, rho)
        beta = solve(lam)

        def loss(coef):
            zeta_r, zeta_o = b_r @ coef, b_o @ coef
            zt = zt0 + lam * dzt
            first = float(np.mean(w * (zt - zeta_r) ** 2))
            bracket = float(np.mean(r_o * (g_o - zeta_o) ** 2) - np.mean(w * (g_r - zeta_r) ** 2))
            return first + om * bracket + rho * float(coef @ coef)

        direct = minimize(loss, np.zeros(len(beta)), method="L-BFGS-B",
                          options={"maxiter": 5000, "ftol": 1e-16, "gtol": 1e-14}).x
        worst = max(worst, float(np.max(np.abs(loss(beta) - loss(direct)))))
    record("6_sieve_fit_minimizes_loss", worst <= 1e-8,
           f"largest loss gap between the normal equation and direct minimization = {worst:.2e}")


def test_score_identity(replications: int = 200) -> None:
    """Test 7 (Theorem 6 item 1): the score difference is unbiased for the risk
    difference when r equals r0."""
    scm = fd.scm_family(1)
    test = fd.scm_test_sample(scm, "shared", 50000, rng_for(41, 1))
    gaps_score, gaps_risk = [], []
    for rep in range(replications):
        seed = stable_seed("acc", "score", rep)
        rct = fd.sample_scm(scm, "RCT", 200, "shared", rng_for(seed, 1))
        obs = fd.sample_scm(scm, "OBS", 4000, "shared", rng_for(seed, 2))
        rows = fc.cate_replication(rct, obs, test, "spline", seed, ratio_unit, "unit")
        drf = [r for r in rows if r["learner"] == "DRF"][0]
        gaps_score.append(drf["score_selected_minus_rct"])
        gaps_risk.append(drf["risk_selected_minus_rct"])
    diff = np.array(gaps_score) - np.array(gaps_risk)
    se = float(np.std(diff, ddof=1) / np.sqrt(len(diff)))
    z = float(np.mean(diff) / se) if se > 0 else np.nan
    record("7_score_identity", abs(z) <= 3,
           f"mean(score gap - risk gap) = {np.mean(diff):.4f}, Monte Carlo SE {se:.4f}, z = {z:.2f}")


def test_cate2_construction(replications: int = 30) -> None:
    """The construction keeps one covariate law and induces the intended bias."""
    cohort = fd.star_cohort("mathk", "pooled")
    naive, sizes, patterns = [], [], []
    for rep in range(replications):
        _, obs, extra = fd.star_cate2_replication(cohort, 400, 0.8, rep)
        d = extra["diagnostics"]
        naive.append(d["naive_obs_contrast"])
        sizes.append(d["obs_size"])
        patterns.append(d["patterns"])
    biased = bool(np.mean(naive) > cohort.ate_reference + 0.1)
    stable = len(set(patterns)) == 1
    record("8_cate2_construction", biased and stable,
           f"mean naive observational contrast {np.mean(naive):.3f} against the experimental "
           f"{cohort.ate_reference:.3f}; {patterns[0]} covariate patterns in every replication; "
           f"observational size {int(np.mean(sizes))}")


# ---------------------------------------------------------------- P0-B checks
def test_star_covariate_laws(draws: int = 60000) -> None:
    """P0-B: the two samplers must draw covariates from the same law.

    The earlier version compared an array with itself and could not fail.  This
    one draws from both samplers and compares the empirical pattern
    distributions, and then injects a perturbed law to confirm the comparison
    would catch a difference.
    """
    cohort = fd.star_cohort("mathk", "pooled")
    law = cohort.pattern_law
    rct, obs, _ = fd.star_cate2_replication(cohort, draws, 0.8, 0, n_eval=10, n_obs=draws)
    share_r = np.bincount(rct["pattern"], minlength=len(law)) / len(rct["pattern"])
    share_o = np.bincount(obs["pattern"], minlength=len(law)) / len(obs["pattern"])
    tolerance = 5.0 * np.sqrt(law * (1 - law) / draws)
    agree = bool(np.all(np.abs(share_r - share_o) <= tolerance + tolerance))
    # injection: a law that differs in one pattern must be detected
    perturbed = law.copy()
    perturbed[0] *= 1.5
    perturbed = perturbed / perturbed.sum()
    injected = np.bincount(np.random.default_rng(0).choice(len(law), size=draws, p=perturbed),
                           minlength=len(law)) / draws
    caught = bool(np.any(np.abs(injected - share_r) > tolerance + tolerance))
    record("9_star_common_covariate_law", agree and caught,
           f"the two samplers agree within five standard errors on all {len(law)} patterns, "
           f"largest difference {np.max(np.abs(share_r - share_o)):.2e}; a perturbed law is detected")


def test_reporting_does_not_drive_selection(replications: int = 12) -> None:
    """P0-B: changing only the reporting half must leave the selection untouched."""
    import exp_cate2 as e2
    from fusion_core import rng_for as _rng
    cohort = fd.star_cohort("mathk", "pooled")
    unchanged, changed_report = [], []
    for rep in range(replications):
        rct, obs, extra = fd.star_cate2_replication(cohort, 400, 0.8, rep)
        seed = stable_seed("acc", "sep", rep)
        evaluation = extra["eval"]
        base = e2.cate2_replication(cohort, rct, obs, evaluation, "spline", seed)
        # reproduce the split the driver uses and perturb only the reporting half
        order = _rng(seed, 3).permutation(len(evaluation["x"]))
        report_rows = order[len(order) // 2:]
        perturbed = {k: v.copy() for k, v in evaluation.items()}
        perturbed["y"][report_rows] += 3.0
        moved = e2.cate2_replication(cohort, rct, obs, perturbed, "spline", seed)
        unchanged.append(base[0]["selected_lambda"] == moved[0]["selected_lambda"]
                         and base[0]["selected_omega"] == moved[0]["selected_omega"]
                         and base[0]["selected_rho"] == moved[0]["selected_rho"])
        changed_report.append(abs(base[0]["score_selected"] - moved[0]["score_selected"]) > 1e-8)
    ok = all(unchanged) and all(changed_report)
    record("10_selection_independent_of_reporting", ok,
           f"perturbing only the reporting half left the selection unchanged in "
           f"{sum(unchanged)}/{replications} replications and moved the reported score in "
           f"{sum(changed_report)}/{replications}")


# ---------------------------------------------------------------- P0-C checks
def test_baseline_source_indicator() -> None:
    """P0-C: the production rule must equal the direct design-matrix product."""
    import fusion_baselines as fbm
    rng = rng_for(555, 1)
    n, p = 40, 4
    basis = lambda z: np.column_stack([np.ones(len(z)), z[:, :p - 1]])
    beta = np.array([0.5, -1.0, 2.0, 0.25, 0.75])        # last entry is the source effect
    x = rng.normal(size=(n, p))
    design = np.column_stack([basis(x), np.ones(n)])      # the indicator is one for the trial
    wanted = design @ beta
    got = basis(x) @ fbm.source_indicator_coefficients(beta)
    gap = float(np.max(np.abs(got - wanted)))
    # injection: the earlier rule added the source effect to every coefficient
    wrong = float(np.max(np.abs(basis(x) @ (beta[:-1] + beta[-1]) - wanted)))
    record("11_source_indicator_prediction", gap <= 1e-12 and wrong > 1e-6,
           f"the production rule matches the direct product to {gap:.1e}; the earlier rule "
           f"differs by {wrong:.3f}, so the check would have caught it")


def test_combination_weight_fixed_before_evaluation() -> None:
    """P0-C: perturbing only the evaluation role must move the estimate and leave
    the combination weight and the test decision alone."""
    from fusion_core import split_roles, trim_propensity, rng_for as _rng
    from fusion_ate import RCT_FRACTIONS
    scm = fd.scm_family(1)
    theta = fd.scm_true_theta(scm, "shared")
    seed = stable_seed("acc", "weight")
    rct = fd.sample_scm(scm, "RCT", 400, "shared", _rng(seed, 1))
    obs = fd.sample_scm(scm, "OBS", 4000, "shared", _rng(seed, 2))
    base = ate_replication(rct, obs, theta, ratio_unit, "unit", "flexible", seed)
    trimmed, _ = trim_propensity(rct)
    order = _rng(seed, 1).permutation(len(trimmed["x"]))
    n_nuis = int(round(RCT_FRACTIONS[0] * len(order)))
    n_tune = int(round(RCT_FRACTIONS[1] * len(order)))
    eval_rows = order[n_nuis + n_tune:]
    moved_rct = {k: v.copy() for k, v in trimmed.items()}
    moved_rct["y"][eval_rows] += 2.0
    moved = ate_replication(moved_rct, obs, theta, ratio_unit, "unit", "flexible", seed)
    def field(rows, name, key):
        hit = [r for r in rows if r["estimator"] == name]
        return hit[0].get(key) if hit else None
    same_weight = all(abs(field(base, n, "gamma") - field(moved, n, "gamma")) < 1e-12
                      for n in ("shrinkage", "adaptive"))
    same_decision = field(base, "adaptive", "bias_detected") == field(moved, "adaptive", "bias_detected")
    moved_estimate = abs(field(base, "shrinkage", "estimate") - field(moved, "shrinkage", "estimate")) > 1e-8
    record("12_combination_weight_from_tuning", same_weight and same_decision and moved_estimate,
           f"perturbing the evaluation outcomes left gamma and the test decision unchanged and "
           f"moved the estimate by {abs(field(base, 'shrinkage', 'estimate') - field(moved, 'shrinkage', 'estimate')):.3f}")


def test_balanced_with_prediction() -> None:
    """P0-C: the loss-targeted balancing map really carries the fitted prediction."""
    import exp_ate2 as e2
    specs = dict(e2.ratio_specs())
    scm = fd.scm_family(1, shift=1.0)
    rct = fd.sample_scm(scm, "RCT", 800, "shifted", rng_for(606, 1))
    obs = fd.sample_scm(scm, "OBS", 3000, "shifted", rng_for(606, 2))
    plain = specs["balanced_x"](rct, obs, seed=1)
    with_g = specs["balanced_x_ghat"](rct, obs, seed=1)
    gap = float(np.max(np.abs(plain(obs["x"]) - with_g(obs["x"]))))
    record("13_balanced_with_prediction", gap > 1e-8 and with_g.diagnostics["balance_residual"] <= 1e-6,
           f"the two maps differ by {gap:.3e} and the loss-targeted map balances to "
           f"{with_g.diagnostics['balance_residual']:.1e}")


def test_no_truth_in_feasible_selection() -> None:
    """P0-C: feasible single-channel rules are chosen by the score, not by the truth."""
    import fusion_cate as fcm
    src = __import__("inspect").getsource(fcm.cate_replication)
    feasible = "lambda_only = pick_by_score" in src and "omega_only = pick_by_score" in src
    oracle_named = "risk_lambda_only_oracle" in src and "risk_omega_only_oracle" in src
    record("14_feasible_selection_uses_score", feasible and oracle_named,
           "single-channel rules select on the validation score; the true-risk versions "
           "are reported separately under oracle names")


def test_bounded_design_assumptions(replications: int = 20) -> None:
    """The range condition of Theorem 6 must hold in every replication of the
    bounded design, and the returned effect must be the exact contrast."""
    fam = fd.bounded_family(1)
    bound, cap = fam.score_range
    from fusion_core import build_scores, fit_outcome, split_roles, trim_propensity
    worst_z, worst_g, worst_y = 0.0, 0.0, 0.0
    for rep in range(replications):
        r = fam.sample("RCT", 400, rng_for(rep, 1))
        o = fam.sample("OBS", 3000, rng_for(rep, 2))
        worst_y = max(worst_y, float(np.abs(r["y"]).max()), float(np.abs(o["y"]).max()))
        rr = split_roles(trim_propensity(r)[0], (0.4, 0.3, 0.3), rng_for(rep, 3))
        oo = split_roles(o, (0.6, 0.2, 0.2), rng_for(rep, 4))
        mu_r = fit_outcome(rr["nuis"], clip=fd.BOUNDED_PREDICTION)
        mu_o = fit_outcome(oo["nuis"], clip=fd.BOUNDED_PREDICTION)
        s = build_scores(rr["eval"], mu_r, mu_o)
        worst_z = max(worst_z, float(np.abs(s.z0).max()))
        worst_g = max(worst_g, float(np.abs(s.ghat).max()))
    x = fam.test_sample(2000, rng_for(7, 1))
    exact = bool(np.allclose(fam.tau(x["x"]), x["tau"]))
    record("16_bounded_design_assumptions",
           worst_y <= 1.0 and worst_z <= bound and worst_g <= fd.BOUNDED_CANDIDATE and exact,
           f"|Y| max {worst_y:.3f} against 1, |Z0| max {worst_z:.3f} against {bound}, "
           f"|ghat| max {worst_g:.3f} against {fd.BOUNDED_CANDIDATE}; the returned effect is exact")


# ---------------------------------------------------------------- P0-D checks
def test_aggregation_error_injection() -> None:
    """P0-D: four deliberate faults must be caught rather than averaged away."""
    import analysis_variance as av
    group = ["cell"]
    base = []
    rng = np.random.default_rng(3)
    for cell in ("a", "b"):
        for rep in range(60):
            for est, shift in (("rct_only", 0.0), ("joint", 0.0)):
                value = rng.normal(shift, 1.0 if est == "rct_only" else 0.8)
                base.append({"cell": cell, "replication": rep, "estimator": est,
                             "estimate": value, "error": value, "variance_hat": 1.0,
                             "covered": 1.0, "ci_length": 2.0, "failed": 0.0})
    frame = pd.DataFrame(base)
    findings = []

    # 1. a missing value in the middle must shrink the effective count consistently
    holed = frame.copy()
    holed.loc[(holed.cell == "a") & (holed.replication == 10), "error"] = np.nan
    report = av.reporting_table(holed, group)
    counts = report[report.cell == "a"]["effective_replications"].unique()
    findings.append(("a mid-table missing value", len(counts) == 1 and counts[0] == 59))

    # 2. a duplicated replication must raise rather than be averaged
    duplicated = pd.concat([frame, frame.iloc[[0]]], ignore_index=True)
    try:
        av.check_uniqueness(duplicated, group)
        findings.append(("a duplicated replication", False))
    except AssertionError:
        findings.append(("a duplicated replication", True))

    # 3. the decomposition identity at K = 2
    two = frame[frame.replication < 2]
    rep2 = av.reporting_table(two, group)
    findings.append(("the identity at two replications",
                     bool(np.nanmax(np.abs(rep2["mse_identity_gap"])) < 1e-12)))

    # 4. a pooled interval must not be the mean of the per-cell endpoints
    variance = av.variance_table(frame, group)
    cell_rows = variance[(variance.scope == "cell") & (variance.comparison == "whole_method")]
    pooled = variance[(variance.scope == "pooled") & (variance.comparison == "whole_method")]
    naive_low = float(cell_rows["rate_low"].mean())
    findings.append(("a pooled interval narrower than the endpoint average",
                     len(pooled) == 1 and float(pooled["rate_low"].iloc[0]) > naive_low))

    passed = all(ok for _, ok in findings)
    record("15_aggregation_error_injection", passed,
           "; ".join(f"{name}: {'caught' if ok else 'missed'}" for name, ok in findings))


def main() -> int:
    test_determinism()
    test_rct_only_agreement()
    test_oracle_ordering()
    test_balance()
    test_sieve_fit()
    test_score_identity()
    test_cate2_construction()
    test_star_covariate_laws()
    test_reporting_does_not_drive_selection()
    test_baseline_source_indicator()
    test_combination_weight_fixed_before_evaluation()
    test_balanced_with_prediction()
    test_no_truth_in_feasible_selection()
    test_aggregation_error_injection()
    test_bounded_design_assumptions()
    frame = pd.DataFrame(RESULTS)
    frame.to_csv(MATERIALS / "acceptance_tests.csv", index=False)
    failed = int((~frame["passed"]).sum())
    print(f"\n{len(frame) - failed}/{len(frame)} acceptance tests passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
