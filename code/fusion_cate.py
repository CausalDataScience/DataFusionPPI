"""One CATE replication: Algorithm 2 on a linear sieve (handoff Section 3.4).

Checks Theorem 6 (honest validation), Theorem 7 (risk of a sieve candidate),
and Corollary 7.1 (oracle sieve coefficients).
"""
from __future__ import annotations

import numpy as np

import fusion_baselines as fb
from fusion_core import (DEFAULT_GRID, RIDGE_AXIS, build_scores, cate_risk, clip_prediction,
                         eval_radius, prediction_bound,
                         fit_outcome, fit_propensity, ghat_on, learner_parts, project,
                         rng_for, sieve_solver, split_roles, stable_seed, trim_propensity,
                         validation_score, var)

RCT_FRACTIONS = (0.4, 0.3, 0.3)
OBS_FRACTIONS = (0.6, 0.2, 0.2)
LEARNERS = ("DRF", "RF")


def _candidates(grid, ridges):
    return [(lam, om, rho) for rho in ridges for (lam, om) in grid]


def cate_replication(rct, obs, test, sieve: str, seed: int, ratio_builder, ratio_name: str,
                     grid=DEFAULT_GRID, oracle_sampler=None, spline_knots: int = 3,
                     spec: str = "flexible", oracle_fn=None, anchor: np.ndarray | None = None,
                     predictions: dict | None = None, score_range: tuple | None = None,
                     candidate_bound: float | None = None, outcome_clip: float | None = None):
    """When ``anchor`` is given, the fitted function of each reported rule is
    evaluated on those fixed covariates and stored in ``predictions``, so that a
    later pass can split the risk into prediction variance and squared bias."""
    rct_t, rct_drop = trim_propensity(rct)
    r_roles = split_roles(rct_t, RCT_FRACTIONS, rng_for(seed, 1))
    o_roles = split_roles(obs, OBS_FRACTIONS, rng_for(seed, 2))

    mu_r = fit_outcome(r_roles["nuis"], spec="flexible", clip=outcome_clip)
    mu_o = fit_outcome(o_roles["nuis"], spec=spec, clip=outcome_clip)
    prop_o = fit_propensity(o_roles["nuis"], seed=seed)
    ratio = ratio_builder(r_roles["nuis"], o_roles["nuis"], seed=seed, oracle_fn=oracle_fn)

    ridges = RIDGE_AXIS if sieve == "neural" else (1e-2,)
    if sieve == "neural":
        from fusion_core import neural_basis_factory
        basis_raw, ghat_fn = neural_basis_factory(o_roles["nuis"], seed=seed)
    else:
        from fusion_core import spline_basis_factory
        cat = np.array([len(np.unique(r_roles["nuis"]["x"][:, j])) <= 3
                        for j in range(r_roles["nuis"]["x"].shape[1])])
        basis_raw = spline_basis_factory(r_roles["nuis"]["x"], spline_knots, categorical=cat)
        ghat_fn = None
    from fusion_core import standardize_factory
    basis = standardize_factory(basis_raw, r_roles["nuis"]["x"])

    s_tune = build_scores(r_roles["tune"], mu_r, mu_o)
    s_eval = build_scores(r_roles["eval"], mu_r, mu_o)
    g_tune_r = s_tune.ghat if ghat_fn is None else ghat_fn(r_roles["tune"]["x"])
    g_eval_r = s_eval.ghat if ghat_fn is None else ghat_fn(r_roles["eval"]["x"])
    g_tune_o = ghat_on(o_roles["tune"]["x"], mu_o) if ghat_fn is None else ghat_fn(o_roles["tune"]["x"])
    g_eval_o = ghat_on(o_roles["eval"]["x"], mu_o) if ghat_fn is None else ghat_fn(o_roles["eval"]["x"])

    b_tune_r, b_eval_r = basis(r_roles["tune"]["x"]), basis(r_roles["eval"]["x"])
    b_tune_o, b_eval_o = basis(o_roles["tune"]["x"]), basis(o_roles["eval"]["x"])
    b_test = basis(test["x"])
    r_tune_o, r_eval_o = ratio(o_roles["tune"]["x"]), ratio(o_roles["eval"]["x"])
    n, n_obs = len(b_eval_r), len(b_eval_o)

    cand = _candidates(grid, ridges)
    m_candidates = len(cand)
    # a prespecified bound where the design supplies one, otherwise read off the
    # tuning sample, which is independent of the evaluation sample
    bound = candidate_bound if candidate_bound is not None else prediction_bound(s_tune.z0)
    g_eval_r = clip_prediction(g_eval_r, bound)
    g_eval_o = clip_prediction(g_eval_o, bound)
    g_tune_r = clip_prediction(g_tune_r, bound)
    g_tune_o = clip_prediction(g_tune_o, bound)
    rows = []
    for learner in LEARNERS:
        w, zt0, dzt = learner_parts(r_roles["tune"], s_tune, learner)
        betas, risks, scores = {}, {}, {}
        for rho in ridges:
            for om in sorted({c[1] for c in cand}):
                solve, _ = sieve_solver(b_tune_r, w, zt0, dzt, g_tune_r,
                                        b_tune_o, r_tune_o, g_tune_o, om, rho)
                for lam in sorted({c[0] for c in cand}):
                    beta = solve(lam)
                    betas[(lam, om, rho)] = beta
                    risks[(lam, om, rho)] = cate_risk(
                        clip_prediction(b_test @ beta, bound), test["tau"])
                    scores[(lam, om, rho)] = validation_score(
                        s_eval.z0,
                        clip_prediction(b_eval_r @ beta, bound), g_eval_r,
                        clip_prediction(b_eval_o @ beta, bound), g_eval_o, r_eval_o, om)

        def pick_by_score(subset):
            return min(subset, key=lambda c: (scores[c], c))

        def pick_by_risk(subset):
            return min(subset, key=lambda c: (risks[c], c))

        selected = pick_by_score(cand)
        grid_oracle = pick_by_risk(cand)
        # the trial-only rule may choose its ridge on the same score, as the joint rule does
        rct_only = pick_by_score([c for c in cand if c[0] == 0.0 and c[1] == 0.0])
        lambda_only = pick_by_score([c for c in cand if c[1] == 0.0])
        omega_only = pick_by_score([c for c in cand if c[0] == 0.0])
        # The range condition of Theorem 6 must hold by design, not by inspection.
        # score_range is supplied by the caller only when the data-generating
        # process bounds the pseudo-outcome and the ratio a priori.
        if score_range is None:
            bound_b, bar_r, radius = np.nan, np.nan, np.nan
        else:
            bound_b, bar_r = score_range
            radius = eval_radius(bound_b, bar_r, n, n_obs, m_candidates)
        regret = risks[selected] - risks[grid_oracle]

        row = {"learner": learner, "sieve": sieve, "ratio": ratio_name,
               "m_candidates": m_candidates,
               "risk_selected": risks[selected], "risk_grid_oracle": risks[grid_oracle],
               "risk_rct_only": risks[rct_only],
               "risk_lambda_only": risks[lambda_only],
               "risk_omega_only": risks[omega_only],
               "risk_lambda_only_oracle": min(risks[c] for c in cand if c[1] == 0.0),
               "risk_omega_only_oracle": min(risks[c] for c in cand if c[0] == 0.0),
               "rct_only_rho": rct_only[2],
               "regret": regret, "radius": radius,
               "regret_within_radius": float(regret <= 2 * radius) if np.isfinite(radius) else np.nan,
               "regret_over_radius": regret / (2 * radius) if np.isfinite(radius) and radius > 0 else np.nan,
               "radius_applicable": float(score_range is not None),
               "selected_lambda": selected[0], "selected_omega": selected[1],
               "selected_rho": selected[2],
               "oracle_grid_lambda": grid_oracle[0], "oracle_grid_omega": grid_oracle[1],
               "bound_b": bound_b, "bar_r": bar_r, "prediction_bound": bound,
               "n_rct_eval": n, "n_obs_eval": n_obs, "rct_trim_fraction": rct_drop,
               "score_selected_minus_rct": scores[selected] - scores[rct_only],
               "risk_selected_minus_rct": risks[selected] - risks[rct_only]}

        for key in ((0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5)):
            c = (key[0], key[1], ridges[0])
            row[f"risk_at_{key[0]}_{key[1]}"] = risks[c]

        if learner == "DRF" and oracle_sampler is not None:
            row.update(_oracle_calculus(oracle_sampler, mu_r, mu_o, ratio, basis, b_test,
                                        test, n_tune_r=len(b_tune_r), n_tune_o=len(b_tune_o),
                                        seed=stable_seed(seed, "cal")))
        if anchor is not None and predictions is not None:
            b_anchor = basis(anchor)
            for tag, c in (("selected", selected), ("rct_only", rct_only),
                           ("grid_oracle", grid_oracle)):
                predictions.setdefault((learner, tag), []).append(
                    clip_prediction(b_anchor @ betas[c], bound))
        rows.append(row)

    beta_eg = fb.experimental_grounding(basis, r_roles["tune"], s_tune, o_roles["tune"], mu_o, prop_o)
    _, beta_dp = fb.domain_indicator_pool(basis, r_roles["tune"], s_tune, o_roles["tune"], mu_o, prop_o)
    beta_dp_rct = fb.source_indicator_coefficients(beta_dp)
    for name, beta in (("experimental_grounding", beta_eg),
                       ("domain_indicator_pool", beta_dp_rct)):
        rows.append({"learner": name, "sieve": sieve, "ratio": ratio_name,
                     "m_candidates": m_candidates,
                     "risk_selected": cate_risk(clip_prediction(b_test @ beta, bound), test["tau"]),
                     "risk_rct_only": rows[0]["risk_rct_only"],
                     "n_rct_eval": n, "n_obs_eval": n_obs})
    return rows


def _oracle_calculus(sampler, mu_r, mu_o, ratio, basis, b_test, test, n_tune_r, n_tune_o, seed):
    """Corollary 7.1 for DRF: A_p, B_p, C_p, D_p and the oracle coefficients."""
    gamma = b_test.T @ b_test / len(b_test)
    beta_p = np.linalg.solve(gamma + 1e-10 * np.eye(len(gamma)),
                             (b_test * test["tau"][:, None]).mean(0))
    zeta_p_test = b_test @ beta_p
    a_p2 = cate_risk(zeta_p_test, test["tau"])

    rct, obs = sampler(20000, seed)
    s = build_scores(rct, mu_r, mu_o)
    b_r, b_o = basis(rct["x"]), basis(obs["x"])
    g_o = ghat_on(obs["x"], mu_o)
    zeta_r, zeta_o = b_r @ beta_p, b_o @ beta_p
    r_o = ratio(obs["x"])
    gi = np.linalg.inv(gamma + 1e-10 * np.eye(len(gamma)))

    def tr_var(mat):
        return float(np.trace(gi @ np.cov(mat.T, ddof=1)))

    def tr_cov(m1, m2):
        c = np.cov(np.column_stack([m1, m2]).T, ddof=1)
        p = m1.shape[1]
        return float(np.trace(gi @ c[:p, p:]))

    u_base = b_r * (s.z0 - zeta_r)[:, None]
    u_delta = b_r * s.delta[:, None]
    u_gap_r = b_r * (s.ghat - zeta_r)[:, None]
    u_gap_o = b_o * (r_o * (g_o - zeta_o))[:, None]
    a_p = tr_var(u_delta)
    b_p = tr_var(u_gap_r) + (n_tune_r / n_tune_o) * tr_var(u_gap_o)
    c_p = tr_cov(u_base, u_delta)
    d_p = tr_cov(u_base, u_gap_r)
    return {"a_p2": a_p2, "A_p": a_p, "B_p": b_p, "C_p": c_p, "D_p": d_p,
            "oracle_sieve_lambda": project(-c_p / a_p) if a_p > 1e-12 else 0.0,
            "oracle_sieve_omega": project(d_p / b_p) if b_p > 1e-12 else 0.0,
            "J_gain_fraction": float((c_p ** 2 / a_p if a_p > 1e-12 else 0.0)
                                     + (d_p ** 2 / b_p if b_p > 1e-12 else 0.0))
            / max(tr_var(u_base), 1e-12)}
