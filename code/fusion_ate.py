"""One ATE replication, shared by experiments ATE-1 and ATE-2.

Implements handoff Section 3.2 (Algorithm 1, the evaluation-sample variance,
the oracle coefficients) and attaches the baselines of Section 4.
"""
from __future__ import annotations

import numpy as np

import fusion_baselines as fb
from fusion_core import (CI_Z, ORACLE_DRAWS, aipw_score, algorithm1, ate_estimate, build_scores,
                         evaluation_variance, fit_outcome, fit_propensity, ghat_on,
                         project, split_roles, stable_seed, theorem1_variance,
                         trim_propensity, tuning_moments, var, rng_for)

RCT_FRACTIONS = (0.4, 0.3, 0.3)
OBS_FRACTIONS = (0.6, 0.2, 0.2)


def oracle_coefficients(sampler, mu_r, mu_o, ratio, n: int, n_obs: int, seed: int):
    """Handoff Section 3.2: fresh draws with this replication's nuisances fixed."""
    rct, obs = sampler(ORACLE_DRAWS, seed)
    s = build_scores(rct, mu_r, mu_o)
    g_o = ghat_on(obs["x"], mu_o)
    r_o = ratio(obs["x"])
    a = var(s.delta)
    b = var(s.ghat) + (n / n_obs) * var(r_o * g_o)
    c = float(np.cov(s.z0, s.delta, ddof=1)[0, 1])
    d = float(np.cov(s.z0, s.ghat, ddof=1)[0, 1])
    lam = project(-c / a) if a > 1e-12 else 0.0
    om = project(d / b) if b > 1e-12 else 0.0
    v_theory = theorem1_variance(s.z0, s.delta, s.ghat, r_o * g_o, lam, om, n, n_obs)
    return {"lambda": lam, "omega": om, "A": a, "B": b, "C": c, "D": d,
            "theorem1_variance": v_theory,
            "gain_vs_rct": float((c ** 2 / a if a > 1e-12 else 0.0)
                                 + (d ** 2 / b if b > 1e-12 else 0.0)) / n}


def ate_replication(rct, obs, theta, ratio_builder, ratio_name, spec, seed,
                    oracle_sampler=None, oracle_fn=None, with_baselines=True,
                    cross_fit_note: str = "honest"):
    # Trim the RCT only.  The trial propensity enters the AIPW score and the RF
    # weight, so extreme values must go; the OBS propensity enters no fusion
    # quantity, and trimming on it would distort the OBS covariate law and so
    # break the transport identity E_O(r0 g) = E_R(g).
    rct_t, rct_drop = trim_propensity(rct)
    obs_t, obs_drop = obs, 0.0
    r_roles = split_roles(rct_t, RCT_FRACTIONS, rng_for(seed, 1))
    o_roles = split_roles(obs_t, OBS_FRACTIONS, rng_for(seed, 2))

    mu_r = fit_outcome(r_roles["nuis"], spec="flexible")
    mu_o = fit_outcome(o_roles["nuis"], spec=spec)
    prop_o = fit_propensity(o_roles["nuis"], seed=seed)
    ratio = ratio_builder(r_roles["nuis"], o_roles["nuis"], seed=seed, oracle_fn=oracle_fn)

    s_tune = build_scores(r_roles["tune"], mu_r, mu_o)
    s_eval = build_scores(r_roles["eval"], mu_r, mu_o)
    g_o_tune = ghat_on(o_roles["tune"]["x"], mu_o)
    g_o_eval = ghat_on(o_roles["eval"]["x"], mu_o)
    r_tune = ratio(o_roles["tune"]["x"])
    r_eval = ratio(o_roles["eval"]["x"])
    n, n_obs = len(s_eval.z0), len(g_o_eval)

    moments = tuning_moments(s_tune, g_o_tune, r_tune, n, n_obs)
    coef = algorithm1(moments)
    rg_eval = r_eval * g_o_eval

    def pack(name, lam, om, extra=None):
        z_lambda = s_eval.z0 + lam * s_eval.delta
        est = ate_estimate(z_lambda, s_eval.ghat, rg_eval, om)
        v = evaluation_variance(z_lambda, s_eval.ghat, rg_eval, om)
        half = CI_Z * np.sqrt(max(v, 0.0))
        row = {"estimator": name, "estimate": est, "error": est - theta,
               "variance_hat": v, "ci_length": 2 * half,
               "covered": float(abs(est - theta) <= half),
               "lambda_used": lam, "omega_used": om,
               "failed": float(not (np.isfinite(est) and np.isfinite(v) and v >= 0))}
        if extra:
            row.update(extra)
        return row

    rows = [pack("rct_only", 0.0, 0.0),
            pack("lambda_only", coef["lambda"], 0.0),
            pack("omega_only", 0.0, coef["omega"]),
            pack("joint", coef["lambda"], coef["omega"],
                 {"lambda_raw": coef["lambda_raw"], "omega_raw": coef["omega_raw"],
                  "lambda_fallback": coef["lambda_fallback"], "omega_fallback": coef["omega_fallback"],
                  "lambda_boundary": coef["lambda_boundary"], "omega_boundary": coef["omega_boundary"]}),
            pack("fixed_half", 0.5, 0.5)]

    if oracle_sampler is not None:
        orc = oracle_coefficients(oracle_sampler, mu_r, mu_o, ratio, n, n_obs, stable_seed(seed, "oracle"))
        rows.append(pack("oracle_joint", orc["lambda"], orc["omega"],
                         {"oracle_lambda": orc["lambda"], "oracle_omega": orc["omega"],
                          "oracle_A": orc["A"], "oracle_B": orc["B"], "oracle_C": orc["C"],
                          "oracle_D": orc["D"], "theorem1_variance": orc["theorem1_variance"],
                          "oracle_gain": orc["gain_vs_rct"]}))
        # the omega-off partner at the same oracle lambda, so that the oracle
        # comparison is paired on lambda exactly as the estimated one is
        rows.append(pack("oracle_lambda_only", orc["lambda"], 0.0,
                         {"oracle_lambda": orc["lambda"], "oracle_omega": orc["omega"]}))
        rows.append(pack("oracle_omega_only", 0.0, orc["omega"],
                         {"oracle_lambda": orc["lambda"], "oracle_omega": orc["omega"]}))
        rows[3]["oracle_lambda"] = orc["lambda"]
        rows[3]["oracle_omega"] = orc["omega"]

    if with_baselines:
        gr_tune = s_tune.ghat * 0 + _rct_contrast(mu_r, r_roles["tune"]["x"])
        gr_eval = _rct_contrast(mu_r, r_roles["eval"]["x"])
        gr_obs_tune = _rct_contrast(mu_r, o_roles["tune"]["x"])
        gr_obs_eval = _rct_contrast(mu_r, o_roles["eval"]["x"])
        ss = fb.semi_supervised(s_tune.z0, gr_tune, gr_obs_tune, r_tune,
                                s_eval.z0, gr_eval, gr_obs_eval, r_eval, n, n_obs)
        rows.append(_baseline_row("semi_supervised", ss, theta, {"eta": ss["eta"]}))

        # Handoff Sections 4.2 to 4.4.  Both component estimates share the RCT
        # covariate sample, so the combination carries their covariance, and the
        # weight and the test decision are fixed on the tuning roles.
        g_r_tune = _rct_side_prediction(mu_o, r_roles["tune"]["x"])
        g_r_eval = _rct_side_prediction(mu_o, r_roles["eval"]["x"])
        u_tune = fb.obs_pseudo(o_roles["tune"], mu_o, prop_o)
        u_eval = fb.obs_pseudo(o_roles["eval"], mu_o, prop_o)
        m_tune = fb.combination_moments(s_tune.z0, g_r_tune, u_tune, r_tune, n, n_obs)
        m_eval = fb.combination_moments(s_eval.z0, g_r_eval, u_eval, r_eval, n, n_obs)
        rows.append(_baseline_row("obs_transported",
                                  {"estimate": m_eval["theta_O"], "variance": m_eval["V_O"]},
                                  theta))
        shrink = fb.shrinkage_weight(m_tune)
        adapt = fb.adaptive_weight(m_tune)
        for name, weights in (("shrinkage", shrink), ("adaptive", adapt)):
            result = fb.combine(m_eval["theta_R"], m_eval["theta_O"], weights["gamma"], m_eval)
            rows.append(_baseline_row(name, result, theta,
                                      {"gamma": weights["gamma"],
                                       "bias_detected": weights.get("bias_detected", np.nan),
                                       "covariance_RO": m_tune["C_RO"]}))

        # Handoff Section 4.4: the literal pooled average with denominator n_R + n_O
        pooled_model = fb.fit_pooled_outcome(r_roles["tune"], o_roles["tune"], spec="flexible")
        z_rct = _pooled_pseudo(r_roles["eval"], pooled_model, r_roles["eval"]["propensity"])
        e_obs = np.clip(prop_o.predict_proba(o_roles["eval"]["x"])[:, 1], 0.02, 0.98)
        z_obs = _pooled_pseudo(o_roles["eval"], pooled_model, e_obs)
        rows.append(_baseline_row("naive_pool",
                                  fb.naive_pool_spec(z_rct, r_eval * z_obs, n, n_obs), theta))

    diagnostics = {"ratio": ratio_name, "outcome_spec": spec, "split": cross_fit_note,
                   "n_rct_eval": n, "n_obs_eval": n_obs,
                   "rct_trim_fraction": rct_drop, "obs_trim_fraction": obs_drop,
                   "tuning_A": moments["A"], "tuning_B": moments["B"],
                   "tuning_C": moments["C"], "tuning_D": moments["D"]}
    for key, value in ratio.diagnostics.items():
        diagnostics[f"ratio_{key}"] = value
    for row in rows:
        row.update(diagnostics)
    return rows


def _rct_contrast(mu_r, x):
    m0, m1 = mu_r.predict(x)
    return m1 - m0


def _rct_side_prediction(mu_o, x):
    """The OBS prediction evaluated on RCT covariates, which is what makes the two
    component estimates of the combination rules share a sample."""
    m0, m1 = mu_o.predict(x)
    return m1 - m0


def _pooled_pseudo(data, model, propensity):
    m0, m1 = model.predict(data["x"])
    return aipw_score({"a": data["a"], "y": data["y"], "propensity": propensity}, m0, m1)


def _baseline_row(name, result, theta, extra=None):
    v = float(result.get("variance", np.nan))
    half = CI_Z * np.sqrt(max(v, 0.0)) if np.isfinite(v) else np.nan
    row = {"estimator": name, "estimate": result["estimate"],
           "error": result["estimate"] - theta, "variance_hat": v,
           "failed": float(not (np.isfinite(result["estimate"]) and np.isfinite(v))),
           "ci_length": 2 * half if np.isfinite(half) else np.nan,
           "covered": float(abs(result["estimate"] - theta) <= half) if np.isfinite(half) else np.nan,
           "lambda_used": np.nan, "omega_used": np.nan}
    if extra:
        row.update(extra)
    return row


def exact_variance_check(sampler, theta, spec, seed, n_rct, n_obs, repeats=200):
    """Acceptance test 3: fix the nuisances, then compare the Monte Carlo
    variance of a fixed-coefficient estimator with the Theorem 1 formula."""
    rct0, obs0 = sampler(max(n_rct, 400), seed)
    _, obs_big = sampler(max(n_obs, 2000), seed + 1)
    mu_r = fit_outcome(rct0, spec="flexible")
    mu_o = fit_outcome(obs_big, spec=spec)
    lam, om = 0.5, 0.5
    estimates = []
    for k in range(repeats):
        rct, obs = sampler(n_rct, stable_seed(seed, "exact", k))
        s = build_scores(rct, mu_r, mu_o)
        g_o = ghat_on(obs["x"], mu_o)
        estimates.append(ate_estimate(s.z0 + lam * s.delta, s.ghat, g_o, om))
    rct_l, obs_l = sampler(100000, stable_seed(seed, "exact", "pop"))
    s_l = build_scores(rct_l, mu_r, mu_o)
    g_l = ghat_on(obs_l["x"], mu_o)
    formula = theorem1_variance(s_l.z0, s_l.delta, s_l.ghat, g_l, lam, om, n_rct, n_obs)
    mc = float(np.var(estimates, ddof=1))
    se = mc * np.sqrt(2.0 / (repeats - 1))
    return {"monte_carlo_variance": mc, "theorem1_variance": formula,
            "standard_error": se, "z_gap": (mc - formula) / se if se > 0 else np.nan,
            "repeats": repeats}
