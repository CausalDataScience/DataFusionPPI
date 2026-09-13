"""CATE-2: STAR with real treatment, real outcomes, and a constructed
confounded observational sample (handoff 5.4).

The conditional effect is unknown here, so the reported CATE metric is the
validation-score difference of Theorem 6 item 1, which estimates the risk
difference without bias because the construction makes r0 identically one.
Selection and reporting use disjoint halves of the evaluation pool, so the
reported difference carries no selection bias.
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd

import fusion_baselines as fb
import fusion_data as fd
from fusion_ate import ate_replication
from fusion_core import (DEFAULT_GRID, RIDGE_AXIS, build_scores, clip_prediction,
                         prediction_bound, fit_outcome, fit_propensity,
                         ghat_on, learner_parts, ratio_unit, rng_for, sieve_solver,
                         split_roles, stable_seed, standardize_factory, trim_propensity,
                         validation_score)
from fusion_io import summarize, write_outputs

OBS_FRACTIONS = (0.6, 0.2, 0.2)
RCT_FRACTIONS = (0.57, 0.43, 0.0)     # the evaluation pool plays the evaluation role


def cate2_replication(cohort, rct, obs, evaluation, sieve: str, seed: int,
                      grid=DEFAULT_GRID, spline_knots: int = 3, spec: str = "linear",
                      anchor=None, predictions: dict | None = None):
    rct_t, _ = trim_propensity(rct)
    r_roles = split_roles(rct_t, RCT_FRACTIONS, rng_for(seed, 1))
    o_roles = split_roles(obs, OBS_FRACTIONS, rng_for(seed, 2))

    half_r = rng_for(seed, 3).permutation(len(evaluation["x"]))
    sel_r = {k: v[half_r[:len(half_r) // 2]] for k, v in evaluation.items()}
    rep_r = {k: v[half_r[len(half_r) // 2:]] for k, v in evaluation.items()}
    half_o = rng_for(seed, 4).permutation(len(o_roles["eval"]["x"]))
    sel_o = {k: v[half_o[:len(half_o) // 2]] for k, v in o_roles["eval"].items()}
    rep_o = {k: v[half_o[len(half_o) // 2:]] for k, v in o_roles["eval"].items()}

    mu_r = fit_outcome(r_roles["nuis"], spec=spec)
    mu_o = fit_outcome(o_roles["nuis"], spec=spec)
    prop_o = fit_propensity(o_roles["nuis"], seed=seed)

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
    basis = standardize_factory(basis_raw, r_roles["nuis"]["x"])

    s_tune = build_scores(r_roles["tune"], mu_r, mu_o)
    s_sel = build_scores(sel_r, mu_r, mu_o)
    s_rep = build_scores(rep_r, mu_r, mu_o)
    g = (lambda x: ghat_on(x, mu_o)) if ghat_fn is None else ghat_fn
    b = basis
    cand = [(lam, om, rho) for rho in ridges for (lam, om) in grid]

    bound = prediction_bound(s_tune.z0)
    g_raw = g
    g = lambda x: clip_prediction(g_raw(x), bound)
    rows = []
    for learner in ("DRF", "RF"):
        w, zt0, dzt = learner_parts(r_roles["tune"], s_tune, learner)
        scores_sel, scores_rep, betas = {}, {}, {}
        for rho in ridges:
            for om in sorted({c[1] for c in cand}):
                solve, _ = sieve_solver(b(r_roles["tune"]["x"]), w, zt0, dzt,
                                        g(r_roles["tune"]["x"]), b(o_roles["tune"]["x"]),
                                        np.ones(len(o_roles["tune"]["x"])),
                                        g(o_roles["tune"]["x"]), om, rho)
                for lam in sorted({c[0] for c in cand}):
                    beta = solve(lam)
                    betas[(lam, om, rho)] = beta
                    scores_sel[(lam, om, rho)] = validation_score(
                        s_sel.z0,
                        clip_prediction(b(sel_r["x"]) @ beta, bound), g(sel_r["x"]),
                        clip_prediction(b(sel_o["x"]) @ beta, bound), g(sel_o["x"]),
                        np.ones(len(sel_o["x"])), om)
                    scores_rep[(lam, om, rho)] = validation_score(
                        s_rep.z0,
                        clip_prediction(b(rep_r["x"]) @ beta, bound), g(rep_r["x"]),
                        clip_prediction(b(rep_o["x"]) @ beta, bound), g(rep_o["x"]),
                        np.ones(len(rep_o["x"])), om)
        selected = min(cand, key=lambda c: (scores_sel[c], c))
        # the trial-only rule chooses its ridge on the selection half, exactly as
        # the joint rule does, so the comparison is between equals
        rct_only = min([c for c in cand if c[0] == 0.0 and c[1] == 0.0],
                       key=lambda c: (scores_sel[c], c))
        if anchor is not None and predictions is not None:
            b_anchor = b(anchor)
            for tag, c in (("selected", selected), ("rct_only", rct_only)):
                predictions.setdefault((learner, tag), []).append(
                    clip_prediction(b_anchor @ betas[c], bound))
        rows.append({"learner": learner, "sieve": sieve,
                     "score_selected": scores_rep[selected],
                     "score_rct_only": scores_rep[rct_only],
                     "score_gap": scores_rep[selected] - scores_rep[rct_only],
                     "score_gap_in_sample": scores_sel[selected] - scores_sel[rct_only],
                     "selected_lambda": selected[0], "selected_omega": selected[1],
                     "selected_rho": selected[2], "rct_only_rho": rct_only[2],
                     "m_candidates": len(cand),
                     "n_select": len(sel_r["x"]), "n_report": len(rep_r["x"]),
                     "prediction_bound": bound})

    beta_eg = fb.experimental_grounding(b, r_roles["tune"], s_tune, o_roles["tune"], mu_o, prop_o)
    _, beta_dp = fb.domain_indicator_pool(b, r_roles["tune"], s_tune, o_roles["tune"], mu_o, prop_o)
    base = rows[0]["score_rct_only"]
    beta_dp_rct = fb.source_indicator_coefficients(beta_dp)
    for name, beta in (("experimental_grounding", beta_eg),
                       ("domain_indicator_pool", beta_dp_rct)):
        sc = validation_score(s_rep.z0,
                              clip_prediction(b(rep_r["x"]) @ beta, bound),
                              g(rep_r["x"]), clip_prediction(b(rep_o["x"]) @ beta, bound),
                              g(rep_o["x"]), np.ones(len(rep_o["x"])), 0.0)
        rows.append({"learner": name, "sieve": sieve, "score_selected": sc,
                     "score_rct_only": base, "score_gap": sc - base})
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replications", type=int, default=100)
    parser.add_argument("--stem", type=str, default="cate2_star_real")
    parser.add_argument("--sieves", type=str, default="spline3,neural")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()

    sieves = tuple(args.sieves.split(","))
    outcomes = ("mathk",) if args.smoke else ("mathk", "readk")
    sizes = (400,) if args.smoke else (200, 400)
    alphas = (0.8,) if args.smoke else (0.4, 0.8)

    rows, ate_rows, variability, started = [], [], [], time.time()
    for outcome in outcomes:
        cohort = fd.star_cohort(outcome, "pooled")
        anchor = cohort.pattern_x            # the retained covariate patterns
        anchor_weight = cohort.pattern_law   # weighted by the trial population law
        for n_rct in sizes:
            for alpha in alphas:
                stores: dict = {}
                for rep in range(args.replications):
                    rct, obs, extra = fd.star_cate2_replication(
                        cohort, n_rct, alpha, rep, n_rct_ate=n_rct + 2000)
                    label = {"outcome": outcome, "n_rct": n_rct, "alpha": alpha,
                             "replication": rep, **extra["diagnostics"]}
                    seed = stable_seed("cate2", outcome, n_rct, alpha, rep)
                    for sieve in sieves:
                        got = cate2_replication(cohort, rct, obs, extra["eval"], sieve, seed,
                                                spline_knots=3, anchor=anchor,
                                                predictions=stores.setdefault(sieve, {}))
                        for row in got:
                            row.update(label)
                            row["seed"] = seed
                        rows.extend(got)
                    ate_trial = extra["ate_rct"]
                    got = ate_replication(ate_trial, obs, cohort.ate_reference, ratio_unit, "unit",
                                          "linear", seed, with_baselines=True)
                    for row in got:
                        row.update(label)
                        row["n_trial_for_ate"] = len(ate_trial["x"])
                    ate_rows.extend(got)
                for sieve, store in stores.items():
                    for (learner, tag), stack in store.items():
                        matrix = np.vstack(stack)
                        variability.append({"outcome": outcome, "n_rct": n_rct, "alpha": alpha,
                                            "sieve": sieve, "learner": learner, "rule": tag,
                                            "replications": len(matrix),
                                            "anchor_points": matrix.shape[1],
                                            "prediction_variance":
                                                float(np.sum(anchor_weight * np.var(matrix, axis=0, ddof=1))),
                                            "prediction_spread":
                                                float(np.sum(anchor_weight * np.std(matrix, axis=0, ddof=1))),
                                            "prediction_variance_unweighted":
                                                float(np.mean(np.var(matrix, axis=0, ddof=1))),
                                            "mean_prediction_range":
                                                float(matrix.mean(axis=0).max() - matrix.mean(axis=0).min())})
                print(f"  {outcome} n_R={n_rct} alpha={alpha}  [{time.time() - started:.0f}s]", flush=True)

    frame = pd.DataFrame(rows)
    ate_frame = pd.DataFrame(ate_rows)
    summary = summarize(frame, ["outcome", "n_rct", "alpha", "sieve", "learner"],
                        {"score_gap": "score_gap", "score_selected": "score_selected",
                         "score_rct_only": "score_rct_only",
                         "selected_lambda": "selected_lambda", "selected_omega": "selected_omega",
                         "obs_size": "obs_size", "naive_obs_contrast": "naive_obs_contrast",
                         "n_trial_total": "n_trial_total", "n_eval_pool": "n_eval_pool",
                         "n_trial_for_ate": "n_trial_for_ate"})
    ate_summary = summarize(ate_frame, ["outcome", "n_rct", "alpha", "estimator"],
                            {"rmse": "error", "bias": "error", "coverage": "covered",
                             "estimate": "estimate", "ci_length": "ci_length"})
    paths = write_outputs(args.stem, frame, summary)
    ate_frame.to_csv(paths["replications"].parent / f"{args.stem}_ate_replications.csv", index=False)
    ate_summary.to_csv(paths["replications"].parent / f"{args.stem}_ate_summary.csv", index=False)
    pd.DataFrame(variability).to_csv(
        paths["replications"].parent / f"{args.stem}_prediction_variability.csv", index=False)
    print(f"wrote {paths['replications'].name}; {len(frame)} CATE rows, {len(ate_frame)} ATE rows "
          f"in {time.time() - started:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
