"""Runner for the manuscript-minimal plan.

One replication draws its randomness once and builds both cells from it, so the
two cells are paired.  Inside a cell every ratio route and both learners see the
same data, the same nuisance fits, the same basis and the same selection sample,
so an omega of zero gives the same prediction on every route.

Risk, squared bias and variance are all averages over one fixed truth grid, and
the plan's own identity holds:

    mean risk = squared bias + (K - 1) / K * variance.

The three quantities the plan needs for every summary and for the paired
bootstrap are G[k, l] = P_T[f_k f_l], v[k] = P_T[f_k tau] and t = P_T[tau^2].
They are sufficient for the mean risk, the squared bias and the variance at any
resampling weights, so the run stores them rather than two hundred prediction
vectors of fifty thousand points each.

Run:
    python3 -m confirm.mm_run --phase preflight
    python3 -m confirm.mm_run --phase smoke
    python3 -m confirm.mm_run --phase pilot
    python3 -m confirm.mm_run --phase main
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time

import numpy as np

from confirm import mm_design as md
from confirm import mm_engine as me
from confirm import provenance as pv

PLAN = "manuscript-minimal-v1"
PHASES = {"smoke": 1, "pilot": 20, "main": 200}
ROUTES = {"shared": ("unit_exact",),
          "shifted": ("exact", "classifier", "wrong_unit")}
RULES = ("rct_only", "lambda_only", "omega_only", "joint")
LEARNERS = ("DRF", "RF")
DIAGNOSTIC_CANDIDATE = (0.5, 0.5)


def seed_for(phase: str, replication: int, source: str, role: str) -> int:
    key = f"{PLAN}|{phase}|{replication}|{source}|{role}"
    return int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)


def group_key(cell: str, learner: str, route: str, rule: str) -> str:
    return f"{cell}|{learner}|{route}|{rule}"


def run_replication(phase: str, replication: int, params, truth, cells) -> tuple[dict, list, list]:
    """Returns truth-grid predictions by group, long rows, and failure records."""
    predictions: dict[str, np.ndarray] = {}
    rows: list[dict] = []
    failures: list[dict] = []
    started = time.time()

    draws = {}
    for source, sizes in (("RCT", md.RCT_ROLE_SIZES), ("OBS", md.OBS_ROLE_SIZES)):
        for role, n in sizes.items():
            rng = np.random.default_rng(seed_for(phase, replication, source, role))
            draws[(source, role)] = md.Draw.make(n, rng)

    for cell in cells:
        rct = {role: md.realise(draws[("RCT", role)], "RCT", params, cell)
               for role in md.RCT_ROLE_SIZES}
        obs = {role: md.realise(draws[("OBS", role)], "OBS", params, cell)
               for role in md.OBS_ROLE_SIZES}
        nuis = me.fit_nuisance(rct["nuisance"], obs["nuisance"], md.TRIAL_PROPENSITY)
        tune = me.trial_parts(rct["tuning"], nuis)
        sel = me.trial_parts(rct["selection"], nuis)
        b_obs_tune, ghat_obs_tune = me.obs_parts(obs["tuning"], nuis)
        b_obs_sel, ghat_obs_sel = me.obs_parts(obs["selection"], nuis)
        b_truth = nuis.basis(truth["x"])

        classifier_converged = True
        ratio_cache = {}
        for route in ROUTES[cell.name]:
            try:
                if route in ("unit_exact", "wrong_unit"):
                    tune_r = np.ones(len(b_obs_tune))
                    sel_r = np.ones(len(b_obs_sel))
                elif route == "exact":
                    tune_r = md.exact_ratio(obs["tuning"]["x"], params, cell)
                    sel_r = md.exact_ratio(obs["selection"]["x"], params, cell)
                else:
                    ratio_fn, classifier_converged = me.classifier_ratio(
                        rct["nuisance"]["x"], obs["nuisance"]["x"])
                    tune_r = ratio_fn(obs["tuning"]["x"])
                    sel_r = ratio_fn(obs["selection"]["x"])
                ratio_cache[route] = (tune_r, sel_r)
            except me.NumericalFailure as exc:
                failures.append({"replication": replication, "cell": cell.name,
                                 "route": route, "learner": "", "rule": "",
                                 "failure_class": "ratio", "detail": str(exc)})
                continue

            for key, value in me.ratio_diagnostics(sel_r).items():
                rows.append({"phase": phase, "replication": replication,
                             "cell": cell.name, "learner": "", "route": route,
                             "rule": "", "metric": f"ratio/{key}", "value": value,
                             "lambda": "", "omega": "", "status": "ok"})
            rows.append({"phase": phase, "replication": replication, "cell": cell.name,
                         "learner": "", "route": route, "rule": "",
                         "metric": "ratio/classifier_converged",
                         "value": float(classifier_converged), "lambda": "",
                         "omega": "", "status": "ok"})

            for learner in LEARNERS:
                try:
                    betas, scores = {}, {}
                    for lam, om in me.GRID:
                        beta = me.solve(learner, tune, b_obs_tune, ghat_obs_tune,
                                        tune_r, lam, om)
                        betas[(lam, om)] = beta
                        scores[(lam, om)] = me.selection_score(
                            beta, sel, b_obs_sel, ghat_obs_sel, sel_r, om)
                except me.NumericalFailure as exc:
                    failures.append({"replication": replication, "cell": cell.name,
                                     "route": route, "learner": learner, "rule": "",
                                     "failure_class": "solve", "detail": str(exc)})
                    continue

                for (lam, om), value in scores.items():
                    rows.append({"phase": phase, "replication": replication,
                                 "cell": cell.name, "learner": learner,
                                 "route": route, "rule": "candidate",
                                 "metric": "selection_score", "value": value,
                                 "lambda": lam, "omega": om, "status": "ok"})

                for rule in RULES:
                    lam, om = me.pick(scores, rule)
                    predictions[group_key(cell.name, learner, route, rule)] = (
                        (b_truth @ betas[(lam, om)]).astype(np.float32))
                    rows.append({"phase": phase, "replication": replication,
                                 "cell": cell.name, "learner": learner,
                                 "route": route, "rule": rule,
                                 "metric": "selected_lambda", "value": lam,
                                 "lambda": lam, "omega": om, "status": "ok"})
                    rows.append({"phase": phase, "replication": replication,
                                 "cell": cell.name, "learner": learner,
                                 "route": route, "rule": rule,
                                 "metric": "selected_omega", "value": om,
                                 "lambda": lam, "omega": om, "status": "ok"})

                # the grid oracle is a label only; it never enters selection
                risks = {key: float(np.mean((b_truth @ beta - truth["tau"]) ** 2))
                         for key, beta in betas.items()}
                best = min(risks, key=lambda k: (risks[k], k[0], k[1]))
                predictions[group_key(cell.name, learner, route, "grid_oracle")] = (
                    (b_truth @ betas[best]).astype(np.float32))

                # fixed-candidate diagnostic, no extra fit
                fixed = betas[DIAGNOSTIC_CANDIDATE]
                base = betas[(0.0, 0.0)]
                rows.append({"phase": phase, "replication": replication,
                             "cell": cell.name, "learner": learner, "route": route,
                             "rule": "diagnostic", "metric": "fixed_candidate_drift",
                             "value": (scores[DIAGNOSTIC_CANDIDATE] - scores[(0.0, 0.0)])
                                      - (risks[DIAGNOSTIC_CANDIDATE] - risks[(0.0, 0.0)]),
                             "lambda": DIAGNOSTIC_CANDIDATE[0],
                             "omega": DIAGNOSTIC_CANDIDATE[1], "status": "ok"})

    rows.append({"phase": phase, "replication": replication, "cell": "", "learner": "",
                 "route": "", "rule": "", "metric": "runtime_seconds",
                 "value": time.time() - started, "lambda": "", "omega": "",
                 "status": "ok"})
    return predictions, rows, failures


def summarise(store: dict[str, list[np.ndarray]], truth) -> dict:
    """G, v and t for every group, which is all the reporting needs."""
    tau = truth["tau"].astype(np.float64)
    t = float(np.mean(tau ** 2))
    out = {"t": t, "groups": {}}
    for key, stack in store.items():
        matrix = np.vstack(stack).astype(np.float64)
        n_grid = matrix.shape[1]
        gram = matrix @ matrix.T / n_grid
        v = matrix @ tau / n_grid
        out["groups"][key] = {"K": len(stack), "G": gram.tolist(), "v": v.tolist()}
    return out


def _parse(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", required=True, choices=["preflight", *PHASES])
    parser.add_argument("--replace", action="store_true")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = _parse(argv)
    code_files = sorted(pv.CODE.glob("*.py")) + sorted((pv.CODE / "confirm").glob("*.py"))
    fingerprint = pv.code_fingerprint(code_files)
    params = md.scm()
    cells = md.cells(params)

    if args.phase == "preflight":
        payload = {"mode": "preflight", "plan": PLAN,
                   "calculation_performed": False, "write_performed": False,
                   "phases": {k: {"replications_per_cell": v,
                                  "cell_records": 2 * v} for k, v in PHASES.items()},
                   "cells": [c.name for c in cells],
                   "routes": {k: list(v) for k, v in ROUTES.items()},
                   "learners": list(LEARNERS), "rules": list(RULES),
                   "grid": [list(g) for g in me.GRID], "ridge": me.RIDGE,
                   "trial_propensity": md.TRIAL_PROPENSITY,
                   "truth_grid": md.TRUTH_GRID,
                   "mahalanobis_shift": float(np.sqrt(
                       md.shift_direction(params)
                       @ np.linalg.inv(params.covariance)
                       @ md.shift_direction(params))),
                   "fingerprint": fingerprint,
                   "targets": {p: str(pv.CONFIRM / f"MM-{p}-{fingerprint[:12]}")
                               for p in PHASES},
                   "collision_state": {
                       p: ("present" if (pv.CONFIRM / f"MM-{p}-{fingerprint[:12]}").exists()
                           else "absent") for p in PHASES},
                   "qualification_status": "ready"}
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    replications = PHASES[args.phase]
    truth = md.truth_grid(params)
    store: dict[str, list[np.ndarray]] = {}
    rows, failures, started = [], [], time.time()
    scheduled = replications

    for replication in range(replications):
        preds, unit_rows, unit_failures = run_replication(
            args.phase, replication, params, truth, cells)
        rows.extend(unit_rows)
        failures.extend(unit_failures)
        for key, vector in preds.items():
            store.setdefault(key, []).append(vector)
        if replications <= 20 or (replication + 1) % 20 == 0:
            print(f"  replication {replication + 1}/{replications} "
                  f"[{time.time() - started:.0f}s]", flush=True)

    summary = summarise(store, truth)
    checks = validate(rows, failures, store, summary, scheduled, cells)

    name = f"MM-{args.phase}-{fingerprint[:12]}"
    stage = pv.Stage(name, fingerprint)
    scratch = stage.open()
    (scratch / "raw").mkdir()
    with (scratch / "raw" / "rows.csv").open("w", newline="") as handle:
        fields = ["phase", "replication", "cell", "learner", "route", "rule",
                  "metric", "lambda", "omega", "value", "status"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    pv.write_json(scratch / "sufficient_statistics.json", summary)
    pv.write_json(scratch / "config.json", {
        "plan": PLAN, "phase": args.phase, "replications_per_cell": replications,
        "trial_propensity": md.TRIAL_PROPENSITY, "ridge": me.RIDGE,
        "grid": [list(g) for g in me.GRID], "truth_grid": md.TRUTH_GRID,
        "rct_role_sizes": md.RCT_ROLE_SIZES, "obs_role_sizes": md.OBS_ROLE_SIZES,
        "routes": {k: list(v) for k, v in ROUTES.items()}, "rules": list(RULES),
        "scm_id": md.SCM_ID})
    pv.write_json(scratch / "code_manifest.json",
                  {"fingerprint": fingerprint,
                   "files": pv.manifest_for(code_files, pv.CODE)})
    pv.write_json(scratch / "environment.json",
                  {**pv.environment(), "git_head": pv.git_head()})
    pv.write_json(scratch / "checks.json", checks)
    pv.write_json(scratch / "failures.json", {"failures": failures})
    pv.write_json(scratch / "runtime.json", {"total_seconds": time.time() - started})

    for row in checks["rows"]:
        print(f"[{row['status']}] {row['check']}: {row['detail']}")
    if checks["failed"]:
        print(f"\n{args.phase} STOP: {checks['failed']} of {checks['total']} checks failed")
        return 1
    final = stage.promote(allow_replace=args.replace)
    print(f"\n{args.phase} PASS: {checks['passed']}/{checks['total']} checks; "
          f"{len(rows)} rows written to {final.relative_to(pv.PROJECT)}")
    return 0


def validate(rows, failures, store, summary, scheduled, cells) -> dict:
    out = []

    def add(name, passed, detail):
        out.append({"check": name, "status": "PASS" if passed else "FAIL",
                    "detail": detail})

    expected = set()
    for cell in cells:
        for route in ROUTES[cell.name]:
            for learner in LEARNERS:
                for rule in (*RULES, "grid_oracle"):
                    expected.add(group_key(cell.name, learner, route, rule))
    add("groups/complete", set(store) == expected,
        f"{len(store)} groups present, {len(expected)} expected; "
        f"missing {sorted(expected - set(store))[:3]}")

    counts = {len(v) for v in store.values()}
    add("accounting/scheduled_equals_recorded",
        counts == {scheduled} and not failures,
        f"every group carries {sorted(counts)} replications against {scheduled} "
        f"scheduled, with {len(failures)} failures")

    # the plan's decomposition identity, on the stored sufficient statistics
    worst = 0.0
    for key, block in summary["groups"].items():
        gram = np.array(block["G"]); v = np.array(block["v"]); k = block["K"]
        a = np.full(k, 1.0 / k)
        mean_risk = float(np.sum(a * (np.diag(gram) - 2 * v + summary["t"])))
        bias2 = float(a @ gram @ a - 2 * a @ v + summary["t"])
        variance = float((np.sum(a * np.diag(gram)) - a @ gram @ a) * k / max(k - 1, 1))
        gap = abs(mean_risk - (bias2 + (k - 1) / k * variance)) / max(abs(mean_risk), 1e-12)
        worst = max(worst, gap)
    add("decomposition/identity_holds", worst <= 1e-10,
        f"largest relative gap between the mean risk and its decomposition {worst:.2e}")

    # omega zero must give the same prediction on every route inside a cell
    gaps = []
    for cell in cells:
        if len(ROUTES[cell.name]) < 2:
            continue
        for learner in LEARNERS:
            base = store[group_key(cell.name, learner, ROUTES[cell.name][0], "rct_only")]
            for route in ROUTES[cell.name][1:]:
                other = store[group_key(cell.name, learner, route, "rct_only")]
                gaps.append(max(float(np.max(np.abs(b - o)))
                                for b, o in zip(base, other)))
    add("routes/omega_zero_agrees", not gaps or max(gaps) == 0.0,
        f"largest difference between routes at omega zero {max(gaps) if gaps else 0.0:.2e}")

    selected = [r for r in rows if r["metric"] == "selected_omega"
                and r["rule"] == "rct_only"]
    add("rules/rct_only_is_the_origin",
        all(r["value"] == 0.0 for r in selected) and bool(selected),
        f"{len(selected)} rct_only rows, all at omega zero")

    ratio_rows = [r for r in rows if r["metric"] == "ratio/mean"]
    shared = [r["value"] for r in ratio_rows if r["cell"] == "shared"]
    add("ratio/shared_is_exactly_one", bool(shared) and max(abs(v - 1.0) for v in shared) == 0.0,
        f"{len(shared)} shared-cell ratio means, largest deviation from one "
        f"{max(abs(v - 1.0) for v in shared) if shared else float('nan'):.2e}")

    return {"rows": out, "total": len(out),
            "passed": sum(1 for r in out if r["status"] == "PASS"),
            "failed": sum(1 for r in out if r["status"] == "FAIL")}


if __name__ == "__main__":
    sys.exit(main())
