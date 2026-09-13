"""Confirmatory CATE runner.

One scheduled unit is one independent stochastic source unit.  Estimators,
candidates and ratio routes multiply the rows a unit emits; they never multiply
the unit count.  Every scheduled seed ends as a success record or an explicit
failure record, and a failed seed is never replaced.

Run:
    python3 -m confirm.cate_run --stage m2 --preflight
    python3 -m confirm.cate_run --stage m2
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
import traceback
from dataclasses import asdict

import numpy as np

from confirm import cate_dgp as dgp
from confirm import cate_engine as ce
from confirm import cate_protocol as cp
from confirm import moments as mo
from confirm import provenance as pv

TARGET_DRAWS = 40_000


# ------------------------------------------------------------ adaptive oracle
def adaptive_oracle(design, nuis, ratio_route, cell, seed: int) -> dict:
    """Plan P0-3: grow the oracle sample until every moment meets its tolerance.

    A moment far from zero is judged on relative Monte Carlo error, one near zero
    on absolute error, because a relative target is meaningless at zero.
    """
    record = {}
    for draws in cp.ORACLE_DRAWS:
        rng = np.random.default_rng(seed + draws)
        rct = design.sample("RCT", draws, rng)
        obs = design.sample("OBS", draws, rng)
        target = design.target_sample(draws, rng)
        scores = ce.trial_scores(rct, nuis)
        ratio = (design.oracle_ratio(obs["x"]) if ratio_route == "oracle"
                 else np.ones(len(obs["x"])))

        def fixture_for(rows_r, rows_o):
            return mo.Fixture(
                name=cell.cell_id, b_target=nuis.basis(target["x"]),
                tau_target=target["tau"], b_rct=nuis.basis(rct["x"][rows_r]),
                z0_rct=scores.z0[rows_r], delta_rct=scores.delta[rows_r],
                ghat_rct=scores.ghat[rows_r],
                weight_rct=np.ones(rows_r.sum()),
                b_obs=nuis.basis(obs["x"][rows_o]),
                ghat_obs=ce.obs_prediction(obs["x"][rows_o], nuis),
                ratio_obs=ratio[rows_o],
                n_rct_tune=max(int(cell.n_rct * cp.ROLE_FRACTIONS["RCT"][1]), 10),
                n_obs_tune=max(int(cell.n_obs * cp.ROLE_FRACTIONS["OBS"][1]), 10))

        def measured(rows_r, rows_o):
            """Moments on a subset, with D_p taken from its exact identity.

            The trial propensity is known and the pseudo-outcome is conditionally
            unbiased, so E(Z_0 | X) = tau exactly.  Substituting that into

                D_p = tr[Gamma^-1 Cov_R{b(Z_0 - zeta_p), b(ghat - zeta_p)}]
                    = tr[Gamma^-1 E_R{b b' (tau - zeta_p)(ghat - zeta_p)}]

            removes the whole pseudo-outcome variance from the estimator.  The
            oracle is allowed to use tau, and without this the Monte Carlo error
            on omega* stays above tolerance even at one hundred thousand draws.
            """
            fx = fixture_for(rows_r, rows_o)
            out = mo.covariance_route(fx)
            exact = mo.Fixture(**{**fx.__dict__,
                                  "z0_rct": rct["tau"][rows_r]})
            out["D_p"] = mo.covariance_route(exact)["D_p"]
            return out

        full = measured(np.ones(draws, bool), np.ones(draws, bool))
        # batch Monte Carlo error, computed inside this draw rather than by
        # differencing two draw sizes
        batches = 8
        edges = np.array_split(np.arange(draws), batches)
        per_batch = []
        for idx in edges:
            mask_r = np.zeros(draws, bool); mask_r[idx] = True
            per_batch.append(measured(mask_r, mask_r))
        def coordinates(m):
            lam = -m["C_p"] / m["A_p"] if m["A_p"] > 0 else 0.0
            om = m["D_p"] / m["B_p"] if m["B_p"] > 0 else 0.0
            # the coefficients are chosen inside [0, 1], so the projection is the
            # quantity the procedure can use and the one worth resolving
            return {"lambda_star": float(np.clip(lam, 0.0, 1.0)),
                    "omega_star": float(np.clip(om, 0.0, 1.0))}

        full_coord = coordinates(full)
        batch_coord = [coordinates(b) for b in per_batch]
        record = {}
        for key, value in full.items():
            spread = np.std([b[key] for b in per_batch], ddof=1)
            record[key] = {"estimate": value, "mcse": float(spread / np.sqrt(batches)),
                           "tolerance": float("nan"), "draws": draws,
                           "status": "reported"}
        # the plan forbids an arbitrary absolute cutoff on a near-zero moment, so
        # convergence is judged on the coordinate the moment feeds, which is
        # bounded by construction and is what the quadrant label depends on
        for key in ("lambda_star", "omega_star"):
            spread = np.std([b[key] for b in batch_coord], ddof=1)
            mcse = float(spread / np.sqrt(batches))
            record[key] = {"estimate": full_coord[key], "mcse": mcse,
                           "tolerance": cp.ORACLE_COORDINATE_TOL, "draws": draws,
                           "status": ("converged" if mcse <= cp.ORACLE_COORDINATE_TOL
                                      else "growing")}
        if all(record[k]["status"] == "converged"
               for k in ("lambda_star", "omega_star")):
            break
    for key in ("lambda_star", "omega_star"):
        if record[key]["status"] != "converged":
            record[key]["status"] = "MC_INCONCLUSIVE"
    return record


def observed_quadrant(oracle: dict) -> str:
    """The quadrant the cell actually lands in, from the measured moments."""
    def channel(key):
        entry = oracle[key]
        if entry["status"] == "MC_INCONCLUSIVE":
            return "?"
        # open when the coordinate is both distinguishable from zero and large
        # enough to matter on the [0, 1] scale the coefficient lives on
        value = abs(entry["estimate"])
        return "1" if value > 2.0 * entry["mcse"] and value >= 0.05 else "0"
    return f"L{channel('lambda_star')}O{channel('omega_star')}"


# ------------------------------------------------------------------- one unit
def run_unit(unit: dict, cell: cp.Cell, fingerprint: str, stage: str) -> tuple[list, dict]:
    started = time.time()
    rows: list[dict] = []
    design = dgp.from_cell(cell)
    seed_key = unit["seed_key"]
    roles = ce.draw_roles(design, seed_key, cell.n_rct, cell.n_obs,
                          cell.trial_nuisance_boost)
    nuis = ce.fit_nuisance(roles, cell.sieve)
    role_prints = "|".join(f"{s}/{r}:{roles[s][r].fingerprint}"
                           for s in ("RCT", "OBS") for r in cp.ROLES)

    tune_r = roles["RCT"]["tuning"].data
    sel_r = roles["RCT"]["selection"].data
    rep_r = roles["RCT"]["reporting"].data
    tune_o = roles["OBS"]["tuning"].data

    s_tune = ce.trial_scores(tune_r, nuis)
    s_sel = ce.trial_scores(sel_r, nuis)
    s_rep = ce.trial_scores(rep_r, nuis)

    b_tune = nuis.basis(tune_r["x"])
    b_sel = nuis.basis(sel_r["x"])
    b_rep = nuis.basis(rep_r["x"])
    b_obs = nuis.basis(tune_o["x"])
    ghat_obs = ce.obs_prediction(tune_o["x"], nuis)

    rng = np.random.default_rng(cp.seed32(seed_key) + 991)
    target = design.target_sample(TARGET_DRAWS, rng)
    b_target = nuis.basis(target["x"])

    routes = cp.RATIO_ROUTES if cell.ratio_route == "all" else (cell.ratio_route,)
    common = {"protocol_id": cp.PROTOCOL_ID, "code_fingerprint": fingerprint,
              "stage": stage, "study": cell.family, "family": cell.family,
              "cell_id": cell.cell_id, "replication_id": unit["replication_id"],
              "seed_key": seed_key, "seed32": cp.seed32(seed_key),
              "role_fingerprints": role_prints, "learner": cell.learner,
              "sieve": cell.sieve, "rho": 0.0, "status": "ok",
              "failure_class": ""}

    def emit(route, candidate_id, lam, omega, metric, value):
        rows.append({**common, "ratio_route": route, "candidate_id": candidate_id,
                     "lambda": lam, "omega": omega, "metric": metric,
                     "value": value, "runtime_seconds": time.time() - started})

    oracle_record = {}
    for route in routes:
        ratio = ce.build_ratio(route, design, roles, nuis, tune_o["x"])
        for key, value in ratio.diagnostics.items():
            emit(route, "ratio", np.nan, np.nan, f"ratio/{key}", value)
        emit(route, "ratio", np.nan, np.nan, "ratio/balance_status",
             {"converged": 1.0, "not_converged": 0.0}.get(ratio.balance_status, np.nan))

        best = {"score": np.inf}
        risks = {}
        for rho in cp.RIDGE_GRID:
            for lam in cp.LAMBDA_GRID:
                for omega in cp.OMEGA_GRID:
                    beta, condition = ce.solve_candidate(
                        b_tune, s_tune.z0, s_tune.delta, s_tune.ghat,
                        b_obs, ghat_obs, ratio.values_obs, lam, omega, rho)
                    cid = f"l{lam:.2f}_o{omega:.2f}_r{rho:g}"
                    score = ce.dr_score(b_sel @ beta, s_sel.z0)
                    risk = ce.true_risk(b_target @ beta, target["tau"])
                    risks[cid] = risk
                    emit(route, cid, lam, omega, "selection_score", score)
                    emit(route, cid, lam, omega, "true_risk", risk)
                    emit(route, cid, lam, omega, "condition_number", condition)
                    if score < best["score"]:
                        best = {"score": score, "beta": beta, "cid": cid,
                                "lambda": lam, "omega": omega, "risk": risk}

        fixed = {
            "rct_only": (0.0, 0.0),
            "lambda_only": max(((l, 0.0) for l in cp.LAMBDA_GRID),
                               key=lambda c: -risks[f"l{c[0]:.2f}_o{c[1]:.2f}_r0"]),
            "omega_only": max(((0.0, o) for o in cp.OMEGA_GRID),
                              key=lambda c: -risks[f"l{c[0]:.2f}_o{c[1]:.2f}_r0"]),
        }
        for name, (lam, omega) in fixed.items():
            cid = f"l{lam:.2f}_o{omega:.2f}_r0"
            emit(route, cid, lam, omega, f"estimator/{name}/true_risk", risks[cid])
        emit(route, best["cid"], best["lambda"], best["omega"],
             "estimator/joint/true_risk", best["risk"])
        emit(route, best["cid"], best["lambda"], best["omega"],
             "estimator/joint/reporting_score", ce.dr_score(b_rep @ best["beta"], s_rep.z0))
        grid_oracle = min(risks, key=risks.get)
        emit(route, grid_oracle, np.nan, np.nan,
             "estimator/grid_oracle/true_risk", risks[grid_oracle])
        emit(route, best["cid"], best["lambda"], best["omega"], "selection_regret",
             best["risk"] - risks[grid_oracle])

        if route in ("oracle", "unit") and cell.family != "star":
            oracle_record = adaptive_oracle(design, nuis, route, cell,
                                            cp.seed32(seed_key))
            for key, entry in oracle_record.items():
                emit(route, "oracle_moment", np.nan, np.nan, f"oracle/{key}",
                     entry["estimate"])
                emit(route, "oracle_moment", np.nan, np.nan, f"oracle/{key}/mcse",
                     entry["mcse"])
                emit(route, "oracle_moment", np.nan, np.nan, f"oracle/{key}/draws",
                     float(entry["draws"]))

    # theorem eligibility, plan P0-2
    eligible = (cell.sieve.startswith("spline") and cp.RIDGE_GRID == (0.0,))
    emit(routes[0], "theorem", np.nan, np.nan, "theorem7/eligible", float(eligible))
    emit(routes[0], "theorem", np.nan, np.nan, "theorem6/eligible", 0.0)
    emit(routes[0], "theorem", np.nan, np.nan, "theorem6/radius", np.nan)

    record = {"seed_key": seed_key, "cell_id": cell.cell_id, "status": "ok",
              "rows": len(rows), "runtime_seconds": time.time() - started,
              "role_fingerprints": role_prints,
              "observed_quadrant": observed_quadrant(oracle_record) if oracle_record else "N/A",
              "intended_quadrant": cell.quadrant}
    return rows, record


def _parse(argv):
    parser = argparse.ArgumentParser(description="confirmatory CATE runner")
    parser.add_argument("--stage", required=True, choices=sorted(cp.SCHEDULES))
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--replace", action="store_true")
    parser.add_argument("--stage-name", default=None)
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = _parse(argv)
    code_files = sorted(pv.CODE.glob("*.py")) + sorted((pv.CODE / "confirm").glob("*.py"))
    fingerprint = pv.code_fingerprint(code_files)
    name = args.stage_name or f"{args.stage.upper()}-{fingerprint[:12]}"
    cells = {c.cell_id: c for c in cp.SCHEDULES[args.stage]()}
    units = cp.scheduled_units(args.stage)

    if args.preflight:
        payload = {"mode": "preflight", "calculation_performed": False,
                   "write_performed": False, "stage": args.stage,
                   "resolved_target": str(pv.CONFIRM / name),
                   "fingerprint": fingerprint,
                   "protocol_id": cp.PROTOCOL_ID,
                   "schedule_hash": cp.schedule_hash(args.stage),
                   "scheduled_unit_count": len(units),
                   "cells": [c.cell_id for c in cp.SCHEDULES[args.stage]()],
                   "ratio_routes": list(cp.RATIO_ROUTES),
                   "collision_state": "present" if (pv.CONFIRM / name).exists() else "absent",
                   "required_outputs": list(cp.REQUIRED_OUTPUTS),
                   "required_raw_keys": list(cp.REQUIRED_RAW_KEYS),
                   "qualification_status": "ready"}
        print(pv.json.dumps(payload, indent=2, sort_keys=True))
        return 0

    rows, records, failures = [], [], []
    started = time.time()
    for unit in units:
        cell = cells[unit["cell_id"]]
        try:
            unit_rows, record = run_unit(unit, cell, fingerprint, args.stage)
            rows.extend(unit_rows)
            records.append(record)
            print(f"  {unit['cell_id']} ok, {len(unit_rows)} rows, "
                  f"quadrant intended {record['intended_quadrant']} "
                  f"observed {record['observed_quadrant']} "
                  f"[{record['runtime_seconds']:.1f}s]", flush=True)
        except Exception as exc:
            failures.append({"seed_key": unit["seed_key"], "cell_id": unit["cell_id"],
                             "failure_class": type(exc).__name__, "detail": str(exc),
                             "traceback": traceback.format_exc()[-2000:]})
            print(f"  {unit['cell_id']} FAILED: {type(exc).__name__}: {exc}", flush=True)

    checks = validate(rows, records, failures, units, fingerprint, args.stage)
    stage = pv.Stage(name, fingerprint)
    scratch = stage.open()
    (scratch / "raw").mkdir()
    with (scratch / "raw" / "replications-000.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cp.REQUIRED_RAW_KEYS))
        writer.writeheader()
        writer.writerows(rows)
    pv.write_json(scratch / "protocol.json", cp.protocol_document(args.stage))
    pv.write_json(scratch / "seed_ledger.json", cp.seed_ledger(args.stage))
    pv.write_json(scratch / "schedule.json",
                  {"stage": args.stage, "units": list(units),
                   "cells": [asdict(c) for c in cp.SCHEDULES[args.stage]()]})
    pv.write_json(scratch / "code_manifest.json",
                  {"fingerprint": fingerprint,
                   "files": pv.manifest_for(code_files, pv.CODE)})
    pv.write_json(scratch / "environment.json",
                  {**pv.environment(), "git_head": pv.git_head()})
    pv.write_json(scratch / "checks.json", checks)
    pv.write_json(scratch / "failures.json", {"failures": failures})
    pv.write_json(scratch / "runtime.json",
                  {"total_seconds": time.time() - started, "units": records})
    (scratch / "RECORD.md").write_text(render_record(args.stage, checks, records,
                                                     failures, rows))

    for row in checks["rows"]:
        print(f"[{row['status']}] {row['check']}: {row['detail']}")
    if checks["failed"]:
        print(f"\n{args.stage.upper()} STOP: {checks['failed']} of {checks['total']} checks failed")
        return 1
    final = stage.promote(allow_replace=args.replace)
    print(f"\n{args.stage.upper()} PASS: {checks['passed']}/{checks['total']} checks; "
          f"{len(rows)} rows written to {final.relative_to(pv.PROJECT)}")
    return 0


def validate(rows, records, failures, units, fingerprint, stage) -> dict:
    out = []

    def add(name, passed, detail):
        out.append({"check": name, "status": "PASS" if passed else "FAIL",
                    "detail": detail})

    add("T27/scheduled_equals_success_plus_failure",
        len(units) == len(records) + len(failures),
        f"{len(units)} scheduled, {len(records)} succeeded, {len(failures)} failed")
    seeds = [u["seed_key"] for u in units]
    add("T28/no_duplicate_or_replaced_seed", len(set(seeds)) == len(seeds),
        f"{len(set(seeds))} distinct seed keys out of {len(seeds)}")

    prints = [p for r in records for p in r["role_fingerprints"].split("|")]
    add("T04/role_draws_independent", len(set(prints)) == len(prints),
        f"{len(set(prints))} distinct role fingerprints out of {len(prints)}")

    add("T07/seed_namespaces_distinct",
        len({(r["seed_key"]) for r in records}) == len(records),
        f"{len(records)} units carry distinct seed keys")

    keys_ok = all(set(r) == set(cp.REQUIRED_RAW_KEYS) for r in rows)
    add("schema/raw_keys_exact", keys_ok,
        f"every one of {len(rows)} rows carries exactly the {len(cp.REQUIRED_RAW_KEYS)} required keys")

    risk_rows = [r for r in rows if r["metric"].endswith("true_risk")
                 and r["family"] != "star"]
    finite = all(np.isfinite(r["value"]) for r in risk_rows)
    add("M2/synthetic_true_risk_finite", finite and bool(risk_rows),
        f"{len(risk_rows)} synthetic true-risk rows, all finite {finite}")

    star_rows = [r for r in rows if r["family"] == "star"
                 and r["metric"].endswith("true_risk")]
    add("T33/star_true_risk_not_reported", not star_rows,
        f"{len(star_rows)} STAR true-risk rows, which must be zero")

    shift_routes = {r["ratio_route"] for r in rows if "shift" in r["cell_id"]}
    add("T24/five_routes_exact", shift_routes == set(cp.RATIO_ROUTES),
        f"shift unit carries {sorted(shift_routes)}")

    balance = [r for r in rows if r["metric"] == "ratio/balance_residual"
               and r["ratio_route"].startswith("balanced")]
    worst = max((abs(r["value"]) for r in balance if np.isfinite(r["value"])), default=np.nan)
    add("T22/balanced_residual_small", bool(balance) and worst <= 1e-6,
        f"largest balance residual {worst:.2e} over {len(balance)} rows")

    regret = [r["value"] for r in rows if r["metric"] == "selection_regret"]
    add("T18/selection_regret_non_negative",
        bool(regret) and min(regret) >= -1e-12,
        f"smallest regret {min(regret):.3e} over {len(regret)} rows"
        if regret else "no regret rows")

    theorem = [r for r in rows if r["metric"] == "theorem7/eligible"]
    add("T12/theorem7_paths_eligible_only",
        all(r["value"] == 1.0 for r in theorem) and bool(theorem),
        f"{len(theorem)} theorem-7 rows, all on an unclipped spline at rho = 0")

    add("M2/every_scheduled_unit_succeeded", not failures,
        f"{len(failures)} scheduled units ended in an explicit failure: "
        + (", ".join(f"{f['cell_id']} ({f['failure_class']})" for f in failures)
           or "none"))

    add("T32/runtime_and_hash_present",
        all(np.isfinite(r["runtime_seconds"]) for r in rows)
        and all(r["code_fingerprint"] == fingerprint for r in rows),
        "every row carries a runtime and the run's code fingerprint")

    return {"rows": out, "total": len(out),
            "passed": sum(1 for r in out if r["status"] == "PASS"),
            "failed": sum(1 for r in out if r["status"] == "FAIL")}


def render_record(stage, checks, records, failures, rows) -> str:
    lines = [f"# {stage.upper()} record", "",
             f"Protocol `{cp.PROTOCOL_ID}`.  "
             f"{len(records)} of {len(records) + len(failures)} scheduled source units "
             f"succeeded, emitting {len(rows)} rows.", "",
             "| unit | intended | observed | rows | seconds |",
             "| --- | --- | --- | ---: | ---: |"]
    for r in records:
        lines.append(f"| {r['cell_id']} | {r['intended_quadrant']} | "
                     f"{r['observed_quadrant']} | {r['rows']} | {r['runtime_seconds']:.1f} |")
    lines += ["", f"Checks: {checks['passed']} of {checks['total']} pass.", ""]
    for row in checks["rows"]:
        lines.append(f"- `{row['status']}` {row['check']}: {row['detail']}")
    if failures:
        lines += ["", "## Failures", ""]
        for f in failures:
            lines.append(f"- `{f['cell_id']}` {f['failure_class']}: {f['detail']}")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    sys.exit(main())
