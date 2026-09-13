"""M1 driver: deterministic mechanism micro fixtures.

Plan section 9, milestone M1.  These are not statistical replications.  Each
check is an algebraic identity on a fixed array, so the tolerances are numerical
rather than statistical.

Go conditions, each executed below:

1. the independent formula and the covariance-route formula agree to 1e-10
2. the oracle coordinates and the fixed-coordinate risk ordering match the
   fixture design
3. deliberate sign, role and ratio mutations are detected

Run: python3 -m confirm.m1_fixtures
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import sys

import numpy as np

from confirm import fixtures as fxs
from confirm import m1_production_tie as tie
from confirm import m1_protocol as protocol
from confirm import moments as mo
from confirm import provenance as pv

FORMULA_TOL = 1e-10
ZERO_TOL = 1e-9          # a fixture zero is exact by construction, not estimated
GRID = np.linspace(0.0, 1.0, 201)


def record_text(check_count: int) -> str:
    return (
        "# M1 deterministic mechanism record\n\n"
        f"Protocol: `{protocol.PROTOCOL_ID}`. Checks: {check_count}/{check_count} PASS.\n\n"
        "The omega span result is restricted to DRF with an exact score, fixed basis, "
        "exact transport, and population first-order moments. RF scalar proportionality "
        "is checked only under constant propensity with all other laws fixed. The L0 "
        "sign flip is null invariance, not fault detection. These deterministic fixtures "
        "are mechanism checks, not performance experiments or manuscript evidence.\n"
    )


class Check:
    def __init__(self) -> None:
        self.rows: list[dict] = []

    def record(self, name: str, passed: bool, detail: str) -> None:
        self.rows.append({"check": name, "status": "PASS" if passed else "FAIL",
                          "detail": detail})

    @property
    def failures(self) -> list[dict]:
        return [r for r in self.rows if r["status"] == "FAIL"]


def check_projection(check: Check, fx: mo.Fixture) -> None:
    """The normal equation E_R{b (tau - zeta_p)} = 0 defines the target."""
    beta_p = fx.projection()
    residual = (fx.b_target * (fx.tau_target - fx.b_target @ beta_p)[:, None]).mean(axis=0)
    worst = float(np.max(np.abs(residual)))
    check.record(f"{fx.name}/normal_equation", worst <= 1e-12,
                 f"largest normal-equation residual {worst:.2e}")


def check_routes(check: Check, fx: mo.Fixture) -> dict[str, float]:
    """Route 1 against route 2, and the quadratic against the direct risk."""
    first = mo.covariance_route(fx)
    second = mo.whitened_route(fx)
    worst = max(abs(first[k] - second[k]) for k in first)
    check.record(f"{fx.name}/two_routes_agree", worst <= FORMULA_TOL,
                 f"largest moment difference {worst:.2e} against {FORMULA_TOL:.0e}")

    gaps = []
    for lam in (0.0, 0.25, 0.6, 1.0):
        for omega in (0.0, 0.4, 1.0):
            gaps.append(abs(mo.direct_risk(fx, lam, omega)
                            - mo.quadratic_risk(first, fx, lam, omega)))
    worst_gap = max(gaps)
    check.record(f"{fx.name}/quadratic_identity", worst_gap <= FORMULA_TOL,
                 f"largest gap between Corollary 7.1 and the risk definition "
                 f"{worst_gap:.2e} over 12 coordinate pairs")
    return first


def check_quadrant(check: Check, name: str, moments: dict[str, float],
                   fx: mo.Fixture) -> None:
    """The moment conditions of plan section 7.2, and the risk ordering."""
    wants_lambda, wants_omega = name[1] == "1", name[3] == "1"
    lam_star, omega_star = mo.oracle_coordinates(moments)

    if wants_lambda:
        ok = moments["C_p"] < -ZERO_TOL and 0.0 < lam_star < 1.0
        check.record(f"{name}/lambda_opportunity", ok,
                     f"C_p {moments['C_p']:.5f} and lambda* {lam_star:.4f} interior")
    else:
        ok = abs(moments["C_p"]) <= ZERO_TOL and lam_star <= ZERO_TOL
        check.record(f"{name}/no_lambda_opportunity", ok,
                     f"C_p {moments['C_p']:.2e} within {ZERO_TOL:.0e} and lambda* {lam_star:.4f}")

    if wants_omega:
        ok = moments["D_p"] > ZERO_TOL and 0.0 < omega_star < 1.0
        check.record(f"{name}/omega_opportunity", ok,
                     f"D_p {moments['D_p']:.5f} and omega* {omega_star:.4f} interior")
    else:
        ok = abs(moments["D_p"]) <= ZERO_TOL and omega_star <= ZERO_TOL
        check.record(f"{name}/no_omega_opportunity", ok,
                     f"D_p {moments['D_p']:.2e} within {ZERO_TOL:.0e} and omega* {omega_star:.4f}")

    # the oracle coordinates must actually minimise the risk surface
    surface = np.array([[mo.direct_risk(fx, float(l), float(o)) for o in GRID] for l in GRID])
    flat = int(np.argmin(surface))
    best_lambda, best_omega = GRID[flat // len(GRID)], GRID[flat % len(GRID)]
    close = abs(best_lambda - lam_star) <= 2 / (len(GRID) - 1) and \
        abs(best_omega - omega_star) <= 2 / (len(GRID) - 1)
    check.record(f"{name}/oracle_minimises_risk", close,
                 f"grid minimum at ({best_lambda:.3f}, {best_omega:.3f}) against "
                 f"the formula's ({lam_star:.3f}, {omega_star:.3f})")

    # fixed-coordinate ablation: each channel alone, at its own oracle value
    base = mo.direct_risk(fx, 0.0, 0.0)
    lam_only = mo.direct_risk(fx, lam_star, 0.0)
    om_only = mo.direct_risk(fx, 0.0, omega_star)
    both = mo.direct_risk(fx, lam_star, omega_star)
    improves = lambda value: (base - value) / base
    order_ok = both <= min(lam_only, om_only) + 1e-15 and min(lam_only, om_only) <= base + 1e-15
    check.record(f"{name}/fixed_coordinate_ordering", order_ok,
                 f"risk falls {improves(lam_only) * 100:.2f}% with lambda alone, "
                 f"{improves(om_only) * 100:.2f}% with omega alone, "
                 f"{improves(both) * 100:.2f}% with both")


def check_mutations(check: Check, name: str, truth: dict[str, float]) -> None:
    """Three deliberate faults, each of which the moments must reveal."""
    base = fxs.quadrant(name)

    flipped = dataclasses.replace(base, delta_rct=-base.delta_rct)
    moved = mo.covariance_route(flipped)["C_p"]
    if abs(truth["C_p"]) > ZERO_TOL:
        detected = abs(moved - truth["C_p"]) > 1e-6 and np.sign(moved) == -np.sign(truth["C_p"])
        check.record(f"{name}/mutation_sign_detected", detected,
                     f"under L1, flipping Delta changed C_p from {truth['C_p']:.5f} "
                     f"to {moved:.5f} with the opposite sign")
    else:
        invariant = abs(moved - truth["C_p"]) <= ZERO_TOL
        check.record(f"{name}/sign_flip_null_invariance", invariant,
                     f"under L0, C_p stayed at its null value ({moved:.2e}); this is "
                     "null invariance, not fault detection")

    # Exchanging the two arguments of the covariance is a no-op, because
    # tr(Gamma^-1 C') equals tr(Gamma^-1 C) for a symmetric Gamma^-1, so a role
    # fault has to be a substitution rather than a swap.
    substituted = dataclasses.replace(base, z0_rct=base.ghat_rct)
    moved_d = mo.covariance_route(substituted)["D_p"]
    check.record(f"{name}/mutation_role_substitution", abs(moved_d - truth["D_p"]) > 1e-6,
                 f"using ghat where the base pseudo-outcome belongs moved D_p from "
                 f"{truth['D_p']:.5f} to {moved_d:.5f}")

    # forgetting to centre the observational prediction on the projection
    beta_p = base.projection()
    uncentred = dataclasses.replace(
        base, ghat_rct=base.ghat_rct + base.b_rct @ beta_p,
        ghat_obs=base.ghat_obs + base.b_obs @ beta_p)
    moved_b = mo.covariance_route(uncentred)["B_p"]
    check.record(f"{name}/mutation_missing_projection", abs(moved_b - truth["B_p"]) > 1e-6,
                 f"leaving zeta_p out of the prediction gap moved B_p from "
                 f"{truth['B_p']:.5f} to {moved_b:.5f}")

    rescaled = dataclasses.replace(base, ratio_obs=base.ratio_obs * 2.0)
    moved_b = mo.covariance_route(rescaled)["B_p"]
    check.record(f"{name}/mutation_ratio", abs(moved_b - truth["B_p"]) > 1e-6,
                 f"doubling the transport ratio moved B_p from "
                 f"{truth['B_p']:.5f} to {moved_b:.5f}")


def check_rf_cross_term(check: Check) -> dict[str, float]:
    """Corollary 7.1 for the R-Fusion learner (4.tex:359).

    With constant propensity and every other law held fixed, the corollary says
    that the cross term vanishes at e = 1/2.  Working the
    weight and the increment through gives more: chi^2 Delta_RF carries the factor
    E{(A-e)^3 | X} = e(1-e)(1-2e), so the cross term should be proportional to
    (1-2e)/[e(1-e)].  Testing the proportionality is a much sharper check than
    testing that one number is zero.
    """
    measured = {}
    for e in (0.5, 0.2, 0.6, 0.8):
        measured[e] = mo.covariance_route(fxs.rf_cross_term(e))["E_p"]
    check.record("RF/cross_term_vanishes_at_half", abs(measured[0.5]) <= ZERO_TOL,
                 f"cross term {measured[0.5]:.2e} at e = 1/2")

    constants = {e: measured[e] / ((1 - 2 * e) / (e * (1 - e)))
                 for e in measured if e != 0.5}
    spread = max(constants.values()) - min(constants.values())
    check.record("RF/cross_term_matches_predicted_factor", spread <= 1e-9,
                 f"cross term over (1-2e)/[e(1-e)] is "
                 + ", ".join(f"{v:.6f} at e={e}" for e, v in sorted(constants.items()))
                 + f"; spread {spread:.2e}")
    return {str(e): v for e, v in measured.items()}


def check_drf_cross_term(check: Check, name: str, truth: dict[str, float]) -> None:
    check.record(f"{name}/drf_cross_term_vanishes", abs(truth["E_p"]) <= ZERO_TOL,
                 f"cross term {truth['E_p']:.2e}, as Corollary 7.1 states for DRF")



def _parse(argv: list[str] | None):
    parser = argparse.ArgumentParser(description="stage options")
    parser.add_argument("--stage", default=None,
                        help="name of the stage directory under materials/confirm")
    parser.add_argument("--replace", action="store_true",
                        help="replace an existing stage that carries this exact "
                             "fingerprint; a stage written under a different "
                             "fingerprint is never replaced, because results from "
                             "different code states are not pooled")
    parser.add_argument("--preflight", "--dry-run", action="store_true", dest="preflight",
                        help="resolve protocol and collision state without calculations or writes")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse(argv)
    code_files = pv.code_files()
    fingerprint = pv.code_fingerprint(code_files)
    required_inputs = ["code/confirm/fixtures.py", "code/confirm/moments.py",
                       "code/confirm/m1_production_tie.py", "code/fusion_cate.py"]
    if args.preflight:
        payload = pv.preflight_payload(
            prefix=protocol.STAGE_ID, explicit_stage=args.stage, fingerprint=fingerprint,
            schedule_hash=protocol.schedule_hash(), required_inputs=required_inputs,
            required_outputs=protocol.REQUIRED_OUTPUTS)
        payload["scheduled_unit_count"] = len(protocol.scheduled_units())
        payload["expected_check_count"] = len(protocol.expected_check_names())
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    pv.require_code_fingerprint(fingerprint)
    stage_name = args.stage or pv.qualified_stage(protocol.STAGE_ID, fingerprint)
    ledger = protocol.SeedLedger()
    check = Check()
    summary = {}
    for name in fxs.QUADRANTS:
        fx = fxs.quadrant(name)
        check_projection(check, fx)
        truth = check_routes(check, fx)
        check_quadrant(check, name, truth, fx)
        check_drf_cross_term(check, name, truth)
        check_mutations(check, name, truth)
        lam, omega = mo.oracle_coordinates(truth)
        summary[name] = {**truth, "lambda_star": lam, "omega_star": omega,
                         "oracle_gain": mo.oracle_gain(truth),
                         "J_00": mo.direct_risk(fx, 0.0, 0.0)}
    rf = check_rf_cross_term(check)

    # M1's first Go condition compares the independent formula against the
    # production moments, not against a second copy of itself
    tied = tie.run(sampler_seed=ledger.seed("production_tie/sampler"),
                   bootstrap_seed=ledger.seed("production_tie/bootstrap"))
    check.record("production/moments_agree", tied["status"] == "PASS",
                 f"largest difference {tied['worst_gap']:.2e} against "
                 f"{tied['tolerance']:.0e} once the exploratory 1e-10 ridge is matched; "
                 f"at ridge zero the difference is "
                 f"{max(r['gap_at_rho0'] for r in tied['moments']):.2e}, which is the ridge")
    check.record("production/oracle_coordinates_interior", tied["oracle"]["interior"],
                 f"tie runs at lambda* {tied['oracle']['confirmatory'][0]:.4f} and "
                 f"omega* {tied['oracle']['confirmatory'][1]:.4f}, both strictly inside")
    check.record("production/drf_cross_term_consistent_with_zero",
                 abs(tied["cross_term"]["z"]) <= 3.0,
                 f"E_p {tied['cross_term']['E_p']:+.2e} against a Monte Carlo standard "
                 f"error of {tied['cross_term']['mcse']:.2e}, z = {tied['cross_term']['z']:+.2f}")

    for row in check.rows:
        print(f"[{row['status']}] {row['check']}: {row['detail']}")
    print()
    if check.failures:
        print(f"M1 STOP and REDESIGN: {len(check.failures)} of {len(check.rows)} checks failed")
        return 1
    actual_names = tuple(row["check"] for row in check.rows)
    expected_names = protocol.expected_check_names()
    if set(actual_names) != set(expected_names) or len(actual_names) != len(expected_names):
        raise RuntimeError("M1 check-name set differs from the frozen protocol")
    pv.require_code_fingerprint(fingerprint)

    record = record_text(len(check.rows))
    stage = pv.Stage(stage_name, fingerprint)
    with stage.scratch() as scratch:
        pv.write_json(scratch / "code_manifest.json", {
            "fingerprint": fingerprint, "files": pv.manifest_for(code_files, pv.CODE)})
        pv.write_json(scratch / "protocol.json", protocol.protocol_payload())
        pv.write_json(scratch / "seed_ledger.json", ledger.payload())
        pv.write_json(scratch / "environment.json", {**pv.environment(), "git_head": pv.git_head()})
        pv.write_json(scratch / "moments.json",
                      {"quadrants": summary, "rf_cross_term": rf, "production_tie": tied})
        pv.write_json(scratch / "checks.json", {"rows": check.rows,
                                                "passed": len(check.rows),
                                                "total": len(check.rows)})
        (scratch / "RECORD.md").write_text(record)
        pv.require_code_fingerprint(fingerprint)
        final = stage.promote(allow_replace=args.replace,
                              required_outputs=protocol.REQUIRED_OUTPUTS)
    print(f"M1 GO: {len(check.rows)}/{len(check.rows)} checks passed; "
          f"written to {final.relative_to(pv.PROJECT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
