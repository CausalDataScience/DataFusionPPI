"""Frozen protocol for the confirmatory CATE stages.

Everything a run needs to know before it starts lives here: the cells, the
estimator set, the ratio routes, the role fractions, the schedule for each
milestone, and the seed ledger.  A run reads this module and never invents a
unit, a seed or a cell of its own.

The data-generating process follows audit plan section 7.1.  Two knobs move the
two channels' opportunities, and M1 showed why each one works.

``confounding`` moves the lambda opportunity.  C_p is the whitened covariance of
the base pseudo-outcome with the increment, and the increment helps exactly when
the observational outcome regression predicts the trial outcome surface better
than the trial regression does.  Confounding is what spoils that, so a small
confounding coefficient opens the lambda channel and a large one closes it.

``roughness`` moves the omega opportunity.  D_p is the whitened covariance of
b(tau - zeta_p) with b(ghat - zeta_p) and is exactly zero when the basis spans
tau.  Putting a component of tau outside the spline span opens the omega channel.

Quadrant labels are names, not claims.  Every run measures A_p, B_p, C_p and D_p
with the adaptive oracle and records which quadrant the cell actually lands in.
"""
from __future__ import annotations

import hashlib
import json
from statistics import NormalDist
from dataclasses import dataclass, asdict

PROTOCOL_ID = "DataFusionPPI-CATE-confirm-v1"
BASE_SEED = 20260912

ROLES = ("nuisance", "tuning", "selection", "reporting")

ESTIMATORS = ("rct_only", "lambda_only", "omega_only", "joint", "grid_oracle")
RATIO_ROUTES = ("oracle", "classifier", "balanced_x", "balanced_x_ghat", "unit")
LEARNERS = ("DRF", "RF")

LAMBDA_GRID = (0.0, 0.25, 0.5, 0.75, 1.0)
OMEGA_GRID = (0.0, 0.25, 0.5, 0.75, 1.0)
RIDGE_GRID = (0.0,)            # P0-2: the Theorem 7 path runs unclipped at rho = 0

# adaptive oracle, plan P0-3
ORACLE_DRAWS = (20_000, 40_000, 80_000, 100_000)
TARGET_TRUTH_DRAWS = 100_000
ORACLE_JACKKNIFE_BLOCKS = 20
ORACLE_SIMULTANEOUS_ALPHA = 0.05
ORACLE_SIMULTANEOUS_M = 6
ORACLE_SIMULTANEOUS_Z = NormalDist().inv_cdf(
    1.0 - ORACLE_SIMULTANEOUS_ALPHA / (2.0 * ORACLE_SIMULTANEOUS_M))
MOMENT_VARIANCE_MCSE_FRACTION = 0.01
MOMENT_COVARIANCE_MCSE_FRACTION = 0.005
COORDINATE_RISK_MCSE_FRACTION = 0.01


@dataclass(frozen=True)
class Cell:
    """One design point.  ``quadrant`` is the intended label, not a claim."""
    cell_id: str
    family: str                # "bounded" or "gaussian_shift" or "star"
    quadrant: str
    n_rct: int
    n_obs: int
    confounding: float         # c0 in plan 7.1, the constant part
    shape_confounding: float   # c1 in plan 7.1, the X-dependent part
    roughness: float           # amplitude of the component outside the spline span
    obs_baseline_shift: float  # mismatch between the two sources' baseline surfaces
    sieve: str
    learner: str
    ratio_route: str
    shift: str = "none"        # "none", "mild", "strong"
    dimension: int = 4
    trial_nuisance_boost: float = 1.0  # enlarges only the trial nuisance draw

    def key(self) -> str:
        return self.cell_id


def _cell(cell_id, quadrant, **kw) -> Cell:
    base = dict(family="bounded", n_rct=200, n_obs=5000, confounding=1.0,
                shape_confounding=0.0, roughness=0.0, obs_baseline_shift=0.0,
                trial_nuisance_boost=1.0, sieve="spline3", learner="DRF",
                ratio_route="unit")
    base.update(kw)
    return Cell(cell_id=cell_id, quadrant=quadrant, **base)


# --------------------------------------------------------------- M2 schedule
# Six independent stochastic source units, exactly as PRD section 9 fixes them.
# Estimators, candidates and routes multiply rows, never source units.
def m2_cells() -> tuple[Cell, ...]:
    # obs_baseline_shift closes the lambda channel, roughness opens the omega
    # channel; M2 measures which quadrant each cell actually lands in
    # trial_nuisance_boost closes the lambda channel by making the trial outcome
    # regression already good, which leaves ghat and so the omega channel alone;
    # roughness opens the omega channel.  M2 measures where each cell lands.
    quadrant_knobs = {
        "L0O0": dict(trial_nuisance_boost=60.0, roughness=0.0),
        "L1O0": dict(trial_nuisance_boost=1.0, roughness=0.0),
        "L0O1": dict(trial_nuisance_boost=60.0, roughness=0.45),
        "L1O1": dict(trial_nuisance_boost=1.0, roughness=0.45),
    }
    cells = [_cell(f"m2/{q}/spline3/DRF", q, **knobs) for q, knobs in quadrant_knobs.items()]
    cells.append(_cell("m2/shift_strong/spline3/DRF", "L1O1", family="gaussian_shift",
                       shift="strong", confounding=0.0, roughness=0.45,
                       ratio_route="all"))
    cells.append(Cell(cell_id="m2/star/mathk/a0.8", family="star", quadrant="N/A",
                      n_rct=200, n_obs=1200, confounding=0.0, shape_confounding=0.0,
                      roughness=0.0, obs_baseline_shift=0.0, sieve="spline3",
                      learner="DRF", ratio_route="unit"))
    return tuple(cells)


def role_sizes(cell: Cell) -> dict[str, int]:
    """Canonical per-role sizes; displayed CATE sizes are not split totals."""
    if cell.family == "star":
        return {"RCT/nuisance": 200, "RCT/tuning": 200,
                "RCT/selection": 1000, "RCT/reporting": 1000,
                "OBS/nuisance": 5000, "OBS/tuning": 5000,
                "OBS/selection": 1000, "OBS/reporting": 1000}
    return {f"RCT/{role}": cell.n_rct for role in ROLES} | {
        f"OBS/{role}": cell.n_obs for role in ROLES}


SCHEDULES = {"m2": m2_cells}


def scheduled_units(stage: str, replications: int | None = None) -> tuple[dict, ...]:
    cells = SCHEDULES[stage]()
    reps = {"m2": 1}[stage] if replications is None else replications
    units = []
    for cell in cells:
        for rep in range(reps):
            units.append({"stage": stage, "cell_id": cell.cell_id,
                          "replication_id": rep,
                          "seed_key": f"{PROTOCOL_ID}|{stage}|{cell.cell_id}|{rep}"})
    return tuple(units)


def seed32(seed_key: str) -> int:
    digest = hashlib.sha256(f"{BASE_SEED}|{seed_key}".encode()).hexdigest()
    return int(digest[:8], 16)


def entropy_hex(seed_key: str) -> str:
    return hashlib.sha256(f"{BASE_SEED}|{seed_key}|entropy".encode()).hexdigest()[:32]


def role_seed(seed_key: str, source: str, role: str) -> int:
    """A separate RNG namespace per source and role, plan P0-1."""
    key = f"{seed_key}|{source}|{role}"
    return int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)


def seed_ledger(stage: str, replications: int | None = None) -> dict:
    units = []
    for unit in scheduled_units(stage, replications):
        units.append({**unit, "seed32": seed32(unit["seed_key"]),
                      "entropy_hex": entropy_hex(unit["seed_key"]),
                      "role_seeds": {f"{s}/{r}": role_seed(unit["seed_key"], s, r)
                                     for s in ("RCT", "OBS") for r in ROLES}})
    return {"protocol_id": PROTOCOL_ID, "stage": stage, "base": BASE_SEED,
            "units": units, "schedule_hash": schedule_hash(stage, replications)}


def schedule_hash(stage: str, replications: int | None = None) -> str:
    payload = json.dumps([u["seed_key"] for u in scheduled_units(stage, replications)],
                         sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()


def protocol_document(stage: str, replications: int | None = None) -> dict:
    return {
        "protocol_id": PROTOCOL_ID, "stage": stage,
        "cells": [asdict(c) for c in SCHEDULES[stage]()],
        "roles": ROLES, "role_sizes": {c.cell_id: role_sizes(c)
                                         for c in SCHEDULES[stage]()},
        "estimators": ESTIMATORS, "ratio_routes": RATIO_ROUTES, "learners": LEARNERS,
        "lambda_grid": LAMBDA_GRID, "omega_grid": OMEGA_GRID, "ridge_grid": RIDGE_GRID,
        "oracle_draws": ORACLE_DRAWS, "target_truth_draws": TARGET_TRUTH_DRAWS,
        "oracle_jackknife_blocks": ORACLE_JACKKNIFE_BLOCKS,
        "oracle_simultaneous": {"alpha": ORACLE_SIMULTANEOUS_ALPHA,
                                "m": ORACLE_SIMULTANEOUS_M,
                                "z": ORACLE_SIMULTANEOUS_Z},
        "moment_mcse_fractions": {"A_p": MOMENT_VARIANCE_MCSE_FRACTION,
                                  "B_p": MOMENT_VARIANCE_MCSE_FRACTION,
                                  "C_p": MOMENT_COVARIANCE_MCSE_FRACTION,
                                  "D_p": MOMENT_COVARIANCE_MCSE_FRACTION},
        "coordinate_risk_mcse_fraction": COORDINATE_RISK_MCSE_FRACTION,
        "schedule_hash": schedule_hash(stage, replications),
        "scheduled_unit_count": len(scheduled_units(stage, replications)),
    }


REQUIRED_OUTPUTS = (
    "FINGERPRINT", "artifact_manifest.json", "code_manifest.json", "protocol.json",
    "seed_ledger.json", "schedule.json", "environment.json", "source_manifest.json",
    "checks.json", "failures.json", "runtime.json", "RECORD.md", "raw",
)
REQUIRED_RAW_KEYS = (
    "protocol_id", "code_fingerprint", "stage", "study", "family", "cell_id",
    "replication_id", "seed_key", "seed32", "role_fingerprints",
    "learner", "sieve", "ratio_route", "candidate_id", "lambda", "omega", "rho",
    "metric", "value", "status", "failure_class", "runtime_seconds",
)
