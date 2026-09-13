"""Frozen M1 protocol, schedule, seeds, checks, and required artifacts."""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

import numpy as np

PROTOCOL_ID = "DataFusionPPI-CATE-confirm-M1-v2"
STAGE_ID = "M1_fixtures"
BASE_SEED = 20260912
QUADRANTS = ("L0O0", "L1O0", "L0O1", "L1O1")
MUTATIONS = ("sign", "role_substitution", "missing_projection", "ratio")
RF_PROPENSITIES = ("0.2", "0.5", "0.6", "0.8")
REQUIRED_OUTPUTS = (
    "FINGERPRINT", "code_manifest.json", "protocol.json", "seed_ledger.json",
    "environment.json", "checks.json", "moments.json", "RECORD.md",
)


def scheduled_units() -> tuple[str, ...]:
    units = [f"quadrant/{q}/base" for q in QUADRANTS]
    units += [f"quadrant/{q}/mutation/{m}" for q in QUADRANTS for m in MUTATIONS]
    units += [f"rf_cross_term/e/{e}" for e in RF_PROPENSITIES]
    units += ["production_tie/sampler", "production_tie/bootstrap"]
    assert len(units) == 26 and len(set(units)) == 26
    return tuple(units)


def expected_check_names() -> tuple[str, ...]:
    names: list[str] = []
    for q in QUADRANTS:
        names += [
            f"{q}/normal_equation", f"{q}/two_routes_agree",
            f"{q}/quadratic_identity",
            f"{q}/{'lambda_opportunity' if q[1] == '1' else 'no_lambda_opportunity'}",
            f"{q}/{'omega_opportunity' if q[3] == '1' else 'no_omega_opportunity'}",
            f"{q}/oracle_minimises_risk", f"{q}/fixed_coordinate_ordering",
            f"{q}/drf_cross_term_vanishes",
            f"{q}/{'mutation_sign_detected' if q[1] == '1' else 'sign_flip_null_invariance'}",
            f"{q}/mutation_role_substitution", f"{q}/mutation_missing_projection",
            f"{q}/mutation_ratio",
        ]
    names += ["RF/cross_term_vanishes_at_half", "RF/cross_term_matches_predicted_factor",
              "production/moments_agree", "production/oracle_coordinates_interior",
              "production/drf_cross_term_consistent_with_zero"]
    assert len(names) == 53 and len(set(names)) == 53
    return tuple(names)


def schedule_hash() -> str:
    body = json.dumps({"protocol_id": PROTOCOL_ID, "stage": STAGE_ID,
                       "units": scheduled_units()}, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(body.encode()).hexdigest()


@dataclass(frozen=True)
class SeedLedger:
    """Prospective ledger whose key contains protocol, stage, base, and unit."""
    protocol_id: str = PROTOCOL_ID
    stage: str = STAGE_ID
    base: int = BASE_SEED

    def row(self, unit: str) -> dict[str, object]:
        if unit not in scheduled_units():
            raise ValueError(f"unscheduled M1 unit: {unit}")
        key = f"{self.protocol_id}|{self.stage}|{self.base}|{unit}"
        entropy = hashlib.sha256(key.encode()).digest()[:16]
        words = np.frombuffer(entropy, dtype=">u4").astype(np.uint32)
        seed32 = int(np.random.SeedSequence(words.tolist()).generate_state(1)[0])
        return {"protocol_id": self.protocol_id, "stage": self.stage,
                "base": self.base, "unit": unit,
                "entropy_hex": entropy.hex(), "seed32": seed32}

    def payload(self) -> dict[str, object]:
        rows = [self.row(unit) for unit in scheduled_units()]
        if len({r["unit"] for r in rows}) != 26:
            raise RuntimeError("M1 schedule is not exact")
        if len({r["entropy_hex"] for r in rows}) != 26:
            raise RuntimeError("M1 entropy collision")
        if len({r["seed32"] for r in rows}) != 26:
            raise RuntimeError("M1 seed32 collision")
        return {"protocol_id": self.protocol_id, "stage": self.stage,
                "base": self.base, "schedule_hash": schedule_hash(), "units": rows}

    def seed(self, unit: str) -> int:
        return int(self.row(unit)["seed32"])


def protocol_payload() -> dict[str, object]:
    return {"protocol_id": PROTOCOL_ID, "stage": STAGE_ID,
            "schedule_hash": schedule_hash(), "scheduled_units": list(scheduled_units()),
            "expected_checks": list(expected_check_names()),
            "required_outputs": list(REQUIRED_OUTPUTS),
            "theory_scope": {
                "omega_span": "DRF only, exact score, fixed basis, exact transport, population first-order",
                "rf_proportionality": "constant propensity with all other laws fixed",
                "L0_sign_flip": "null invariance, not fault detection",
            }}
