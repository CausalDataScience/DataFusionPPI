"""M0 driver: freeze the existing evidence and lock the protocol.

Writes ``materials/confirm/M0_freeze/`` containing the artifact manifest, the
code and document fingerprints, the seed ledger for M1, and the label map that
records which existing artifact is exploratory and which is confirmatory.

The two Go conditions of M0 are checked here rather than asserted in prose:
promotion is atomic, and a stage written under one fingerprint cannot be
replaced by a run carrying another.

Run: python3 -m confirm.m0_freeze
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

from confirm import provenance as pv
from confirm import m1_protocol as protocol

# Every existing CATE artifact is exploratory.  The audit judged the current
# CATE results FAIL against the manuscript-ready bar, so none of them may be
# promoted to a confirmatory claim without being regenerated under this plan.
EXPLORATORY_GLOBS = [
    "cate1_*.csv", "cate2_*.csv", "cate3_*.csv",
    "figure0_*.png", "figure2_*.png", "figure4_*.png", "figure5_*.png",
    "figure6_*.png", "figureA_*.png",
    "confounding_shape.csv",
]
# ATE artifacts are frozen as context.  They are not the subject of this plan.
CONTEXT_GLOBS = [
    "ate1_*.csv", "ate2_*.csv", "conditional_variance_*.csv",
    "figure1_*.png", "figure3_*.png", "paper_tables.csv", "judgment_axes.csv",
    "acceptance_tests.csv",
]
DOCUMENT_PATHS = [
    "materials/2026-09-12-cate-experiment-audit-and-confirmatory-plan.md",
    "materials/2026-09-09-implementation-handoff-JA.md",
    "materials/2026-09-12-implementation-monitoring-report.md",
    "materials/2026-09-12-traceability-map.md",
    "manuscript/main/4.tex",
    "manuscript/main/2.tex",
    "manuscript/main/3.tex",
]
QUADRANTS = ("L0O0", "L1O0", "L0O1", "L1O1")


def collect(globs: list[str]) -> list[Path]:
    found: list[Path] = []
    for pattern in globs:
        found.extend(sorted(pv.MATERIALS.glob(pattern)))
    return found


def self_test_promotion() -> list[str]:
    """M0 Go conditions, executed rather than asserted."""
    notes = []
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        first = pv.Stage("probe", "fingerprint-A", root=root)
        with first.scratch() as scratch:
            (scratch / "payload.txt").write_text("first run\n")
            partial_visible = (root / "probe").exists()
            first.promote()
        notes.append(f"a stage is invisible until it is promoted: {not partial_visible}")
        notes.append(f"promotion produced {(root / 'probe' / 'payload.txt').exists()}")

        second = pv.Stage("probe", "fingerprint-B", root=root)
        try:
            with second.scratch():
                second.promote()
            notes.append("FAILED: a different fingerprint was allowed to replace the stage")
        except pv.StagePromotionError as exc:
            notes.append(f"a different fingerprint is refused: {'fingerprint' in str(exc)}")

        third = pv.Stage("probe", "fingerprint-A", root=root)
        try:
            with third.scratch():
                third.promote()
            notes.append("FAILED: an identical fingerprint replaced the stage without consent")
        except pv.StagePromotionError as exc:
            notes.append(f"an identical fingerprint still needs consent: {'allow_replace' in str(exc)}")
        kept = (root / "probe" / "payload.txt").read_text().strip()
        notes.append(f"the original payload survived both attempts: {kept == 'first run'}")
    return notes



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
                        help="resolve provenance and collision state without calculations or writes")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse(argv)
    code_files = pv.code_files()
    fingerprint = pv.code_fingerprint(code_files)
    required_outputs = ("FINGERPRINT", "artifact_manifest.json", "code_manifest.json",
                        "document_manifest.json", "seed_ledger.json", "environment.json",
                        "label_map.json", "promotion_self_test.txt")
    if args.preflight:
        payload = pv.preflight_payload(
            prefix="M0_freeze", explicit_stage=args.stage, fingerprint=fingerprint,
            schedule_hash=protocol.schedule_hash(), required_inputs=DOCUMENT_PATHS,
            required_outputs=required_outputs)
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    stage_name = args.stage or pv.qualified_stage("M0_freeze", fingerprint)

    exploratory = collect(EXPLORATORY_GLOBS)
    context = collect(CONTEXT_GLOBS)
    documents = [pv.PROJECT / p for p in DOCUMENT_PATHS]

    ledger = protocol.SeedLedger()
    schedule = ledger.payload()

    notes = self_test_promotion()
    failures = [n for n in notes if n.startswith("FAILED") or n.endswith("False")]

    pv.require_code_fingerprint(fingerprint)
    stage = pv.Stage(stage_name, fingerprint)
    with stage.scratch() as scratch:
        pv.write_json(scratch / "artifact_manifest.json", {
            "exploratory": pv.manifest_for(exploratory, pv.MATERIALS),
            "context": pv.manifest_for(context, pv.MATERIALS),
        })
        pv.write_json(scratch / "code_manifest.json", {
            "fingerprint": fingerprint, "files": pv.manifest_for(code_files, pv.CODE)})
        pv.write_json(scratch / "document_manifest.json", pv.manifest_for(documents, pv.PROJECT))
        pv.write_json(scratch / "seed_ledger.json", schedule)
        pv.write_json(scratch / "environment.json", {**pv.environment(), "git_head": pv.git_head()})
        pv.write_json(scratch / "label_map.json", {
            "rule": "every CATE artifact that existed before this plan is exploratory",
            "exploratory": sorted(p.name for p in exploratory),
            "context_not_in_scope": sorted(p.name for p in context), "confirmatory": []})
        (scratch / "promotion_self_test.txt").write_text("\n".join(notes) + "\n")
        if failures:
            raise RuntimeError("M0 atomic promotion self-test failed: " + "; ".join(failures))
        pv.require_code_fingerprint(fingerprint)
        final = stage.promote(allow_replace=args.replace, required_outputs=required_outputs)
    print(f"M0 freeze written to {final.relative_to(pv.PROJECT)}")
    print(f"  code fingerprint   {fingerprint}")
    print(f"  exploratory files  {len(exploratory)}")
    print(f"  context files      {len(context)}")
    print(f"  documents          {sum(1 for d in documents if d.is_file())} of {len(documents)}")
    print(f"  M1 seed units      {len(schedule['units'])}")
    for line in notes:
        print("  self-test: " + line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
