"""M0: evidence freeze, seed ledger and atomic stage promotion.

Three jobs.

1. Hash every artifact that already exists, so a later run cannot silently
   overwrite or contradict it.
2. Fingerprint the code and the data-generating process, so a result computed
   under one fingerprint is never pooled with a result computed under another.
3. Promote a finished stage directory atomically, so a partial run never
   appears as a complete one.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import secrets
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
CODE = PROJECT / "code"
MATERIALS = PROJECT / "materials"
CONFIRM = MATERIALS / "confirm"
CHUNK = 1 << 20


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(CHUNK)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def manifest_for(paths: list[Path], root: Path) -> dict[str, dict]:
    """Path relative to ``root`` mapped to its digest, size and modification time."""
    out = {}
    for path in sorted(paths):
        if not path.is_file():
            continue
        stat = path.stat()
        out[str(path.relative_to(root))] = {
            "sha256": sha256_file(path),
            "bytes": stat.st_size,
            "mtime": int(stat.st_mtime),
        }
    return out


def code_fingerprint(files: list[Path]) -> str:
    """One digest over the sorted (relative path, content digest) pairs.

    Two runs may be pooled only when this matches.  It deliberately covers the
    file names as well as the contents, so adding or removing a module changes
    the fingerprint even when no existing file is edited.
    """
    parts = [f"{p.relative_to(CODE)}:{sha256_file(p)}" for p in sorted(files) if p.is_file()]
    return sha256_text("\n".join(parts))


def environment() -> dict[str, str]:
    import numpy
    import pandas
    import sklearn
    versions = {"python": sys.version.split()[0], "numpy": numpy.__version__,
                "pandas": pandas.__version__, "sklearn": sklearn.__version__,
                "platform": platform.platform()}
    try:
        import torch
        versions["torch"] = torch.__version__
    except Exception:
        versions["torch"] = "absent"
    return versions


def git_head() -> str:
    try:
        out = subprocess.run(["git", "-C", str(PROJECT), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=20)
        return out.stdout.strip() or "unavailable"
    except Exception:
        return "unavailable"


# --------------------------------------------------- atomic stage promotion
class StagePromotionError(RuntimeError):
    pass


@dataclass
class Stage:
    """A stage writes into a scratch directory and is promoted in one rename.

    ``fingerprint`` is stored beside the payload.  ``promote`` refuses to replace
    an existing stage whose fingerprint differs, which is the prefix
    incompatibility check the plan asks for: an M3 prefix may be reused by M5
    only when protocol, code path and seeds all match.
    """
    name: str
    fingerprint: str
    root: Path = CONFIRM
    _scratch: Path | None = field(default=None, init=False, repr=False)
    _owner: str | None = field(default=None, init=False, repr=False)
    _inode: int | None = field(default=None, init=False, repr=False)

    MARKER = ".stage-owner.json"

    @property
    def final(self) -> Path:
        return self.root / self.name

    def open(self) -> Path:
        if self._scratch is not None:
            raise StagePromotionError("this Stage instance already owns an open scratch directory")
        self.root.mkdir(parents=True, exist_ok=True)
        self._scratch = Path(tempfile.mkdtemp(prefix=f".{self.name}.staging.", dir=self.root))
        self._owner = secrets.token_hex(16)
        self._inode = self._scratch.stat().st_ino
        write_json(self._scratch / self.MARKER,
                   {"owner": self._owner, "inode": self._inode,
                    "name": self.name, "parent": str(self.root.resolve())})
        return self._scratch

    def _validate_owned_scratch(self) -> Path:
        if self._scratch is None or self._owner is None or self._inode is None:
            raise StagePromotionError("this Stage instance has no open scratch directory")
        scratch = self._scratch
        if scratch.parent.resolve() != self.root.resolve():
            raise StagePromotionError("scratch parent changed")
        if not scratch.name.startswith(f".{self.name}.staging."):
            raise StagePromotionError("scratch prefix changed")
        if not scratch.is_dir() or scratch.stat().st_ino != self._inode:
            raise StagePromotionError("scratch inode changed")
        marker = scratch / self.MARKER
        try:
            payload = json.loads(marker.read_text())
        except Exception as exc:
            raise StagePromotionError("scratch ownership marker is missing or invalid") from exc
        expected = {"owner": self._owner, "inode": self._inode,
                    "name": self.name, "parent": str(self.root.resolve())}
        if payload != expected:
            raise StagePromotionError("scratch ownership marker does not match this Stage instance")
        return scratch

    def cleanup(self) -> None:
        scratch = self._validate_owned_scratch()
        shutil.rmtree(scratch)
        self._scratch = None
        self._owner = None
        self._inode = None

    @contextmanager
    def scratch(self):
        path = self.open()
        try:
            yield path
        finally:
            if self._scratch is not None:
                self.cleanup()

    def promote(self, *, allow_replace: bool = False,
                required_outputs: tuple[str, ...] | None = None) -> Path:
        scratch = self._validate_owned_scratch()
        backup: Path | None = None
        try:
            (scratch / "FINGERPRINT").write_text(self.fingerprint + "\n")
            if required_outputs is not None:
                found = {p.name for p in scratch.iterdir() if p.name != self.MARKER}
                if found != set(required_outputs):
                    raise StagePromotionError(
                        f"required output mismatch: missing={sorted(set(required_outputs)-found)}, "
                        f"extra={sorted(found-set(required_outputs))}")
            if self.final.exists():
                existing = self.final / "FINGERPRINT"
                previous = existing.read_text().strip() if existing.exists() else "missing"
                if previous != self.fingerprint:
                    raise StagePromotionError(
                        f"stage {self.name} exists under fingerprint {previous[:12]} and this run "
                        f"has {self.fingerprint[:12]}; use a fingerprint-qualified stage")
                if not allow_replace:
                    raise StagePromotionError(
                        f"stage {self.name} already exists with the same fingerprint; "
                        "pass allow_replace=True only to redo an identical run")
                backup = self.final.with_suffix(".superseded")
                if backup.exists():
                    raise StagePromotionError(f"{backup.name} already exists; remove it by hand")
                os.rename(self.final, backup)
            (scratch / self.MARKER).unlink()
            try:
                os.rename(scratch, self.final)
            except Exception:
                if backup is not None and backup.exists() and not self.final.exists():
                    os.rename(backup, self.final)
                raise
            self._scratch = None
            self._owner = None
            self._inode = None
            return self.final
        except Exception:
            if backup is not None and backup.exists() and not self.final.exists():
                os.rename(backup, self.final)
            if self._scratch is not None:
                # The marker is removed only immediately before the final rename.
                # Recreate it if that rename failed, then delete only this instance's
                # exact scratch directory after the full ownership check.
                marker = self._scratch / self.MARKER
                if not marker.exists() and self._scratch.exists():
                    write_json(marker, {"owner": self._owner, "inode": self._inode,
                                        "name": self.name,
                                        "parent": str(self.root.resolve())})
                self.cleanup()
            raise


def code_files() -> list[Path]:
    return sorted(CODE.glob("*.py")) + sorted((CODE / "confirm").glob("*.py"))


def require_code_fingerprint(expected: str) -> None:
    current = code_fingerprint(code_files())
    if current != expected:
        raise RuntimeError(f"source fingerprint drift: expected {expected}, found {current}")


def qualified_stage(prefix: str, fingerprint: str) -> str:
    return f"{prefix}-{fingerprint[:12]}"


def collision_state(path: Path, fingerprint: str) -> str:
    if not path.exists():
        return "absent"
    marker = path / "FINGERPRINT"
    if not path.is_dir() or not marker.is_file():
        return "incomplete"
    return "same_fingerprint" if marker.read_text().strip() == fingerprint else "different_fingerprint"


def preflight_payload(*, prefix: str, explicit_stage: str | None, fingerprint: str,
                      schedule_hash: str, required_inputs: list[str],
                      required_outputs: tuple[str, ...], root: Path = CONFIRM) -> dict[str, object]:
    stage_name = explicit_stage or qualified_stage(prefix, fingerprint)
    target = root / stage_name
    state = collision_state(target, fingerprint)
    inputs = [{"path": item, "exists": (PROJECT / item).exists()} for item in required_inputs]
    missing = [item["path"] for item in inputs if not item["exists"]]
    status = "ready" if not missing and state == "absent" else (
        "blocked_missing_inputs" if missing else f"collision_{state}")
    return {"mode": "preflight", "calculation_performed": False,
            "write_performed": False, "resolved_target": str(target),
            "fingerprint": fingerprint, "schedule_hash": schedule_hash,
            "collision_state": state, "fingerprint_qualified": explicit_stage is None,
            "qualification_status": status,
            "required_inputs": inputs, "required_outputs": list(required_outputs)}


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")
