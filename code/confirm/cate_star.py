"""Frozen Tennessee STAR input and source-law construction for CATE-2.

The scientific runner never downloads data.  A caller supplies ``--star-csv``
or ``DATAFUSIONPPI_STAR_CSV``; the host cache is only a final candidate and is
accepted only after the canonical hash and schema checks pass.
"""
from __future__ import annotations

import os
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.preprocessing import OneHotEncoder


STAR_URL = "https://vincentarelbundock.github.io/Rdatasets/csv/AER/STAR.csv"
STAR_SHA256 = "0e8b179ea3d883730b25008ca1293c4c8f886aa8d5d299c03ab6ff05d0123ae6"
STAR_CACHE_CANDIDATE = Path("/private/tmp/claude-502/datafusionppi-data/STAR.csv")
STAR_ROWS = 11_598
STAR_REQUIRED = ("rownames", "stark", "gender", "ethnicity", "lunchk",
                 "schoolk", "schoolidk", "mathk", "readk")


def resolve_star_csv(cli_path: str | None) -> Path:
    """Resolve an existing local input without downloading or copying it."""
    candidate = cli_path or os.environ.get("DATAFUSIONPPI_STAR_CSV")
    path = Path(candidate).expanduser() if candidate else STAR_CACHE_CANDIDATE
    return path.resolve()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_star_csv(path: Path) -> tuple[pd.DataFrame, dict]:
    if not path.is_file():
        raise FileNotFoundError(f"resolved STAR input does not exist: {path}")
    digest = _sha256(path)
    if digest != STAR_SHA256:
        raise RuntimeError(f"STAR input hash mismatch: {digest}")
    frame = pd.read_csv(path)
    missing = sorted(set(STAR_REQUIRED) - set(frame.columns))
    if len(frame) != STAR_ROWS or missing:
        raise RuntimeError(f"STAR schema mismatch: rows={len(frame)}, missing={missing}")
    return frame, {"resolved_path": str(path), "source_url": STAR_URL,
                   "source_sha256": digest, "source_rows": len(frame)}


def _hash_payload(value) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"),
                     ensure_ascii=True).encode()
    return hashlib.sha256(raw).hexdigest()


@dataclass(frozen=True)
class FrozenStarSource:
    frame: pd.DataFrame
    x: np.ndarray
    learner_columns: tuple[str, ...]
    cell_keys: tuple[tuple, ...]
    cell_probability: np.ndarray
    propensity: np.ndarray
    members: tuple[tuple[np.ndarray, np.ndarray], ...]
    source_weight: np.ndarray
    outcome: np.ndarray
    oof_training_excludes_row: np.ndarray
    metadata: dict

    def source_marginal_error(self) -> float:
        # Reconstruct P_O^X from the arm and within-cell sampling weights rather
        # than comparing the frozen vector with itself.
        po = np.empty_like(self.cell_probability)
        for k, (control, treated) in enumerate(self.members):
            control_mass = float(np.sum(self.source_weight[control]
                                        / self.source_weight[control].sum()))
            treated_mass = float(np.sum(self.source_weight[treated]
                                        / self.source_weight[treated].sum()))
            arm_mass = ((1.0 - self.propensity[k]) * control_mass
                        + self.propensity[k] * treated_mass)
            po[k] = self.cell_probability[k] * arm_mass
        return float(np.max(np.abs(po - self.cell_probability)))

    def propensity_identity_error(self) -> float:
        empirical = np.array([
            self.frame.iloc[np.r_[pair[0], pair[1]]]["A"].mean()
            for pair in self.members
        ])
        return float(np.max(np.abs(self.propensity - empirical)))

    def draw(self, source: str, n: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
        if source not in {"RCT", "OBS"}:
            raise ValueError(f"unknown STAR source {source}")
        chosen = rng.choice(len(self.cell_keys), size=n, p=self.cell_probability)
        arms = (rng.random(n) < self.propensity[chosen]).astype(int)
        rows = np.empty(n, dtype=int)
        for j, (cell, arm) in enumerate(zip(chosen, arms)):
            pool = self.members[cell][arm]
            probability = None
            if source == "OBS":
                weight = self.source_weight[pool]
                probability = weight / weight.sum()
            rows[j] = int(rng.choice(pool, p=probability))
        return {"x": self.x[rows], "a": arms.astype(float),
                "y": self.outcome[rows], "propensity": self.propensity[chosen],
                "ratio": np.ones(n),
                "row_id": self.frame.iloc[rows]["rownames"].to_numpy(),
                "x_index": chosen}


def freeze_star_source(path: Path, *, outcome: str, alpha: float,
                       fold_seed: int) -> FrozenStarSource:
    """Build the empirical source law once, before any performance unit."""
    if outcome not in {"mathk", "readk"} or not 0.0 <= alpha < 1.0:
        raise ValueError("STAR outcome/alpha is outside the frozen protocol")
    raw, metadata = validate_star_csv(path)
    keep = (raw["stark"].isin(["small", "regular", "regular+aide"])
            & raw["mathk"].notna() & raw["readk"].notna())
    frame = raw.loc[keep].copy()
    frame["A"] = (frame["stark"] == "small").astype(int)
    raw_columns = ("gender", "ethnicity", "lunchk", "schoolk", "schoolidk")
    missing_columns = []
    imputation = {}
    for column in raw_columns:
        missing = f"{column}__missing"
        frame[missing] = frame[column].isna().astype(int)
        missing_columns.append(missing)
        mode = frame[column].mode(dropna=True)
        if mode.empty:
            raise RuntimeError(f"STAR column has no observed mode: {column}")
        imputation[column] = str(mode.iloc[0])
        frame[column] = frame[column].fillna(mode.iloc[0]).astype(str)
    frame["ethnicity3"] = np.where(
        frame["ethnicity"].eq("afam"), "afam",
        np.where(frame["ethnicity"].eq("cauc"), "cauc", "other"))
    learner_columns = ("gender", "ethnicity3", "lunchk", "schoolk", "schoolidk",
                       *missing_columns)
    frame["x_key"] = list(map(tuple, frame[list(learner_columns)].to_numpy(object)))
    share = frame.groupby("x_key", sort=True)["A"].mean()
    eligible = set(share[(share >= 0.15) & (share <= 0.85)].index)
    frame = frame[frame["x_key"].isin(eligible)].copy()
    frame = frame.sort_values("rownames", kind="stable").reset_index(drop=True)
    if frame.empty:
        raise RuntimeError("STAR exact-X eligibility produced an empty cohort")

    encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    x = encoder.fit_transform(frame[list(learner_columns)].astype(str))
    a = frame["A"].to_numpy(int)
    y_raw = frame[outcome].to_numpy(float)
    fold_rng = np.random.default_rng(fold_seed)
    fold = np.empty(len(frame), dtype=int)
    fold[fold_rng.permutation(len(frame))] = np.arange(len(frame)) % 5
    prediction = np.empty(len(frame))
    design = np.column_stack([x, a])
    excluded = np.ones(len(frame), dtype=bool)
    for k in range(5):
        train, test = fold != k, fold == k
        model = Ridge(alpha=5.0, solver="lsqr", tol=1e-10).fit(design[train], y_raw[train])
        prediction[test] = model.predict(design[test])
        excluded[test] = ~train[test]
    residual = y_raw - prediction
    rank = np.empty(len(frame))
    row_id = frame["rownames"].to_numpy()
    for arm in (0, 1):
        idx = np.flatnonzero(a == arm)
        order = np.lexsort((row_id[idx], residual[idx]))
        ordinal = np.empty(len(idx), dtype=int)
        ordinal[order] = np.arange(len(idx))
        rank[idx] = (ordinal + 0.5) / len(idx)
    source_weight = np.where(a == 1, 1.0 + alpha * (2.0 * rank - 1.0),
                             1.0 - alpha * (2.0 * rank - 1.0))
    mean, sd = float(y_raw.mean()), float(y_raw.std(ddof=0))
    outcome_values = (y_raw - mean) / sd

    grouped = frame.groupby("x_key", sort=True).indices
    keys = tuple(sorted(grouped, key=repr))
    counts = np.array([len(grouped[key]) for key in keys], dtype=float)
    probability = counts / counts.sum()
    propensity = np.array([a[grouped[key]].mean() for key in keys])
    members = tuple((np.asarray(grouped[key])[a[grouped[key]] == 0],
                     np.asarray(grouped[key])[a[grouped[key]] == 1]) for key in keys)
    ordered_ids = [str(value) for value in frame["rownames"]]
    probability_payload = [[repr(k), format(float(p), ".17g")]
                           for k, p in zip(keys, probability)]
    propensity_payload = [[repr(k), format(float(p), ".17g")]
                          for k, p in zip(keys, propensity)]
    metadata.update({
        "outcome": outcome, "alpha": alpha, "eligible_rows": len(frame),
        "exact_cells": len(keys), "learner_columns": list(learner_columns),
        "missingness_indicators": missing_columns, "imputation_modes": imputation,
        "ethnicity3_mapping": "afam->afam;cauc->cauc;other observed->other",
        "oof": "five-fold deterministic Ridge(alpha=5) on one-hot X plus A",
        "outcome_standardization": {"mean": mean, "sd": sd},
        "ordered_row_id_sha256": _hash_payload(ordered_ids),
        "fold_assignment_sha256": hashlib.sha256(
            np.ascontiguousarray(fold, dtype="<i8").tobytes()).hexdigest(),
        "probability_table_sha256": _hash_payload(probability_payload),
        "propensity_table_sha256": _hash_payload(propensity_payload),
    })
    metadata["frozen_support_sha256"] = _hash_payload({
        key: metadata[key] for key in ("ordered_row_id_sha256",
                                       "fold_assignment_sha256",
                                       "probability_table_sha256",
                                       "propensity_table_sha256")})
    return FrozenStarSource(frame=frame, x=x, learner_columns=learner_columns,
                            cell_keys=keys, cell_probability=probability,
                            propensity=propensity, members=members,
                            source_weight=source_weight, outcome=outcome_values,
                            oof_training_excludes_row=excluded, metadata=metadata)
