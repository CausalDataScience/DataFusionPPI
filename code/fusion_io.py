"""Output conventions shared by the four experiments (handoff Section 2)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[1]
MATERIALS = PROJECT_DIR / "materials"


def summarize(frame: pd.DataFrame, group: list[str], metrics: dict[str, str]) -> pd.DataFrame:
    rows = []
    for keys, block in frame.groupby(group, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        base = dict(zip(group, keys))
        for metric, column in metrics.items():
            if column not in block or block[column].isna().all():
                continue
            values = block[column].to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            if len(values) == 0:
                continue
            if metric == "rmse":
                point = float(np.sqrt(np.mean(values ** 2)))
                se = float(np.std(values ** 2, ddof=1) / (2 * point * np.sqrt(len(values)))) if point > 0 else 0.0
            else:
                point = float(np.mean(values))
                se = float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else 0.0
            rows.append({**base, "metric": metric, "value": point,
                         "monte_carlo_se": se, "replications": len(values)})
    return pd.DataFrame(rows)


def write_outputs(stem: str, replications: pd.DataFrame, summary: pd.DataFrame) -> dict[str, Path]:
    MATERIALS.mkdir(exist_ok=True)
    paths = {"replications": MATERIALS / f"{stem}_replications.csv",
             "summary": MATERIALS / f"{stem}_summary.csv"}
    replications.to_csv(paths["replications"], index=False)
    summary.to_csv(paths["summary"], index=False)
    return paths
