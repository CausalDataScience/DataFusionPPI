"""Command line: CSV files in, JSON out.

    python -m dfppi ate  --trial trial.csv --obs obs.csv --trial-propensity 0.5
    python -m dfppi cate --trial trial.csv --obs obs.csv --trial-propensity 0.5 --learner DR \\
                         --predict points.csv --out cate.csv

Each CSV has a header.  The treatment column (default A) holds 0 or 1, the outcome column (default Y) is numeric,
and the covariates are --covariates, or else every other column.  --no-shift uses r_0 = 1 (the OBS covariates have
the trial law).  The result is printed as JSON; on an error the JSON is {"status": "error", "error": ...} and the
exit code is 1.
"""
import argparse
import csv
import json
import sys

import numpy as np

from .api import fuse_ate, fuse_cate


def read_csv(path, covariates=None, treatment=None, outcome=None):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"{path} has no rows")
    columns = list(rows[0])
    if covariates is None:
        covariates = [c for c in columns if c not in (treatment, outcome)]
    missing = [c for c in covariates + [c for c in (treatment, outcome) if c] if c not in columns]
    if missing:
        raise ValueError(f"{path} lacks the columns {missing}")
    x = np.array([[float(r[c]) for c in covariates] for r in rows])
    if treatment is None:
        return {"x": x}, covariates
    return {"x": x, "a": np.array([float(r[treatment]) for r in rows]),
            "y": np.array([float(r[outcome]) for r in rows])}, covariates


def main(argv=None):
    parser = argparse.ArgumentParser(prog="dfppi", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("estimand", choices=("ate", "cate"))
    parser.add_argument("--trial", required=True)
    parser.add_argument("--obs", required=True)
    parser.add_argument("--trial-propensity", type=float, required=True,
                        help="the known probability of treatment in the trial")
    parser.add_argument("--treatment", default="A")
    parser.add_argument("--outcome", default="Y")
    parser.add_argument("--covariates", help="comma-separated covariate columns (default: every other column)")
    parser.add_argument("--no-shift", action="store_true", help="r_0 = 1: the OBS covariates have the trial law")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--level", type=float, default=0.95)
    parser.add_argument("--learner", choices=("DR", "R"), default="DR")
    parser.add_argument("--predict", help="CSV of covariate rows at which to evaluate the CATE (default: the trial)")
    parser.add_argument("--out", help="CSV for the CATE predictions (default: in the JSON, up to 10,000 rows)")
    a = parser.parse_args(argv)
    try:
        covariates = a.covariates.split(",") if a.covariates else None
        trial, covariates = read_csv(a.trial, covariates, a.treatment, a.outcome)
        obs, _ = read_csv(a.obs, covariates, a.treatment, a.outcome)
        if a.estimand == "ate":
            result = fuse_ate(trial, obs, a.trial_propensity, not a.no_shift, a.seed, level=a.level)
        else:
            fit = fuse_cate(trial, obs, a.trial_propensity, a.learner, not a.no_shift, a.seed)
            points = read_csv(a.predict, covariates)[0]["x"] if a.predict else trial["x"]
            tau, tau_0 = fit.predict(points), fit.predict_trial_only(points)
            result = {**fit.summary(), "covariates": covariates, "n_predictions": len(tau)}
            if a.out:
                with open(a.out, "w", newline="") as handle:
                    writer = csv.writer(handle)
                    writer.writerow(["cate", "cate_trial_only"])
                    writer.writerows(zip(tau, tau_0))
                result["predictions_file"] = a.out
            elif len(tau) <= 10_000:
                result["predictions"] = {"cate": tau.tolist(), "cate_trial_only": tau_0.tolist()}
            else:
                result["warnings"].append("more than 10,000 predictions: pass --out to write them")
        print(json.dumps({"status": "ok", **result}, indent=2))
        return 0
    except (ValueError, RuntimeError, OSError) as error:
        print(json.dumps({"status": "error", "error": str(error)}, indent=2))
        return 1


if __name__ == "__main__":
    sys.exit(main())
