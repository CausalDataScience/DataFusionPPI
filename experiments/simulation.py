"""Runs the simulation of Sec. 5 and saves the results, one CSV file per cell (source, trial size).

    python3 simulation.py                                   every source of experiment.SOURCES, 300 replications;
                                                            this includes scm3 and lung, which the paper does not
                                                            use (lung needs optional files).  README.md lists the
                                                            commands of the paper (1,000 replications).
    python3 simulation.py --sources scm1 ihdp acic --reps 1000 --workers 16
    python3 simulation.py --sources scm1 --m 100 --shifts strong_shift --reps 20       one small cell
    python3 simulation.py --sources scm1_bias0 scm1_bias0.25 scm1_bias0.5 scm1_bias0.75 scm1_bias1 scm1_bias1.25
                                                            the bias sweep: SCM 1 with a smaller confounding bias
    python3 simulation.py --sources scm1 ihdp --m 100 300 --obs 1500 5000 45000 --shifts strong_shift
                                                            the OBS-size sweep (files <source>_n<3m>_N<N>.csv)

The covariate laws: "same" has r_0 = 1; "weak_shift" and "strong_shift" have r_0 != 1, and there every
estimator estimates its balancing ratio.  "misspec_shift" (the misspecified shift) runs only when asked for
with --shifts, and "native" is the law of the lung design.

A cell that is already saved is skipped, so an interrupted run continues where it stopped.  Results made
with other settings belong in another folder (--out).  visualization.py draws the figures from the folder.
"""
import argparse
import csv
from functools import partial
from multiprocessing import Pool
from pathlib import Path

from experiment import COLUMNS, EXTRA_SHIFTS, N_OBS, SHIFTS, SOURCES, SWEEP, design_for, replicate


def one(job, shifts):
    source, m, rep, N = job
    return replicate(source, m, rep, shifts=shifts, N=N)


def saved(path):
    """(number of replications, covariate laws) of a saved cell."""
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    return len({r["rep"] for r in rows}), {r["shift"] for r in rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sources", nargs="+", default=list(SOURCES), choices=SOURCES + tuple(SWEEP))
    parser.add_argument("--m", nargs="+", type=int, default=[100, 200, 300, 500, 1000],
                        help="rows per trial role; the trial has 3m rows")
    parser.add_argument("--shifts", nargs="+", default=list(SHIFTS), choices=SHIFTS + EXTRA_SHIFTS)
    parser.add_argument("--obs", nargs="+", type=int, default=[N_OBS], help="OBS sample sizes N")
    parser.add_argument("--reps", type=int, default=300)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--out", default="results")
    a = parser.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    for source in a.sources:
        design_for(source, 0)                        # downloads the benchmark files once, before any worker starts
        for m, N in ((m, N) for m in a.m for N in a.obs):
            path = out / (f"{source}_n{3 * m}.csv" if N == N_OBS else f"{source}_n{3 * m}_N{N}.csv")
            if path.exists():
                laws = {law for law in a.shifts if law in design_for(source, 0).laws}
                if saved(path) != (a.reps, laws):
                    raise SystemExit(f"{path} was made with other settings; choose another --out")
                print(f"{path}: already saved")
                continue
            jobs = [(source, m, rep, N) for rep in range(a.reps)]
            if a.workers > 1:
                with Pool(a.workers) as pool:
                    results = pool.map(partial(one, shifts=a.shifts), jobs)
            else:
                results = [one(job, a.shifts) for job in jobs]
            with open(path, "w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=COLUMNS, restval="")
                writer.writeheader()
                writer.writerows(row for result in results for row in result)
            print(f"{path}: {a.reps} replications saved", flush=True)
