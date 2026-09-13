"""M1 tie: the new moment engine against the exploratory production engine.

M1's first Go condition compares an independent formula against the production
``A_p, B_p, C_p, D_p``.  The production implementation is
``fusion_cate._oracle_calculus``, which draws its own sample internally.  This
module hands it a sampler that returns one fixed array whatever it is asked for,
then rebuilds the same arrays through the public score builders and feeds them to
``confirm.moments``.  Both engines therefore see byte-identical data and any
difference is a formula difference rather than a sampling difference.

Run: python3 -m confirm.m1_production_tie
"""
from __future__ import annotations

import sys

import numpy as np

import fusion_cate as fc
from fusion_core import build_scores, ghat_on
from confirm import moments as mo
from confirm import m1_protocol as protocol

TOL = 1e-10


def _outcome_pair(kind: str, gap: float = 0.30):
    """Two deterministic outcome regressions with a ``predict`` method.

    ``fusion_cate._oracle_calculus`` only ever calls ``predict``, so a
    deterministic stand-in keeps the comparison free of any fitting step.  The
    observational pair carries the true contrast, including the part of tau the
    basis cannot span, which is what makes D_p positive.  The trial pair is the
    worse of the two by a factor ``gap``, which makes Delta variance-reducing and
    C_p negative; ``gap`` is tuned so that lambda* lands strictly inside [0, 1]
    rather than clipping at a boundary, where a sign error could hide.
    """
    class Fixed:
        def __init__(self, baseline_error: float, contrast_error: float) -> None:
            self.baseline_error = baseline_error
            self.contrast_error = contrast_error

        def predict(self, x: np.ndarray):
            base = x[:, 0]
            mu0 = 0.4 * base + self.baseline_error * np.cos(2.0 * base)
            mu1 = mu0 + true_tau(x) + self.contrast_error * base
            return mu0, mu1

    # The two regressions must err in opposite directions for the optimal blend
    # to be interior.  The pseudo-outcome is affine in the regression, so
    # Z_lambda is the score of (1-lambda) mu_R + lambda mu_O; when both err the
    # same way, the better one wins outright and lambda* clips at a boundary.
    if kind == "O":
        return Fixed(-gap * 0.45, -gap * 0.30)
    return Fixed(gap, gap * 0.65)


def true_tau(x: np.ndarray) -> np.ndarray:
    """A cubic the basis spans, plus a cosine it does not.

    The off-span term is what makes D_p non-zero, because D_p is the whitened
    covariance of b(tau - zeta_p) with b(ghat - zeta_p) and vanishes identically
    when the basis spans tau.
    """
    base = x[:, 0]
    return 0.3 + 0.45 * base - 0.2 * base ** 2 + 0.30 * np.cos(6.0 * base)


def build_case(seed: int, n_rct: int = 20000, n_obs: int = 30000):
    rng = np.random.default_rng(seed)
    def draw(n: int) -> dict[str, np.ndarray]:
        x = rng.uniform(-1.0, 1.0, size=(n, 2))
        e = np.full(n, 0.5)
        a = (rng.uniform(size=n) < e).astype(float)
        y = 0.4 * x[:, 0] + a * true_tau(x) + rng.normal(scale=0.5, size=n)
        return {"x": x, "a": a, "y": y, "propensity": e}

    rct, obs = draw(n_rct), draw(n_obs)
    test_x = rng.uniform(-1.0, 1.0, size=(40000, 2))
    test = {"x": test_x, "tau": true_tau(test_x)}
    basis = lambda x: np.vander(x[:, 0], 4, increasing=True)
    ratio = lambda x: np.ones(len(x))
    return rct, obs, test, basis, ratio


def cross_term_with_error(fixture: mo.Fixture, resamples: int = 400,
                          *, seed: int) -> tuple[float, float]:
    """E_p and a bootstrap standard error over the trial rows."""
    rng = np.random.default_rng(seed)
    import dataclasses
    point = mo.covariance_route(fixture)["E_p"]
    n = len(fixture.b_rct)
    draws = []
    for _ in range(resamples):
        idx = rng.integers(0, n, size=n)
        sub = dataclasses.replace(
            fixture, b_rct=fixture.b_rct[idx], z0_rct=fixture.z0_rct[idx],
            delta_rct=fixture.delta_rct[idx], ghat_rct=fixture.ghat_rct[idx],
            weight_rct=fixture.weight_rct[idx])
        draws.append(mo.covariance_route(sub)["E_p"])
    return point, float(np.std(draws, ddof=1))


def run(*, sampler_seed: int, bootstrap_seed: int) -> dict:
    """The tie, as data.  ``main`` prints it; the M1 driver records it."""
    rct, obs, test, basis, ratio = build_case(sampler_seed)
    mu_r, mu_o = _outcome_pair("R"), _outcome_pair("O")
    # both oracle coordinates must be strictly interior for this tie to be able
    # to catch a sign error; the check below asserts that rather than trusting it
    b_test = basis(test["x"])
    n_tune_r, n_tune_o = 120, 900

    # the production engine, given a sampler that always returns the same arrays
    frozen = lambda n, seed: (rct, obs)
    produced = fc._oracle_calculus(frozen, mu_r, mu_o, ratio, basis, b_test, test,
                                   n_tune_r=n_tune_r, n_tune_o=n_tune_o,
                                   seed=sampler_seed)

    # the same arrays, routed through the confirmatory engine
    scores = build_scores(rct, mu_r, mu_o)
    fixture = mo.Fixture(
        name="production_tie", b_target=b_test, tau_target=test["tau"],
        b_rct=basis(rct["x"]), z0_rct=scores.z0, delta_rct=scores.delta,
        ghat_rct=scores.ghat, weight_rct=np.ones(len(rct["x"])),
        b_obs=basis(obs["x"]), ghat_obs=ghat_on(obs["x"], mu_o),
        ratio_obs=ratio(obs["x"]), n_rct_tune=n_tune_r, n_obs_tune=n_tune_o)

    import dataclasses
    mine = mo.covariance_route(fixture)
    other = mo.whitened_route(fixture)
    # the exploratory engine stabilises Gamma with a 1e-10 ridge; matching it is
    # the only way to tell a formula difference from a conditioning choice
    matched = mo.covariance_route(dataclasses.replace(fixture, ridge=1e-10))

    rows, worst = [], 0.0
    for key in ("A_p", "B_p", "C_p", "D_p"):
        gap_matched = abs(produced[key] - matched[key])
        gap_routes = abs(mine[key] - other[key])
        worst = max(worst, gap_matched, gap_routes)
        rows.append((key, produced[key], mine[key], matched[key],
                     abs(produced[key] - mine[key]), gap_matched, gap_routes))

    lam_p, om_p = produced["oracle_sieve_lambda"], produced["oracle_sieve_omega"]
    lam_m, om_m = mo.oracle_coordinates(matched)
    coord_gap = max(abs(lam_p - lam_m), abs(om_p - om_m))
    worst = max(worst, coord_gap)
    interior = 1e-6 < lam_m < 1 - 1e-6 and 1e-6 < om_m < 1 - 1e-6

    e_p, mcse = cross_term_with_error(fixture, seed=bootstrap_seed)
    return {
        "moments": [{"moment": k, "production": prod, "confirmatory_rho0": mynum,
                     "confirmatory_rho1e10": match, "gap_at_rho0": g0,
                     "gap_matched": gm, "gap_between_routes": gr}
                    for k, prod, mynum, match, g0, gm, gr in rows],
        "oracle": {"production": [lam_p, om_p], "confirmatory": [lam_m, om_m],
                   "gap": coord_gap, "interior": bool(interior)},
        "cross_term": {"E_p": e_p, "mcse": mcse, "z": e_p / mcse},
        "worst_gap": worst, "tolerance": TOL,
        "status": "PASS" if (worst <= TOL and interior) else "FAIL",
    }


def main() -> int:
    ledger = protocol.SeedLedger()
    result = run(sampler_seed=ledger.seed("production_tie/sampler"),
                 bootstrap_seed=ledger.seed("production_tie/bootstrap"))
    print(f"{'moment':6s}{'production':>15s}{'rho=0':>15s}{'rho=1e-10':>15s}"
          f"{'gap at rho=0':>14s}{'gap matched':>13s}{'route gap':>12s}")
    for row in result["moments"]:
        print(f"{row['moment']:6s}{row['production']:15.10f}"
              f"{row['confirmatory_rho0']:15.10f}{row['confirmatory_rho1e10']:15.10f}"
              f"{row['gap_at_rho0']:14.2e}{row['gap_matched']:13.2e}"
              f"{row['gap_between_routes']:12.2e}")
    oracle = result["oracle"]
    print(f"oracle coordinates: production ({oracle['production'][0]:.6f}, "
          f"{oracle['production'][1]:.6f}) confirmatory ({oracle['confirmatory'][0]:.6f}, "
          f"{oracle['confirmatory'][1]:.6f}), gap {oracle['gap']:.2e}, "
          f"both interior {oracle['interior']}")
    cross = result["cross_term"]
    print(f"cross term E_p {cross['E_p']:+.3e} with Monte Carlo standard error "
          f"{cross['mcse']:.3e}, z = {cross['z']:+.2f}; the DRF corollary puts it "
          "at zero in population")
    print()
    print(f"M1 production tie {result['status']}: largest difference "
          f"{result['worst_gap']:.2e} against {TOL:.0e}")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
