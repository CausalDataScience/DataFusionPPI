"""Verification for the manuscript-minimal CATE plan, PRD section 8.

Each test names the manuscript object it checks.  Mutation tests inject the
specific error they are meant to catch and assert that the check fires, so a
test cannot pass by being vacuous.
"""
from __future__ import annotations

import unittest

import numpy as np
from scipy.optimize import minimize

import ssem_ate_pilot as legacy
from confirm import mm_design as md
from confirm import mm_engine as me
from confirm import mm_run as mr


def fixture(n_rct=300, n_obs=900, seed=3, propensity=None):
    params = md.scm()
    shared, shifted = md.cells(params)
    rng = np.random.default_rng(seed)
    if propensity is not None:
        md.TRIAL_PROPENSITY = propensity
    rct_draw, obs_draw = md.Draw.make(n_rct, rng), md.Draw.make(n_obs, rng)
    rct = md.realise(rct_draw, "RCT", params, shared)
    obs = md.realise(obs_draw, "OBS", params, shared)
    nuis = me.fit_nuisance(rct, obs, md.TRIAL_PROPENSITY)
    tune = me.trial_parts(rct, nuis)
    b_obs, ghat_obs = me.obs_parts(obs, nuis)
    return params, shared, shifted, rct, obs, nuis, tune, b_obs, ghat_obs


class PseudoOutcome(unittest.TestCase):
    def test_aipw_matches_the_direct_formula(self):
        """2.tex def:fusion-score, written out term by term."""
        _, _, _, rct, _, nuis, tune, _, _ = fixture()
        e = np.full(len(rct["x"]), nuis.trial_propensity)
        q0, q1 = nuis.mu_r.predict(rct["x"])
        direct = (q1 - q0
                  + rct["a"] * (rct["y"] - q1) / e
                  - (1 - rct["a"]) * (rct["y"] - q0) / (1 - e))
        self.assertLess(float(np.max(np.abs(direct - tune.z0))), 1e-12)


class NormalEquations(unittest.TestCase):
    """4.tex eq:dr-fusion-loss and eq:rf-fusion-loss against their derivatives."""

    def objective(self, learner, tune, b_obs, ghat_obs, ratio, lam, omega, beta):
        zeta_trial = tune.b @ beta
        zeta_obs = b_obs @ beta
        correction = (np.mean(ratio * (ghat_obs - zeta_obs) ** 2)
                      - np.mean((np.ones_like(tune.chi) if learner == "DRF" else tune.chi)
                                * (tune.ghat - zeta_trial) ** 2))
        if learner == "DRF":
            z_lambda = (1 - lam) * tune.z0 + lam * tune.z1
            main = np.mean((z_lambda - zeta_trial) ** 2)
        else:
            m_lambda = (1 - lam) * tune.m_r + lam * tune.m_o
            residual = tune.y - m_lambda - tune.centred_arm * zeta_trial
            main = np.mean(residual ** 2) * tune.weight_scale
        return main + omega * correction + me.RIDGE * float(beta @ beta)

    def gradient(self, learner, tune, b_obs, ghat_obs, ratio, lam, omega, beta,
                 ridge=None):
        """The loss gradient, assembled from per-observation residuals.

        The solver builds Gram matrices and solves once.  This route never forms
        a Gram matrix and instead weights each observation's basis row by its own
        residual, so a mistake in the Gram assembly cannot cancel out of both.
        """
        ridge = me.RIDGE if ridge is None else ridge
        zeta_trial = tune.b @ beta
        zeta_obs = b_obs @ beta
        weight = np.ones_like(tune.chi) if learner == "DRF" else tune.chi
        if learner == "DRF":
            z_lambda = (1 - lam) * tune.z0 + lam * tune.z1
            main = -2.0 * (tune.b * (z_lambda - zeta_trial)[:, None]).mean(0)
        else:
            m_lambda = (1 - lam) * tune.m_r + lam * tune.m_o
            residual = tune.y - m_lambda - tune.centred_arm * zeta_trial
            main = -2.0 * tune.weight_scale * (
                tune.b * (tune.centred_arm * residual)[:, None]).mean(0)
        obs_term = -2.0 * (b_obs * (ratio * (ghat_obs - zeta_obs))[:, None]).mean(0)
        trial_term = 2.0 * (tune.b * (weight * (tune.ghat - zeta_trial))[:, None]).mean(0)
        return main + omega * (obs_term + trial_term) + 2.0 * ridge * beta

    def test_solution_solves_the_direct_derivative(self):
        """PRD section 8: the direct derivative and the normal equation agree."""
        _, _, _, _, _, _, tune, b_obs, ghat_obs = fixture(propensity=0.35)
        ratio = np.ones(len(b_obs))
        worst = 0.0
        for learner in ("DRF", "RF"):
            for lam, omega in me.GRID:
                beta = me.solve(learner, tune, b_obs, ghat_obs, ratio, lam, omega)
                grad = self.gradient(learner, tune, b_obs, ghat_obs, ratio,
                                     lam, omega, beta)
                # scale by the size of the gradient's own terms at the origin, so
                # the criterion does not depend on the units of the outcome
                at_origin = self.gradient(learner, tune, b_obs, ghat_obs, ratio,
                                          lam, omega, np.zeros_like(beta))
                relative = float(np.max(np.abs(grad))
                                 / max(float(np.max(np.abs(at_origin))), 1e-12))
                worst = max(worst, relative)
                self.assertLess(relative, 1e-10,
                                f"{learner} ({lam}, {omega}) relative gradient "
                                f"{relative:.2e}")
        print(f"\n    largest relative gradient at the solution {worst:.2e}")

    def test_a_wrong_ridge_is_detected(self):
        """The check fails when the solver uses a penalty the loss does not."""
        _, _, _, _, _, _, tune, b_obs, ghat_obs = fixture(propensity=0.35)
        ratio = np.ones(len(b_obs))
        beta = me.solve("DRF", tune, b_obs, ghat_obs, ratio, 0.5, 0.5,
                        ridge=me.RIDGE * 10)
        grad = self.gradient("DRF", tune, b_obs, ghat_obs, ratio, 0.5, 0.5, beta)
        at_origin = self.gradient("DRF", tune, b_obs, ghat_obs, ratio, 0.5, 0.5,
                                  np.zeros_like(beta))
        self.assertGreater(float(np.max(np.abs(grad))
                                 / np.max(np.abs(at_origin))), 1e-6)


class RFWeight(unittest.TestCase):
    def test_learners_coincide_at_one_half(self):
        """4.tex: at e = 1/2 chi is one and the two pseudo-outcomes agree."""
        try:
            _, _, _, _, _, _, tune, b_obs, ghat_obs = fixture(propensity=0.5)
            ratio = np.ones(len(b_obs))
            worst = max(
                float(np.max(np.abs(me.solve("DRF", tune, b_obs, ghat_obs, ratio, l, o)
                                    - me.solve("RF", tune, b_obs, ghat_obs, ratio, l, o))))
                for l, o in me.GRID)
        finally:
            md.TRIAL_PROPENSITY = 0.35
        self.assertLess(worst, 1e-10, f"largest coefficient gap {worst:.2e}")

    def test_dropping_chi_is_detected_off_one_half(self):
        _, _, _, _, _, _, tune, b_obs, ghat_obs = fixture(propensity=0.35)
        ratio = np.ones(len(b_obs))
        beta = me.solve("RF", tune, b_obs, ghat_obs, ratio, 0.5, 0.5)
        mutated = me.Parts(**{**tune.__dict__, "chi": np.ones_like(tune.chi)})
        broken = me.solve("RF", mutated, b_obs, ghat_obs, ratio, 0.5, 0.5)
        self.assertGreater(float(np.max(np.abs(beta - broken))), 1e-6)


class SelectionScore(unittest.TestCase):
    """4.tex eq:common-reference-validation."""

    def setUp(self):
        (_, _, _, _, _, self.nuis, self.tune, self.b_obs,
         self.ghat_obs) = fixture()
        self.ratio = np.ones(len(self.b_obs))
        self.beta = me.solve("DRF", self.tune, self.b_obs, self.ghat_obs,
                             self.ratio, 0.5, 0.5)

    def score(self, omega, ratio=None):
        return me.selection_score(self.beta, self.tune, self.b_obs, self.ghat_obs,
                                  self.ratio if ratio is None else ratio, omega)

    def test_observational_correction_is_present(self):
        """Removing the observational term must move the score."""
        truth = self.score(0.5)
        zeta = self.tune.b @ self.beta
        without = float(np.mean((self.tune.z0 - zeta) ** 2
                                - 0.5 * (self.tune.ghat - zeta) ** 2))
        self.assertGreater(abs(truth - without), 1e-8)

    def test_selection_omega_is_the_candidate_omega(self):
        """A separate validation omega would change the score at fixed beta."""
        self.assertGreater(abs(self.score(0.5) - self.score(1.0)), 1e-8)

    def test_selection_trial_term_carries_no_chi(self):
        zeta = self.tune.b @ self.beta
        stated = self.score(0.5)
        with_chi = float(np.mean((self.tune.z0 - zeta) ** 2
                                 - 0.5 * self.tune.chi * (self.tune.ghat - zeta) ** 2))
        with_chi += float(np.mean(0.5 * self.ratio * (self.ghat_obs
                                                      - self.b_obs @ self.beta) ** 2))
        self.assertGreater(abs(stated - with_chi), 1e-8)


class TargetAndRatio(unittest.TestCase):
    def test_exact_ratio_matches_the_gaussian_density_difference(self):
        params = md.scm()
        _, shifted = md.cells(params)
        rng = np.random.default_rng(11)
        x = md.realise(md.Draw.make(4000, rng), "OBS", params, shifted)["x"]
        cov = params.covariance
        log_r = (legacy.gaussian_log_density(x, params.covariate_mean, cov)
                 - legacy.gaussian_log_density(
                     x, params.covariate_mean + shifted.obs_offset, cov))
        analytic = md.exact_ratio(x, params, shifted)
        gap = float(np.max(np.abs(np.log(analytic) - log_r)))
        self.assertLess(gap, 1e-10, f"largest log-ratio gap {gap:.2e}")

    def test_zero_shift_gives_exactly_one(self):
        params = md.scm()
        shared, _ = md.cells(params)
        rng = np.random.default_rng(12)
        x = md.realise(md.Draw.make(500, rng), "OBS", params, shared)["x"]
        self.assertTrue(np.all(md.exact_ratio(x, params, shared) == 1.0))

    def test_trial_arrays_are_identical_across_cells(self):
        params = md.scm()
        shared, shifted = md.cells(params)
        draw = md.Draw.make(400, np.random.default_rng(13))
        a = md.realise(draw, "RCT", params, shared)
        b = md.realise(draw, "RCT", params, shifted)
        for key in ("x", "a", "y", "tau"):
            self.assertTrue(np.array_equal(a[key], b[key]), key)


class FairCandidates(unittest.TestCase):
    def test_every_rule_reads_the_same_nine_candidates(self):
        scores = {key: float(i) for i, key in enumerate(me.GRID)}
        self.assertEqual(me.pick(scores, "rct_only"), (0.0, 0.0))
        self.assertTrue(all(k[1] == 0.0 for k in
                            [me.pick({k: v for k, v in scores.items()
                                      if me.RULES["lambda_only"](*k)}, "lambda_only")]))
        self.assertEqual(len([k for k in scores if me.RULES["joint"](*k)]), 9)

    def test_ties_break_lexicographically(self):
        scores = {key: 1.0 for key in me.GRID}
        self.assertEqual(me.pick(scores, "joint"), (0.0, 0.0))
        scores[(0.0, 0.0)] = 2.0
        self.assertEqual(me.pick(scores, "joint"), (0.0, 0.5))

    def test_truth_never_reaches_selection(self):
        """The selection score takes no tau argument, by signature."""
        import inspect
        names = set(inspect.signature(me.selection_score).parameters)
        self.assertNotIn("tau", names)
        self.assertNotIn("truth", names)


class Decomposition(unittest.TestCase):
    def test_identity_holds_on_random_predictions(self):
        rng = np.random.default_rng(17)
        k, grid = 12, 800
        tau = rng.normal(size=grid)
        matrix = rng.normal(size=(k, grid))
        gram = matrix @ matrix.T / grid
        v = matrix @ tau / grid
        t = float(np.mean(tau ** 2))
        a = np.full(k, 1.0 / k)
        mean_risk = float(np.sum(a * (np.diag(gram) - 2 * v + t)))
        bias2 = float(a @ gram @ a - 2 * a @ v + t)
        variance = float((np.sum(a * np.diag(gram)) - a @ gram @ a) * k / (k - 1))
        direct = float(np.mean(np.mean((matrix - tau) ** 2, axis=1)))
        self.assertLess(abs(mean_risk - direct) / abs(direct), 1e-12)
        self.assertLess(abs(mean_risk - (bias2 + (k - 1) / k * variance))
                        / abs(mean_risk), 1e-12)

    def test_unpaired_resampling_is_detected(self):
        """Resampling groups separately breaks the pairing the plan requires."""
        rng = np.random.default_rng(19)
        k, grid = 40, 400
        tau = rng.normal(size=grid)
        a_mat = rng.normal(size=(k, grid)) + tau
        b_mat = a_mat + 0.10 * rng.normal(size=(k, grid))
        risk = lambda m, idx: float(np.mean(np.mean((m[idx] - tau) ** 2, axis=1)))
        paired, split = [], []
        for _ in range(400):
            idx = rng.integers(0, k, size=k)
            other = rng.integers(0, k, size=k)
            paired.append(risk(a_mat, idx) - risk(b_mat, idx))
            split.append(risk(a_mat, idx) - risk(b_mat, other))
        self.assertLess(np.std(paired), 0.5 * np.std(split))


if __name__ == "__main__":
    unittest.main()
