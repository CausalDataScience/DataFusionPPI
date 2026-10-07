"""Checks that aipwf.py, drf.py, rf.py and method.py compute what the paper defines, and that data.py
generates the designs it describes.  Run in this folder:  python3 -m unittest test_paper -v"""
import unittest

import numpy as np
import pandas as pd

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "experiments")]

from dfppi import aipwf as AIPWF, drf as DRF, rf as RF, method as M  # noqa: E402
import experiment as E  # noqa: E402
from data import (E_TRIAL, MEAN_BIAS_IN_SD_Y, NU0, RMS_BIAS_IN_SD_TAU, SIGMA, Benchmark, LungCancer, SyntheticSCM,
                  UnrelatedSCM, confounding_bias, dataset)

LEARNERS = (DRF, RF)


def split(sample, sizes):
    cuts = np.cumsum([0] + list(sizes.values()))
    return {role: {k: v[cuts[i]:cuts[i + 1]] for k, v in sample.items()} for i, role in enumerate(sizes)}


def proxy_risk(learner, beta, lam, omega, nuis, r, R, O):
    """Eq. (9) for the candidate t = zeta'beta, written directly from its definition."""
    kappa, Z = learner.pseudo_outcome(R, M.blended_regression(nuis, lam), nuis.e)
    t_R, t_O = DRF.sieve(nuis, R["x"]) @ beta, DRF.sieve(nuis, O["x"]) @ beta
    return (np.mean(kappa * (Z - t_R) ** 2) + omega * np.mean(r(O["x"]) * (nuis.g(O["x"]) - t_O) ** 2)
            - omega * np.mean(kappa * (nuis.g(R["x"]) - t_R) ** 2))


class PaperChecks(unittest.TestCase):
    """The estimators under a covariate shift (r_0 != 1), where each of them estimates its balancing ratio."""

    @classmethod
    def setUpClass(cls):
        cls.scm, rng = SyntheticSCM(), np.random.default_rng(0)
        cls.R = split(cls.scm.trial(cls.scm.draw(600, rng)), {"nuis": 200, "tune": 200, "eval": 200})
        cls.O = split(cls.scm.obs(cls.scm.draw(5000, rng), "strong_shift"), {"nuis": 3000, "tune": 1000, "eval": 1000})
        cls.nuis = M.fit_nuisances(cls.R["nuis"], cls.O["nuis"], E_TRIAL, same_domain=False)
        cls.r_ate = staticmethod(AIPWF.density_ratio(cls.nuis))           # staticmethod: a function kept on a class
        cls.r_cate = staticmethod(DRF.density_ratio(cls.nuis))
        cls.x_R, cls.x_O = cls.R["nuis"]["x"], cls.O["nuis"]["x"]

    def test_def5_ate_ratio_balances_g_exactly(self):
        g, r = self.nuis.g, self.r_ate
        self.assertAlmostEqual(np.mean(r(self.x_O) * g(self.x_O)), g(self.x_R).mean(), places=10)
        self.assertAlmostEqual(r(self.x_O).mean(), 1.0, places=12)

    def test_prop7_cate_ratio_balances_its_features_within_their_standard_errors(self):
        n = self.nuis
        F_O = DRF.balance_features(n.components(self.x_O), n.g(self.x_O))
        F_R = DRF.balance_features(n.components(self.x_R), n.g(self.x_R))
        q = n.components(self.x_R).shape[1] - 1                           # number of principal components
        self.assertEqual(F_O.shape[1], 2 + 2 * q + q * (q + 1) // 2)      # g, g^2, Z_k, Z_k g, Z_k Z_l
        # these are the distinct products of Prop. 6(ii) for zeta = (1, Z, g): g^2, zeta g, zeta_k zeta_l
        zeta = DRF.sieve(n, self.x_O)
        products = [n.g(self.x_O) ** 2] + [zeta[:, k] * n.g(self.x_O) for k in range(q + 2)] + \
                   [zeta[:, k] * zeta[:, l] for k in range(q + 2) for l in range(k, q + 2)]
        for column in products:
            if np.ptp(column) > 0:                                        # skip the constant 1 * 1
                self.assertLess(np.min(np.abs(F_O - column[:, None]).max(0)), 1e-12)
        r_base = n.r_base(self.x_O)
        delta = np.sqrt(np.var(r_base[:, None] * F_O, axis=0, ddof=1) / len(F_O)
                        + np.var(F_R, axis=0, ddof=1) / len(F_R))
        gap = np.abs((self.r_cate(self.x_O)[:, None] * F_O).mean(0) - F_R.mean(0))
        self.assertTrue(np.all(gap <= delta + 1e-6 * F_O.std(0)))
        self.assertAlmostEqual(self.r_cate(self.x_O).mean(), 1.0, places=12)

    def test_the_two_ratios_start_from_the_same_classifier_and_differ(self):
        base, ate, cate = self.nuis.r_base(self.x_O), self.r_ate(self.x_O), self.r_cate(self.x_O)
        self.assertTrue(np.all(base > 0) and np.all(ate > 0) and np.all(cate > 0))
        self.assertGreater(np.abs(ate - cate).max(), 1e-6)                # each estimator balances its own functions
        self.assertGreater(np.corrcoef(np.log(ate), np.log(base))[0, 1], 0.9)
        self.assertGreater(np.corrcoef(np.log(cate), np.log(base))[0, 1], 0.9)

    def test_assumption3i_a_known_ratio_is_one_and_adds_no_variance(self):
        O = split(self.scm.obs(self.scm.draw(5000, np.random.default_rng(3)), "same"),
                  {"nuis": 3000, "tune": 1000, "eval": 1000})
        nuis = M.fit_nuisances(self.R["nuis"], O["nuis"], E_TRIAL, same_domain=True)
        r_ate, r_cate = AIPWF.density_ratio(nuis), DRF.density_ratio(nuis)
        self.assertTrue(np.all(r_ate(O["eval"]["x"]) == 1.0) and np.all(r_cate(O["eval"]["x"]) == 1.0))
        out = AIPWF.aipwf(nuis, r_ate, self.R["tune"], O["tune"], self.R["eval"], O["eval"])
        g, omega = nuis.g, out["omega"]
        psi_R = M.fused_score(self.R["eval"], nuis, out["lambda"]) - omega * g(self.R["eval"]["x"])
        V = np.var(psi_R, ddof=1) / 200 + np.var(omega * g(O["eval"]["x"]), ddof=1) / 1000
        self.assertAlmostEqual(out["se"], np.sqrt(V), places=12)          # no calibration variance

    def test_thm3ii_an_estimated_ratio_adds_the_calibration_variance(self):
        out = AIPWF.aipwf(self.nuis, self.r_ate, self.R["tune"], self.O["tune"], self.R["eval"], self.O["eval"])
        g, r, omega = self.nuis.g, self.r_ate, out["omega"]
        self.assertGreater(omega, 0.0)
        psi_R = M.fused_score(self.R["eval"], self.nuis, out["lambda"]) - omega * g(self.R["eval"]["x"])
        psi_O = omega * r(self.O["eval"]["x"]) * g(self.O["eval"]["x"])
        V = np.var(psi_R, ddof=1) / 200 + np.var(psi_O, ddof=1) / 1000
        V_cal = omega ** 2 * (np.var(g(self.x_R), ddof=1) / 200 + np.var(r(self.x_O) * g(self.x_O), ddof=1) / 3000)
        self.assertAlmostEqual(out["se"], np.sqrt(V + V_cal), places=12)

    def test_prop5_normal_equation_minimises_eq8_plus_ridge(self):
        rng, ridge = np.random.default_rng(1), 0.01
        for learner in LEARNERS:
            for lam, omega in DRF.GRID:
                beta = DRF.fit_candidate(lam, omega, self.nuis, self.r_cate, self.R["tune"], self.O["tune"], ridge,
                                         learner.pseudo_outcome)
                J = lambda b: proxy_risk(learner, b, lam, omega, self.nuis, self.r_cate, self.R["tune"],
                                         self.O["tune"]) + ridge * b @ b
                for _ in range(5):
                    step = 1e-2 * rng.normal(size=len(beta))
                    self.assertGreater(J(beta + step), J(beta))
                    self.assertAlmostEqual((J(beta + step) - J(beta - step)) / 2, 0.0, places=9)   # zero gradient

    def test_def10_score_is_the_drf_risk_at_lambda_zero_on_the_evaluation_samples(self):
        beta = DRF.fit_candidate(0.5, 0.5, self.nuis, self.r_cate, self.R["tune"], self.O["tune"],
                                 pseudo_outcome=RF.pseudo_outcome)
        for omega in (0.0, 0.5, 1.0):
            self.assertAlmostEqual(
                DRF.validation_score(beta, omega, self.nuis, self.r_cate, self.R["eval"], self.O["eval"]),
                proxy_risk(DRF, beta, 0.0, omega, self.nuis, self.r_cate, self.R["eval"], self.O["eval"]), places=10)

    def test_rf_is_drf_with_the_pseudo_outcome_of_def8(self):
        samples = (self.R["tune"], self.O["tune"], self.R["eval"], self.O["eval"])
        lam, omega, beta = RF.fit(self.nuis, self.r_cate, *samples)
        lam_d, omega_d, beta_d = DRF.fit(self.nuis, self.r_cate, *samples, pseudo_outcome=RF.pseudo_outcome)
        self.assertEqual((lam, omega), (lam_d, omega_d))
        self.assertTrue(np.array_equal(beta, beta_d))
        kappa, Z = RF.pseudo_outcome(self.R["tune"], self.nuis.mu_R, E_TRIAL)
        a, y = self.R["tune"]["a"], self.R["tune"]["y"]
        mu0, mu1 = self.nuis.mu_R(self.R["tune"]["x"])
        self.assertTrue(np.allclose(kappa * Z ** 2, (y - E_TRIAL * mu1 - (1 - E_TRIAL) * mu0) ** 2
                                    / (E_TRIAL * (1 - E_TRIAL))))        # kappa Z^2 is the squared residual
        self.assertTrue(np.all(DRF.pseudo_outcome(self.R["tune"], self.nuis.mu_R, E_TRIAL)[0] == 1.0))
        self.assertIs(RF.density_ratio, DRF.density_ratio)
        self.assertIs(RF.sieve, DRF.sieve)

    def test_shrinkage_and_pretest_pool_as_far_as_the_gap_allows(self):
        self.assertEqual(AIPWF.shrinkage(1.0, 0.5, 1.0), 1.0)                  # no gap: the pooled estimate
        self.assertAlmostEqual(AIPWF.shrinkage(0.0, 1.0, 1.0), 0.5)            # gap of one standard error: halfway
        self.assertLess(abs(AIPWF.shrinkage(0.0, 0.1, 10.0)), 0.01)            # large gap: the trial estimate
        self.assertEqual(AIPWF.pretest(0.0, 1.0, 1.9), 1.9)
        self.assertEqual(AIPWF.pretest(0.0, 1.0, 2.0), 0.0)

    def test_two_step_corrects_g_by_a_ridge_fit_on_the_trial(self):
        rng, ridge = np.random.default_rng(2), 0.01
        basis, g = self.nuis.components(self.R["tune"]["x"]), self.nuis.g(self.R["tune"]["x"])
        for learner in LEARNERS:
            beta = learner.two_step(self.nuis, self.R["tune"], ridge)
            kappa, Z = learner.pseudo_outcome(self.R["tune"], self.nuis.mu_R, E_TRIAL)
            J = lambda b: np.mean(kappa * (Z - g - basis @ b) ** 2) + ridge * b @ b
            for _ in range(5):
                step = 1e-2 * rng.normal(size=len(beta))
                self.assertGreater(J(beta + step), J(beta))
                self.assertAlmostEqual((J(beta + step) - J(beta - step)) / 2, 0.0, places=9)   # zero gradient

    def test_integrative_r_learner_reduces_to_the_trial_r_learner_when_c_is_free(self):
        e_O = M.fit_propensity(self.O["nuis"])
        p = e_O(self.O["eval"]["x"])
        self.assertTrue(np.all((p >= 0.01) & (p <= 0.99)))
        samples = (self.R["tune"], self.O["tune"], self.R["eval"], self.O["eval"])
        beta = RF.integrative(self.nuis, e_O, *samples)
        self.assertEqual(len(beta), self.nuis.components(self.x_R).shape[1])
        # with no penalty on tau and almost none on c, c absorbs the OBS rows and tau is the trial R-learner,
        # whether c lies in the basis of tau or also uses g
        trial_only = RF.rct_only(self.nuis, self.R["tune"], ridge=1e-10)
        for c_basis in (None, lambda x: DRF.sieve(self.nuis, x)):
            free = RF.integrative(self.nuis, e_O, *samples, grid=((0.0, 1e-10),), c_basis=c_basis)
            self.assertTrue(np.allclose(free, trial_only, atol=1e-6))

    def test_ir_confounding_basis_takes_its_components_from_the_obs(self):
        basis = RF.c_basis(self.nuis, self.x_O)
        F = basis(self.x_O)
        self.assertEqual(F.shape[1], 1 + 5 + 1)                            # 1, five OBS components, g
        self.assertTrue(np.all(F[:, 0] == 1.0))
        self.assertTrue(np.allclose(F[:, 1:6].std(0), 1.0))                 # unit sd on the OBS nuisance rows
        self.assertTrue(np.array_equal(F[:, -1], self.nuis.g(self.x_O)))

    def test_sieve_is_the_principal_components_and_g(self):
        x = self.R["tune"]["x"]
        zeta = DRF.sieve(self.nuis, x)
        self.assertTrue(np.array_equal(zeta[:, :-1], self.nuis.components(x)))
        self.assertTrue(np.array_equal(zeta[:, -1], self.nuis.g(x)))
        self.assertTrue(np.all(zeta[:, 0] == 1.0))

    def test_outcome_regressions_do_not_predict_beyond_their_own_range(self):
        far = 50.0 * self.x_R                                             # covariates far outside both samples
        for mu, own in ((self.nuis.mu_R, self.x_R), (self.nuis.mu_O, self.x_O)):
            for fitted, extrapolated in zip(mu(own), mu(far)):            # the two arms
                self.assertGreaterEqual(extrapolated.min(), fitted.min())
                self.assertLessEqual(extrapolated.max(), fitted.max())

    def test_rct_only_ate_is_the_aipw_mean_on_the_evaluation_sample(self):
        out = AIPWF.aipwf(self.nuis, self.r_ate, self.R["tune"], self.O["tune"], self.R["eval"], self.O["eval"],
                          False, False)
        Z0 = M.aipw_score(self.R["eval"], self.nuis.mu_R, E_TRIAL)
        self.assertEqual((out["lambda"], out["omega"]), (0.0, 0.0))
        self.assertAlmostEqual(out["estimate"], Z0.mean(), places=12)
        self.assertAlmostEqual(out["se"], np.sqrt(np.var(Z0, ddof=1) / len(Z0)), places=12)

    def test_rct_only_cate_does_not_depend_on_the_obs(self):
        other = split(self.scm.obs(self.scm.draw(5000, np.random.default_rng(9)), "weak_shift"),
                      {"nuis": 3000, "tune": 1000, "eval": 1000})
        nuis_other = M.fit_nuisances(self.R["nuis"], other["nuis"], E_TRIAL, same_domain=False)
        for learner in LEARNERS:
            a = learner.rct_only(self.nuis, self.R["tune"])
            b = learner.rct_only(nuis_other, self.R["tune"])
            self.assertTrue(np.array_equal(a, b))
            self.assertEqual(len(a), self.nuis.components(self.x_R).shape[1])              # no coefficient for g

    def test_coefficients_stay_in_the_unit_square_and_the_grid_contains_the_origin(self):
        out = AIPWF.aipwf(self.nuis, self.r_ate, self.R["tune"], self.O["tune"], self.R["eval"], self.O["eval"])
        self.assertTrue(0.0 <= out["lambda"] <= 1.0 and 0.0 <= out["omega"] <= 1.0)
        self.assertIn((0.0, 0.0), DRF.GRID)


class DesignChecks(unittest.TestCase):
    """The designs of data.py against their description at the top of that file."""

    def test_synthetic_scm(self):
        scm = SyntheticSCM()
        self.assertAlmostEqual(1.0 - scm.load @ scm.precision @ scm.load, 0.62, places=2)      # Var(U | X)
        self.assertAlmostEqual(scm.truth_tau.var() / SIGMA ** 2, 3.0, delta=0.1)               # Var_R tau = 3 sigma^2
        x = scm.truth_x[:20_000]
        delta = scm.delta(x)
        self.assertAlmostEqual(np.sqrt(np.mean(delta ** 2)) / scm.tau(x).std(), RMS_BIAS_IN_SD_TAU, delta=0.04)
        self.assertGreater(delta.std(), 0.5 * abs(delta.mean()))                               # delta varies with x
        # a shift moves the mean of X by a d, with a^2 = -log(effective sample fraction) and d'Sigma^{-1}d = 1
        numbers = scm.draw(200_000, np.random.default_rng(0))
        moved = scm.covariates(numbers, "strong_shift").mean(0) - scm.covariates(numbers).mean(0)
        self.assertAlmostEqual(np.exp(-moved @ scm.precision @ moved), 0.40, places=10)

    def test_obs_sizes_keep_the_split_of_the_paper(self):
        self.assertEqual(E.obs_sizes(15_000), {"nuis": 9000, "tune": 3000, "eval": 3000})
        for N in (1500, 5000, 45_000):
            sizes = E.obs_sizes(N)
            self.assertEqual(sum(sizes.values()), N)
            self.assertEqual(sizes["nuis"], 3 * sizes["tune"])

    def test_bias_sweep_changes_only_the_obs_treatment_and_hits_its_bias_size(self):
        base, numbers = SyntheticSCM(), SyntheticSCM.draw(4000, np.random.default_rng(0))
        for bias in NU0:
            scm = SyntheticSCM(bias)
            x = scm.truth_x[:20_000]
            self.assertAlmostEqual(np.sqrt(np.mean(scm.delta(x) ** 2)) / scm.tau(x).std(), bias, delta=0.04)
            for key in ("x", "a", "y"):                                    # the trial does not change
                self.assertTrue(np.array_equal(scm.trial(numbers)[key], base.trial(numbers)[key]))
            self.assertTrue(np.array_equal(scm.obs(numbers)["x"], base.obs(numbers)["x"]))
        self.assertLess(np.abs(SyntheticSCM(0.0).delta(x)).max(), 1e-12)   # bias 0: no unmeasured confounding
        self.assertEqual(SyntheticSCM().nu0, NU0[RMS_BIAS_IN_SD_TAU])      # SCM 1 is the largest size of the sweep

    def test_scm3_shares_the_trial_of_scm1_and_uses_other_covariates_in_the_obs(self):
        one, three = SyntheticSCM(), UnrelatedSCM()
        numbers = SyntheticSCM.draw(20_000, np.random.default_rng(0))
        for key in ("x", "a", "y"):
            self.assertTrue(np.array_equal(one.trial(numbers)[key], three.trial(numbers)[key]))
        obs1, obs3 = one.obs(numbers), three.obs(numbers)
        self.assertTrue(np.array_equal(obs1["x"], obs3["x"]))                              # same covariates
        self.assertEqual((three.theta, three.truth_tau.shape), (one.theta, one.truth_tau.shape))   # same truth
        x = obs3["x"]
        mirrored_effect = three.tau(x[:, ::-1])                                            # the OBS effect
        self.assertLess(abs(np.corrcoef(mirrored_effect, three.tau(x))[0, 1]), 0.1)        # unrelated to tau_0
        self.assertGreater(np.mean(obs1["a"] != obs3["a"]), 0.2)                           # other treatment rule

    def test_lung_cancer_pair(self):
        lung = LungCancer()
        self.assertEqual((len(lung.x_R), len(lung.x_O), lung.x_R.shape[1]), (335, 16217, 5))
        self.assertEqual(lung.laws, ("native",))
        self.assertAlmostEqual(lung.theta, np.mean(lung.t1_R - lung.t0_R), places=12)
        numbers = lung.draw(3000, np.random.default_rng(0))
        trial, obs = lung.trial(numbers), lung.obs(numbers)
        rows = (numbers["row_uniform"] * 335).astype(int)
        self.assertTrue(np.array_equal(trial["y"], np.where(trial["a"] == 1, lung.t1_R[rows], lung.t0_R[rows])))
        self.assertAlmostEqual(trial["a"].mean(), E_TRIAL, delta=0.03)                     # randomised here
        self.assertTrue(set(np.unique(obs["a"])) == {0.0, 1.0} and 0.2 < obs["a"].mean() < 0.35)   # as in the EHR
        self.assertGreater(abs(lung.x_O[:, 1].mean() - lung.x_R[:, 1].mean()), 5.0)         # age differs: r_0 != 1

    def test_benchmarks(self):
        for name, index, units, covariates in (("ihdp", 0, 747, 25), ("acic", 1, 4802, 58)):
            b = Benchmark(name, index)
            self.assertEqual(len(b.x), units)
            raw_covariates = b.x.shape[1] if name == "ihdp" else pd.read_csv(dataset("acic2016/x.csv")).shape[1]
            self.assertEqual(raw_covariates, covariates)
            e = E_TRIAL                                                    # Var_R of the noiseless outcome is one
            noiseless = np.mean(e * b.mu1 ** 2 + (1 - e) * b.mu0 ** 2) - np.mean(e * b.mu1 + (1 - e) * b.mu0) ** 2
            self.assertAlmostEqual(noiseless, 1.0, places=10)
            delta = confounding_bias(b.eta, b.nu, np.zeros(units), 1.0, 1.0)
            self.assertAlmostEqual(abs(delta.mean()) / np.sqrt(3.0), MEAN_BIAS_IN_SD_Y, places=8)
            self.assertGreater(delta.std(), 0.5 * abs(delta.mean()))                           # delta varies with x
            for shift, fraction in (("weak_shift", 0.70), ("strong_shift", 0.40)):
                r0 = b.row_law["same"] / b.row_law[shift]                                      # density ratio dP_R / dP_O
                self.assertAlmostEqual(1.0 / np.sum(b.row_law[shift] * r0 ** 2), fraction, places=8)


if __name__ == "__main__":
    unittest.main()
