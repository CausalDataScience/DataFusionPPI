"""The public interface (dfppi.api, dfppi.cli) against the experiment code of the paper: on the same draw, with the
trial and OBS rows in their drawn order (shuffle=False), fuse_ate and fuse_cate return what experiment.replicate
records for the cross-fitted Algorithms 1 and 2."""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "experiments")]

import dfppi  # noqa: E402
import experiment as E  # noqa: E402
from data import E_TRIAL  # noqa: E402
from dfppi.cli import main  # noqa: E402


def draw(source, m, rep, law):
    """The trial and OBS of experiment.replicate(source, m, rep) under one covariate law."""
    design = E.design_for(source, rep)
    rng = np.random.default_rng([E.SEED, E.SOURCES.index(source), m, rep])
    trial = design.trial(design.draw(3 * m, rng))
    obs = design.obs(design.draw(E.N_OBS, rng), law)
    return design, trial, obs


class InterfaceMatchesTheExperiment(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m, cls.rep = 100, 2
        cls.rows = {(r["shift"], r["estimand"], r.get("learner", ""), r["method"]): r
                    for r in E.replicate("scm1", cls.m, cls.rep, shifts=("same", "strong_shift"))}

    def test_ate_under_both_covariate_laws(self):
        for law in ("same", "strong_shift"):
            design, trial, obs = draw("scm1", self.m, self.rep, law)
            out = dfppi.fuse_ate(trial, obs, E_TRIAL, covariate_shift=law != "same", shuffle=False)
            fusion, rct = self.rows[(law, "ATE", "", "fusion")], self.rows[(law, "ATE", "", "rct_only")]
            self.assertAlmostEqual(out["estimate"] - design.theta, fusion["value"], places=10)
            self.assertAlmostEqual(out["se"], fusion["se"], places=10)
            self.assertAlmostEqual(out["trial_only"]["estimate"] - design.theta, rct["value"], places=10)
            json.dumps(out)                                                # serialisable as it is

    def test_cate_of_both_learners(self):
        design, trial, obs = draw("scm1", self.m, self.rep, "strong_shift")
        for learner, name in (("DR", "DRF"), ("R", "RF")):
            fit = dfppi.fuse_cate(trial, obs, E_TRIAL, learner=learner, shuffle=False)
            risk = np.mean((fit.predict(design.truth_x) - design.truth_tau) ** 2)
            risk_0 = np.mean((fit.predict_trial_only(design.truth_x) - design.truth_tau) ** 2)
            self.assertAlmostEqual(risk, self.rows[("strong_shift", "CATE", name, "fusion")]["value"], places=9)
            self.assertAlmostEqual(risk_0, self.rows[("strong_shift", "CATE", name, "rct_only")]["value"], places=9)
            json.dumps(fit.summary())


class InputsAndCommandLine(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _, cls.trial, cls.obs = draw("scm1", 100, 0, "strong_shift")

    def test_bad_inputs_are_refused(self):
        with self.assertRaises(ValueError):
            dfppi.fuse_ate({"x": self.trial["x"], "a": self.trial["a"]}, self.obs, E_TRIAL)
        with self.assertRaises(ValueError):
            dfppi.fuse_ate(self.trial, self.obs, 1.0)
        with self.assertRaises(ValueError):
            dfppi.fuse_ate(self.trial, {**self.obs, "x": self.obs["x"][:, :3]}, E_TRIAL)

    def test_a_wrong_trial_propensity_is_flagged(self):
        out = dfppi.fuse_ate(self.trial, self.obs, 0.7)
        self.assertTrue(any("trial_propensity" in w for w in out["warnings"]))

    def test_command_line_writes_json(self):
        with tempfile.TemporaryDirectory() as folder:
            paths = {}
            for name, sample in (("trial", self.trial), ("obs", self.obs)):
                paths[name] = Path(folder) / f"{name}.csv"
                p = sample["x"].shape[1]
                rows = np.column_stack([sample["x"], sample["a"], sample["y"]])
                np.savetxt(paths[name], rows, delimiter=",", comments="",
                           header=",".join([f"x{j}" for j in range(p)] + ["A", "Y"]))
            for args in (["ate"], ["cate", "--learner", "R", "--out", str(Path(folder) / "cate.csv")]):
                buffer = io.StringIO()
                with contextlib.redirect_stdout(buffer):
                    code = main(args + ["--trial", str(paths["trial"]), "--obs", str(paths["obs"]),
                                        "--trial-propensity", str(E_TRIAL)])
                result = json.loads(buffer.getvalue())
                self.assertEqual((code, result["status"]), (0, "ok"))
            self.assertEqual(len((Path(folder) / "cate.csv").read_text().splitlines()), len(self.trial["y"]) + 1)


if __name__ == "__main__":
    unittest.main()
