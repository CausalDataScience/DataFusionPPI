#!/usr/bin/env python3
"""Unit tests for the deterministic DataFusionPPI v3 Agile gate."""

from __future__ import annotations

import contextlib
import io
import json
import inspect
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd

import v3_gate as gate
from confirm import m0_freeze as confirm_m0
from confirm import m1_fixtures as confirm_m1
from confirm import m1_production_tie as production_tie
from confirm import m1_protocol as protocol
from confirm import provenance as confirm_pv


def valid_rows(stage: str = "g0") -> list[dict[str, object]]:
    rows=[]
    for study,cell,rep,estimator,ratio,metric in gate.expected_schedule(stage):
        not_applicable=metric=="balance_residual" and ratio not in ("BAL-X","BAL-X+g")
        gate.row(rows,stage,study,cell,rep,estimator,ratio,metric,
                 "N/A" if not_applicable else 0.0,"not_applicable" if not_applicable else "ok")
    return rows


def toy_star_frame() -> pd.DataFrame:
    rows=[]
    row_id=1
    for school in (1,2,3):
        for i in range(20):
            treated=i<8
            rows.append({"rownames":row_id,"stark":"small" if treated else "regular",
                         "gender":"female" if school%2 else "male","ethnicity":"afam" if school==1 else "cauc",
                         "lunchk":"free","schoolk":"rural","schoolidk":school,
                         "mathk":450+school*10+i+5*treated,"readk":440+school*8+i})
            row_id+=1
    return pd.DataFrame(rows)


class GateV3Tests(unittest.TestCase):
    def test_named_seed_namespaces_are_stable_and_distinct(self) -> None:
        first = [gate.seed32("T", "C", 0, role) for role in gate.ROLE_NAMES]
        second = [gate.seed32("T", "C", 0, role) for role in gate.ROLE_NAMES]
        self.assertEqual(first, second)
        self.assertEqual(len(first), len(set(first)))
        draws=[gate.sample_scmb(20,"RCT",gate.rng("ROLE","CELL",0,role),1.) for role in ("R_nuis","R_tune","R_select","R_report")]
        self.assertEqual(len({gate.fingerprint(d) for d in draws}),len(draws))

    def test_fixed_algebra_assertions(self) -> None:
        assertions = gate.fixed_assertions()
        hard = [a for a in assertions if a["hard"]]
        self.assertTrue(hard)
        self.assertTrue(all(a["status"] == "pass" for a in hard))
        self.assertTrue(all(a["status"] == "not_evaluable_at_g0" for a in assertions if a["study"] == "AGGREGATE"))

    def test_scmb_exact_contract(self) -> None:
        r = gate.sample_scmb(500, "RCT", gate.rng("T", "B", 0, "R_nuis"), 1.0)
        o0 = gate.sample_scmb(500, "OBS", gate.rng("T", "B", 0, "O_nuis"), 0.0)
        self.assertEqual(r["x"].shape, (500, 5))
        self.assertGreaterEqual(float(r["x"].min()), -1.0)
        self.assertLessEqual(float(r["x"].max()), 1.0)
        np.testing.assert_allclose(o0["e"], 0.5 + 0.2 * o0["x"][:, 0], atol=1e-14)
        self.assertTrue(np.all(np.abs(r["y"]) < 1.0))

    def test_algorithm_one_projection_and_threshold(self) -> None:
        z0 = np.array([0.0, 1.0, 2.0, 3.0])
        z1 = z0 + np.array([1.0, -1.0, 1.0, -1.0])
        g = np.array([0.0, 0.5, 1.0, 1.5])
        go = np.array([0.2, -0.1, 0.4, 0.3, 0.1])
        lam, omega, moments = gate.algo1(z0, z1, g, go, 4, 5)
        self.assertTrue(0 <= lam <= 1 and 0 <= omega <= 1)
        self.assertEqual(set(moments), {"A", "B", "C", "D"})
        lam2, _, _ = gate.algo1(z0, z0.copy(), g, go, 4, 5)
        self.assertEqual(lam2, 0.0)

    def test_spline_normal_equation(self) -> None:
        rr = gate.sample_scmb(80, "RCT", gate.rng("T", "S", 0, "R_tune"), 1.0)
        oo = gate.sample_scmb(120, "OBS", gate.rng("T", "S", 0, "O_tune"), 1.0)
        mur, muo = gate.fit_outcome(rr), gate.fit_outcome(oo)
        z0, gr = gate.pseudo(rr, mur, True)
        z1, _ = gate.pseudo(rr, muo, True)
        _, go = gate.pseudo(oo, muo, True)
        basis = gate.spline_basis(rr["x"])
        beta, residual, _ = gate.solve_drf(basis, rr, oo, z0, z1, gr, go, .5, .5, .1)
        self.assertLessEqual(residual, 1e-8)
        br, bo = basis.transform(rr["x"]), basis.transform(oo["x"])
        h = .5 * br.T @ br / len(br) + .5 * bo.T @ bo / len(bo)
        rhs = br.T @ (z0 + .5 * (z1-z0) - .5*gr) / len(br) + .5 * bo.T @ go / len(bo)
        np.testing.assert_allclose(beta, np.linalg.solve(h + .1*np.eye(h.shape[0]), rhs), rtol=1e-8, atol=1e-10)

    def test_fusion_g_comes_from_obs_fit_on_both_sources(self) -> None:
        data_r = {"x": np.array([[0.], [1.]]), "a": np.array([0., 1.]),
                  "y": np.array([0., 1.]), "e": np.array([.5, .5])}
        data_o = {"x": np.array([[2.], [3.]]), "a": np.array([1., 0.]),
                  "y": np.array([2., 1.]), "e": np.array([.5, .5])}
        mur = gate.NeuralOutcomeFit(lambda x: np.column_stack((np.zeros(len(x)), np.ones(len(x)))))
        muo = gate.NeuralOutcomeFit(lambda x: np.column_stack((2*x[:, 0], 3+4*x[:, 0])))
        _, _, gr, go = gate.fusion_scores(data_r, data_o, mur, muo)
        np.testing.assert_allclose(gr, 3 + 2*data_r["x"][:, 0])
        np.testing.assert_allclose(go, 3 + 2*data_o["x"][:, 0])

    def test_artifacts_are_complete_and_protocol_locked(self) -> None:
        assertions = gate.fixed_assertions()
        rows = valid_rows()
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            gate.write_artifacts(out, rows, assertions, "g0", 0.0)
            gate.validate(out, "g0")
            meta = json.loads((out / "gate_v3_assertions.json").read_text())
            self.assertEqual(meta["decision"], "INCONCLUSIVE")
            meta["protocol"] = "alien"
            (out / "gate_v3_assertions.json").write_text(json.dumps(meta))
            with self.assertRaises(RuntimeError):
                gate.write_artifacts(out, rows, assertions, "g0", 0.0)

    def test_g1_is_explicitly_not_implemented(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaisesRegex(RuntimeError,"G1_NOT_IMPLEMENTED"):
                gate.run("g1",Path(td))
            with self.assertRaisesRegex(RuntimeError,"G1_NOT_IMPLEMENTED"):
                gate.validate(Path(td),"g1")

    def test_validator_rejects_code_hash_schedule_and_failure(self) -> None:
        assertions=gate.fixed_assertions()
        with tempfile.TemporaryDirectory() as td:
            out=Path(td);gate.write_artifacts(out,valid_rows(),assertions,"g0",0.)
            meta_path=out/"gate_v3_assertions.json";meta=json.loads(meta_path.read_text());meta["code_sha256"]="0"*64;meta_path.write_text(json.dumps(meta))
            with self.assertRaisesRegex(RuntimeError,"code hash"):gate.validate(out,"g0")
        with tempfile.TemporaryDirectory() as td:
            out=Path(td);gate.write_artifacts(out,valid_rows(),assertions,"g0",0.)
            frame=pd.read_csv(out/"gate_v3_replications.csv",keep_default_na=False);frame.iloc[1:].to_csv(out/"gate_v3_replications.csv",index=False)
            with self.assertRaisesRegex(RuntimeError,"schedule/cardinality"):gate.validate(out,"g0")
        with tempfile.TemporaryDirectory() as td:
            out=Path(td);gate.write_artifacts(out,valid_rows(),assertions,"g0",0.)
            frame=pd.read_csv(out/"gate_v3_replications.csv",keep_default_na=False);frame.loc[0,"status"]="failure";frame.to_csv(out/"gate_v3_replications.csv",index=False)
            with self.assertRaisesRegex(RuntimeError,"REDESIGN"):gate.validate(out,"g0")

    def test_lbfgs_success_and_failure_status(self) -> None:
        result=gate.minimize_lbfgs(lambda x:(float(np.sum((x-2)**2)),2*(x-2)),np.array([8.,-3.]))
        self.assertTrue(result.success);np.testing.assert_allclose(result.x,[2.,2.],atol=1e-8)
        failed=gate.minimize_lbfgs(lambda x:(1.0,np.ones_like(x)),np.zeros(2),max_iter=3)
        self.assertFalse(failed.success);self.assertEqual(failed.status,"line_search_failed")

    def test_balance_protocol_status_is_separate_from_optimizer_status(self) -> None:
        scm=gate.generate_scm_parameters(1)
        cov=np.diag(scm.covariate_noise_sd**2)+np.outer(scm.latent_loading,scm.latent_loading)
        inv=np.linalg.inv(cov);direction=scm.rct_mean_shift;d2=np.log(2.)
        delta=direction*np.sqrt(d2/float(direction@inv@direction))
        rn=gate.sample_scm(scm,40,"RCT",gate.rng("ATE2","mild",0,"R_nuis"),1.,delta)
        on=gate.sample_scm(scm,3000,"OBS",gate.rng("ATE2","mild",0,"O_nuis"),1.,delta)
        ratios=gate.fit_ratios(rn,on,gate.fit_outcome(on),scm,delta,"ATE2","mild",0)
        balanced=[r for r in ratios if r.name.startswith("BAL-")]
        self.assertTrue(all(r.protocol_status=="converged_to_protocol_tolerance" for r in balanced))
        self.assertTrue(all(r.residual is not None and r.residual<=gate.BALANCE_TOLERANCE for r in balanced))
        self.assertTrue(any(not r.optimizer_success for r in balanced))
        self.assertTrue(all(r.optimizer_status!="N/A" for r in balanced))
        self.assertTrue(all(r.residual is None and r.protocol_status=="N/A" for r in ratios[:3]))

    def test_normalization_is_diagnostic_and_nonbalance_status_is_na(self) -> None:
        assertions=gate.fixed_assertions();rows=valid_rows()
        for item in rows:
            if item["metric"]=="normalization_error":item["value"]=1e6
        with tempfile.TemporaryDirectory() as td:
            out=Path(td);gate.write_artifacts(out,rows,assertions,"g0",0.);gate.validate(out,"g0")
            frame=pd.read_csv(out/"gate_v3_replications.csv",keep_default_na=False)
            nonbal=frame[(frame.metric=="balance_residual") & ~frame.ratio.isin(("BAL-X","BAL-X+g"))]
            self.assertTrue(((nonbal.value=="N/A") & (nonbal.status=="not_applicable")).all())

    def test_cate2_freeze_fingerprints_and_role_streams(self) -> None:
        h,meta=gate.prepare_star_frame(toy_star_frame(),{"sha256":"fixture"})
        for key in ("ordered_row_id_sha256","fold_assignment_sha256","probability_table_sha256","propensity_table_sha256","frozen_support_sha256"):
            self.assertEqual(len(meta[key]),64)
        roles=[gate.draw_star_role(h,"RCT",50,gate.rng("CATE2_SMOKE","toy",0,role)) for role in ("R_nuis","R_tune","R_select","R_report")]
        self.assertEqual(len({gate.fingerprint(x) for x in roles}),4)
        source=inspect.getsource(gate.run)
        self.assertLess(source.index("prepare_star_source"),source.index("for rep in reps"))

    def test_import_is_plotting_free_and_legacy_parity_is_frozen(self) -> None:
        self.assertNotIn("matplotlib",sys.modules)
        self.assertNotIn("ssem_ate_pilot",sys.modules)
        self.assertNotIn("total_budget_nested_benchmark",sys.modules)
        parity=[a for a in gate.fixed_assertions() if a["name"]=="internalized_legacy_parity"]
        self.assertEqual(len(parity),1);self.assertEqual(parity[0]["status"],"pass")

    def test_lexicographic_tie_break(self) -> None:
        scores = np.array([2.0, 1.0, 1.0, 3.0])
        selected = int(np.lexsort((np.arange(len(scores)), scores))[0])
        self.assertEqual(selected, 1)


class ConfirmInfrastructureTests(unittest.TestCase):
    """Infrastructure-only checks.  No M0/M1 normal run or performance run."""

    def test_protocol_schedule_and_seed_ledger_are_exact(self) -> None:
        units = protocol.scheduled_units()
        self.assertEqual(len(units), 26)
        self.assertEqual(len(set(units)), 26)
        self.assertEqual(sum("/base" in unit for unit in units), 4)
        self.assertEqual(sum("/mutation/" in unit for unit in units), 16)
        self.assertEqual(sum(unit.startswith("rf_cross_term/") for unit in units), 4)
        self.assertIn("quadrant/L0O0/mutation/missing_projection", units)
        ledger = protocol.SeedLedger().payload()
        self.assertEqual([row["unit"] for row in ledger["units"]], list(units))
        self.assertEqual(len({row["entropy_hex"] for row in ledger["units"]}), 26)
        self.assertEqual(len({row["seed32"] for row in ledger["units"]}), 26)
        self.assertTrue(all(len(row["entropy_hex"]) == 32 for row in ledger["units"]))

    def test_production_tie_requires_injected_sampler_and_bootstrap_seeds(self) -> None:
        sig = inspect.signature(production_tie.run)
        self.assertEqual(list(sig.parameters), ["sampler_seed", "bootstrap_seed"])
        self.assertTrue(all(p.default is inspect.Parameter.empty for p in sig.parameters.values()))
        build = inspect.signature(production_tie.build_case)
        cross = inspect.signature(production_tie.cross_term_with_error)
        self.assertIs(build.parameters["seed"].default, inspect.Parameter.empty)
        self.assertIs(cross.parameters["seed"].default, inspect.Parameter.empty)
        source = inspect.getsource(production_tie.run)
        self.assertIn("build_case(sampler_seed)", source)
        self.assertIn("seed=bootstrap_seed", source)
        self.assertNotIn("seed=0", source)

    def test_preflight_is_calculation_free_and_write_free(self) -> None:
        before = sorted(p.name for p in confirm_pv.CONFIRM.iterdir())
        for module in (confirm_m0, confirm_m1):
            with mock.patch.object(
                    module, "self_test_promotion" if module is confirm_m0 else "check_projection",
                    side_effect=AssertionError("scientific calculation reached")):
                out = io.StringIO()
                with contextlib.redirect_stdout(out):
                    self.assertEqual(module.main(["--preflight"]), 0)
                payload = json.loads(out.getvalue())
                self.assertFalse(payload["calculation_performed"])
                self.assertFalse(payload["write_performed"])
                self.assertEqual(payload["schedule_hash"], protocol.schedule_hash())
                self.assertIn(payload["collision_state"],
                              {"absent", "same_fingerprint", "different_fingerprint", "incomplete"})
        self.assertEqual(before, sorted(p.name for p in confirm_pv.CONFIRM.iterdir()))

    def test_collision_states_cover_absent_incomplete_same_and_different(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            target = root / "stage"
            self.assertEqual(confirm_pv.collision_state(target, "abc"), "absent")
            target.mkdir()
            self.assertEqual(confirm_pv.collision_state(target, "abc"), "incomplete")
            (target / "FINGERPRINT").write_text("abc\n")
            self.assertEqual(confirm_pv.collision_state(target, "abc"), "same_fingerprint")
            self.assertEqual(confirm_pv.collision_state(target, "xyz"), "different_fingerprint")

    def test_scratch_lifecycle_cleans_only_owned_directory(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            sentinel = root / ".probe.staging.preexisting"
            sentinel.mkdir()
            (sentinel / "keep").write_text("keep")
            stage = confirm_pv.Stage("probe", "fp", root=root)
            with stage.scratch() as scratch:
                owned = scratch
                with self.assertRaises(confirm_pv.StagePromotionError):
                    stage.open()
            self.assertFalse(owned.exists())
            self.assertTrue((sentinel / "keep").is_file())

            stage = confirm_pv.Stage("probe", "fp", root=root)
            with self.assertRaisesRegex(RuntimeError, "writer"):
                with stage.scratch() as scratch:
                    owned = scratch
                    raise RuntimeError("writer failure")
            self.assertFalse(owned.exists())
            self.assertTrue(sentinel.exists())

            final = root / "probe"
            final.mkdir()
            (final / "FINGERPRINT").write_text("other\n")
            stage = confirm_pv.Stage("probe", "fp", root=root)
            with self.assertRaises(confirm_pv.StagePromotionError):
                with stage.scratch() as scratch:
                    owned = scratch
                    stage.promote()
            self.assertFalse(owned.exists())
            self.assertTrue(sentinel.exists())

    def test_successful_promotion_and_required_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            stage = confirm_pv.Stage("probe", "fp", root=root)
            with stage.scratch() as scratch:
                for name in protocol.REQUIRED_OUTPUTS:
                    if name != "FINGERPRINT":
                        (scratch / name).write_text("fixture\n")
                final = stage.promote(required_outputs=protocol.REQUIRED_OUTPUTS)
            self.assertEqual(set(p.name for p in final.iterdir()), set(protocol.REQUIRED_OUTPUTS))
            self.assertEqual((final / "FINGERPRINT").read_text(), "fp\n")

    def test_required_output_mismatch_and_fingerprint_drift_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            stage = confirm_pv.Stage("probe", "fp", root=root)
            with self.assertRaisesRegex(confirm_pv.StagePromotionError, "required output mismatch"):
                with stage.scratch() as scratch:
                    owned = scratch
                    (scratch / "checks.json").write_text("{}\n")
                    stage.promote(required_outputs=protocol.REQUIRED_OUTPUTS)
            self.assertFalse(owned.exists())
        with self.assertRaisesRegex(RuntimeError, "source fingerprint drift"):
            confirm_pv.require_code_fingerprint("0" * 64)

    def test_check_names_and_record_scope_are_frozen(self) -> None:
        names = protocol.expected_check_names()
        self.assertEqual(len(names), 53)
        self.assertEqual(len(set(names)), 53)
        self.assertIn("L0O0/sign_flip_null_invariance", names)
        self.assertIn("L1O0/mutation_sign_detected", names)
        record = confirm_m1.record_text(53)
        for phrase in ("DRF with an exact score", "fixed basis", "exact transport",
                       "population first-order", "constant propensity with all other laws fixed",
                       "null invariance, not fault detection"):
            self.assertIn(phrase, record)
        self.assertNotIn("fault detected under L0", record)


if __name__ == "__main__":
    unittest.main()
