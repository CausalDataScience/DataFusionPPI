"""ATE-only acceptance checks; no CATE execution or artifacts."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from fusion_core import ate_estimate, evaluation_variance, stable_seed, var
from exp_ate_star_real import one_replication
from make_ate_figures import _shift_bootstrap

M=Path(__file__).resolve().parents[1]/"materials"


def test_fixed_ate_and_variance_formula():
    z=np.array([1.,2.,4.,5.]); gr=np.array([.2,.5,.7,.9]); go=np.array([.1,.3,.6,1.1]); om=.4
    direct=z.mean()+om*(go.mean()-gr.mean())
    v=var(z-om*gr)/len(z)+om**2*var(go)/len(go)
    assert abs(ate_estimate(z,gr,go,om)-direct)<1e-12
    assert abs(evaluation_variance(z,gr,go,om)-v)<1e-12


def test_existing_artifact_cardinalities_and_failures():
    a=pd.read_csv(M/"ate1_when_fusion_helps_replications.csv")
    keys=["family","axis","n_rct","n_obs","confounding","spec"]
    assert len(a[keys].drop_duplicates())==44 and a.replication.nunique()==100
    assert int(a.failed.fillna(0).sum())==0
    s=pd.read_csv(M/"ate2_shift_replications.csv")
    assert len(s)==42000 and int(s.failed.fillna(0).sum())==0
    c=pd.read_csv(M/"conditional_variance_400_replications.csv")
    assert len(c)==15*400*7 and int(c.failed.fillna(0).sum())==0
    z=pd.read_csv(M/"conditional_variance_400_summary.csv")
    assert len(z)==105 and float(z.z_gap.abs().max())<=3
    q=pd.read_csv(M/"cate2_star_real_ate_replications.csv")
    assert len(q)==8*100*10 and int(q.failed.fillna(0).sum())==0


def test_reporting_identity():
    d=pd.read_csv(M/"ate1_when_fusion_helps_reporting.csv")
    assert float(d.mse_identity_gap.abs().max())<1e-12
    assert (d.effective_replications==100).all()


def test_shift_statistic_is_mean_absolute_cell_bias_and_deterministic():
    d=pd.read_csv(M/"ate2_shift_replications.csv"); d=d[(d.panel=="synthetic")&(d.estimator=="joint")]
    routes=["oracle","classifier","balanced_x","balanced_x_ghat","unit"]
    for shift in ["mild","strong"]:
        p1,c1=_shift_bootstrap(d,shift,routes,100,7); p2,c2=_shift_bootstrap(d,shift,routes,100,7)
        assert np.array_equal(p1,p2) and np.array_equal(c1,c2)
        expected=[]
        for route in routes:
            cell_bias=(d[(d["shift"]==shift)&(d.ratio==route)]
                       .groupby(["family","n_rct"],sort=True).error.mean())
            expected.append(cell_bias.abs().mean())
        assert np.allclose(p1,np.asarray(expected),rtol=0,atol=1e-15)


def test_star_driver_seed_estimand_schema_parity():
    import fusion_data as fd
    import exp_ate_star_real as driver
    cohort=fd.star_cohort("mathk","pooled")
    captured={}
    production=driver.ate_replication
    def capture(*args, **kwargs):
        captured["seed"]=args[6]
        return production(*args, **kwargs)
    driver.ate_replication=capture
    try:
        rows=one_replication("mathk",200,.4,0,cohort)
    finally:
        driver.ate_replication=production
    old=pd.read_csv(M/"cate2_star_real_ate_replications.csv")
    old=old[(old.outcome=="mathk")&(old.n_rct==200)&(old.alpha==.4)&(old.replication==0)]
    assert len(rows)==len(old)==10
    assert captured["seed"]==stable_seed("cate2","mathk",200,.4,0)
    assert all(abs(r["law_weighted_reference"]-cohort.ate_reference)<1e-15 for r in rows)
    assert set(old.estimator)=={r["estimator"] for r in rows}
    common=["estimate","error","variance_hat","covered","lambda_used","omega_used"]
    new=pd.DataFrame(rows).sort_values("estimator").reset_index(drop=True)
    old=old.sort_values("estimator").reset_index(drop=True)
    assert np.allclose(new[common],old[common],equal_nan=True,atol=1e-12,rtol=0)
