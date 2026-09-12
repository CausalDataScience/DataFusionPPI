"""Generate only the publication ATE figures from frozen result CSVs."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
import numpy as np
import pandas as pd

from fusion_io import MATERIALS

NAVY, TEAL, RUST, GREY = "#173B7A", "#0E7490", "#9A3412", "#53657A"
COLORS = {"rct_only": "#11233F", "joint": NAVY, "lambda_only": TEAL,
          "omega_only": GREY, "oracle_joint": RUST, "shrinkage": "#777777",
          "obs_transported": "#BBBBBB"}
LABELS = {"rct_only": "trial only", "joint": "joint", "lambda_only": "lambda only",
          "omega_only": "omega only", "oracle_joint": "oracle coefficients",
          "shrinkage": "shrinkage", "obs_transported": "observational AIPW"}


def _save(fig, path: Path):
    fig.tight_layout()
    fig.savefig(path, dpi=200, metadata={"Software": "DataFusionPPI"})
    plt.close(fig)


def rmse_figure(out: Path):
    s = pd.read_csv(MATERIALS / "ate1_when_fusion_helps_summary.csv")
    r = s[s.metric.eq("rmse_ratio_to_rct")]
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.1), sharey=True)
    order = ["rct_only", "joint", "lambda_only", "omega_only", "oracle_joint",
             "shrinkage", "obs_transported"]
    for est in order:
        b = r[(r.axis == "rct_size") & (r.estimator == est)].groupby("n_rct").value.mean()
        ax[0].plot(b.index, b.values, marker="o", lw=1.3, ms=3,
                   color=COLORS[est], label=LABELS[est])
    quality = [(0.0, "flexible"), (2.0, "flexible"), (0.0, "linear"),
               (1.0, "linear"), (2.0, "linear")]
    qlabels = ["c=0\nflex", "c=2\nflex", "c=0\nlinear", "c=1\nlinear", "c=2\nlinear"]
    for est in order[1:]:
        vals = [r[(r.axis == "quality") & (r.estimator == est) &
                  (r.confounding == c) & (r.spec == spec)].value.mean()
                for c, spec in quality]
        ax[1].plot(range(5), vals, marker="s", lw=1.3, ms=3, color=COLORS[est])
    ax[0].set_xscale("log")
    ax[0].xaxis.set_major_locator(FixedLocator([50,100,200,400]))
    ax[0].xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x)}"))
    ax[0].xaxis.set_minor_locator(NullLocator())
    ax[0].set_xlabel("trial size n_R"); ax[0].set_ylabel("mean cell RMSE ratio to trial only")
    ax[0].set_title("RCT-size axis, correct OBS model")
    ax[1].set_xticks(range(5)); ax[1].set_xticklabels(qlabels, fontsize=7)
    ax[1].set_xlabel("all five quality additions at n_R=100"); ax[1].set_title("prediction quality axis")
    for a in ax: a.axhline(1, color="#999", ls=":", lw=.8); a.spines[["top","right"]].set_visible(False)
    fig.legend(*ax[0].get_legend_handles_labels(), frameon=False, ncol=4, fontsize=7, loc="upper center")
    fig.tight_layout(rect=(0,0,1,.87)); fig.savefig(out/"figure1_ate_rmse.png",dpi=200,metadata={"Software":"DataFusionPPI"}); plt.close(fig)


def variance_figure(out: Path):
    fig, axes = plt.subplots(1,2,figsize=(7.2,3.0),sharex=True)
    for a,title,stem in zip(axes,["conditional on frozen fits","whole procedure refitted"],
                            ["conditional_variance_400","ate1_when_fusion_helps"]):
        d=pd.read_csv(MATERIALS/f"{stem}_variance.csv")
        for comp,color,label,off in [("omega_given_lambda",NAVY,"estimated omega",-.12),
                                     ("omega_given_lambda_oracle",RUST,"oracle omega",.12)]:
            b=d[(d.comparison==comp)&d.scope.astype(str).str.startswith("pooled:")].copy()
            b["family"]=b.scope.str.split(":").str[1]
            p=d[(d.comparison==comp)&(d.scope=="pooled")]
            b=pd.concat([b,p.assign(family="pooled")])
            y=np.arange(len(b))+off; x=100*b.reduction_rate.to_numpy(); lo=100*b.rate_low.to_numpy(); hi=100*b.rate_high.to_numpy()
            a.errorbar(x,y,xerr=[x-lo,hi-x],fmt="o",capsize=2,color=color,label=label)
        labs=list(b.family); a.set_yticks(range(len(labs))); a.set_yticklabels(labs); a.axvline(0,color="#999",ls=":")
        a.set_title(title); a.set_xlabel("variance reduction (%)\npaired bootstrap 95% CI"); a.spines[["top","right"]].set_visible(False)
    fig.legend(*axes[0].get_legend_handles_labels(),frameon=False,ncol=2,loc="upper center")
    fig.tight_layout(rect=(0,0,1,.9)); fig.savefig(out/"figure1_omega_variance_reduction.png",dpi=200,metadata={"Software":"DataFusionPPI"}); plt.close(fig)


def decomposition_figure(out: Path):
    d=pd.read_csv(MATERIALS/"ate1_when_fusion_helps_reporting.csv"); d=d[d.axis=="rct_size"].copy(); d["bias2"]=d.bias**2
    order=["rct_only","omega_only","lambda_only","joint","shrinkage","obs_transported"]
    fig,axes=plt.subplots(1,3,figsize=(7.8,2.8))
    for a,col,title in zip(axes,["empirical_variance","bias2","rmse"],["empirical variance","squared bias","RMSE"]):
        x=np.arange(len(order)); vals=d.groupby("estimator")[col].mean().reindex(order)
        a.bar(x,vals,color=[COLORS.get(e,GREY) for e in order]); a.set_xticks(x); a.set_xticklabels([LABELS[e] for e in order],rotation=35,ha="right",fontsize=7); a.set_yscale("log"); a.set_title(title); a.spines[["top","right"]].set_visible(False)
    _save(fig,out/"figure2_variance_bias_rmse.png")


def inference_figure(out: Path):
    d=pd.read_csv(MATERIALS/"ate1_when_fusion_helps_reporting.csv"); d=d[d.axis=="rct_size"]
    order=["rct_only","omega_only","lambda_only","joint","shrinkage","adaptive","naive_pool","obs_transported"]
    pooled=pd.read_csv(MATERIALS/"ate1_when_fusion_helps_pooled.csv")
    fig,a=plt.subplots(1,2,figsize=(7.2,3.0)); y=np.arange(len(order))
    vr=d.groupby("estimator").variance_ratio_hat_to_empirical.mean().reindex(order)
    a[0].barh(y,vr,color=NAVY); a[0].axvline(1,color=RUST); a[0].set_yticks(y); a[0].set_yticklabels(order,fontsize=7); a[0].set_xlabel("mean estimated / empirical variance")
    c=pooled[pooled.quantity=="covered"].set_index("estimator").reindex(order)
    a[1].barh(y,c.value,xerr=[c.value-c.low,c.high-c.value],color=TEAL); a[1].axvline(.95,color=RUST); a[1].set_yticks(y); a[1].set_yticklabels([]); a[1].set_xlabel("95% interval coverage (bootstrap 95% CI)")
    for z in a: z.spines[["top","right"]].set_visible(False)
    _save(fig,out/"figure3_variance_and_coverage.png")


def coefficient_figure(out: Path):
    d=pd.read_csv(MATERIALS/"ate1_when_fusion_helps_replications.csv"); d=d[(d.estimator=="joint")&d.oracle_lambda.notna()]
    fig,a=plt.subplots(1,2,figsize=(5.6,2.6))
    for z,x,y,label in zip(a,["oracle_lambda","oracle_omega"],["lambda_used","omega_used"],["lambda","omega"]):
        z.scatter(d[x],d[y],s=2,alpha=.12,color=NAVY); z.plot([0,1],[0,1],color=RUST); z.set_xlabel(f"oracle {label}"); z.set_ylabel(f"selected {label}"); z.spines[["top","right"]].set_visible(False)
    _save(fig,out/"figureA_coefficient_recovery.png")


def _shift_bootstrap(d, shift, routes, draws=2000, seed=20260912):
    cells=sorted(d[d["shift"]==shift][["family","n_rct"]].drop_duplicates().itertuples(index=False,name=None))
    point=[]; boot=np.empty((draws,len(routes))); rng=np.random.default_rng(seed+(shift=="strong"))
    arrays={(c,r):d[(d["shift"]==shift)&(d.family==c[0])&(d.n_rct==c[1])&(d.ratio==r)].sort_values("replication").error.to_numpy() for c in cells for r in routes}
    for r in routes: point.append(np.mean([abs(arrays[(c,r)].mean()) for c in cells]))
    for b in range(draws):
        vals={r:[] for r in routes}
        for c in cells:
            k=len(arrays[(c,routes[0])]); ix=rng.integers(0,k,k)
            for r in routes: vals[r].append(abs(arrays[(c,r)][ix].mean()))
        boot[b]=[np.mean(vals[r]) for r in routes]
    return np.array(point),np.percentile(boot,[2.5,97.5],axis=0)


def shift_figure(out: Path):
    d=pd.read_csv(MATERIALS/"ate2_shift_replications.csv"); d=d[(d.panel=="synthetic")&(d.estimator=="joint")]
    routes=["oracle","classifier","balanced_x","balanced_x_ghat","unit"]; labels=["oracle","classifier","BAL-X","BAL-X+g","ignore r=1"]
    fig,a=plt.subplots(figsize=(6.4,3.4)); x=np.arange(len(routes)); width=.36
    for j,shift in enumerate(["mild","strong"]):
        p,ci=_shift_bootstrap(d,shift,routes); pos=x+(j-.5)*width
        a.bar(pos,p,width,label=shift,color=[TEAL,RUST][j],alpha=[.55,.9][j],yerr=[p-ci[0],ci[1]-p],capsize=2)
    a.set_xticks(x); a.set_xticklabels(labels,rotation=20,ha="right"); a.set_ylabel("mean absolute cell bias (bootstrap 95% CI)"); a.set_title("ATE under covariate shift, joint estimator"); a.legend(frameon=False); a.spines[["top","right"]].set_visible(False)
    _save(fig,out/"figure4_ate_shift_routes.png")


def main():
    p=argparse.ArgumentParser(); p.add_argument("--output-dir",type=Path,default=MATERIALS); args=p.parse_args(); args.output_dir.mkdir(parents=True,exist_ok=True)
    rmse_figure(args.output_dir); variance_figure(args.output_dir); decomposition_figure(args.output_dir); inference_figure(args.output_dir); coefficient_figure(args.output_dir); shift_figure(args.output_dir)

if __name__=="__main__": main()
