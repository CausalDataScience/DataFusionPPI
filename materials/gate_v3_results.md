# DataFusionPPI v3 Agile Gate Results

- Protocol: `DataFusionPPI-v3-agile-gate-2026-09-10`
- Stage: `G0`
- Decision: `INCONCLUSIVE`
- Scheduled performance replications: 8
- Hard assertions: 32 total, 0 failed
- Runtime seconds: 86.772687
- Code SHA-256: `4d0cd7c87d3a23cf0b4449e5f69d94fcad5f43094f4805461a35698b96756e55`
- Handoff SHA-256: `e1b2207beac31094f18856ed74391e1dc1f3bb665efbbf986f3f98231e950b39`
- Git commit: `7d7cdcdaff2287bf14c8c5a7c6d6fb200713649f`

## Interpretation

G0 is an engineering and algebraic smoke test. Aggregate population variance, oracle ordering, coverage, and stochastic three-MCSE checks are not evaluable at G0; G1/G2 remain unimplemented and locked. Gaussian ratios are analytically common-support and unbounded; their Theorem 6 radius is N/A. The neural path is an empirical finite-grid candidate, not a Theorem 7 check.

## SCM-B contract

`X in [-1,1]^5`, and `U, epsilon` are independent Uniform[-1,1]. The RCT propensity and outcome are exactly those in v3. OBS treatment uses `expit(logit(e_R(X)) + c U)` with `c=0,1,2`; potential outcomes are shared.
