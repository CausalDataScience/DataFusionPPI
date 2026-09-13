# AISTATS Experiment Protocol: One-Page Overview

Superseded as an implementation specification on 2026-09-10.  The single
authority for the experiments is
`materials/2026-09-09-implementation-handoff-JA.md` (version 2), which carries
the formulas, the data with pinned hashes, the sized designs, the compute
budget, and the acceptance tests.  This file is kept as the short overview.

Four studies support Sections 3 and 4 of the manuscript.  Each uses 100
replications per design cell, the honest three-way sample split of
Assumption 1, and the RCT-only estimator as the reference.

| study | question it answers | results checked | data | cells |
| --- | --- | --- | --- | ---: |
| ATE-1 | When does fusion help, and does Algorithm 1 find the coefficients that help? | Theorems 1, 2, 2.1, 3, 4 | three SCM families and STAR real covariates, common covariate marginal | 44 |
| ATE-2 | What does covariate shift cost, and which density ratio controls it? | Theorem 1 drift, Definition 12 with Proposition 6, Theorem 5 with Corollary 5.1 | shifted SCM families, plus the NSW trial against the CPS and PSID comparison groups | 18 |
| CATE-1 | Does honest grid validation deliver its guarantee, and does the sieve oracle explain what it selects? | Theorems 6, 7, Corollary 7.1 | SCM families and STAR real covariates with a known conditional effect | 24 |
| CATE-2 | Does the pipeline hold up with real treatment, real outcomes, and a representation learned from the observational sample? | Theorem 6 item 1 | STAR kindergarten cohort with a constructed confounded observational sample | 8 |

Total workload is 9,400 replications, about eight hours on one core and close
to one hour on eight cores.

Decisions fixed with the author: ATE-2 includes the NSW and CPS panel, CATE-1
includes a frozen neural basis alongside additive cubic B-splines, and every
study uses 100 replications.

Two figures and two tables from these studies go into the body of the paper.
Figure 1 is the ATE root-mean-square-error ratio against the trial size,
Figure 2 is the CATE risk against the trial size, Table 1 collects interval
coverage and validation regret, and Table 2 collects the two real-data panels.
Everything else goes to the appendix.
