# M2 status: pipeline connected, two blockers open

2026-09-12.  Execution PRD stages 0 and 1 pass.  Stage 2 and stage 3 are
partially built and M2 does not yet promote.

## Stage 0, read-only preflight: PASS

Infrastructure tests 8 of 8.  Both preflights report
`calculation_performed=false`, `write_performed=false`, every required input
present, collision absent, qualification ready.

## Stage 1, current-code M0 and M1 reproduction: PASS

Ten criteria checked mechanically:

| criterion | result |
| --- | --- |
| 26-unit ledger exact | 26 units, names equal the frozen protocol |
| seed collisions | 0, all 26 seed32 values distinct |
| M1 checks | 53 of 53 pass |
| production against independent formula | 7.99e-15 against a 1e-10 tolerance |
| M0 required outputs | exactly 8, none missing, none extra |
| M1 required outputs | exactly 8, none missing, none extra |
| fingerprint stable across the run | both stages carry 253e787cc034 |

Promoted to `M0_freeze-253e787cc034` and `M1_fixtures-253e787cc034`.

## Stage 2 and 3: 12 of 13 checks pass, one blocker

Six modules now exist under `code/confirm/`: the frozen protocol, the
data-generating process, the four-role engine, the runner, and the two M1
modules.  The M2 run connects data generation, four independent role draws,
nuisance fitting, the sieve candidate grid, honest selection, independent
reporting, five ratio routes, the adaptive oracle, theorem eligibility and
atomic promotion.

Passing: scheduled equals success plus failure, no duplicate seed, forty distinct
role fingerprints, exact raw schema over 951 rows, finite synthetic true risk,
no STAR true-risk row, five ratio routes exact, balance residual 7.3e-10,
non-negative selection regret, theorem-7 rows all unclipped at rho zero, runtime
and code hash on every row.

Failing: `m2/star/mathk/a0.8`.  The frozen STAR source law of blocker P0-4 is not
implemented.  The design raises rather than falling through to the synthetic
sampler, which would have emitted uniform covariates under a STAR label.

## Two findings that changed the design

**Confounding is not the lambda knob.**  The first design used the confounding
coefficient to close the lambda channel.  It does not.  The trial propensity is
known, so a biased observational outcome regression leaves the pseudo-outcome
unbiased and costs only variance, and with an observational nuisance sample of
2,750 rows against a trial one of 80 it still beats the trial regression.  At a
confounding coefficient of 2.5 the measured `lambda*` was 0.96, wide open.  The
knob is now the size of the trial nuisance draw, which decides how much room the
observational increment has left to improve on.

**The omega channel's Monte Carlo error needed the exact identity, not more
draws.**  `omega*` would not resolve at one hundred thousand draws, with a
standard error of 0.12 against a tolerance of 0.02.  The cause is that `D_p` was
estimated through the pseudo-outcome, whose variance is large.  Because the trial
propensity is known, `E(Z_0 | X) = tau` exactly, so

    D_p = tr[Gamma^-1 E_R{b b' (tau - zeta_p)(ghat - zeta_p)}]

and the oracle may use tau.  Substituting it removes the pseudo-outcome variance
entirely and `omega*` now resolves at twenty thousand draws in two of the four
cells.

The plan's rule against an arbitrary absolute cutoff also forced a change.
Convergence is now judged on `lambda*` and `omega*`, which live on the unit
interval, rather than on an absolute tolerance that a near-zero moment can never
meet.

## Open, in the order they block

1. **STAR source law, blocker P0-4.**  Needs the exact learner tuple, five-fold
   deterministic out-of-fold pooling, exact-cell propensity inside [0.15, 0.85],
   a law frozen before any replication, and the four hashes.  This is the only
   check standing between the current state and an M2 promotion.
2. **Quadrant calibration.**  Two of the four cells land where they should.
   `L1O0` measures `lambda* = 1.000` and `omega* = 0.000` as intended.  `L0O0`
   and `L0O1` put `lambda*` close enough to the boundary that one hundred
   thousand draws cannot separate it, so they report `MC_INCONCLUSIVE` rather
   than a wrong label, which is the behaviour the plan asks for.  The knob needs
   a setting that puts these cells away from the boundary before M3.
3. **Analysis and plotting.**  `cate_analyze.py` and `cate_plot.py` belong to
   stages 9 and 10 and are not built.

Nothing has been promoted for M2.  No manuscript, mirror, commit or push.
