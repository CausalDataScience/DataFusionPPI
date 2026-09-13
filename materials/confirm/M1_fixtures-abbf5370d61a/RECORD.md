# M1 record: deterministic mechanism micro fixtures

Stage `M1_fixtures-abbf5370d61a`, 2026-09-12.  Plan section 9, milestone M1.  These are
algebraic identities on fixed arrays, not statistical replications, so the
tolerances below are numerical rather than statistical.

## Verdict

`GO`.  53 of 53 checks pass.  Nothing here required a change to the manuscript.

## What the fixtures do

The audit's root cause 5.2 is that the two channels' opportunities were never
moved independently.  The fixtures move them with two knobs that turn out to be
algebraically independent, which the numbers below show rather than assert.

Writing `zeta_p` for the basis projection of the true conditional effect, the
two moments that decide the quadrant are

    C_p = tr[Gamma^-1 E{b b' Cov(Z_0, Delta | X)}]
    D_p = tr[Gamma^-1 E{b b' (tau - zeta_p)(ghat - zeta_p)}]

so the conditional correlation between the base pseudo-outcome and the increment
is the lambda knob, and the part of the effect that the basis cannot span is the
omega knob.  Measured on the four fixtures:

| quadrant | A_p | B_p | C_p | D_p | lambda* | omega* |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| L0O0 | 1.0042 | 0.5709 | -9.5e-18 | -3.1e-17 | 0.0000 | 0.0000 |
| L1O0 | 1.0042 | 0.5709 | -0.96402 | -3.1e-17 | 0.9600 | 0.0000 |
| L0O1 | 1.0042 | 2.1981 | -2.4e-17 | 0.65200 | 0.0000 | 0.2966 |
| L1O1 | 1.0042 | 2.1981 | -0.96402 | 0.65200 | 0.9600 | 0.2966 |

`A_p` is identical across all four, `C_p` moves with the lambda knob alone and
`D_p` with the omega knob alone.  The separation is exact, not approximate.

## Three findings worth recording

**The omega channel has no first-order opportunity when the sieve spans the
effect.**  `D_p` is the whitened covariance of `b(tau - zeta_p)` with
`b(ghat - zeta_p)`.  When the basis spans `tau`, `zeta_p` equals `tau`, the first
factor is pure noise with conditional mean zero, and `D_p` is exactly zero.  The
omega channel pays only where the sieve is misspecified.  This is visible in the
table as the exact zeros in the `O0` rows.

**The RF cross term is proportional to a known factor, not merely zero at one
half.**  Corollary 7.1 states that the lambda-omega cross term vanishes when the
trial propensity is one half.  Multiplying the weight through gives more:

    chi * Delta_RF   = -(A - e) delta_m / [e(1-e)]
    chi^2 * Delta_RF = -(A - e)^3 delta_m / [e(1-e)]^2

and `E{(A-e)^3 | X} = e(1-e)(1-2e)`, so the cross term should carry the factor
`(1-2e)/[e(1-e)]`.  Dividing the measured cross term by that factor gives
0.222910 at each of e = 0.2, 0.6 and 0.8, with a spread of 4.1e-15.  The
corollary's statement is correct and the sharper proportionality holds.

**Exchanging the two arguments of the covariance is not a role mutation.**  The
first mutation test written for `D_p` swapped the base pseudo-outcome with the
prediction and detected nothing.  It could not: `tr(Gamma^-1 C')` equals
`tr(Gamma^-1 C)` for symmetric `Gamma^-1`, so `D_p` is invariant under the swap.
The test was replaced by two substitutions, using the prediction where the base
pseudo-outcome belongs and leaving the projection out of the prediction gap, both
of which move the moments by more than 1e-6.

## Go conditions

| condition | result |
| --- | --- |
| independent formula against production, difference at most 1e-10 | 9.2e-15 once the exploratory 1e-10 ridge is matched; 2.2e-9 at ridge zero, which is the ridge itself |
| oracle coordinates and fixed-coordinate risk ordering match the design | grid minimum agrees with the closed form to the grid resolution in all four quadrants |
| deliberate sign, role and ratio mutations detected | all four quadrants detect all four mutations |

The production comparison runs at `lambda* = 0.6840` and `omega* = 0.7603`, both
strictly interior, so a sign error could not pass by clipping at a boundary.

## Two identities checked along the way

The Corollary 7.1 quadratic was checked against the definition of the first-order
risk at twelve coordinate pairs per fixture.  The largest gap is 8.3e-17, so the
quadratic is exact rather than approximate on these arrays.  The normal equation
that defines the projection holds to 3.3e-17.

## Files

`moments.json` holds the four quadrants, the RF cross-term sweep and the
production tie.  `checks.json` holds every check with its status and detail.
`FINGERPRINT` ties the stage to the code state that produced it.
