# Confirmatory CATE stages

Every stage here is written into a scratch directory and promoted in one atomic
rename, so a partial run never appears as a complete one.  Each promoted stage
carries a `FINGERPRINT` file holding the digest of the code that produced it.

A stage is never replaced by a run under a different fingerprint.  Results from
different code states are not pooled, so a changed code state gets a new stage
rather than overwriting an old one.  Stage directories therefore carry the first
twelve characters of the fingerprint in their name.

## Historical pre-fix stages

| milestone | stage | status |
| --- | --- | --- |
| M0 evidence freeze | `M0_freeze` | historical, produced before the current promotion and schedule fixes |
| M1 deterministic fixtures | `M1_fixtures-abbf5370d61a` | historical, 53 of 53 checks passed under pre-fix code |

## Superseded

`M1_fixtures` and `M1_fixtures.superseded` hold two earlier M1 runs made while
the fixtures were still being written.  They pass 50 of 50 checks each, under a
code state that did not yet contain the production tie, the corrected role
mutation or the R-Fusion proportionality test.  They are kept rather than deleted
because the promotion rule forbids replacing a stage across fingerprints, and
they are listed here so that no one mistakes them for a current-code M1 result.

There is no current-code M0 or M1 authority yet.  A future normal run must use
the fingerprint-qualified default stage name and produce the complete protocol,
seed, code, environment, check, moment, record, and fingerprint artifact set.

## Abandoned scratch directories

The following eight directories are retained as historical debris.  They are
not promoted stages and must not be used as evidence:

- `.M1_fixtures.staging.czr6hl_4`
- `.M1_fixtures.staging.ev0c1ggf`
- `.M1_fixtures.staging.hqy6flnm`
- `.M1_fixtures.staging.jo2p8m87`
- `.M1_fixtures.staging.lb5p8t00`
- `.M1_fixtures.staging.xytiltpj`
- `.M1_fixtures.staging.y2boez92`
- `.M1_fixtures.staging.ym493i1z`

They predate instance-owned scratch cleanup.  This documentation does not
authorize deleting them or rewriting their historical contents.
