# SAS PROC MIXED oracle fixtures (NGM check, part 1)

Copied unchanged from `reintervals-simstudy` at commit **d3d8a13** ("SAS NGM oracle
fixtures: part 1 results, part 2 program and data"), directory
`oracles/sas_ngm/part1/`. The simstudy repository is the source of record. Do not
edit these copies; re-copy them if the source changes.

| File | Content |
|---|---|
| `ngm_sas_check_satterth.results.txt` | The `##` result lines of the SAS log, with the provenance header |
| `oneway_bal.csv` | One-way, 6 groups x 3 (seed 20261105) |
| `oneway_unbal.csv` | One-way, n_i = 3, 3, 2, 4, 1, 3 |
| `nested_bal.csv` | Two-fold nested, 4 A x 3 B x 2 |
| `stab_bal.csv` | Stability, 8 lots, months 0, 3, 6, 9, 12, 18, 24 |
| `stab_unbal.csv` | Same lots, staggered |

Program: `ngm_sas_check_satterth.sas` (same source directory; not copied, since the
data are inline and identical to the CSVs). Software: SAS Viya 4 (SAS V.04.00),
PROC MIXED, METHOD=REML, DDFM=SATTERTH and DDFM=KR. Run 2026-11-05.

Models, as PROC MIXED statements, with the lme4 equivalents used in
`test-oracle-sas-ci.R`:

| Tag | PROC MIXED | lmer |
|---|---|---|
| OW_BAL, OW_UNBAL | `class grp; model y = ; random grp;` | `y ~ 1 + (1 \| grp)` |
| NEST_BAL | `class a b; model y = ; random a b(a);` | `y ~ 1 + (1 \| a) + (1 \| a:b)` |
| FS_BAL, FS_UNBAL | `class lot; model y = month; random intercept / subject=lot;` | `y ~ month + (1 \| lot)` |

`##MEAN` lines are OUTPM= rows: the marginal mean x0'b with its SE and DDFM df. They
are the oracle for `ci_lmm()`.

## Documented difference: SAS DDFM=KR is not the Kenward-Roger oracle

In reintervals, `ci_df = "kenward-roger"` means Kenward & Roger (1997) as
implemented in pbkrtest: the SE comes from the adjusted covariance
`pbkrtest::vcovAdj()` and the df from `pbkrtest::Lb_ddf()`. On these fixtures, SAS
DDFM=KR differs:

- it reports the **Satterthwaite df** unchanged in every `##MEAN` row (the KR df
  equals the SATTERTH df to all printed digits) and only adjusts the SE;
- on OW_UNBAL its adjusted SE also differs from pbkrtest's:

| OW_UNBAL, intercept | reintervals = pbkrtest | SAS DDFM=KR |
|---|---|---|
| SE | 0.5217893 | 0.5210774 (1.4e-3 relative) |
| df | 4.783 | 5.150 (= SAS SATTERTH) |

The cause of the SE difference has not been established. `test-oracle-sas-ci.R`
therefore compares Kenward-Roger to SAS on the SE only, and only on the balanced
fits (OW_BAL, NEST_BAL, FS_BAL), where the KR and Satterthwaite df coincide. It
checks Kenward-Roger against pbkrtest itself on all five fits.

SHA-256 of the copies, at the time of copying:

```
99280663f260a140cf68d9635a41be9edfa8b2a061ce5edddd22418ebe0fa2a0  nested_bal.csv
ae8e856657bd2859e25d985524292487a7038caef4b4e6029ce71e773563d1f6  ngm_sas_check_satterth.results.txt
2cc9ba2d61ce6d0404c391a53e884f8bc399237f3fa65dd1152b6ec444ab6fb1  oneway_bal.csv
ac0c3e15760b3c283756ce66e431837dc5d4e47719195f4f3ba432a13135d0c2  oneway_unbal.csv
e72570899e1d6ad7f1275a3111d8abe3920d022f26ba2a44270787c48be491a4  stab_bal.csv
1d014a0ba182a647bc1a32a26e550d1b24a1e1ee27f55dfa6948d287acee2b16  stab_unbal.csv
```
