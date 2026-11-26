# sim_designs

Simulated data sets (no reference values), used as inputs for regression pins
on the package's own values and for the cross-checks against pbkrtest and the
simstudy `lotline()` reference. Copied unchanged from the `reintervals-simstudy`
repository (the source of record) at commit **d3d8a13**; copied on 2026-11-07
(`stab_unbal_v2_a.csv`: 2026-11-08) and the hashes re-checked against the
simstudy copies on 2026-11-09. Do not edit these copies; re-copy them if the
source changes.

| File | Content |
|---|---|
| `oneway_bal.csv` | One-way, 6 groups x 3 (seed 20261105) |
| `oneway_unbal.csv` | One-way, n_i = 3, 3, 2, 4, 1, 3 |
| `nested_bal.csv` | Two-fold nested, 4 A x 3 B x 2 |
| `stab_bal.csv` | Stability, 8 lots, months 0, 3, 6, 9, 12, 18, 24 |
| `stab_unbal.csv` | Same lots, staggered |
| `stab_unbal_v2_a.csv` | Stability, 8 lots, months 0-24, staggered, n = 50 |

Used by `test-ci-sim-designs.R` (random-intercept CI: package pins and the
pbkrtest Kenward-Roger cross-check, oracle `pbkrtest-kr-sim-fits`) and
`test-ngm-arms.R` (new-group-mean pins; `stab_unbal_v2_a.csv` for the
`simstudy-lotline-f4ce68e` lot-line pins).

Models (lmer):

| Tag | lmer | Data |
|---|---|---|
| OW_BAL, OW_UNBAL | `y ~ 1 + (1 \| grp)` | `oneway_bal.csv`, `oneway_unbal.csv` |
| NEST_BAL | `y ~ 1 + (1 \| a) + (1 \| a:b)` | `nested_bal.csv` |
| FS_BAL, FS_UNBAL | `y ~ month + (1 \| lot)` | `stab_bal.csv`, `stab_unbal.csv` |
| RS_UN_BAL | `y ~ month + (1 + month \| lot)` | `stab_bal.csv` |
| RS_VC_BAL | `y ~ month + (1 \| lot) + (0 + month \| lot)` | `stab_bal.csv` |

Evaluation points: none (intercept) for the one-way and nested fits; months 0,
12, 24, 36, 48 for the stability fits.

SHA-256 (also in `../SHA256SUMS`):

```
99280663f260a140cf68d9635a41be9edfa8b2a061ce5edddd22418ebe0fa2a0  nested_bal.csv
2cc9ba2d61ce6d0404c391a53e884f8bc399237f3fa65dd1152b6ec444ab6fb1  oneway_bal.csv
ac0c3e15760b3c283756ce66e431837dc5d4e47719195f4f3ba432a13135d0c2  oneway_unbal.csv
e72570899e1d6ad7f1275a3111d8abe3920d022f26ba2a44270787c48be491a4  stab_bal.csv
1d014a0ba182a647bc1a32a26e550d1b24a1e1ee27f55dfa6948d287acee2b16  stab_unbal.csv
f8c30a75eead51bd3e470691c4d64693de508cc2e5fd1ae02183946e09b4a7aa  stab_unbal_v2_a.csv
```
