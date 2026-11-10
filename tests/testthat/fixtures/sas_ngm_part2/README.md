# Unbalanced stability fixture for the fv-anova reference pins (B2)

`stab_unbal_v2_a.csv` is copied unchanged from `reintervals-simstudy` at commit
**d3d8a13**, `oracles/sas_ngm/part2/stab_unbal_v2_a.csv` (8 lots, months 0-24,
staggered, n = 50). The simstudy repository is the source of record; re-copy if
it changes. Copied on 2026-11-08 (commit `8808cbc`); hash re-checked against the
simstudy copy on 2026-11-09 and also in `../SHA256SUMS`. The SAS part-2 results do not exist yet, so the file is used only as
data for the lot-line pins in `test-ngm-arms.R` (reference: `lotline()` of
simstudy `pilot/methods.R`, commit f4ce68e).

SHA-256 at the time of copying:

```
f8c30a75eead51bd3e470691c4d64693de508cc2e5fd1ae02183946e09b4a7aa  stab_unbal_v2_a.csv
```
