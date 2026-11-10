# slope_oracle.rds

**What it is.** Committed `vfun_components` containers for two random-slope fits
(`comps$M2c`: correlated random slope, t0 = 5, sleepstudy-sized; `comps$M2i`:
independent random slope, t0 = 12) and the CI / PI / new-group-mean limits and df
the method-pure constructor `.vfun_interval()` produced from them
(`oracle`). `test-slope-oracle.R` re-runs the constructor on the committed
containers and compares at 1e-9. It is a **regression snapshot of the interval
arithmetic**, not an external oracle, and is therefore not listed in
`oracles.csv`.

**Provenance.** Added in the initial commit `4676ab4` (2026-07-07), captured on
`aarch64-apple-darwin20` (field `captured_on`). No generator script exists in
this repository or in reintervals-simstudy (`grep -rl slope_oracle` finds only the
test), so the fits it came from cannot be re-derived; the source comment in
`R/slope-intervals.R` calling it an external oracle is not verifiable.

**Stripped at step T (2026-11-09, data unchanged otherwise).** Fields no longer
read after A1 (the closed-form random-slope TI `reml-mls` is defunct) were
removed with `readRDS()` / `saveRDS(version = 2)`:

- `comps$*$n_E`, `comps$*$dfs["e"]`, `attr(comps$*, "target")`;
- the four `oracle` rows with `type == "TI"`.

Kept unchanged (checked `identical()` against the original): `type`, `phi`, `cov`,
`mean`, `var_mean`, `V_G`, `V_T`, `dfs[c("F", "T", "G")]`, `t0`, `singular`,
`boundary`, the class and `attr(, "vg_floored")`, every non-TI `oracle` row, and
`captured_on`. `phi` and `cov` are not read by the constructors either; they are
kept as the record of the fit. The 18 expectations of `test-slope-oracle.R` give
`identical()` constructor outputs and reference values before and after (only the
data-frame row names of `oracle` changed, from dropping the TI rows).

SHA-256:

```
9cec845acab69fbea2dffcad9bd025b18ed9e90f6831867ff6644131d4810488  slope_oracle.rds (4676ab4, original)
4a741794f91328cbebbb0df3088c1fad3cd5b2c843abab01b6aa836b4344d7b5  slope_oracle.rds (stripped, step T)
```
