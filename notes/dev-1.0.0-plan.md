# reintervals 1.0.0 — development plan (branch `dev`)

Status: planning only. No engine code has been changed. Line numbers are at HEAD of
`dev` = `d94d441` (R code identical to `a73e8b7`, the 0.2.0 CRAN submission).
`main` stays frozen at the CRAN submission.

Housekeeping before step 1 of the implementation: `notes/` is not in `.Rbuildignore`
(current entries: `.github`, `README.Rmd`, `paper.*`, `CONTRIBUTING.md`,
`CODE_OF_CONDUCT.md`, `cran-comments.md`, `*.html`, `.DS_Store`, `data-raw`,
`LICENSE.md`, `.claude`, `CRAN-SUBMISSION`). Add `^notes$`, or `R CMD check` will
NOTE a non-standard top-level directory.

---

## 0. Map of the engines (what the items touch)

| Path | Entry | Variance machinery | df | Method tag |
|---|---|---|---|---|
| RI CI / PI | `.coverage_verb` intervals-api.R:366 → `.interval_lmm` slope-intervals.R:192 → `intervals_lmm` intervals-lmm.R:45 | `.fixed` vc-extract.R:164 (`stats::vcov`), `.theta_cov` vc-extract.R:73 | CI `.ci_ddf` vc-extract.R:185; PI `.theta_cov` vc-extract.R:111 | `ems-mls` |
| RI TI (EMS) | `ti_lmm` dispatch.R:452 → `.ti_ems` ti-lmm.R:39 → `ti_francq` ti-francq.R:58 | `.design_components_fit` vc-extract.R:284 (k_j, EMS_j) | stratum df vc-extract.R:254 | `ems-mls` |
| RI TI (HK / Montes) | dispatch.R:422-429, 444-447 → `.montes_core` ti-montes.R:41 | unweighted lot-mean ANOVA | s, r | `anova-mls` |
| RI NGM | slope-intervals.R:227-229 → `.ngm_ri` slope-intervals.R:105 | `.montes_core` attrs `ngm_W`/`ngm_df` (2a/2b) or EMS (2c-2e) | see `.ngm_ri` | `anova-mls` / `ems-mls` |
| RS TI | `ti_lmm` dispatch.R:384-396 → `.ti_slope_path` dispatch.R:184 | `vfun_extract` vfun-extract.R:266 + `ti_vfun` ti-vfun.R:42, or `.gpq_core` gpq.R:52 | `vfun_extract` resolve_df vfun-extract.R:313-323 | `reml-mls` / `gpq` |
| RS CI / PI / NGM | slope-intervals.R:203-224 → `vfun_extract` → `.vfun_interval` slope-intervals.R:30 | `vfun_extract` | dfs F / T / G | `reml-mls` |

Note that the **tag `reml-mls` currently names two different things**: the
random-slope closed-form TI (`ti_vfun`, ti-vfun.R:96) and the random-slope
CI / PI / NGM t-intervals (`.vfun_interval`, slope-intervals.R:70). Only the first is
an MLS construction. Decided: the CI / PI / NGM t-intervals become `"reml"`, and
`"reml-mls"` remains only as the deprecated/defunct TI name (A1.3b). A1 does not
depend on A2: phase 0 validated one-sided GPQ at the 0.2.0 machinery with M = 20000
(A1.2).

---

## A1. `reml-mls` (random-slope closed-form TI)

### A1.1 Every code path that computes it or routes to it

Computation:
- `ti_vfun()` ti-vfun.R:42-101: the whole function. Graybill-Wang bound
  ti-vfun.R:76-79; two-sided half-width ti-vfun.R:86; one-sided Hoffman form
  ti-vfun.R:88; tag ti-vfun.R:96.
- Its inputs from `vfun_extract()` vfun-extract.R:266-337 (see A1.3 for which
  ones only the TI uses).

Routing:
- `ti_lmm()` method mapping `` `reml-mls` = "closed" `` dispatch.R:385-391;
  `match.arg` list dispatch.R:372; random-intercept rejection of `reml-mls`
  dispatch.R:414-418.
- `.choose_slope_engine()` dispatch.R:148-162: `closed` branch dispatch.R:157;
  `auto` → `closed` for future-observation, non-singular fits dispatch.R:160;
  closed + `true_value` error dispatch.R:149-154.
- **One-sided override** dispatch.R:212-214: `one_sided_closed` forces `closed` on
  every fit, singular included. This is the 0.2.0 one-sided random-slope TI default.
- Per-t0 df-unavailable rerouting dispatch.R:220-233 (it builds `comps` only when
  `chosen == "closed"`, and reads `df_fallback[c("G","e")]`).
- Closed-form loop dispatch.R:247-249.
- `route_reason` values dispatch.R:229-233 (`"none"` means the closed form ran).

Docs and user-facing text that name it:
- dispatch.R:5-14, 204-211, 279-292, 329-337; intervals-api.R:199, 211-213,
  226, 333-334, 365, 368; reintervals-package.R:30, 51-54, 72;
  reintervals-models.R:58-62; README.Rmd:101, 113-114, 121;
  vignettes/choosing-an-interval.Rmd:39, 47, 64-65, 138, 157-160, 173;
  vignettes/getting-started.Rmd:140-142; ti-vfun.R roxygen (`@noRd`, internal).

Tests that exercise it:
- test-dispatch.R:55 (routing matrix), :65 (non-singular → reml-mls), :86 (method
  override, refuse reml-mls + batch-mean), :163 (df-unavailable → gpq), :193
  (route_reason), :223 (indefinite Cov(phi) → gpq).
- test-stratum-sides.R:94-120 (one-sided → reml-mls on every fit, singular
  included).
- test-slope-intervals.R:140 (Graybill-Wang re-derivation), :167 (reml vs GPQ
  cross-method anchor).
- test-slope-oracle.R:18-41 (`TI = ti_vfun(...)` in the `build` list,
  slope-oracle.R:29).
- test-intervals-api.R:103 (core column contract on every engine), :130.

### A1.2 One-sided default: resolved by phase 0 (no longer blocked by A2)

The 0.2.0 one-sided random-slope TI default is `reml-mls` on every fit
(dispatch.R:212-214). The reason given was that one-sided GPQ under-covered near the
boundary in W3 (dispatch.R:206-209: worst cell 0.928).

The phase-0 study (`~/dev/reintervals-simstudy/phase0/REPORT.md`, simstudy commit
d20b954, pre-registration `phase0/PREREG.md`) re-ran one-sided GPQ in the W3 regime
against reintervals 0.2.0 @ a73e8b7. Its pre-registered outcome is **A: GPQ
passes**:
- GPQ at M = 20000 passes all 6 side × target combinations (lower / two.sided /
  upper × obs / bm, 65 cells each). The worst coverage is 0.9405.
- The W3 0.928 cell (t0 = 12, P = 0.90, lower) reproduced exactly on the W3 seeds
  (928/1000; the W3 gate passed 8/8). On 5000 fresh replicates it covered
  0.9512-0.9520 at every M (10000 / 20000 / 50000).
- Seed coupling, M and singularity were each ruled out as causes; the 0.928 does
  not reproduce on fresh replicates.

Consequence: under A1, one-sided random-slope TIs default to GPQ (`auto` → `gpq`
for every random-slope TI, both targets, all sides). They run at the 0.2.0 GPQ
machinery with M = 20000, which is what phase 0 validated. A2 is a hardening pass
(CRN, fixed seed, MC error, classed errors). It is not a precondition: A1 can land
before or after A2. If A2 lands first, re-check that its draw-order change does not
move phase-0 coverage, by running the phase-0 harness against the A2 build. Only the
default `M` matters for A1: `lmm_interval_control(M = )` is 10000 in 0.2.0
(intervals-api.R:34), so A1 must raise it to 20000 (shared with A2 step 4), or else
document that the validated configuration is M = 20000.

### A1.3 What `vfun_extract` must keep for the random-slope CI / PI / NGM

`.vfun_interval` (slope-intervals.R:45-56) reads `mean`, `var_mean`, `V_G`, `V_T`,
`dfs[c("F","T","G")]`, `t0`, `type`, `singular`, `boundary`, and the
`vg_floored` attribute.

| `vfun_extract` piece | Lines | CI/PI/NGM need it? | TI-only? |
|---|---|---|---|
| `.vfun_parse` (classify, Z, phi, make_Sigma/make_V, `reml_nll`, `vF_fun`) | vfun-extract.R:65-201 | yes (all of it) | — |
| `.vfun_vc_cov` (Hessian Cov(phi), PD check, boundary flag) | vfun-extract.R:210-222 | yes (all df come from it) | — |
| `mu`, `V_F` | :289-290 | yes | — |
| `V_G` + `.floor_vg` | :298-299, :344-353 | yes (PI, NGM) | — |
| `V_T` | :300 | yes (PI) | — |
| `content_var`, `n_E` | :301-302 | **no** | yes (ti-vfun.R:74, 86) |
| `target` argument and `attr(,"target")` | :267, :270, :332 | **no** (CI/PI/NGM always pass `"observable"`, slope-intervals.R:219) | yes (ti-vfun.R:56-61; dispatch.R:222) |
| `rF` (nu_F) | :319 | yes (CI) | — |
| `rT` (nu_T) | :320 | yes (PI) | — |
| `rG` (nu_G) | :321 | yes (NGM) | also TI |
| `re` (nu_e) | :322 | **no** | yes (ti-vfun.R:73) |
| `attr(,"df_fallback")` | :334-335 | not read today (see A5: silent) | read by dispatch.R:226 |
| `.vfun_df_containment` / between / residual | :49-53 | yes (fallback of F/T/G) | `.vfun_df_residual` only matters for nu_e and the w<0.5 branch of nu_T |

`vfun_components()` (vfun-components.R:96-104) requires `n_E` and `dfs` in its
validator. If the TI is removed, `n_E` and `dfs["e"]` become dead fields, but
dropping them changes the container contract and every fixture in
`slope_oracle.rds`. **Recommendation: keep both fields under either option**;
removing them saves little and costs a fixture regeneration.

### A1.3b Tag rename (decided; applies under either option)

The random-slope CI / PI / NGM t-intervals (`.vfun_interval`, tag at
slope-intervals.R:70) become **method `"reml"`**. `"reml-mls"` survives only as the
name of the deprecated/defunct random-slope closed-form TI (`ti_vfun`,
ti-vfun.R:96).

Sites to change:
- the tag at slope-intervals.R:70;
- the coverage-verb `method` check at intervals-api.R:368 (`ok <- "reml"` for
  slope models; `method = "reml-mls"` on a slope CI/PI/NGM becomes an error naming
  `"reml"`);
- docs at intervals-api.R:199, 211-213, 226, 333-334, 365; reintervals-package.R:27-31,
  51-54; README.Rmd:101, 113-114; vignettes/choosing-an-interval.Rmd:39, 47,
  64-65, 173; vignettes/getting-started.Rmd:140-142;
- tests at test-intervals-api.R:68 (provenance) and :103 (column contract),
  test-slope-intervals.R:77, and any `expect_identical(.., "reml-mls")` on
  CI/PI/NGM rows.

This is a user-visible output change (the `method` column), so it needs a NEWS
entry.

### A1.4 Option (a): defunct

Change:
- dispatch.R:385-391: `` `reml-mls` `` → a classed `reintervals_error_defunct`
  saying why (the failure regimes; see A1.6) and pointing to `method = "gpq"`.
- dispatch.R:148-162: delete the `closed` engine; `auto` → `gpq` for every
  random-slope TI.
- dispatch.R:212-214: delete the `one_sided_closed` override (A1.2: one-sided goes
  to GPQ).
- dispatch.R:220-233, 245-249: delete the per-t0 closed/df_bad split; the slope
  path becomes "GPQ for all t0". `route_reason` collapses to
  `group_mean_target` / `method` / a new `default_gpq`, or is dropped.
- Delete ti-vfun.R (101 lines). Delete `.gpq_rows`/closed-form column alignment
  dispatch.R:250-252 if only GPQ rows remain.
- Tests: rewrite test-dispatch.R:55/65/86/163/193/223 and
  test-stratum-sides.R:94-120; delete test-slope-intervals.R:140-183; drop
  `TI` from test-slope-oracle.R:25-30 (and the TI rows of `slope_oracle.rds`).
- Docs: every site listed in A1.1.

What it leaves in the code: all of `vfun_extract` except that `n_E`/`nu_e`/`target`
lose their only consumer; `.vfun_interval` with its new `"reml"` tag (A1.3b).
`"reml-mls"` remains only as the string recognised by `match.arg` (dispatch.R:372)
so that the defunct error can fire.

Cost: about -150 lines of R, about 10 test blocks rewritten or removed, all doc
tables. Risk: **medium-high**. It changes the default numbers for every
random-slope TI (two-sided non-singular and all one-sided) to GPQ, which phase 0
validated (A1.2). It also removes the only way to reproduce a 0.2.0 random-slope
closed-form TI.

### A1.5 Option (b): deprecated

Change:
- dispatch.R:385-391: `reml-mls` still maps to `closed`, but first emits a
  `reintervals_warning_deprecated` naming the failure regimes.
- dispatch.R:158-160: `auto` never returns `closed`.
- dispatch.R:212-214: delete the `auto` one-sided override (one-sided goes to GPQ,
  A1.2); the forced `reml-mls` one-sided limit (ti-vfun.R:88) stays reachable.
- dispatch.R:224-233: the df_bad rerouting stays (it only applies when the user
  forces `reml-mls`; at present it requires `engine == "auto"`, so it would never
  fire. Decide whether a forced `reml-mls` should still drop to GPQ on df failure;
  recommend **no**: a forced method runs that method, and the warning plus a
  `df_fallback` diagnostic say what happened).
- ti-vfun.R is kept unchanged; add the deprecation note to its header.

What it leaves in the code: everything in A1.1, with `ti_vfun` reachable only
through `method = "reml-mls"`. The `.choose_slope_engine` `closed` branch, the
`closed` + `true_value` error, and the `df_fallback["e"]` read all stay. Tests
move from "auto routes here" to "forced, warns, still computes"; the
test-slope-intervals.R:140 Graybill-Wang check and the slope-oracle TI rows stay
as they are.

Cost: about 30 changed lines plus docs. Risk: **medium**. Default numbers change
exactly as in (a); the old interval stays reproducible for audit (useful for GxP
users who validated against 0.2.0). Removal moves to 1.1/2.0.

**Decision: PENDING (Percy).** Options (a) and (b) have the same default
behaviour; (b) also keeps a reproducible path for 0.2.0 results. The earlier
recommendation for (b) stands as input, not as a decision.

### A1.6 Failure regimes for the error/warning text

Source: the decision memo (v4), random-slope closed-form TI, observation target.
- **One-sided:** coverage 0.916-0.934 in the residual-dominant regime at an
  extrapolated t0.
- **Two-sided, when forced** (`method = "reml-mls"`): coverage 0.84-0.93.

Message sketch (either class):
"method = \"reml-mls\" (the REML closed-form MLS tolerance interval for
random-slope models) under-covers: one-sided 0.916-0.934 when the residual
variance dominates and t0 is extrapolated; two-sided 0.84-0.93 where it was
forced. Use method = \"gpq\" (the default)."

Already known from the code, for completeness: the batch-mean target was never
served by the closed form (non-singular breakdown about 0.89, dispatch.R:10-12).
Before writing the final message text, cite the memo by path and version in the
NEWS entry.

---

## A2. GPQ hardening (gpq.R)

Current state:
- Draws happen **inside** the per-t0 `lapply` (gpq.R:97, draws at gpq.R:113-115),
  so each t0 gets a different pivot set, and the draw stream depends on how many
  t0 values come before it.
- One call = one `content`, one `sides`, one `c_resid` (gpq.R:52-54).
- `seed = NULL` = ambient RNG (gpq.R:56-61; intervals-api.R:15-16, 34).
- `M` default 10000 in `lmm_interval_control` (intervals-api.R:17, 34) and
  `.gpq_core` (gpq.R:54).
- No Monte-Carlo error is reported. The audit fields are `tau2_med`, `s2e_med`,
  `frac_boundary` (gpq.R:154-156).
- The < 2-distinct-time-points check exists (gpq.R:83-84), but it is an unclassed
  `stop()` without `call. = FALSE` that fires inside the batch loop, and only via
  `qr(ZtZ)$rank < 2`. Same for the residual-df stop at gpq.R:70-71.

Plan:
1. **Common random numbers.** Move `Ustar`, `Qstar`, `Zstar` (gpq.R:113-115) out of
   the t0 loop: draw them once per call, in a fixed order (U, Q, Z), right after the
   data checks (after gpq.R:71) and before any t0-dependent work. Every t0 then
   reuses the same pivots. `Rs2e` (gpq.R:116) is t0-free and moves out too.
2. **One pivot set for all P, sides, targets.** None of the draws depend on
   `content`, `sides` or `c_resid`, so with (1) plus a fixed seed, two calls that
   differ only in P/sides/target already share draws. To make that true *inside*
   one call, split `.gpq_core` into
   - `.gpq_pivots(y, time, batch, M, seed)` → per-batch OLS pieces (gpq.R:73-91) +
     `Ustar/Qstar/Zstar/Rs2e`;
   - `.gpq_eval(piv, t0, content, sides, c_resid, conf)` → the gpq.R:98-159 body.
   The dispatcher and `ti_gpq_raw` call pivots once, eval many times. The
   B-arm RS NGM `gpq` (section B) reuses the same `.gpq_pivots`.
3. **Fixed documented default seed.** `lmm_interval_control(seed = <constant>)`,
   with the constant defined once (e.g. `.gpq_default_seed`) and stated in
   `?lmm_interval_control` and NEWS. Keep an explicit opt-out for the ambient
   stream (`seed = NA`, since `NULL` currently means ambient; or keep `NULL` =
   ambient and change only the default; recommend the latter, it is the smaller
   API change). `withr::with_seed` (gpq.R:58-61) already restores the caller's
   stream and pins the default `RNGkind`; document that.
4. **Default M ≥ 20000** in intervals-api.R:34 and gpq.R:54; update
   intervals-api.R:17. Memory: the M×K matrix `g` (gpq.R:124) is 20000×K doubles
   (K = 30 → 4.8 MB), so fine.
5. **Monte-Carlo error.** For each reported quantile limit, report a
   distribution-free MC interval from the order statistics: the
   `conf ± z_{0.975} sqrt(conf(1-conf)/M)` quantiles of the same draws (two-sided:
   of `Rzeta`, gpq.R:136-138; one-sided: of `Rdelta`, gpq.R:142-148), converted to
   a half-width on the limit scale. Add `mc_se_lower`, `mc_se_upper` (or
   `mc_halfwidth`) to the row and to `.value_cols`/the diagnostics list
   (intervals-api.R:114-118, 160-163). Optional: a `reintervals_warning_mc` when
   MC half-width / interval half-width exceeds a documented threshold.
6. **Classed error for any lot with < 2 distinct time points**, checked once,
   up-front, before any draw: `length(unique(time[batch == k])) < 2` for any k →
   `reintervals_error_data` subclass `reintervals_error_gpq_time_points`, naming
   the lots. Replaces gpq.R:83-84. Also give the dispatcher a pre-check, so that
   `auto` on such data errors with the same class *before* computing, rather than
   only when a t0 happens to route to GPQ (dispatch.R:236-244).

Affected functions: `.gpq_core`, `ti_gpq_raw` (gpq.R:262-272), `.ti_slope_path`
(dispatch.R:236-244), `.gpq_rows` (dispatch.R:167-179), `lmm_interval_control`.

Tests affected: test-gpq.R:38, :46 (ambient-seed reproducibility, which changes
meaning under a fixed default seed); test-dispatch.R:102 (seed honoured);
test-gpq.R:21 (uses explicit seeds, OK); test-slope-intervals.R:167 (explicit seed
and M; numbers move by MC noise only, 2% bound should hold); every GPQ number in
the vignettes (choosing-an-interval.Rmd:183-186; getting-started.Rmd:95-122)
changes because the draw order changes. New tests: CRN (the same draws give
identical `Rs2e` across t0; monotone/consistent limits across P); the
time-point error class; MC fields present and shrinking like 1/sqrt(M).

Risk: **medium**. All GPQ outputs change numerically (expected); no
method change.

---

## A3. Kenward-Roger df paired with the unadjusted vcov

The bug, with evidence:
- `.ci_ddf()` vc-extract.R:189-201 takes the KR **df** from
  `pbkrtest::get_Lb_ddf(model, Lrow)` (vc-extract.R:195).
- The **variance** the CI uses is `var_fix` from `.fixed()` vc-extract.R:168-172,
  i.e. `diag(l %*% stats::vcov(model) %*% t(l))`: the unadjusted REML vcov.
- `intervals_lmm` puts that `var_fix` in the container (intervals-lmm.R:64-68:
  `var_mean = fx$var_fix[i]`) with the KR df (intervals-lmm.R:72-75).
- `ci_mean` combines them: ci-mean.R:52-53, `qt(.., df) * sqrt(comp$var_mean)`.

So the KR CI uses a KR df with the naive SE. KR proper pairs the KR df with the
**adjusted** covariance `Phi_A`. The SAS oracle in hand shows the gap (DDFM=KR
reports the adjusted SE):
`OW_UNBAL` intercept SE 0.51495 (SATTERTH) vs **0.52108** (KR), with identical df
5.1503; `FS_UNBAL` month = 12 mean SE 0.27586 vs **0.27619**
(`~/dev/reintervals-simstudy/oracles/sas_ngm/part1/ngm_sas_check_satterth.results.txt`,
`##FIXD`/`##MEAN` lines).

Plan:
1. In `intervals_lmm` (intervals-lmm.R:53-75), when `ci_method == "kenward-roger"`:
   compute `Vadj <- pbkrtest::vcovAdj(model)` **once** per call (today
   `get_Lb_ddf` recomputes it per row inside the `lapply`, intervals-lmm.R:72-73,
   which is O(rows) KR fits), `V0 <- stats::vcov(model)`; per row,
   `ddf = pbkrtest::Lb_ddf(L_i, V0, Vadj)` and `var_mean = L_i Vadj L_i'`.
2. Move the KR branch of `.ci_ddf` (vc-extract.R:189-202) into a helper that
   returns `list(df, var_mean, method)`, keeping the `requireNamespace` guard
   (vc-extract.R:190-194, mocked by test-ci-df.R:75-88 through the
   `requireNamespace <- NULL` binding vc-extract.R:213) and the no-silent-fallback
   error (vc-extract.R:197-200).
3. `make_comp` (intervals-lmm.R:64-68) takes the per-row `var_mean` for CI rows.
   PI rows keep `fx$var_fix` (KR applies to the CI only, unless a B-arm says
   otherwise).
4. `ci_mean` (ci-mean.R:36-66) needs no code change; update its doc (ci-mean.R:6-9)
   to say the variance comes with the df.
5. Also report `se` in the diagnostics for the RI CI (not emitted today) so the
   adjusted SE can be audited.

Tests: test-ci-df.R:62-73 checks the df only (`expect_equal(diag_of(k)$df,
pbkrtest::get_Lb_ddf(fm, 1))`), so it passes both before and after the fix, which
is why the bug went unnoticed. (Original plan: assert the KR SE and limits against
the SAS `##MEAN ... KR` lines on all five fits. Superseded by decision (a) below.)

**Implemented on `dev`, with a finding and a decision (Percy: option a).**
Fix: `.kr_parts()` computes `vcovAdj()` once per call; `.ci_ddf()` returns the KR
df from `pbkrtest::Lb_ddf(l, vcov, vcovAdj)` (= `get_Lb_ddf`, which cannot take a
precomputed Vadj) together with `var_mean = l Vadj l'`; `intervals_lmm` uses that
`var_mean` for the KR CI rows.

Finding: against the SAS `##MEAN` lines (fixtures `tests/testthat/fixtures/sas_ngm/`),
Satterthwaite matches on all 13 rows (worst SE 1.9e-6 relative, df 3.9e-5 absolute).
SAS DDFM=KR, however, reports the **Satterthwaite df unchanged** in every row, with
only the SE adjusted. On OW_UNBAL its adjusted SE also differs from pbkrtest's:
SE 0.5217893 (reintervals = pbkrtest) vs 0.5210774 (SAS), 1.4e-3 relative; df 4.783
vs 5.150. The cause of the SE difference has not been established.

Decision (a): `"kenward-roger"` in reintervals means Kenward & Roger (1997) as
implemented in pbkrtest (`vcovAdj()` SE, `Lb_ddf()` df). SAS DDFM=KR is not the KR
oracle. Tests (`test-oracle-sas-ci.R`): Satterthwaite vs SAS on every row; KR SE vs
SAS on the balanced fits only (OW_BAL, NEST_BAL, FS_BAL); KR SE and df vs pbkrtest
on all five fits at 1e-8 (the pairing guard); the OW_UNBAL difference is pinned as a
documented SAS difference. Documented in roxygen, the vignettes, and NEWS (Bug
fixes).

Dependency: pbkrtest stays in Suggests (DESCRIPTION:37). `vcovAdj`/`Lb_ddf`
are exported by pbkrtest; check the minimum version that exports `Lb_ddf` and
pin it in Suggests. Risk: **low**. Changes only the opt-in KR numbers; the default
path is untouched.

---

## A4. Graybill-Wang non-negativity guard

**Status: DONE on `dev`.** `R/mls.R` `.mls_upper()`: Graybill-Wang when every
k_j >= 0, otherwise Ting-Burdick-Graybill (H_q on the positive terms, G_r on the
negative ones, H_qr over every (positive, negative) pair; for one term of each sign
it reproduces simstudy `pilot/methods.R` `mls_upper()` exactly). `ti_francq()` uses
it; `mls_bound` is reported in the diagnostics for every MLS TI (`ems-mls`
computed; `anova-mls` and `reml-mls` are always `"graybill-wang"`, since their
coefficients are non-negative by construction). The crossed NGM Gv df falls back to
`.ddf_containment()` with `df_rule = "containment"` when pA or the df is undefined.

Deviations / findings:
- **No negative k_j found from the synthesis.** Searches: all 625 2×2 crossed tables
  with n_ij in 1..5; 3000+ random crossed tables (A, B in 2..6, small and large
  counts); 3000 random nested designs. The smallest k was **-2.8e-16 round-off**
  (the exact value is 0, e.g. counts [[3,1],[4,1]] give k_AB = -3.6e-17). The
  analysis agrees for nested designs: k_A = 1/e_AA > 0 always. So the guard uses a
  scale-relative tolerance, `k < -1e-12 * max|k|`; a bare `k < 0` would have switched
  those designs to TBG on round-off.
- Test (ii) therefore builds the negative coefficient at the **container** level
  (`re_components` with `k_AB = -0.3`), not from a design. The crossed-NGM fallback
  is tested by mocking `.design_components_fit` to zero k_A and k_B (pA = 0/0).
- With several negative terms the TBG form is the pairwise generalisation; it has no
  oracle (the reference has one negative term). It is unreachable today.
- The design-synthesis.R:57-69 code itself is unchanged: the guard lives where the
  k_j meet the bound (`ti_francq`). The synthesis still only enforces sum(k) = 1.


Where the k_j enter the MLS bound:
- `ti_francq()` ti-francq.R:88 (`sigma2_t <- sum(ems$k * ems$ms)`), ti-francq.R:96-97
  (`h <- df/qchisq(a, df) - 1`; `inner <- sum(h^2 * k^2 * ms^2)`), used in
  ti-francq.R:103-109 for both sides. `h` is the Graybill-Wang upper factor H_j,
  which is valid only for **non-negative** coefficients. With k_j < 0 the
  `k_j^2` silently treats a subtracted term as an added one.
- Sources of k_j:
  - design-oneway.R:23: `1/n0`, `1 - 1/n0` (≥ 0 because n0 ≥ 1);
  - design-nested.R:28: `1/(nB)`, `1/n - 1/(nB)`, `1 - 1/n` (≥ 0 for B ≥ 1, n ≥ 1);
  - design-crossed.R:32-35: `k_AB = (1 - 1/A - 1/B)/n`, ≥ 0 iff A, B ≥ 2 (true
    for any fitted crossed design, = 0 at A = B = 2);
  - **design-synthesis.R:57-69 `.solve_k_observable`**: unbalanced nested/crossed,
    `solve(t(E), 1)`. Nothing constrains the sign. This is where k_j < 0 can occur.
    Only `sum(k) = 1` is enforced (design-synthesis.R:63-67).
- Not affected: `.montes_core` ti-montes.R:110 (its coefficient `1 - 1/JH` ≥ 0 as
  JH ≥ 1); `ti_vfun` ti-vfun.R:79 (coefficients 1 on V_G and s2e).
- Same k used outside the TI: `.ngm_ri` Gv df slope-intervals.R:169-175. With a
  negative k_1 or k_2, `pA` (slope-intervals.R:169) can leave [0, 1] or divide by
  zero (NaN df, NaN limits; see A5). The Gm df (slope-intervals.R:178) does not
  use k.

Plan:
1. New pure helper `.mls_upper(ms, k, df, conf)` (new file R/mls.R, method-pure)
   that returns `U` (the upper confidence bound on `sum(k * ms)`) and the rule used:
   - all `k ≥ 0` → Graybill-Wang (current arithmetic, bit-for-bit);
   - any `k < 0` → the Ting-Burdick-Graybill (1990) MLS bound for a difference of
     positive and negative combinations (H_q on the positive terms, G_r on the
     negative terms, plus the H_qr cross terms). Pin the formula against
     Burdick & Graybill (1992) before coding; it has no oracle in the repo yet.
2. `ti_francq()` calls `.mls_upper` at ti-francq.R:95-109, so `inner` becomes
   `U - sigma2_t`. Carry the rule as a diagnostic (`mls_rule = "graybill-wang" |
   "tbg"`) through `.ti_ems` (ti-lmm.R:64-79) into `.unify_intervals`
   (intervals-api.R:160-163).
3. `.ngm_ri` Gv: guard `pA` (slope-intervals.R:169) against non-finite values or
   values outside [0, 1]; classed error, or the B-arm Satterthwaite with signed
   coefficients.
4. Tests: one synthesis case with a negative k_j (build a crossed count table with
   strong imbalance; check `any(k < 0)`); bit-identity of all current
   Graybill-Wang cases (Francq 5.2 / 5.3 / 4.4); a TBG oracle (none in hand; add
   to the oracle placeholders, section T).

Risk: **low-medium**. Bit-identical on every design with k ≥ 0. Only changes
unbalanced designs that today give a wrong (mis-signed) bound.

---

## A5. Classed conditions

**Status: DONE on `dev`.** `R/conditions.R`: `.abort()` / `.warn()` (base
`errorCondition` / `warningCondition`, no call recorded) and the registered
hierarchy. All 111 `stop()` / `warning()` sites in R/ are converted with unchanged
message text. The only bare `stop()` / `warning()` left are inside the two helpers
(guarded by a grep test in `test-conditions.R`). Names follow the step-2 brief:
`reintervals_warning_limitation` (the plan said `reintervals_known_limitation`).
Reserved, defined but not yet emitted: `reintervals_error_defunct`,
`reintervals_warning_deprecated`, `reintervals_warning_mc`.

Silent paths surfaced (inventory numbers from A5.1):
- #5 `pi_df_fallback` is now carried into `attr(., "diagnostics")`.
- #2 slope CI/PI/NGM `df_rule` = `"satterthwaite"` / `"containment"` (NA for a
  hand-built container without `df_fallback`).
- #1 the 1e6 substitution is replaced by the containment df **K - 1** (as
  specified). `vfun_extract` now attaches `attr(, "n_groups")`; a container
  without it is a `reintervals_error_numerical`. **No test value moved**: no test
  reaches this branch (the committed `slope_oracle.rds` dfs are finite, and
  `vfun_extract` always returns a finite df). Note that K - 1 is a different rule
  from the slope engine's own between-containment `B - k_re` (= K - 2 with a random
  slope, vfun-extract.R:49).
- #11 `.solve_tau2` marks draws without a solution as NA (before, an NA would have
  errored inside the indexed assignments; it was not actually let through as a
  number); `.gpq_core` drops them, reports `n_na_draws`, and raises
  `reintervals_error_numerical` above 1%.
- #12 crossed NGM NaN path: see A4.
- #18 the few-groups PI limitation: `reintervals_warning_limitation` from `pi_lmm()`
  on random-intercept models whose coarsest random factor has <= 5 levels, using the
  documented text and no ICC threshold. #17 (the single-obs fixed-slope CI/PI
  warning) is now the same class.
- Not done, as not in the step-2 brief: the fallback *warnings* proposed in A5.2
  (`reintervals_warning_fallback` on df fallbacks; the class exists and is used
  only by the truncation/floor warnings); #13, #14, #15 (newdata check on the
  slope coverage path), #19 (extrapolation).
- Found while testing: `ems_synthesis(matrix(1, 1, 1), "crossed")` fails inside
  `stats::model.matrix` with an unclassed base error (no package check runs first).
  It is unreachable from a fitted model (classify_design needs >= 2 levels), but
  it is a gap for the internal function.
- Test churn: the new limitation warning hit 4 existing tests (top factor 3-5
  levels). They use a `quiet_limitation()` helper that muffles only that class
  (helper-conditions.R), or, in covariate-routing, `expect_warning(class = )`. The
  two stray outputs are silenced at the call site (test-vc-extract.R:52
  `suppressMessages`; test-vfun-components.R:145 `expect_output`).


### A5.1 Inventory: places that can return NaN/NA limits or fall back silently

| # | Site | What happens | Visible to the user today? |
|---|---|---|---|
| 1 | slope-intervals.R:58-60 | non-finite or ≤ 0 slope df → `1e6` (normal quantile) | no (the `df` column keeps the bad value, the limits use 1e6) |
| 2 | vfun-extract.R:313-318 | per-df Satterthwaite → containment fallback | TI: via `route_reason` only; **CI/PI/NGM: silent** (`.vfun_interval` never reads `df_fallback`) |
| 3 | vfun-extract.R:214-219 | Hessian/solve/PD failure → `cov = NULL` → all df fall back | no (only through #2) |
| 4 | vfun-extract.R:49-50 | containment df floored at 1 by `max(1, .)` | no |
| 5 | vc-extract.R:98-104, 111 | singular / Hessian fail / non-PD → PI df = containment | `attr(,"pi_df_fallback")` is set (intervals-lmm.R:102), **but `.unify_intervals` does not carry it** (intervals-api.R:160-167), so it is silent through the verbs |
| 6 | vc-extract.R:203-206 | CI Satterthwaite NA/<1 → containment | `df_method = "containment"` in the diagnostics; no condition |
| 7 | vc-extract.R:274 | containment df floored at 1 | no |
| 8 | re-components.R:85-90 | negative VC truncated to 0 | unclassed `warning()` |
| 9 | vfun-extract.R:298, 344-353 | negative V_G floored to 0 | unclassed `warning()` + `vg_floored` flag |
| 10 | gpq.R:107, 125, 132 | eigenvalues, `g`, content variance clamped at 0 | no (`frac_boundary` diagnostic only) |
| 11 | gpq.R:186-188, 193, 207-214 | `.solve_tau2`: `denom ≤ 0 → NA` gives `h_lo = NA`, the draw goes to `interior`; the bracket loop gives up silently after 200 iterations | no. **Possible NA/NaN in `Rtau2`**, which then reaches the quantiles |
| 12 | slope-intervals.R:169-175 | crossed Gv `pA` with signed or zero `k_A EMS_A + k_B EMS_B` → NaN df → NaN limits | no; only `W ≤ 0` is checked (slope-intervals.R:113-117) |
| 13 | dispatch.R:201; vfun-extract.R:278 | `isSingular` error → `NA` → treated as non-singular | no |
| 14 | intervals-api.R:102 | design classification error → `grouping = NA` | no |
| 15 | slope-intervals.R:216-217 | random-slope CI/PI/NGM with `newdata` lacking the time column → `t0s = NULL` → empty `rbind` → obscure downstream error | no clear message (the TI path checks this at dispatch.R:195-199; the coverage path does not, and there is no `.check_newdata` on the slope path) |
| 16 | ti-lmm.R:84-87 | singular EMS TI | unclassed `warning()` |
| 17 | intervals-api.R:386-394 | fixed-slope single-observation lots, CI/PI | unclassed `warning()` (a known limitation) |
| 18 | intervals-api.R:278-283 | PI under-coverage with ≤ 5 top-level groups and ICC ~0.9 | **documented only, no condition** |
| 19 | dispatch.R:303-307 | slope t0 extrapolation outside the observed range | documented "without a guard or warning" |
| 20 | gpq.R:70-71, 83-84 | `stop()` without `call. = FALSE`, unclassed | error, but not classed |
| 21 | ~109 `stop()` sites across R/ (`grep -c "stop("`: vfun-components 18, re-components 16, dispatch 10, vfun-extract 9, intervals-api 8, slope-intervals 8, vc-extract 7, ti-montes 6, design-synthesis 5, ti-francq 5, ti-vfun 4, gpq 3, ci-mean 3, pi-newobs 3, ti-lmm 3, intervals-lmm 1) | all unclassed `simpleError` | message only |

### A5.2 Proposed hierarchy (base R `errorCondition`/`warningCondition`; R ≥ 4.1 per DESCRIPTION:35, so no rlang)

```
reintervals_condition
├─ reintervals_error                         (+ "error")
│  ├─ reintervals_error_input                bad args, newdata, alpha/level/P, control
│  ├─ reintervals_error_unsupported          design/method/target combination not available
│  │    ├─ reintervals_error_design          off-catalog (crossed-no-interaction, 4+ factors)
│  │    └─ reintervals_error_method          forced method incompatible with model/interval
│  ├─ reintervals_error_data                 data cannot support the interval
│  │    ├─ reintervals_error_groups          < 2 lots/groups (gpq.R:67, ti-montes.R:68)
│  │    ├─ reintervals_error_time_points     lot with < 2 distinct times (A2)
│  │    └─ reintervals_error_residual_df     (gpq.R:70, ti-montes.R:96)
│  ├─ reintervals_error_numerical            W ≤ 0, sum(k·EMS) ≤ 0, singular synthesis, NaN guard
│  ├─ reintervals_error_dependency           pbkrtest missing / KR not computable
│  └─ reintervals_error_defunct              (A1 option a)
└─ reintervals_warning                       (+ "warning")
   ├─ reintervals_warning_fallback           a df/variance rule fell back; carries `which` + `from`/`to`
   │    ├─ ..._df_fallback                   #2, #5, #6 (and #1 replaced by it)
   │    ├─ ..._truncation                    #8, #9
   │    └─ ..._clamp                         #10 (only above a frac_boundary threshold)
   ├─ reintervals_warning_singular           #16 (and, uniformly, on every verb)
   ├─ reintervals_warning_deprecated         (A1 option b)
   ├─ reintervals_warning_mc                 MC error above threshold (A2)
   └─ reintervals_known_limitation           a validated-shortfall regime, interval still returned
        ├─ ..._single_obs_fixed_slope        #17
        ├─ ..._few_groups_pi                 #18 (new: emit it)
        └─ ..._extrapolation                 #19 (new: t0 outside the time range)
```

Implementation: two internal helpers, `.abort(msg, class, ...)` and `.warn(msg,
class, ...)`, in new R/conditions.R. Each condition carries structured fields
(e.g. `which_df`, `t0`, `rule`). Every fallback is also written to the
diagnostics (`df_fallback`, `pi_df_fallback`, `mls_rule`), so callers that suppress
the warning keep the audit trail. Fix #1 (delete the 1e6 substitution, which is
unreachable once #2 guarantees a finite df) and #11/#12 (convert to a
`reintervals_error_numerical`, never a NaN limit). Add a final NaN guard in
`.unify_intervals` (intervals-api.R:147-155): any non-finite `lower`/`upper` that
is not a deliberate ±Inf → `reintervals_error_numerical`.

Tests: a `test-conditions.R` that triggers each class with `expect_error(class =
)` / `expect_warning(class = )`; update every `expect_error(..., "regex")` that
the reworded messages break.

Risk: **low** for the numbers (no interval value changes), **medium** for churn
(touches about 110 call sites; message regexes in tests). Doing it early lets every
later item use the helpers.

---

## B. Candidate arms as `lmm_interval_control()` options (plan only)

Current control: intervals-api.R:34-45 (`seed`, `M`, `ci_df`). Only `ci_df`
reaches the engines (intervals-api.R:395-397 → slope-intervals.R:192-201 →
intervals-lmm.R:45-47). **Plumbing change shared by all arms:** pass the whole
`control` object down `.coverage_verb → .interval_lmm → intervals_lmm / .ngm_ri`
rather than single fields, and write the arm that ran into the diagnostics
(`df_method`, `df_rule`, `ngm_rule`), so the step-13 simulations can read it off
every row.

Names: keep the 0.2.0 exported `ci_df` (it is on CRAN; same values), add the rest
with a scope prefix. All defaults reproduce 0.2.0 exactly until step 13 decides.

| Option | Values | Default (= 0.2.0) | Hook |
|---|---|---|---|
| `ci_df` (RI CI) | `"satterthwaite"`, `"kenward-roger"`, `"containment"` | `"satterthwaite"` (with its containment fallback) | `.ci_ddf` vc-extract.R:185-207: add a `containment` branch returning `.ddf_containment` (vc-extract.R:272) directly; `match.arg` sites intervals-api.R:35-36, slope-intervals.R:195/201, intervals-lmm.R:47/52. KR depends on A3. |
| `slope_ci_df` (RS CI) | `"hessian-satterthwaite"`, `"kenward-roger"`, `"k-minus-1"` | `"hessian-satterthwaite"` (`rF`, vfun-extract.R:319, with containment fallback) | `vfun_extract` gains a CI-df argument (vfun-extract.R:266-269, 313-319). `k-minus-1` = `parsed$B - 1` (note: today's containment uses `B - k_re`, vfun-extract.R:49; it is a different rule). KR needs `pbkrtest::vcovAdj` and the adjusted `V_F` (vfun-extract.R:290 uses `parsed$vcov_beta`, vfun-extract.R:197), so A3's helper is reused. Remove the KR rejection slope-intervals.R:204-208. |
| `pi_df` (RI PI) | `"francq-satterthwaite"`, `"containment"`, `"kenward-roger"` | `"francq-satterthwaite"` (with the containment fallback, vc-extract.R:104-111) | `.theta_cov` vc-extract.R:73-119 (or a thin wrapper choosing `df_pi`); consumed at intervals-lmm.R:65, 93. `kenward-roger` is **not defined** for a new-observation PI (KR gives df for fixed-effect contrasts); keep it as a reserved value that raises `reintervals_error_unsupported` until a definition is pinned. |
| `ngm` (RI NGM) | `"current"` (2a-2e as in `.ngm_ri`), `"satterthwaite-broad"` (SAS DDFM=SATTERTH broad-inference new-group prediction), `"hts"` (Higgins-Thompson-Spiegelhalter), `"classical-pi"` | `"current"` | `.ngm_ri` slope-intervals.R:105-186: route before the 2a/2b branch (slope-intervals.R:126) and before the EMS branch (slope-intervals.R:150). `value_rows` (slope-intervals.R:112-124) is shared by all arms. `satterthwaite-broad` reuses `.theta_cov` (Hessian, vc-extract.R:73) plus a gradient of the new-group variance (pattern: `.ci_ddf_satterthwaite` vc-extract.R:219-235). Oracle: SAS part1 `##NGM ... SATTERTH`. HTS and classical-pi formulas to be pinned in their step spec. |
| `slope_ngm` (RS NGM) | `"current"` (`mu ± t(nu_G) sqrt(V_F + V_G)`), `"gpq"` (generalized prediction interval), `"fv-anova"` (lot-line PI on per-lot fitted values) | `"current"` | slope-intervals.R:218-221 (kind `new_group_mean`). `gpq` reuses A2's `.gpq_pivots` (the theta/tau² GPQs, gpq.R:116-133), drawing a new-lot value `Rtheta + Z' sqrt(Rtau2)`; it **depends on A2**. `fv-anova` reuses the per-batch OLS fitted values `NY` (gpq.R:99-100), so factor `.batch_ols()` out of gpq.R:73-91 first. Oracle: SAS part1 `RS_*` NGM lines (RS_UN_BAL, RS_VC_BAL; RS_UN_UNBAL is singular, boundary case only, per the oracle README) and part2 once it has been run (`ngm_sas_check_part2.results.txt` does not exist yet). |

`lmm_interval_control` validation (intervals-api.R:36-44) gets one `match.arg` per
option. Applying an option to a model it does not concern (e.g. `slope_ngm` on a
random-intercept fit) is a no-op, recorded in the diagnostics; it is not an
error, so one control object can be reused across a simulation grid.

---

## C. One-sided NGM for random-intercept designs

**Already supported** in 0.2.0 for every current RI arm: `.ngm_ri` takes `sides`
(slope-intervals.R:105-106), `tcrit` switches to `t_{1-α}` (slope-intervals.R:110-111),
`value_rows` opens the right side (slope-intervals.R:120-121), and
`new_group_mean_lmm` passes `sides` through (intervals-api.R:351-360).
test-stratum-sides.R:72-92 checks one-sided 0.95 = two-sided 0.90 for the HK (2b),
fixed-slope (2a) and nested (2c) NGM.

What is still needed:
1. Fix the stale comment at slope-intervals.R:104 ("Two-sided only.").
2. Extend test-stratum-sides.R:72 to the crossed Gv (2d) and the one-way-with-covariates
   (2e) arms; neither is in that test today.
3. Every new B arm must implement `sides`: the t-based arms trivially through
   `value_rows`; the RS `gpq` NGM arm needs the lower/upper quantile pattern of
   gpq.R:141-149.
4. Coverage evidence: the one-sided NGM is the same pivot with `t_{1-α}`. Whether
   the simstudy validated it one-sided per arm is not recorded in this repo; add
   it to the step-13 grid.

---

## T. Oracle-test scaffolding

### T.1 Current self-referential tests (code checked against its own formula or output)

| Test | Lines | What is being compared |
|---|---|---|
| test-ti-francq.R | 3-23 | `ti_by_hand` restates Eq. 26 exactly as coded |
| test-ti-francq.R / test-stratum-sides.R | 102-117 / 160-180 | one-sided Hoffman form on the Francq bound, restated (labelled "provisional") |
| test-ci-mean.R | 10-22 | Eq. 12 half-width restated |
| test-pi-newobs.R | 10-32 | Eq. 21 half-width restated |
| test-slope-intervals.R | 16-42 | `.vfun_interval` vs the same formula on the same container |
| test-slope-intervals.R | 140-165 | Graybill-Wang TI restated; its dfs come from `vfun_extract` (not independently anchored) |
| test-slope-oracle.R | 18-41 | committed `slope_oracle.rds`, described as an external oracle (slope-intervals.R:17-18), but **no generator script is in the repo** (`grep -rl slope_oracle` finds only the test); it was added in the initial commit 4676ab4. Provenance unverifiable, so treat it as a regression snapshot |
| test-ngm.R | 1-200 | references "written here from its defining formula": an independent re-implementation of the package's own constructions (2b-vs-Ca on balanced data is the only external identity) |
| test-hk.R | 105-131 | HK vs ems-mls coincidence, and the core vs an in-test HK re-implementation |
| test-design-synthesis.R | 33-92 | balanced-limit, trace, and n0 consistency checks (internal) |
| test-dispatch.R | 27-44, 113-127 | ti_lmm vs the "frozen EMS engine", and M1 core vs EMS (regression/cross-engine) |
| test-vfun-extract.R | 42-55 | balanced M1: REML route = ANOVA route (cross-method) |
| test-slope-intervals.R | 167-183 | reml-mls vs GPQ agreement (cross-method) |
| test-ci-df.R | 62-73 | KR df vs `pbkrtest::get_Lb_ddf` (checks the df only, see A3) |

Tests that already hit published numbers, which move into the oracle files:
Francq 5.2: test-ci-mean.R:3, test-pi-newobs.R:3, test-ti-francq.R:25,
test-vc-extract.R:3/19, test-intervals-lmm.R:3, test-ti-lmm.R:3, test-design-oneway.R:19.
Francq 5.3: test-hk.R:90, test-vc-extract.R:56. Francq 4.4: test-design-synthesis.R:4,
test-stratum-sides.R:144, test-ti-lmm-unbalanced.R:11. MBL: test-ti-montes.R:5, :71,
test-vfun-extract.R:56. HK: test-hk.R:38-82. Oliva-Aviles & Hauser: test-gpq.R:13
(closed-form tau² = 8.439 only; there is no published GPQ interval to compare
against).

### T.2 Proposed layout

```
inst/extdata/oracles/
  README.md                       provenance table: source, table/eq, transcription date, who checked
  montes2019/                     mbl_appendix_a (already data/), tables_I_II.csv (TI, n_E, U, C1, C2, one-sided H)
  hoffman_kringle2007/            hk_assay (already data/), cases.csv (O1-O4 limits + intermediates)
  sas_ngm/part1/                  COPY of ~/dev/reintervals-simstudy/oracles/sas_ngm/part1/*.csv
                                  + results parsed to ngm_part1.csv (case, ddfm, row, pred, se, df, lower, upper)
                                  + the .results.txt kept verbatim, plus the .sas program for provenance
  sas_ngm/part2/                  placeholder until the part2 results exist
  nist_recipe_ex4/                placeholder (README only)
  francq2019_table2_corrected/    placeholder
  chow_shao1991/                  placeholder
  oliva_aviles_hauser2025/        dgp (already data/oliva_hauser_dgp) + any printed interval, if one is recoverable
tests/testthat/
  helper-oracles.R                read_oracle(name) via system.file("extdata/oracles", ...); parser for SAS ## lines
  test-oracle-montes2019.R        MBL Appendix A two-sided + Tables I-II one-sided (moved from test-ti-montes.R)
  test-oracle-hk2007.R            O1-O6 (moved from test-hk.R)
  test-oracle-francq2019.R        5.2 / 5.3 / 4.4 (moved), + corrected Table 2 placeholder (skip())
  test-oracle-sas-ngm.R           part1: ##MEAN SATTERTH vs ci_lmm (df, SE, limits);
                                  ##MEAN KR vs ci_lmm(ci_df = "kenward-roger") (A3 acceptance: SE!);
                                  ##NGM SATTERTH vs ngm = "satterthwaite-broad"; ##NGM RS_* vs slope_ngm arms;
                                  ##COVP vs lme4 VarCorr (fit agreement gate, tolerance ~1e-4)
  test-oracle-nist-recipe.R       skip("oracle not yet transcribed")
  test-oracle-chow-shao1991.R     skip(...)
  test-oracle-oliva-hauser2025.R  tau² exact (moved from test-gpq.R:13) + skip() for any interval values
```

Rules: each oracle file reads only `inst/extdata/oracles`, calls only exported
verbs (plus the method-pure constructors, where the paper prints the components),
and states its tolerance next to the printed precision. Copy the SAS fixtures (the
simstudy tree is read-only and must not be referenced by path from tests). Cross-
method and self-consistency tests stay where they are, renamed `test-consistency-*`
so the distinction is visible. Size: all CSV/txt, well under 100 KB. Check the
`.results.txt` header carries no host details (per the oracle README it does not).

---

## 9. Dependencies and proposed order

```
A5 (conditions helpers) ──┬──> A2 (GPQ: classed time-point error) ───> B slope_ngm = "gpq"
                          ├──> A1 (deprecated/defunct class; NOT blocked by A2: phase 0 outcome A)
                          ├──> A3 (KR vcovAdj; dependency error class) ──> B ci_df / slope_ci_df KR
                          ├──> A4 (mls_rule; numerical error class) ──> B ngm Gv guard
T (oracle scaffolding) ───┴── needed as the acceptance gate for A3, A4, B (SAS fixtures)
B plumbing (control object passed down) ──> all B arms ──> C (sides for new arms)
```

Proposed implementation order (each step is one reviewable commit series, with
defaults unchanged unless stated):
1. **T scaffolding + `.Rbuildignore` `^notes$`**: move the published-number tests,
   add the SAS part1 fixtures and parser. No R/ changes. (Low risk.)
2. **A5 helpers + inventory fixes**. **Done on `dev`** (see A5 status; #13-#15 and
   #19 remain).
3. **A3** KR vcovAdj. **Done on `dev`** (decision a). Gates: KR vs pbkrtest
   (SE + df, 1e-8) on all five SAS fits; KR SE vs SAS on the balanced fits only.
4. **A4** `.mls_upper` + TBG. **Done on `dev`** (no negative k_j reachable from the
   synthesis; guard with a round-off tolerance).
5. **A2** GPQ hardening (CRN, pivots/eval split, default seed, M = 20000, MC error,
   time-point class). **Changes the GPQ numbers** (NEWS).
6. **B plumbing**: control passed down, diagnostics carry the arm; then arms one at
   a time (`ci_df = "containment"`, `pi_df`, `slope_ci_df`, `ngm`, `slope_ngm`),
   each with its oracle or a reference test and a one-sided check (C).
7. ~~Simstudy one-sided GPQ re-run~~: done (phase 0, outcome A; A1.2).
8. **A1** (defunct vs deprecated: PENDING, Percy): `auto` routes every random-slope
   TI to GPQ with M = 20000, the `reml-mls` TI is defunct or deprecated, and the
   RS CI/PI/NGM tag becomes `"reml"`. **Changes the default RS TI numbers** (NEWS,
   vignettes, README tables). Can be scheduled independently of step 5 (A2); if A2
   lands first, re-run the phase-0 harness against the A2 build.
9. Step 13 (simstudy) chooses the defaults for the B arms; flip the defaults in one
   commit with NEWS.

Risk summary: A1 medium-high (default change to a validated GPQ; defunct vs
deprecated pending); A2 medium (all GPQ numbers move); A5 medium churn / zero
numeric change; A4 low-medium; A3 low (done); B low while the defaults are
unchanged; T low.
