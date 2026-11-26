# reintervals (development version)

* Development towards 1.0.0.

## Breaking changes

* **`method = "reml-mls"` is defunct** for `ti_lmm()`: it signals an error of
  class `reintervals_error_defunct`. The REML closed-form MLS tolerance
  interval for random-slope models under-covers: one-sided 0.916-0.934 when the
  residual variance dominates and `t0` is extrapolated, two-sided 0.84-0.93
  where forced. Use `method = "gpq"` (the default; Oliva-Aviles & Hauser
  2025). The internal closed-form constructor (`ti_vfun()`) is removed.
* **Every random-slope tolerance interval is now GPQ**: both targets
  (`over = "observation"` and `"group_mean"`), two- and one-sided, singular
  fits included. In 0.2.0 the future-observation target used the closed form
  on non-singular fits and for every one-sided interval, so **those numbers
  change**. A simulation study (one-sided GPQ in the regime that had motivated
  the closed-form one-sided default) validated GPQ at `M = 20000` for all six
  side x target combinations (worst coverage 0.9405 for a nominal 0.95).
  `attr(., "diagnostics")$route_reason` now takes the values `"auto"` and
  `"method"` (the 0.2.0 values `"none"`, `"singular"`, `"df_unavailable"`
  and `"group_mean_target"` described the closed-form routing, which no
  longer exists).
* **Tag rename:** the random-slope confidence, prediction and new-group-mean
  intervals report `method = "reml"` (was `"reml-mls"`, the name of the
  removed tolerance interval). Their numbers are unchanged. Passing
  `method = "reml-mls"` to `ci_lmm()`, `pi_lmm()` or `new_group_mean_lmm()`
  on a random-slope model errors and names `"reml"`.
* **Default `M = 20000`** in `lmm_interval_control()` (was 10000), the GPQ
  size validated by simulation. GPQ results at the default therefore change
  (and take about twice as long).
* **GPQ numbers change (common random numbers).** The GPQ pivot draws are now
  generated once per call and reused for every time point in `newdata` (and
  both limits); before, each time point drew its own set, in turn. For a
  seeded call the first time point is unchanged and every later time point
  moves by Monte-Carlo noise (sleepstudy regression pin, `P = 0.99`,
  `M = 2000`: by up to about 5 units). With `seed = NULL`, every GPQ number changes,
  because the draws now run under a seed taken from the ambient stream (see
  below).
* **A random-slope tolerance interval now errors when any lot has fewer than
  two distinct time points** (class `reintervals_error_time_points`): GPQ needs
  a per-lot slope. In 0.2.0 the removed `reml-mls` closed form handled such data
  for the future-observation target (the batch-mean target already errored).
  The error message suggests refitting without those lots or, if the
  between-lot slope variance is negligible, a random-intercept model with a
  fixed slope. The random-slope CI, PI and new-group mean are unaffected.
* **`lmm_interval_control(M = )` below 2000 is an error** (class
  `reintervals_error_input`), as is a non-integer `M`. The Monte-Carlo
  standard error uses 20 MC sections, and below 2000 a section has fewer than
  100 draws.
* **A character time column in `newdata` is an error** for the fixed-slope
  `new_group_mean_lmm()` (class `reintervals_error_input`). It was coerced to
  numeric before (`"5"` taken as 5).

## New features

* **`ti_lmm(method = "gpq")` on random-intercept models** (step 13b; under
  validation: every call signals `reintervals_warning_unvalidated`, and the
  default `"auto"` route is unchanged). The GPQ engine of Oliva-Aviles &
  Hauser (2025, Technometrics 67(2):193-202, doi:10.1080/00401706.2024.2407324)
  now also runs on the random-intercept members of its model class (their
  sections 4.1 and 4.3): one grouping factor with a random intercept only, any
  fixed effects with `W = (X0, Z)` of full column rank, single-observation
  lots allowed, residual df `n - q - K` > 0. Two- and one-sided, for both
  targets; `over = "group_mean"` (the lot-mean distribution) is new for
  random-intercept models and exists only on this route. On the one-way model
  the pivots are those of Liao, Lin & Iyer (2005, Technometrics 47(3):323-335,
  doi:10.1198/004017005000000102), whose sulfur example is reproduced within
  Monte-Carlo error. Random-slope GPQ results are unchanged (bit-identical).
  Two or more random factors (nested, crossed) still signal
  `reintervals_error_method`. The diagnostics carry `M` on this route.
* **`lmm_interval_control(gpq_two_sided = "content")`** (validation arm,
  package-original; under validation, step 13b; signals
  `reintervals_warning_unvalidated`). A two-sided GPQ tolerance interval that
  solves the content criterion of Liao, Lin & Iyer (2005, eq. 1) draw by draw,
  in place of their Wald-Wolfowitz / Howe approximation
  `T +/- z_{(1+P)/2} R_zeta` (Oliva-Aviles & Hauser 2025, eq. 1), which stays
  the default (`"approx"`, bit-identical). Every GPQ tolerance-interval route
  (random slope, random intercept, `ti_gpq_raw()`), two-sided only; the option
  is recorded per row in `attr(., "diagnostics")$gpq_two_sided`. In a
  1000-replicate check on the 16 two-sided step-13 random-intercept cells
  below nominal (extrapolated, small between-lot variance), its coverage was
  at least nominal - 0.01 in every cell, where the default GPQ reached
  0.85-0.97; widths grew by at most 22% there and 8% in passing control cells.
* **Vectorised `P`, `level` / `alpha` and `sides`.** `ti_lmm()` accepts vectors
  of `P`, `level` (or `alpha`) and `sides`; `ci_lmm()`, `pi_lmm()`,
  `new_group_mean_lmm()` and `lmm_predict()` accept vectors of `level` (or
  `alpha`) and `sides`. Values must be unique, with no `NA`; each is validated
  as before. The result has one row per combination of evaluation point x
  `sides` x `P` x `level`, the full cross product, with the evaluation point
  varying fastest, then `sides`, then `P`, then `level` (each in the order
  given); the `P`, `level` and `sides` columns identify each row. `over` stays
  a single value. Scalar inputs give exactly the previous output. The
  default of `sides` is now the single value "two.sided" (it was the
  match.arg choice vector, which meant the same); a vector of sides gives one
  block of rows per side.
* **One GPQ draw set per call.** A vectorised GPQ call (the random-slope
  tolerance interval, both targets, and the `slope_ngm_method = "gpq"` arm)
  draws the pivots once and solves for the between-lot variance once per time
  point, then evaluates every content, confidence and side on those same
  draws; each row keeps its own Monte-Carlo SE. With an explicit `seed`, every
  row is identical to the corresponding scalar call. With `seed = NULL` a
  vectorised call draws one seed (one draw set), while separate scalar calls
  each draw their own, so the equality then holds only with an explicit seed.
  `gpq_seed`, `rng_kind` and `mc_sections` stay per call. The other engines
  compute each combination with their scalar code. Warnings (extrapolation,
  df fallback, unvalidated arm, limitation, singular fit) are signalled once
  per call. Per-call diagnostics of the scalar anova-mls engine that depend on
  the combination or the time point (`n_E`, `U`, `C1`, `C2`) are row-aligned
  vectors in a vectorised call.

* **Degrees-of-freedom validation arms** in `lmm_interval_control()`, for the
  1.0.0 validation study. The defaults are unchanged and remain the validated
  rules; every non-default value below is a validation arm whose coverage is
  under evaluation, and applying one signals `reintervals_warning_unvalidated`
  on every call (muffle it by class).
  - `ci_df = "containment"` (random-intercept CI): the containment df alone
    (the rank rule of the existing fallback) with the unadjusted standard
    error. `ci_df = "kenward-roger"` is unchanged and does not warn (a 0.2.0
    option).
  - `slope_ci_df` (random-slope CI), new: `"hessian-satterthwaite"` (default,
    the current df with its `K - 2` fallback), `"kenward-roger"` (pbkrtest
    `vcovAdj()` + `Lb_ddf()` on the random-slope fit, with the KR-adjusted
    standard error) and `"k-minus-1"` (`K - 1`, `K` lots, with the unadjusted
    standard error).
  - `pi_df` (random-intercept PI), new: `"francq-satterthwaite"` (default) and
    `"containment"` (the same variance with the containment df). There is no
    Kenward-Roger prediction interval.

  Each option applies to one model structure and verb and is otherwise a
  no-op, so one control object can be reused across models. The diagnostics of
  every CI, PI and new-group mean record `ci_df_rule`, `slope_ci_df_rule` and
  `pi_df_rule`: the arm applied, or `"not_applicable"`. Consequently
  `ci_df = "kenward-roger"` on a random-slope model is now ignored (it was an
  error); the random-slope Kenward-Roger arm is `slope_ci_df = "kenward-roger"`.
  A failed Kenward-Roger df or variance under `slope_ci_df = "kenward-roger"`
  is an error naming the arm (`reintervals_error_numerical`); a validation arm
  never falls back to another df.
* **New-group-mean validation arms** in `lmm_interval_control()`, under the
  same rules (defaults unchanged; non-default values warn with
  `reintervals_warning_unvalidated`; recorded as `ngm_arm` / `slope_ngm_arm`
  in the diagnostics; an arm not defined for the model or design is a no-op).
  - `ngm_method` (random intercept): `"package"` (default, the 0.2.0
    constructions); `"broad-inference"` (the broad-inference prediction of a
    new group, McLean, Sanders & Stroup 1991: `W = V_F` plus every
    between-group variance component, Satterthwaite df on the variance scale;
    NA limits and
    `ngm_df_failed = TRUE`, no fallback, only when the df cannot be computed
    or is non-finite or not positive; a positive df below 1 is used as is);
    `"hts"` (Higgins,
    Thompson & Spiegelhalter 2009, with the REML between-group variance;
    intercept-only one-way, `I >= 3`); `"classical"` (Hahn & Meeker 1991 on the
    group means; intercept-only one-way, balanced, where it equals the default).
  - `slope_ngm_method` (random slope): `"reml"` (default);
    `"broad-inference"` (the broad-inference prediction of a new lot: the
    default's center and variance `V_F + V_G`, with the Satterthwaite df of
    the whole `V_F + V_G` instead of `V_G` alone; NA limits and
    `ngm_df_failed = TRUE`, no fallback, only when the df cannot be computed
    or is non-finite or not positive; a positive df below 1 is used as is);
    `"gpq"` (a
    generalized prediction interval on the GPQ tolerance-interval pivots, with
    the GPQ seed contract, common random numbers across `t0` and Monte-Carlo
    SE; the tolerance-interval draws are unchanged); `"fv-anova"`
    (package-original lot-line interval on the per-lot OLS fitted values,
    `Ybar_fit +/- t_{K-1} S_L sqrt(1 + 1/K)`, conservative for the new-lot mean).
* **`slope_pi_method`** in `lmm_interval_control()`: the random-slope
  prediction interval. `"reml"` (default) is the existing interval, unchanged
  bit for bit. `"fv-anova"` (validation arm, package-original) is the
  `slope_ngm_method = "fv-anova"` construction plus the pooled within-lot
  residual variance, `Ybar_fit +/- t_df sqrt(S_L^2 (1 + 1/K) + S2_E)`, with
  the Satterthwaite df of the two mean squares (`K - 1` and `n - 2K` df).
  `"hybrid"` (under validation, step 13b) uses `"fv-anova"` on a singular fit
  (or a failed singularity check) and `"reml"` otherwise. Non-default values
  warn (`reintervals_warning_unvalidated`); the option is recorded as
  `slope_pi_arm`, and the random-slope prediction rows carry the construction
  actually used in `method` (`"reml"` or `"fv-anova"`).
  The `"fv-anova"` df uses the effective (Satterthwaite) df of `S_L^2`,
  `tr(AV)^2 / tr((AV)^2)` with the per-lot variances
  `v_k = max(VG_hat, 0) + a_k S2_E` (`a_k` the lot's leverage at `t0`),
  reported as `fv_df_lots`: `K - 1` in balanced designs (identical results),
  smaller when the lots' leverages differ. The `slope_ngm_method = "fv-anova"`
  arm keeps `K - 1`. `"hybrid"` never errors where `"reml"` computes: when the
  fv-anova interval is infeasible on a singular fit it uses `"reml"` with a
  warning of class `reintervals_warning_hybrid_fallback` (a
  `reintervals_warning_fallback`), and records why in `slope_pi_route_reason`.
* Random-slope prediction-interval rows now report `df_fallback_reason` in the
  diagnostics, per row, with the random-intercept semantics: why the default
  df fell back to containment (e.g. `"singular fit"`), `NA` where it did not.
* **GPQ reproducibility contract.** Every GPQ interval (`ti_lmm()` random-slope
  route and `ti_gpq_raw()`) records `gpq_seed`, the integer seed actually used,
  and `rng_kind`, the `RNGkind()` triple in effect, in
  `attr(., "diagnostics")`. `lmm_interval_control(seed = <gpq_seed>)` with the
  same `M` under the same `rng_kind` returns an identical interval. With the
  default `seed = NULL`, one integer seed is drawn from the caller's stream (so
  `set.seed()` before the call still controls the result; the stream advances
  by exactly one draw) and the draws run under it, restoring the caller's RNG
  state afterwards.
* **Common random numbers across time points.** Within one call, a GPQ limit
  at a time point no longer depends on which other time points were requested
  or their order.
* **Monte-Carlo standard error of each GPQ limit**, `mc_se_lower` and
  `mc_se_upper` in the diagnostics (`NA` on the infinite side of a one-sided
  interval), by 20 MC sections of consecutive draws (`mc_sections = 20`). It
  is the run-to-run jitter of the Monte-Carlo computation, not statistical
  uncertainty: a limit within approximately 2 MC SE of a specification limit
  has a seed-dependent pass/fail. At `M = 2000` the sectioning estimate ran
  about 8% below the observed run-to-run SD in a sleepstudy check (sections of
  100 draws); the gap is expected to shrink at the default `M = 20000`. It is
  not printed.
* **Per-lot time-point guard.** A random-slope GPQ tolerance interval
  (`ti_lmm()` and `ti_gpq_raw()`) signals `reintervals_error_time_points`
  (a `reintervals_error_data`), before any draw, when a lot has fewer than two
  distinct time points, naming the lots. Previously this was reported only for
  the first such lot found. lme4 fits such data, so the check is reachable
  through `ti_lmm()`.

## Bug fixes

* The advisory `boundary` flag (random-slope confidence, prediction and
  new-group-mean rows) depended on the time units: it compared the slope
  variance (y^2 per squared time unit) with the other variances. It is now
  judged on the response scale, the slope variance entering as its
  contribution `s2_1 * var(t)` (`var(t)` the population variance of the
  observed times), so it is unchanged by rescaling the response or time. Its
  value can change on fits near the `1e-4` threshold; nothing routes on it.

* **Dfs and intervals could depend on the measurement units** when variance
  components were below about 1.8e-5 (numerical derivatives used an absolute
  step there): the random-intercept PI / CI Satterthwaite dfs could fall back
  to containment and the random-slope verbs could stop with a `chol()` error on
  data in small units (e.g. an impurity in %, or the same data in g/mL instead
  of mg/mL). Derivatives are now taken on a scale-free parameterization (the
  variance parameters divided by the residual variance, back-transformed
  exactly), and every interval and df is scale-equivariant (tested from 1e-4
  to 1e4). In the scale-free parameters every non-zero parameter gets a
  relative step (numDeriv's `zero.tol` is set to 1e-12, so only an exactly
  zero, singular component gets the absolute step), so results are also
  unchanged when time is rescaled (e.g. months to years; tested from 1e-2 to
  1e2); and the positive-definiteness check of the covariance is made on its
  correlation form, which rescaling a parameter does not change. Where a
  parameter is below 1.8e-5 of the residual variance (a variance near the
  boundary), the covariance is now accepted only if the Hessian at half the
  step agrees: both positive definite, and every Satterthwaite df built on
  them (random slope: the mean-line, future-observation and between-group dfs
  at the call's time points; random intercept: the PI df and the CI df at the
  call's rows) within 1e-3 relative. Otherwise it is unavailable (random slope:
  `vc_cov_route = "unavailable"`), the default dfs fall back to containment
  with the reason "unstable numerical Hessian" in the df-fallback warning, and
  the `"broad-inference"` arms report `NA` limits. Fits with every parameter
  above 1.8e-5 of the residual variance are unaffected (their numbers move at
  floating-point level, about 1e-10 relative in the dfs). A failure of the
  numerical gradient of the mean-line variance now gives the containment df
  with the reason "numerical derivative failed" (and the df-fallback warning)
  instead of an error.

* The random-slope variance-parameter covariance could be unavailable on
  non-singular fits with a strong intercept-slope correlation (`|rho|` above
  about 0.86), because the numerical Hessian's first step (10% of each
  parameter) left the valid region, where the REML objective returns a
  sentinel value. The default random-slope CI / PI / new-group-mean dfs then
  fell back to containment and the `"broad-inference"` new-group-mean arm
  returned `NA` limits. When the default Hessian gives no valid covariance, it
  is now recomputed with a 0.1% first step, and the retry is accepted only when
  every evaluation stays in the valid region and the inverse is positive
  definite. The route is reported as `vc_cov_route` (`"default"`,
  `"retry-d1e-3"` or `"unavailable"`) in `attr(., "diagnostics")` of the
  random-slope `ci_lmm()`, `pi_lmm()`, `new_group_mean_lmm()` and
  `lmm_predict()`. Fits whose default Hessian was valid are unchanged.
  Singular fits never get a variance-parameter covariance, on either engine:
  when `lme4::isSingular()` is `TRUE` (or the check fails) the random-slope
  path, like the random-intercept path, computes no Hessian (the REML optimum
  is on the boundary, where the inverse observed information is not a
  covariance of the estimates). The route is `"unavailable"`, the default dfs
  fall back to containment, and the `"broad-inference"` arm returns `NA`
  limits; before, a singular random-slope fit got a covariance whenever the
  numerical Hessian happened to be positive definite. Why a covariance or a
  df is unavailable is now reported per row in `attr(., "diagnostics")`:
  `vc_cov_reason` on the random-slope rows (`"singular fit"`, `"singularity
  check failed; treated as singular"`, `"sentinel"`, `"unstable numerical
  Hessian"` or `"Hessian not PD or not invertible"`; `NA` when a covariance
  was obtained) and `df_fallback_reason` on the random-intercept CI, PI and
  new-group-mean rows (`NA` where the df did not fall back).
* The default random-slope new-group-mean interval (and the random-slope
  confidence and prediction intervals) could return `(-Inf, Inf)` when the
  Hessian-Satterthwaite df was a tiny positive number (for example 0.004 on a
  three-lot fit with an independent intercept/slope structure). Every default
  df rule now applies one rule: a Satterthwaite df that is non-finite **or
  below 1** falls back to the containment df, reported as
  `df_rule = "containment"` (`df_method = "containment"` for the
  random-intercept CI) with the reason `"df below 1"` in the
  `reintervals_warning_df_fallback` warning. This was already the rule for the
  random-intercept confidence interval; it is now also the rule for the
  random-slope CI / PI / new-group-mean dfs (which fell back only on a
  non-finite or non-positive df), the random-intercept prediction interval's
  Francq generalized-Satterthwaite df (which fell back only when the Hessian
  route was unusable) and the crossed new-group-mean Satterthwaite df (which
  fell back only on a non-finite or non-positive df). No fit in the package's
  tests or fixtures has a df below 1, so no pinned result changes.
* The `"broad-inference"` new-group-mean arms (`ngm_method`,
  `slope_ngm_method`) returned `NA` limits whenever the df was below 1. They
  now return `NA` (with `ngm_df_failed = TRUE`) only when the df cannot be
  computed or is non-finite or not positive, and otherwise use the fractional
  df in `qt()`.

* `ti_gpq_raw(sides = "lower")` / `"upper"` reported `sides = "two.sided"` in
  its output (the limits themselves were one-sided and correct). It now uses
  the same row schema as the `ti_lmm()` GPQ route and so reports the side, and
  its diagnostics now match that route's fields (adding `df` and
  `route_reason = "method"`).

* `ci_lmm(control = lmm_interval_control(ci_df = "kenward-roger"))` paired the
  Kenward-Roger df with the *unadjusted* standard error (from `vcov()`). The
  confidence interval now uses the KR-adjusted covariance
  (`pbkrtest::vcovAdj()`) for the standard error together with the KR df
  (`pbkrtest::Lb_ddf()`), computed once per call. The default Satterthwaite CI
  is unchanged.
* A random-slope confidence, prediction or new-group-mean interval built from a
  component container with a non-finite df used a normal quantile (the df was
  replaced by 1e6). It now uses the random-slope engine's own containment df,
  the rule its Satterthwaite fallback already uses (`B - 2` groups-based df for
  a random-slope model, or the residual df where the residual dominates the
  prediction interval), and reports `df_rule = "containment"`. Containers built
  by the package always carry a finite df, so no interval computed from a fitted
  model changes.
* The tolerance interval's MLS bound (`ems-mls`) now uses the
  Ting-Burdick-Graybill (1990) bound when an EMS coefficient is negative;
  Graybill-Wang is valid only for non-negative coefficients. The bound used is
  reported in `attr(., "diagnostics")$mls_bound`. No design in the package's
  checks produced a negative coefficient, so no current result changes.
* The crossed new-group-mean interval can no longer return `NaN` limits: when
  its Satterthwaite df is undefined it falls back to the containment df
  (`df_rule = "containment"`).
* The GPQ engine drops Monte-Carlo draws without a valid between-batch variance
  solution (reported as `n_na_draws` in the diagnostics) and errors if they
  exceed 1% of the draws.
* Fixed-slope `ci_lmm()`, `pi_lmm()` and `ti_lmm()` with an `NA` time value in
  `newdata` dropped that row and misaligned the remaining rows (the interval
  of one time point was reported against another). They now signal a
  `reintervals_error_input` naming the time variable.
* A factor time column in `newdata` was used as its level codes (so
  `factor(5)` was evaluated at time 1) by the random-slope `ti_lmm()` and the
  fixed-slope `new_group_mean_lmm()`. It is now a `reintervals_error_input`.

## Other changes

* **One-sided GPQ tolerance limits: documented limitation.** One-sided GPQ
  tolerance-interval rows (`ti_lmm()` random-slope route and
  `method = "gpq"` random-intercept route, `ti_gpq_raw()`) signal a
  `reintervals_warning_limitation`, once per call, (a) at a confidence level
  below 0.95, both targets (in validation, step 13, random-slope designs,
  one-sided GPQ limits at 90% confidence covered 0.879 to 0.90, median 0.894
  to 0.903), and (b) at 0.95 or above with `over = "group_mean"` (one-sided
  batch-mean limits at 95% covered 0.936 to 0.95, median 0.948). Two-sided GPQ
  intervals are not affected at these levels. Each GPQ tolerance row records
  `attr(., "diagnostics")$limitation_flag` (`"gpq_one_sided_lowconf"`,
  `"gpq_one_sided_group_mean"` or `NA`). The intervals are unchanged; the
  random-intercept GPQ route keeps its `reintervals_warning_unvalidated`.
* **`df_fallback_reason` on every random-slope coverage row.** The random-slope
  confidence-interval rows (every `slope_ci_df` arm) and new-group-mean rows
  (every `slope_ngm_method` arm) now carry `df_fallback_reason` in the
  diagnostics, as the random-slope prediction rows already did: why the
  default df fell back to containment (`"singular fit"`, `"singularity check
  failed; treated as singular"`, `"df below 1"`, `"non-finite df"`, or the
  `vc_cov_reason`), `NA` where it did not, and `NA` on the arms with no
  fallback. Every other output is unchanged.
* **`M` on every GPQ row.** The number of Monte-Carlo draws `M` is now in the
  diagnostics of the random-slope GPQ tolerance interval, the GPQ
  new-group-mean arm (`slope_ngm_method = "gpq"`) and `ti_gpq_raw()`, as on
  the random-intercept GPQ route. Every other output is unchanged.
* **Classed errors for invalid `lmm_interval_control()` choices.** An invalid
  value of `ci_df`, `slope_ci_df`, `pi_df`, `ngm_method`, `slope_ngm_method`,
  `slope_pi_method` or `gpq_two_sided` now signals a
  `reintervals_error_input` whose message names the argument and lists its
  allowed values (it was an unclassed `match.arg()` error). Valid values,
  partial matching and the defaults are unchanged.
* **Tests: no results from external commercial software.** The comparisons
  against SAS PROC MIXED runs (registry rows `sas-part1-*`, the pending
  `sas-part2-*` / `sas-part3-*`, the results file and its parser, and the SAS
  program generator in `data-raw/`) are removed. The simulated data sets they
  used are kept in `tests/testthat/fixtures/sim_designs/`; the tests that lost
  their only reference value now pin the package's own values (regression
  pins), and the Kenward-Roger check against pbkrtest is kept
  (`pbkrtest-kr-sim-fits`). The documentation no longer quotes SAS results.
* **Tests: independent references restored for the former pins.** The
  default Satterthwaite CI df and SE are compared with lmerTest
  (`contest1D()`, Suggests only) on ten non-singular random-intercept fits;
  the Kenward-Roger check against pbkrtest now also covers unbalanced one-way
  and fixed-slope fits; the random-intercept broad-inference and random-slope
  `reml` / `broad-inference` new-group-mean intervals are compared with an
  implementation written in the test file from McLean, Sanders & Stroup
  (1991) and Satterthwaite (1946). Registry rows `lmertest-satterthwaite-ci`,
  `pbkrtest-kr-unbalanced`, `indep-ngm-ri-broad`, `indep-ngm-rs-reml-broad`;
  the regression pins are kept. No package code changed.
* The GPQ tolerance interval warns (`reintervals_warning_limitation`) when the
  pooled within-lot residual sum of squares is zero (every lot's data on a
  straight line): the residual-variance pivot is then degenerate and the
  interval reflects only the between-lot variance.

* Every error and warning is now a classed condition: errors inherit from
  `reintervals_error` (`_input`, `_unsupported`, `_design`, `_method`, `_data`,
  `_groups`, `_time_points`, `_residual_df`, `_numerical`, `_dependency`,
  `_internal`) and warnings from `reintervals_warning` (`_fallback`,
  `_truncation`, `_singular`, `_limitation`). Message texts are unchanged.
* `pi_lmm()` now emits a `reintervals_warning_limitation` for a
  random-intercept model whose top random factor has five or fewer levels (the
  known limitation documented in `?pi_lmm` since 0.2.0), at most once per R
  session. The interval is unchanged.
* New diagnostics: `pi_df_fallback` (prediction interval), `df_rule` for the
  random-slope confidence, prediction and new-group-mean intervals, `mls_bound`
  (tolerance intervals), `n_na_draws` (GPQ).
* New warning class `reintervals_warning_extrapolation`: on fixed-slope and
  random-slope models every verb (`ci_lmm()`, `pi_lmm()`, `ti_lmm()`,
  `new_group_mean_lmm()`, `lmm_predict()`), and the raw engines
  `ti_gpq_raw()` and `ti_anova_raw()`, warn once per call when a requested
  `t0` lies outside the observed time range. The message names the
  out-of-range `t0` values and the range. The interval is unchanged.
* New warning class `reintervals_warning_df_fallback` (a child of
  `reintervals_warning_fallback`): emitted once per call
  when a default df rule falls back to the containment df (random-intercept CI
  Satterthwaite, random-intercept PI Francq df, crossed new-group-mean
  Satterthwaite, random-slope Hessian-Satterthwaite df of the CI, PI or
  new-group mean). It names the verb, the df used and the reason when known
  (e.g. singular fit). Validation arms and user-chosen containment arms do not
  emit it. The interval is unchanged.
* For fixed-slope and random-slope models, `newdata` must contain the time
  variable as a numeric, finite, non-missing column; otherwise every verb
  signals a `reintervals_error_input` naming the variable. Before, a missing
  column on the random-slope CI / PI / new-group mean failed with an unclassed
  base error, and a non-numeric, `NA` or `Inf` time value otherwise failed
  with an unclassed base error or a message that did not name the variable. (The silently wrong cases are under "Bug fixes"; the character
  column on the fixed-slope new-group mean under "Breaking changes".)

# reintervals 0.2.0

## New features
* `ti_lmm()` computes the Hoffman & Kringle (2005) two-sided tolerance
  interval for an intercept-only model with one random-intercept factor
  (`y ~ 1 + (1 | g)`), balanced or unbalanced. It is the slope-free branch of
  the `anova-mls` engine and is now the default for that design (previously
  `ems-mls`); `method = "ems-mls"` still gives the Francq et al. (2019)
  interval. The two coincide on balanced data with a non-singular fit and
  differ on singular or unbalanced data. `ti_anova_raw()` gives the same
  interval from raw data with `time = NULL`.
* New dataset `hk_assay`: the Hoffman & Kringle (2007) Table III assay data
  (6 runs x 3 replicates). Hoffman & Kringle (2005) Table 6 misprints run 5 /
  replicate 3 as 1.020 (correct: 0.998); Francq et al. (2019) Sections 5.2-5.3
  used the misprinted value.
* `lmm_interval_control()` gains `ci_df`, the confidence-interval df rule for
  random-intercept models: `"satterthwaite"` (default) or `"kenward-roger"`
  (requires pbkrtest and errors if it is not installed; random-slope models
  keep their own df and reject `"kenward-roger"`). The df source is reported
  in `attr(., "diagnostics")$df_method`.
* Random-slope `ti_lmm()` reports its routing reason per row in
  `attr(., "diagnostics")$route_reason` (`"none"`, `"singular"`,
  `"df_unavailable"`, `"group_mean_target"`, or `"method"`).
* One-sided intervals (`sides = "lower"` / `"upper"`) for `ci_lmm()`,
  `pi_lmm()`, and `new_group_mean_lmm()` on every engine: the two-sided
  interval's center, variance, and df with `t_{1 - alpha}`.
* One-sided tolerance limits: `anova-mls` (Hoffman-Kringle one-way and
  Montes fixed slope) uses Hoffman's (2010) one-sided limit; `ems-mls` and
  random-slope models use the package's extension of Hoffman's limit to the
  Francq and REML closed-form bounds, validated by simulation. For random
  slopes the one-sided default is `reml-mls` on every fit, singular ones
  included; `method = "gpq"` remains available but under-covered near the
  variance-component boundary in validation. `lmm_predict()` passes `sides`
  through.

## Known limitations
* With five or fewer levels of the top random factor and a dominant
  between-group variance (intraclass correlation near 0.9), the Francq et al.
  (2019) prediction interval for random-intercept models can under-cover (worst
  coverage observed in validation about 0.93 for a nominal 0.95), consistent
  with Francq et al. (2019, Fig. 2).

## Bug fixes
* EMS stratum degrees of freedom now account for the fixed effects: each
  stratum loses one df per fixed-effect column estimated in it (a column
  belongs to the coarsest random factor within whose levels it is constant,
  otherwise to the residual), so the residual df is N - rank([X Z]). This
  affects the EMS tolerance interval, the containment df, and the new-group
  mean (e.g. a between-group covariate gives a one-way top df of A - 2, not
  A - 1). Intercept-only designs are unchanged.
* When the observed-information prediction-interval df is not usable
  (singular fit, a covariance that is not positive definite, or a
  non-positive variance of the total variance), the random-intercept
  prediction interval now uses the stratum containment df (e.g. A - 1, or
  min(A, B) - 1 for crossed designs). The former bounded-Satterthwaite
  fallback gave df far above the stratum df on singular fits (8 to 70 against
  1 to 4) and under-covered; the N - p - q fallback is removed as well.
* The one-sided `anova-mls` tolerance limit used the two-sided formula with
  z_P, which is anti-conservative for P < 0.99 (Montes et al. 2019). It now
  uses Hoffman's (2010) one-sided limit.
* The random-intercept confidence-interval df no longer depends on which
  packages are installed. It used Kenward-Roger when pbkrtest was installed
  and otherwise the residual df N - p, because the Satterthwaite df was never
  given the variance-component covariance. It now uses Satterthwaite with that
  covariance computed.
* When the Satterthwaite df is unavailable (typically a singular fit), the
  confidence interval falls back to the SAS-style containment df (A - 1 for
  one-way and nested designs, min(A, B) - 1 for crossed ones) instead of the
  anti-conservative N - p.
* The random-intercept new-group-mean interval (`new_group_mean_lmm()`) was
  anti-conservative, with coverage as low as 0.80: it used the
  observable-variance df and collapsed when a variance component was truncated
  to zero. It now targets the true mean of a new innermost cell (a new group,
  B-within-A cell, or (a, b) cell; a new lot at `t0` for the fixed-slope
  design) with constructions validated by simulation: the unweighted lot-mean
  (Montes) construction for fixed-slope and intercept-only one-way designs
  (`anova-mls`), and an EMS construction with the top factor's df (nested) or
  a Satterthwaite df with `V_F` attached to the main effects (crossed)
  (`ems-mls`). `W`, `df`, and `df_rule` are reported in
  `attr(., "diagnostics")`. Designs without a closed-form decomposition now
  signal an error, and `method =` accepts only the design's own engine for this
  interval. On single-observation fixed-slope designs the
  new-group-mean interval no longer warns (the warning remains for the
  confidence and prediction intervals).
* A model whose single numeric covariate is constant within every group (a
  between-group covariate, e.g. `y ~ z + (1 | g)`) was misclassified as the
  fixed-slope stability design: `new_group_mean_lmm()` and `ti_lmm()` failed
  ("no within-lot time variation") and `ci_lmm()` / `pi_lmm()` warned. A
  covariate is now treated as the fixed-slope time only if it varies within
  at least one group; otherwise all four verbs follow the `ems-mls` path.
* The anti-conservative warning of `ci_lmm()` and `pi_lmm()` on fixed-slope
  designs now fires only when at least one lot has a single observation.
  Staggered designs whose lots all have two or more observations showed no
  shortfall in validation and no longer warn. The tolerance interval still
  uses `anova-mls` on every fixed-slope design that is not clean balanced.
* A variance-component covariance matrix (the inverse observed information)
  that is not positive definite is now treated as unavailable: the
  containment df is used and the random-slope two-sided tolerance interval
  routes to GPQ. Previously an indefinite matrix could yield arbitrarily large
  Satterthwaite df (anti-conservative), and the routing could differ between
  evaluation times of the same fit.
* Two-sided random-slope tolerance intervals switch to GPQ not only on a
  singular fit but also, at a given `t0`, when the Satterthwaite df of the
  closed form's between-group or residual variance is unavailable (its
  containment fallback fired). Such fits are non-singular but sit near the
  variance-component boundary.
* A model whose single non-intercept fixed term is a factor or logical
  covariate (e.g. `y ~ grp + (1 | lot)` with a 2-level `grp`) was misdetected as
  the fixed-slope stability design. When every lot had at least two
  observations, all four verbs failed ("arguments must have same length"). With
  single-observation lots, `ci_lmm()`, `pi_lmm()`, and `new_group_mean_lmm()`
  gave a spurious fixed-slope single-observation warning, and `ti_lmm()` was
  misrouted to `anova-mls` and failed. The fixed-slope design now requires
  that term to be a numeric model-frame column (mapped through the model's
  terms, not the coefficient name). Any other model follows the `ems-mls`
  path. Numeric-covariate routing is unchanged.
* `ci_lmm()` now reports `singular` from `lme4::isSingular()`, as `pi_lmm()`,
  `new_group_mean_lmm()`, and `ti_lmm()` already did. Previously the
  confidence interval always reported `singular = FALSE`.
* Internal documentation corrected to match the code: unbalanced nested and
  crossed designs are supported through EMS synthesis; off-catalog designs make
  `ti_lmm()` signal an error (they do not return `NA` bounds); the
  random-slope `Cov(phi)` is `NULL` only when the Hessian or its inverse fails;
  the `boundary` flag is an advisory flag on the point estimate (time-scale
  dependent) and does not drive the df fallback; `n_E` uses the target's content
  variance; the `true_value` tolerance-interval documentation and errors (in
  both the closed-form constructor and the EMS dispatch path) no longer claim
  to apply only to the fixed-slope design; the `Cov(theta)` and random-slope
  `Cov(phi)` Hessians (and the prediction-interval df built on them) are
  described as the REML log-likelihood on the variance-component scale, with
  beta at its GLS value, not a profiled REML; and printing a
  `vfun_components` with a `NULL` `Cov(phi)` now says "Hessian or its inverse
  unavailable" rather than "singular".
* The getting-started vignette's engine-routing description has been
  corrected (only the tolerance interval is rerouted on singular random-slope
  fits and on single-observation fixed-slope designs; `method =` overrides are
  limited to each design's allowed engines).

# reintervals 0.1.0

Initial release: confidence, prediction, tolerance, and
new-group-mean intervals for linear mixed models, through one uniform interface.

## Interval verbs
* `ci_lmm()` — confidence interval for the mean.
* `pi_lmm()` — prediction interval for a future observation.
* `ti_lmm()` — tolerance interval (content `P`); `over = "group_mean"` gives the
  batch-mean tolerance interval.
* `new_group_mean_lmm()` — the CInew / new-batch-mean interval (prediction for
  one new group's mean).
* `lmm_predict()` — a `predict()`-style wrapper over the four, by `interval =`.
* `lmm_interval_control()` — GPQ tuning (RNG `seed`, Monte-Carlo size `M`).
* `ti_anova_raw()`, `ti_gpq_raw()` — tier-2 direct engines from raw
  `(y, time, group)` vectors (no fitted model required).

## Engines (selected automatically; `method` overrides)
* `ems-mls` — expected-mean-square / modified large-sample, for random-intercept
  models (Francq, Lin & Hoyer 2019).
* `reml-mls` — REML closed form, for random-slope models.
* `anova-mls` — unweighted lot-mean ANOVA closed form, for unbalanced
  single-observation fixed-slope stability designs (Montes, Burdick & Leblond
  2019).
* `gpq` — generalized-pivotal-quantity Monte Carlo (Oliva-Aviles & Hauser 2025);
  used for the batch-mean target and for singular random-slope fits.

## Output
All verbs return a single tidy `lmm_interval` data frame with a uniform set of
columns (`type`, `estimate`, `lower`, `upper`, `level`, `P`, `sides`, `method`,
route-independent `design`, and `singular` / `boundary` / `vg_floored` flags),
the eval-point column(s) named by the model predictor(s), and engine-specific
diagnostics in `attr(., "diagnostics")`.

## Scope (v1)
General fixed effects for random-intercept models; a single continuous covariate for
random-slope models (the stability design). On single-observation / unbalanced
fixed-slope data the EMS confidence/prediction/new-group-mean intervals warn
that they may be anti-conservative; the tolerance interval is routed to the
validated ANOVA closed form there.

## Validation
Every shipped interval is checked against a published or external-source oracle
on continuous integration: the Francq, Lin & Hoyer (2019) worked numbers; the
Montes-Burdick-Leblond (2019) Appendix A tolerance interval (57.25, 75.41); the
random-slope intervals against committed reference values to ~1e-9; and the GPQ
content variance against the Oliva-Aviles & Hauser (2025) closed form
(tau^2 = 8.439).
