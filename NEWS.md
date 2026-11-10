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

## New features

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
    constructions); `"broad-inference"` (SAS PROC MIXED broad-inference
    prediction of a new group: `W = V_F` plus every between-group variance
    component, Satterthwaite df on the variance scale; it reproduces SAS
    DDFM=SATTERTH on the package's SAS fixtures; NA limits and
    `ngm_df_failed = TRUE`, no fallback, when the df fails); `"hts"` (Higgins,
    Thompson & Spiegelhalter 2009, with the REML between-group variance;
    intercept-only one-way, `I >= 3`); `"classical"` (Hahn & Meeker 1991 on the
    group means; intercept-only one-way, balanced, where it equals the default).
  - `slope_ngm_method` (random slope): `"reml"` (default);
    `"broad-inference"` (SAS PROC MIXED broad-inference prediction of a new
    lot: the default's center and variance `V_F + V_G`, with the Satterthwaite
    df of the whole `V_F + V_G` instead of `V_G` alone; it reproduces SAS
    DDFM=SATTERTH on the package's SAS fixtures; NA limits and
    `ngm_df_failed = TRUE`, no fallback, when the df fails); `"gpq"` (a
    generalized prediction interval on the GPQ tolerance-interval pivots, with
    the GPQ seed contract, common random numbers across `t0` and Monte-Carlo
    SE; the tolerance-interval draws are unchanged); `"fv-anova"`
    (package-original lot-line interval on the per-lot OLS fitted values,
    `Ybar_fit +/- t_{K-1} S_L sqrt(1 + 1/K)`, conservative for the new-lot mean).
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
  is unchanged. SAS PROC MIXED `DDFM=KR` can differ on unbalanced data: in the
  package's checks against SAS it reported the Satterthwaite df with an adjusted
  standard error.
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

## Other changes

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
