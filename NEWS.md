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
  (requires pbkrtest and errors if it is not installed). The df source is
  reported in `attr(., "diagnostics")$df_method`.
* Random-slope `ti_lmm()` reports why GPQ was used, per row, in
  `attr(., "diagnostics")$route_reason`.

## Bug fixes
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
  signal an error. On single-observation fixed-slope designs the
  new-group-mean interval no longer warns (the warning remains for the
  confidence and prediction intervals).
* Random-slope tolerance intervals switch to GPQ not only on a singular fit
  but also, at a given `t0`, when the Satterthwaite df of the closed form's
  between-group or residual variance is unavailable (its containment fallback
  fired). Such fits are non-singular but sit near the variance-component
  boundary.
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
