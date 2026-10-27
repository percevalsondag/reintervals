# reintervals 0.1.1

## Bug fixes
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

## Pending
<!-- Placeholder for the two outstanding fixes (new-group-mean construction;
     CI df method), pending design decisions. -->

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
