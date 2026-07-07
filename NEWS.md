# reintervals 0.0.0.9000 (development)

First development version: confidence, prediction, tolerance, and
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
General fixed effects for random-intercept models; a single time covariate for
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
