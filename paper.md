---
title: "reintervals: Confidence, prediction, tolerance, and new-group-mean intervals for linear mixed models"
tags:
  - R
  - linear mixed models
  - tolerance intervals
  - prediction intervals
  - stability analysis
authors:
  - name: Perceval Sondag
    orcid: 0009-0001-9183-6991
date: 30 June 2026
bibliography: paper.bib
---

# Summary

`reintervals` computes four kinds of statistical interval for a linear mixed
model fitted with `lme4` [@bates2015]: a confidence interval for the mean; a
prediction interval for a future observation; a tolerance interval covering a
proportion of the population; and a new-group-mean interval for the mean of one new, 
as-yet-unobserved group. The four interval types are produced by four validated 
computational engines that the package selects automatically from the fitted model's 
random-effects structure, and every interval is returned in one uniform tidy data 
frame so that confidence, prediction, tolerance, and new-group-mean results 
share a single column contract. A `predict()`-style wrapper and two tier-2 
functions that operate on raw `(response, time, group)` vectors complete the interface.

# Statement of need

Several widely used tools each cover part of the linear-mixed-model interval
problem well. `predict.merMod` in `lme4` [@bates2015], with `emmeans`
(estimated marginal means) and `ggeffects` (adjusted predictions / marginal
effects), give point predictions and confidence intervals for the mean or for
marginal effects. `merTools` provides
simulation-based prediction intervals for `merMod` fits. `ciTools` provides
confidence and prediction intervals for random-intercept models. The `tolerance`
package provides tolerance intervals for many distributions and some regression
settings, for simpler designs and without the mixed-model degrees-of-freedom
treatment.

What none of these provides is tolerance and new-group-mean intervals across
both random-intercept and random-slope (and fixed-slope stability) designs,
computed with the observed-information degrees of freedom of @francq2019, the
unweighted-lot-mean analysis-of-variance closed form of @montes2019, and the
generalized-pivotal-quantity method of @oliva2025 behind a single, uniform,
automatically dispatched interface. `reintervals` fills this gap. The
methods it implements are published; the contribution of this package is the
validated, unified implementation of those methods, not new methodology.

# Methods and validation

The package exposes four interval verbs (`ci_lmm`, `pi_lmm`, `ti_lmm`,
`new_group_mean_lmm`) over four engines, tagged by `method`: an
expected-mean-square / modified-large-sample construction for random-intercept
models with observed-information degrees of freedom (`ems-mls`, @francq2019); a
restricted-maximum-likelihood closed form for random-slope models (`reml-mls`);
an unweighted lot-mean analysis-of-variance closed form for unbalanced
single-observation fixed-slope stability designs (`anova-mls`, @montes2019); and
a generalized-pivotal-quantity Monte-Carlo engine for the batch-mean target and
for singular random-slope fits (`gpq`, @oliva2025). The engine is chosen from the
fitted model's structure, and may be overridden. Where more than one engine can 
produce a tolerance interval for the same design, the results are close but not 
identical: on balanced fixed-slope data the reml-mls and anova-mls engines share 
the variance-component decomposition exactly but apply different tolerance 
factors (a Graybill--Wang two-component construction [@graybill1980confidence] versus the Montes factor), 
so their tolerance intervals differ slightly. Both are validated against their 
respective published sources. The Satterthwaite-style degrees-of-freedom behavior 
near variance-component boundaries follows the analysis of @karl2026.

Validation is layered: First, the interval arithmetic is verified against 
committed reference values to a tight tolerance (~`1e-9`) on continuous 
integration; because re-fitting a model incurs `lme4`'s platform-dependent numerical 
jitter, these platform-portable checks run on committed, fixed extracted variance 
components rather than on re-fitted models. Second, the published point numbers 
are reproduced: the @montes2019 Appendix A tolerance interval `(57.25, 75.41)`, and the @oliva2025
closed-form content variance `tau^2 = 8.439`. Third, the full fit-to-interval
pipeline is checked exactly on a single reference platform; the only
cross-platform variation is the model-fitting step itself, which inherits
`lme4`'s numerical behavior and is not this package's code. 

For the GPQ engine we measured tolerance-interval coverage directly. Using the
published @oliva2025 data-generating process at `t0 = 20`, sweeping the
between-batch variance ratio and replicating `N_sim = 3000` datasets per
condition, the GPQ batch-mean `(P, gamma)`-tolerance interval attains its nominal
confidence at the paper's Monte-Carlo resolution (`M = 20000`): empirical
coverage is 0.944--0.955 at `gamma = 0.95`, with no condition more than about 1.3
binomial standard errors below nominal. Sub-0.94 coverage seen at low
Monte-Carlo resolution (`M ~ 4000`) is a result of too few generalized-pivotal
draws and resolves at high `M`; the slightest residual softness (~0.944) occurs
at the smallest variance ratio. We make no claim of uniform or conservative
coverage beyond this measurement.

# Scope

Version 1 supports general fixed effects for random-intercept models, and a
single time covariate for random-slope models (the stability design). The planned
extension is multi-covariate random-slope intervals.

# References
